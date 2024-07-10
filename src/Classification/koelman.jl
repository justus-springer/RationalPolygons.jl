function height_one_points(P :: RationalPolygon{T,N}) where {N,T <: Integer}
    V = vertex_matrix(P)
    Hs = affine_halfplanes(P)
    lines = [line(Hs[i] - 1) for i = 1 : N]
    res = LatticePoint{T}[]
    for i = 1 : N
        p = intersection_point(lines[mod(i-1,1:N)], lines[i])
        q = intersection_point(lines[mod(i+1,1:N)], lines[i])
        nv = normal_vector(lines[i])
        _, a, b = gcdx(-nv[1],-nv[2])
        x0 = LatticePoint{T}(V[1,i] + a, V[2,i] + b)
        append!(res, integral_points_on_line_segment_with_given_integral_point(p,q,x0))
    end
    unique!(res)
    filter!(p -> minimum([distance(p,H) for H ∈ Hs]) ≥ -1, res)
    return res
end

function single_point_extensions(Ps :: Vector{<:RationalPolygon{T}}) where {T <: Integer}

    out_dicts = Vector{Dict{Int, Set{<:RationalPolygon{T}}}}(undef, Threads.nthreads())
    for i = 1 : Threads.nthreads()
        out_dicts[i] = Dict{Int, Set{<:RationalPolygon{T}}}()
    end

    Threads.@threads for P ∈ Ps
        tid = Threads.threadid()
        N = number_of_vertices(P)
        bs = height_one_points(P)
        V = vertex_matrix(P)
        Hs = affine_halfplanes(P)
        for b ∈ bs
            contained_in_last_halfplane = contains_in_interior(b, Hs[N])
            local entry_point :: Int
            local exit_point :: Int
            for i = 1 : N
                contained_in_this_halfplane = contains_in_interior(b, Hs[i])
                if contained_in_last_halfplane && !contained_in_this_halfplane
                    exit_point = i
                elseif contained_in_this_halfplane && !contained_in_last_halfplane
                    entry_point = i
                end
                contained_in_last_halfplane = contained_in_this_halfplane
            end
            exit_point < entry_point && (exit_point += N)
            new_vertices = LatticePoint{T}[b]
            for i = entry_point : exit_point
                push!(new_vertices, V[:,mod(i,1:N)])
            end
            Q = RationalPolygon(new_vertices, one(T))
            (Q,is) = affine_normal_form_with_special_indices(Q)
            (1,0) ∈ is || (1,1) ∈ is || continue
            
            M = number_of_vertices(Q)
            !haskey(out_dicts[tid], M) && (out_dicts[tid][M] = Set{RationalPolygon{T,M,2M}}())
            push!(out_dicts[tid][M], Q)
        end
    end

    for i = 2 : length(out_dicts)
        mergewith!(union!, out_dicts[1], out_dicts[i])
        out_dicts[i] = Dict{Int,Set{<:RationalPolygon{T}}}()
    end

    return out_dicts[1]
end



abstract type KoelmanStorage{T <: Integer} end

mutable struct InMemoryKoelmanStorage{T <: Integer} <: KoelmanStorage{T}
    polygons :: Vector{Vector{RationalPolygon{T}}}
    total_count :: Int

    function InMemoryKoelmanStorage{T}() where {T <: Integer}
        polygons = [RationalPolygon{T}[], RationalPolygon{T}[], [RationalPolygon(SMatrix{2,3,T,6}(0,0,1,0,0,1), one(T))]]
        total_count = 1
        return new{T}(polygons, total_count)
    end
end

function classify_next_number_of_lattice_points(st :: InMemoryKoelmanStorage{T}) where {T <: Integer}
    Ps = last(st.polygons)
    new_Ps = collect(union(values(single_point_extensions(Ps))...))
    push!(st.polygons, new_Ps)
    new_count = length(new_Ps)
    st.total_count += new_count
    return new_count
end

function classify_polygons_by_number_of_lattice_points(st :: InMemoryKoelmanStorage{T}, max_number_of_lattice_points :: Int; logging :: Bool = false) where {T <: Integer}
    for l = length(st.polygons) + 1 : max_number_of_lattice_points
        new_count = classify_next_number_of_lattice_points(st)
        @info "[l = $l]. New polygons: $new_count. Total: $(st.total_count)"
    end
    return st.polygons
end

classify_polygons_by_number_of_lattice_points(max_number_of_lattice_points :: Int, T :: Type{<:Integer} = Int; logging :: Bool = false) =
classify_polygons_by_number_of_lattice_points(InMemoryKoelmanStorage{T}(), max_number_of_lattice_points; logging)


struct HDFKoelmanStoragePreferences{T <: Integer}
    swmr :: Bool
    maximum_number_of_vertices :: Int
    maximum_number_of_lattice_points :: Int
    block_size :: Int
    
    HDFKoelmanStoragePreferences{T}(
        swmr :: Bool = false,
        maximum_number_of_vertices :: Int = 100,
        maximum_number_of_lattice_points :: Int = 200,
        block_size :: Int = 10^6) where {T <: Integer} =
    new{T}(swmr, maximum_number_of_vertices, maximum_number_of_lattice_points, block_size)

end

mutable struct HDFKoelmanStorage{T <: Integer} <: KoelmanStorage{T}
    preferences :: HDFKoelmanStoragePreferences{T}
    file_path :: String
    group_path :: String
    last_completed_number_of_lattice_points :: Int
    total_count :: Int

    HDFKoelmanStorage{T}(preferences :: HDFKoelmanStoragePreferences{T},
                         file_path :: String,
                         group_path :: String = "/") where {T <: Integer} =
    new{T}(preferences, file_path, group_path, 0, 0)

    HDFKoelmanStorage{T}(file_path :: String,
                         group_path :: String = "/";
                         swmr :: Bool = false,
                         maximum_number_of_vertices :: Int = 100,
                         maximum_number_of_lattice_points :: Int = 200,
                         block_size :: Int = 10^6) where {T <: Integer} =
    HDFKoelmanStorage{T}(HDFKoelmanStoragePreferences{T}(swmr, maximum_number_of_vertices, maximum_number_of_lattice_points, block_size), file_path, group_path)

end

function initialize_koelman_storage(st :: HDFKoelmanStorage{T}) where {T <: Integer}

    f = h5open(st.file_path, "cw"; swmr = st.preferences.swmr)
    if !haskey(f, st.group_path)
        g = create_group(f, st.group_path)
    else
        g = f[st.group_path]
    end

    create_group(g, "l3")

    P = RationalPolygon(SMatrix{2,3,T,6}(0,0,1,0,0,1), 1)
    write_polygon_dataset(g, "l3/n3", [P])
    st.last_completed_number_of_lattice_points = 3

    g["numbers_of_polygons"] = zeros(Int, st.preferences.maximum_number_of_lattice_points, st.preferences.maximum_number_of_vertices) 
    numbers_of_polygons = g["numbers_of_polygons"]
    numbers_of_polygons[3,3] = 1

    close(f)

end

function classify_next_number_of_lattice_points(st :: HDFKoelmanStorage{T}; logging :: Bool = false) where {T <: Integer}
    f = h5open(st.file_path, "r+"; swmr = st.preferences.swmr)
    g = f[st.group_path]
    l = st.last_completed_number_of_lattice_points
    last_group = g["l$l"]
    current_group = create_group(g, "l$(l+1)")
    block_size = st.preferences.block_size

    for n_string ∈ keys(last_group)
        N = length(dataspace(last_group[n_string]))
        number_of_blocks = N ÷ block_size + 1

        for b = 1 : number_of_blocks
            elapsed_time = @elapsed begin

                I = (b-1)*block_size + 1 : min(b*block_size, N)
                new_count = 0
                n = parse(Int, n_string[2:end])
                Ps = read_polygon_dataset(one(T), last_group, n_string, I) 

                logging && @info "[l = $l, n = $n, block $b/$number_of_blocks]. Polygons to extend: $(length(Ps))"

                new_Ps = single_point_extensions(Ps)

                logging && @info "[l = $l, n = $n, block $b/$number_of_blocks]. Extension complete. New polygons: $(sum(length.(values(new_Ps))))"

                for (m, Qs) ∈ new_Ps
                    write_polygon_dataset(current_group, "n$m", collect(Qs))
                    new_count += length(Qs)
                    g["numbers_of_polygons"][l+1, m] += length(Qs)
                    st.total_count += length(Qs)
                end
                logging && @info "[l = $l, n = $n, block $b/$number_of_blocks]. Writeout complete. Running total: $(st.total_count)"
            end
            pps = round(new_count / elapsed_time, digits=2)
            elapsed_time = round(elapsed_time, digits=2)
            logging && @info "[l = $l, n = $n, block $b/$number_of_blocks]. Last step took $(elapsed_time) seconds. Averaging $pps polygons/s"
        end
    end

    st.last_completed_number_of_lattice_points = l+1

    close(f)

end

is_finished(st :: HDFKoelmanStorage{T}) where {T <: Integer} =
st.last_completed_number_of_lattice_points == st.preferences.maximum_number_of_lattice_points


function classify_polygons_by_number_of_lattice_points(st :: HDFKoelmanStorage{T}; logging :: Bool = false) where {T <: Integer}
    while !is_finished(st)
        classify_next_number_of_lattice_points(st; logging)
    end
    return st
end

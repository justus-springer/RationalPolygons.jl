function create_polygon_dataset(
        f :: Union{HDF5.File, HDF5.Group},
        path :: String,
        k :: Integer,
        N :: Int;
        T :: Type{<:Integer} = Int,
        chunk_size :: Int = 1000,
        deflate :: Int = 3)

    dspace = dataspace((0,),(-1,))
    dtype = datatype(NTuple{2N, T})
    dset = create_dataset(f, path, dtype, dspace; chunk = (chunk_size,), deflate)
    write_attribute(dset, "count", 0)
    write_attribute(dset, "rationality", k)
    return dset

end

function read_polygon_dataset(
        f :: Union{HDF5.File, HDF5.Group},
        path :: String)

    dset = HDF5.open_dataset(f, path)
    n = read_attribute(dset, "count")
    k = read_attribute(dset, "rationality")
    dtype = HDF5.datatype(dset)
    tuple_type = HDF5.get_jl_type(dtype)
    N = length(fieldnames(tuple_type)) ÷ 2
    T = eltype(tuple_type)
    Ms = HDF5.generic_read(dset, dtype, SMatrix{2,N,T,2N}, 1:n)
    close(dset)

    Ps = RationalPolygon{T,N,2N}[]
    for i = 1 : n
        push!(Ps, RationalPolygon(Ms[i], T(k)))
    end
    return Ps
end

function write_polygon_dataset(
        f :: Union{HDF5.File, HDF5.Group},
        path :: String,
        Ps :: Vector{RationalPolygon{T,N,M}}) where {N,M,T <: Integer}

    isempty(Ps) && return
    k = rationality(first(Ps))
    if !haskey(f, path)
        create_polygon_dataset(f, path, k, N; T)
    end

    dset = open_dataset(f, path)
    n = read_attribute(dset, "count")
    HDF5.set_extent_dims(dset, (n+length(Ps),))
    data = [vertex_matrix(P).data for P ∈ Ps]
    dset[n+1 : n+length(Ps)] = data
    attrs(dset)["count"] += length(Ps)
    close(dset)

end

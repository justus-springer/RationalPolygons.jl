
function export_polygons_to_h5(
        parent :: Union{HDF5.File, HDF5.Group},
        path :: AbstractString,
        Ps :: Vector{RationalPolygon{T,N,M}};
        chunk_size :: Int = 1000,
        deflate :: Int = 3) where {N, M, T <: Integer}

    n = length(Ps)
    As = vertex_matrix.(Ps)
    data = [A.data for A ∈ As]
    dspace = HDF5.dataspace((n,),(-1,))
    dtype = HDF5.datatype(data)
    dset = HDF5.create_dataset(parent, path, dtype, dspace; chunk = (chunk_size,), deflate)
    HDF5.write_dataset(dset, dtype, data)
    close(dset)
end

function append_polygons_to_h5(
        parent :: Union{HDF5.File, HDF5.Group},
        path :: AbstractString,
        Ps :: Vector{RationalPolygon{T,N,M}}) where {N, M, T <: Integer}

    dset = HDF5.open_dataset(parent, path)
    n = length(dset)
    HDF5.set_extent_dims(dset, (n+length(Ps),))
    As = vertex_matrix.(Ps)
    data = [A.data for A ∈ As]
    dset[n+1 : n+length(Ps)] = data
    close(dset)

end


function import_polygons_from_h5(
        k :: Integer,
        parent :: Union{HDF5.File, HDF5.Group},
        path :: AbstractString)

    dset = HDF5.open_dataset(parent, path)
    n = length(dset)
    dtype = HDF5.datatype(dset)
    tuple_type = HDF5.get_jl_type(dtype)
    N = length(fieldnames(tuple_type)) ÷ 2
    T = eltype(tuple_type)
    buf = Vector{SMatrix{2,N,T,2N}}(undef, n)
    HDF5.read_dataset(dset, dtype, buf)
    close(dset)

    Ps = RationalPolygon{T,N,2N}[]
    for i = 1 : n
        push!(Ps, RationalPolygon(buf[i], T(k)))
    end
    return Ps
end




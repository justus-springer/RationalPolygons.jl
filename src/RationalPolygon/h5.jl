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
    return dset

end


function read_polygon_dataset(k :: T,
        f :: Union{HDF5.File, HDF5.Group},
        path :: String, I...) where {T <: Integer}

    dset = HDF5.open_dataset(f, path)
    dtype = HDF5.datatype(dset)
    tuple_type = HDF5.get_jl_type(dtype)
    N = length(fieldnames(tuple_type)) ÷ 2
    Ps = [RationalPolygon(M, k) for M ∈ HDF5.generic_read(dset, dtype, SMatrix{2,N,T,2N}, I...)]
    close(dset)

    return Ps
end

function write_polygon_dataset(
        f :: Union{HDF5.File, HDF5.Group},
        path :: String,
        Ps :: Vector{RationalPolygon{T,N,M}}) where {N,M,T <: Integer}

    isempty(Ps) && return
    k = rationality(first(Ps))
    if !haskey(f, path)
        dset = create_polygon_dataset(f, path, k, N; T)
    else
        dset = open_dataset(f, path)
    end

    n = length(dataspace(dset))
    HDF5.set_extent_dims(dset, (n+length(Ps),))
    data = [vertex_matrix(P).data for P ∈ Ps]
    dset[n+1 : n+length(Ps)] = data
    close(dset)

end

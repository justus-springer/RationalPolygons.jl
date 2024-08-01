@doc raw"""
    create_polygon_dataset(f :: Union{HDF5.File, HDF5.Group}, path :: String, N :: Int)

Create an HDF5 dataset named `path` for storing rational polygons with `N`
vertices. The dataset will have an HDF5 compound datatype with `2*N` integers,
which are the entries of the vertex matrices stored in a column major layout.
This function takes three keyword arguments:

- `T :: Type{<:Integer}`: The integer type to be used, e.g. `Int64`, `Int32`, ... The default is `Int`.
- `chunk_size :: Int`: Chunk size for HDF5. The default is `1000`.
- `deflate :: Int`: Deflate parameter for HDF5. The default is `3`.

"""
function create_polygon_dataset(
        f :: Union{HDF5.File, HDF5.Group},
        path :: String,
        N :: Int;
        T :: Type{<:Integer} = Int,
        chunk_size :: Int = 1000,
        deflate :: Int = 3)

    dspace = dataspace((0,),(-1,))
    dtype = datatype(NTuple{2N, T})
    dset = create_dataset(f, path, dtype, dspace; chunk = (chunk_size,), deflate)
    return dset

end


@doc raw"""
    read_polygon_dataset(k :: T, f :: Union{HDF5.File, HDF5.Group}, path :: String, I...) where {T <: Integer}

Read from an HDF5 dataset containing `k`-rational polygons.

"""
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

@doc raw"""
    write_polygon_dataset(f :: Union{HDF5.File, HDF5.Group}, path :: String, Ps :: Vector{RationalPolygon{T,N,M}}) where {N,M,T <: Integer}

Write the polygons `Ps` to an HDF5 dataset named `path`. Creates the dataset if
it does not exist already. If it does exist, the data will be appended to it.

# Example

Write the reflexive lattice triangles to an HDF5 file:

```julia
julia> using RationalPolygons, StaticArrays, HDF5

julia> Vs = SMatrix{2,3}[[1 0 -2 ; 0 1 -1], [1 1 -3 ; 0 2 -4], [1 0 -1 ; 0 1 -1],[1 0 -3 ; 0 1 -2], [1 1 -2 ; 0 3 -3]];

julia> Ps = [RationalPolygon(V,1) for V ∈ Vs];

julia> f = h5open("/tmp/reflexive_triangles.h5", "cw")
🗂️ HDF5.File: (read-write) /tmp/reflexive_triangles.h5

julia> write_polygon_dataset(f, "my_triangles", Ps)

julia> close(f)
```

The resulting HDF5 file can be inspected using a command line tool such as
`h5dump`:

```shell
$ h5dump /tmp/reflexive_triangles.h5 
HDF5 "/tmp/reflexive_triangles.h5" {
GROUP "/" {
   DATASET "my_triangles" {
      DATATYPE  H5T_COMPOUND {
         H5T_STD_I64LE "1";
         H5T_STD_I64LE "2";
         H5T_STD_I64LE "3";
         H5T_STD_I64LE "4";
         H5T_STD_I64LE "5";
         H5T_STD_I64LE "6";
      }
      DATASPACE  SIMPLE { ( 5 ) / ( H5S_UNLIMITED ) }
      DATA {
      (0): {
            1,
            0,
            0,
            1,
            -2,
            -1
         },
      (1): {
            1,
            0,
            1,
            2,
            -3,
            -4
         },
      (2): {
            1,
            0,
            0,
            1,
            -1,
            -1
         },
      (3): {
            1,
            0,
            0,
            1,
            -3,
            -2
         },
      (4): {
            1,
            0,
            1,
            3,
            -2,
            -3
         }
      }
   }
}
}
```

We can read them back into RationalPolygons.jl at any time using
[`read_polygon_dataset`](@ref):

```julia
julia> using RationalPolygons, HDF5

julia> f = h5open("/tmp/reflexive_triangles.h5", "r")
🗂️ HDF5.File: (read-only) /tmp/reflexive_triangles.h5
└─ 🔢 my_triangles

julia> Ps = read_polygon_dataset(1, f, "my_triangles")
5-element Vector{RationalPolygon{Int64, 3, 6}}:
 Rational polygon of rationality 1 with 3 vertices.
 Rational polygon of rationality 1 with 3 vertices.
 Rational polygon of rationality 1 with 3 vertices.
 Rational polygon of rationality 1 with 3 vertices.
 Rational polygon of rationality 1 with 3 vertices.
```

"""
function write_polygon_dataset(
        f :: Union{HDF5.File, HDF5.Group},
        path :: String,
        Ps :: Vector{RationalPolygon{T,N,M}}) where {N,M,T <: Integer}

    isempty(Ps) && return
    if !haskey(f, path)
        dset = create_polygon_dataset(f, path, N; T)
    else
        dset = open_dataset(f, path)
    end

    n = length(dataspace(dset))
    HDF5.set_extent_dims(dset, (n+length(Ps),))
    data = [vertex_matrix(P).data for P ∈ Ps]
    dset[n+1 : n+length(Ps)] = data
    close(dset)

end

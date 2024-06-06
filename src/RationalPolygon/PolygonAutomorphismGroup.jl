@doc raw"""
    PolygonAutomorphismGroup

A struct holding information about the automorphism group of a rational
polygon. It has two fields `is_cyclic :: Bool` and `n :: Int`.

"""
struct PolygonAutomorphismGroup
    is_cyclic :: Bool
    n :: Int
end


@doc raw"""
    CyclicGroup(n :: Int)

Return the cyclic group of order `n`, as an `PolygonAutomorphismGroup`.

"""
CyclicGroup(n :: Int) = PolygonAutomorphismGroup(true, n)


@doc raw"""
    DihedralGroup(n :: Int)

Return the dihedral group of order `2n`, as an `PolygonAutomorphismGroup`.

"""
DihedralGroup(n :: Int) = PolygonAutomorphismGroup(false, n)


@doc raw"""
    is_cyclic(G :: PolygonAutomorphismGroup)

Check whether `G` is a cyclic group.

"""
is_cyclic(G :: PolygonAutomorphismGroup) = G.is_cyclic


@doc raw"""
    order(G :: PolygonAutomorphismGroup)

Return the order of the group `G`.

"""
order(G :: PolygonAutomorphismGroup) = is_cyclic(G) ? G.n : 2 * G.n

Base.show(io :: IO, G :: PolygonAutomorphismGroup) =
Base.print(io, (is_cyclic(G) ? "Z" : "D") * "$(G.n)")


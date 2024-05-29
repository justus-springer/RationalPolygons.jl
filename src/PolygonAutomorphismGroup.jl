struct PolygonAutomorphismGroup
    is_cyclic :: Bool
    n :: Int
end

CyclicGroup(n :: Int) = PolygonAutomorphismGroup(true, n)
DihedralGroup(n :: Int) = PolygonAutomorphismGroup(false, n)

is_cyclic(G :: PolygonAutomorphismGroup) = G.is_cyclic

order(G :: PolygonAutomorphismGroup) = is_cyclic(G) ? G.n : 2 * G.n

Base.show(io :: IO, G :: PolygonAutomorphismGroup) =
Base.print(io, (is_cyclic(G) ? "Z" : "D") * "$(G.n)")


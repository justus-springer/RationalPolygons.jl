using Test, RationalPolygons

@testset "Koelman's classification" begin
    
    max_number_of_lattice_points = 30

    st = InMemoryKoelmanStorage{Int}()
    classify_polygons_by_number_of_lattice_points(st, max_number_of_lattice_points)
    Pss = st.polygons

    a371917 = [0, 0, 1, 3, 6, 13, 21, 41, 67, 111, 175, 286, 419, 643, 938, 1370, 1939, 2779, 3819, 5293, 7191, 9752, 12991, 17321, 22641, 29687, 38533, 49796, 63621, 81300]

    for l = 1 : max_number_of_lattice_points
        @test length(Pss[l]) == a371917[l]
        @test all(P -> number_of_lattice_points(P) == l, Pss[l])
    end

end

@testset "Castryck's classification" begin

    max_genus = 20

    st = InMemoryCastryckStorage{Int}()
    classify_lattice_polygons_by_genus(st, max_genus)
    Pss = st.all_polygons

    a322343 = [16, 45, 120, 211, 403, 714, 1023, 1830, 2700, 3659, 6125, 8101, 11027, 17280, 21499, 28689, 43012, 52736, 68557, 97733]

    for i = 1 : max_genus
        @test length(Pss[i]) == a322343[i]
        @test all(P -> number_of_interior_lattice_points(P) == i, Pss[i])
    end

end

@testset "Brown & Kasprzyk's classification" begin

    max_side_length = 6

    square(m) = convex_hull(LatticePoint{Int}[(0,0),(m,0),(0,m),(m,m)])

    Pss = [subpolygons(square(m); use_affine_normal_form = true, only_equal_number_of_interior_lattice_points = false) for m = 1 : max_side_length];

    a374975 = [2, 15, 131, 1369, 13842, 129185]
    max_number_of_vertices = [4,6,8,9,10,12]
    number_of_vertex_maximizers = [1,1,1,1,15,2]

    @test all(m -> length(Pss[m]) - length(Pss[m-1]) == a374975[m], 2 : max_side_length)
    @test all(m -> maximum(number_of_vertices, Pss[m]) == max_number_of_vertices[m], 1 : max_side_length)
    @test all(m -> length(filter(P -> number_of_vertices(P) == max_number_of_vertices[m], Pss[m])) == number_of_vertex_maximizers[m], 2 : max_side_length)

end

@testset "Maximal polygons in R x [-1,1]" begin
    expected_numbers = [1 2 4 5 6 7 8 9 10 11 12;
                        4 9 13 18 22 26 30 34 38 42 46;
                        12 26 41 54 68 81 94 107 120 133 146;
                        24 57 86 117 145 174 203 231 259 288 316]

    for k = 1 : 4, i = 0 : 10
        Ps = classify_maximal_polygons_m1p1(k,i)
        @test length(Ps) == expected_numbers[k,i+1]
        @test all(P -> minimal_number_of_interior_integral_lines(P) == 1, Ps)
        @test all(is_maximal, Ps)
    end
end

@testset "Maximal polygons with no interior lattice points" begin
    expected_numbers = [1,4,14,39,134,299]
    for k = 1 : 6
        Ps = classify_maximal_lattice_free_polygons(k)
        @test length(Ps) == expected_numbers[k]
        @test all(P -> number_of_interior_lattice_points(P) == 0, Ps)
        @test all(is_maximal, Ps)
    end
end

@testset "Maximal polygons with one interior latticie point" begin
    expected_numbers = [3,10,39]
    for k = 1 : 3
        Ps = classify_maximal_polygons_genus_one(k)
        @test length(Ps) == expected_numbers[k]
        @test all(P -> number_of_interior_lattice_points(P) == 1, Ps)
        @test all(is_maximal, Ps)
    end
end

@testset "Polygons with one interior lattice point" begin
    expected_numbers = [16,5145,924042]
    for k = 1 : 3
        Ps = classify_polygons_genus_one(k)
        @test length(Ps) == expected_numbers[k]
        @test all(P -> number_of_interior_lattice_points(P) == 1, Ps)
    end
end

@testset "LDP Polygons with one interior lattice point" begin
    expected_numbers = [16,505,48032]
    for k = 1 : 3
        Ps = classify_polygons_genus_one(k; primitive=true)
        @test length(Ps) == expected_numbers[k]
        @test all(P -> number_of_interior_lattice_points(P) == 1, Ps)
    end
end


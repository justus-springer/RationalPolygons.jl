export modified_unit_fraction_solutions
export classify_quadrilaterals_by_gorenstein_index
export get_degree_matrix_solutions

@doc raw"""
    modified_unit_fraction_solutions(r :: T, s :: T, c :: T, d :: T) where {T <: Integer}
    
Return all integral solutions $(x,y) ∈ \mathbb{Z}^2_{\geq 1}$ to the
equation $r/s = 1/x + 1/y - c/(d*y)$. Returns an error if there are
infinitely many solutions, which happens if and only if $c = d$.

"""
function modified_unit_fraction_solutions(r :: T, s :: T, c :: T, d :: T) where {T <: Integer}
    c == d && error("infinitely many solutions")
    res = Tuple{T,T}[]

    # x ≤ y
    lower = 1
    upper = floor(T, 2*s // r)
    for x = lower : upper
        r*x != s || continue
        (s*x*(d-c)) % (d*(r*x-s)) == 0 || continue
        y = (s*x*(d-c)) ÷ (d*(r*x-s))
        x ≤ y || continue
        push!(res, (x,y))
    end

    # y < x
    for y = lower : upper
        d*r*y+c*s != d*s || continue
        (s*d*y) % (d*r*y+c*s-d*s) == 0 || continue
        x = (s*d*y) ÷ (d*r*y+c*s-d*s)
        y < x || continue
        push!(res, (x,y))
    end

    return res
end

modified_unit_fraction_solutions(q :: Rational{T}, c :: T, d :: T) where {T <: Integer} =
modified_unit_fraction_solutions(numerator(q), denominator(q), c, d)

function is_almost_free_grading_matrix(Q :: SMatrix{2,4,T}) where {T <: Integer}
    gcd(det(Q[:,1],Q[:,2]), det(Q[:,2],Q[:,3]), det(Q[:,3],Q[:,1])) == 1 || return false
    gcd(det(Q[:,1],Q[:,2]), det(Q[:,2],Q[:,4]), det(Q[:,4],Q[:,1])) == 1 || return false
    gcd(det(Q[:,1],Q[:,3]), det(Q[:,3],Q[:,4]), det(Q[:,4],Q[:,1])) == 1 || return false
    gcd(det(Q[:,2],Q[:,3]), det(Q[:,3],Q[:,4]), det(Q[:,4],Q[:,2])) == 1 || return false
    return true
end

function get_degree_matrix_solutions(ι :: T) where {T <: Integer}

    res = SMatrix{2,4,T}[]

    for a11 = ι+1 : 3ι, 
        a12 = 1 : a11-1,
        (a32, a21) in modified_unit_fraction_solutions(1 // ι - 1 // a11, a12, a11),
        a31 = 1 : a32-1

        a11*a32-ι*a11-ι*a32 ≠ 0 || continue
        ι*a11*(a32-a31) % (a11*a32-ι*a11-ι*a32) == 0 || continue
        a22 = ι*a11*(a32-a31) ÷ (a11*a32-ι*a11-ι*a32)

        a12*a22-ι*a12-ι*a22+ι*a21 ≠ 0 || continue
        ι*a12*a22 % (a12*a22-ι*a12-ι*a22+ι*a21) == 0 || continue
        a41 = ι*a12*a22 ÷ (a12*a22-ι*a12-ι*a22+ι*a21)
        a41 > 0 || continue

        (ι*a22*a41+ι*a11*a41+ι*a11*a22-a11*a22*a41) % (ι*a22) == 0 || continue
        a42 = (ι*a22*a41+ι*a11*a41+ι*a11*a22-a11*a22*a41) ÷ (ι*a22)
        a41 > 0 || continue

        G = [ι-a11 ι-a12 ι ι ; ι ι-a21 ι-a22 ι ; ι ι ι-a31 ι-a32 ; ι-a42 ι ι ι-a41]
        Q = SMatrix{2,4,T}(Matrix(kernel(matrix(ZZ, G))))
        is_almost_free_grading_matrix(Q) || continue
        push!(res, Q)

    end

    for a11 = ι+1 : 2ι
        a12 = a11
        (ι*a11) % (a11-ι) == 0 || continue
        a31 = (ι*a11) ÷ (a11-ι)
        a32 = a31
        for a42 = 1 : a31, a21 = 1 : a31

            a11*a21+ι*a12-ι*a11-ι*a21 ≠ 0 || continue
            (ι*a11*a21+ι*a12*a42-ι*a21*a42) % (a11*a21+ι*a12-ι*a11-ι*a21) == 0 || continue 
            a41 = (ι*a11*a21+ι*a12*a42-ι*a21*a42) ÷ (a11*a21+ι*a12-ι*a11-ι*a21)
            a41 > 0 || continue

            a11*a41+ι*a42-ι*a11-ι*a41 ≠ 0 || continue
            (ι*a11*a41) % (a11*a41+ι*a42-ι*a11-ι*a41) == 0 || continue 
            a22 = (ι*a11*a41) ÷ (a11*a41+ι*a42-ι*a11-ι*a41)
            a22 > 0 || continue

            G = [ι-a11 ι-a12 ι ι ; ι ι-a21 ι-a22 ι ; ι ι ι-a31 ι-a32 ; ι-a42 ι ι ι-a41]
            Q = SMatrix{2,4,T}(Matrix(kernel(matrix(ZZ, G))))
            is_almost_free_grading_matrix(Q) || continue
            push!(res, Q)

        end
    end

    return res

end

divisors(n :: Integer) = filter(k -> n % k == 0, 1 : abs(n))

# Returns all toric surface with given prescribed gorenstein index
# and free part of the degree matrix.
function classify_quadrilaterals_by_gorenstein_index(ι :: T) where {T <: Integer}

    result = Set{RationalPolygon{T,4,8}}()

    for Q in get_degree_matrix_solutions(ι)

        mu12 = det(Q[:,3],Q[:,4])
        mu23 = det(Q[:,1],Q[:,4])
        mu34 = det(Q[:,1],Q[:,2])
        mu14 = det(Q[:,2],Q[:,3])
        mu13 = det(Q[:,2],Q[:,4])
        mu24 = det(Q[:,1],Q[:,3])

        for ιr in divisors(ι)
            ιr*(mu12 - mu24 - mu14) % mu14 == 0 || continue
            q = ιr*(mu12 - mu24 - mu14) ÷ mu14

            for c in divisors(q), ι1 in divisors(ι), α in 0 : (ι1*c)-1
                β = ι1*c
                gcd(β, 1-α) == c || continue
                gcd(β, α) == 1 || continue

                v1 = LatticePoint{T}(1,0)
                v2 = LatticePoint{T}(α,β)
                gorenstein_index(v1, v2) == ι1 || continue

                (mu23 + α*mu13) % mu12 == 0 || continue
                (β*mu13) % mu12 == 0 || continue
                (mu24 + α*mu14) % mu12 == 0 || continue
                (β*mu14) % mu12 == 0 || continue

                v3 = LatticePoint{T}(-(mu23 + α*mu13) ÷ mu12, -(β*mu13) ÷ mu12)
                gcd(v3[1], v3[2]) == 1 || continue
                v4 = LatticePoint{T}((mu24 + α*mu14) ÷ mu12, (β*mu14) ÷ mu12)
                gcd(v4[1], v4[2]) == 1 || continue
                gorenstein_index(v4, v1) == ιr || continue

                P = convex_hull([v1, v2, v3, v4])
                gorenstein_index(P) == ι || continue

                push!(result, unimodular_normal_form(P))

            end
        end
    end

    return result
end

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

function gorenstein_coefficients_to_degree_matrix_minors(ι :: T,
        A :: SMatrix{4,2,T,8}) where {T <: Integer}

    # The gorenstein matrix
    G = T[ι-A[1,1] ι-A[1,2] ι ι ; ι ι-A[2,1] ι-A[2,2] ι ; ι ι ι-A[3,1] ι-A[3,2] ; ι-A[4,2] ι ι ι-A[4,1]]
    U = hnfc(G).U'

    # The lower two rows of U are now a lattice basis of the kernel of G,
    # i.e. can be viewed as the free part of the degree matrix Q.
    m12 = U[3,3] * U[4,4] - U[4,3] * U[3,4]
    m13 = U[3,2] * U[4,4] - U[4,2] * U[3,4]
    m14 = U[3,2] * U[4,3] - U[4,2] * U[3,3]
    m23 = U[3,1] * U[4,4] - U[4,1] * U[3,4]
    m24 = U[3,1] * U[4,3] - U[4,1] * U[3,3]
    m34 = U[3,1] * U[4,2] - U[4,1] * U[3,2]

    return (m12, m13, m14, m23, m24, m34)

end


function get_gorenstein_coefficients(ι :: T) where {T <: Integer}

    res = SMatrix{4,2,T}[]

    # Case a11 < a12
    for a12 = ι+1 : 3ι-1, a11 = 1 : a12-1
        for (a31,a42) in modified_unit_fraction_solutions(1//ι - 1//a12, a11, a12)
            a12 ≤ a31 || continue
            for a32 = 1 : a31-1
                a41, a41rem = divrem(a42*a12*(a32-a31), a31*(a11-a12))
                a41rem == 0 || continue
                a12 ≤ max(a41,a42) || continue

                a11*a41 + a32*a42 - a11*a32 ≠ 0 || continue

                a21, a21rem = divrem(a12*a32*a42, a11*a41 + a32*a42 - a11*a32)
                a21rem == 0 || continue

                a22, a22rem = divrem(a11*a21*(a31-a32), a32*(a12-a11))
                a22rem == 0 || continue

                a12 ≤ max(a21,a22) || continue

                A = SMatrix{4,2,T}(a11,a21,a31,a41,a12,a22,a32,a42)
                all(a -> a > 0, A) || continue
                push!(res, A)
                
            end
        end
    end

    # Case a11 = a12 (hence also a31 = a32)
    for a11 = ι+1 : 3ι-1
        a31, a31rem = divrem(ι*a11, a11-ι)
        a31rem == 0 || continue
        a11 ≤ a31 || continue

        # Subcase a22 ≤ a21.
        # Then we must have a22 ≤ a31 and a41 ≤ a31.
        for a22 = 1 : a31, a41 = 1 : a31
            # Now we can compute a21 and a42
            a21, a21rem = divrem(a11*a22*(a41-a31), a31*a41)
            a21 = a11 - a21
            a21rem == 0 || continue
            a22 ≤ a21 || continue

            a42, a42rem = divrem(a21*a41, a22)
            a42rem == 0 || continue
            a41 ≤ a42 || continue

            A = SMatrix{4,2,T}(a11,a21,a31,a41,a11,a22,a31,a42)
            all(a -> a > 0, A) || continue
            push!(res, A)

        end

        # Subcase a21 < a22.
        # Then min(a21,a42) < 2ι. By symmetry, we may assume
        # a21 = min(a21,a42) < 2ι.
        for a21 = 1 : 2ι-1
            # Subsubcase a22 ≤ a31.
            for a22 = 1 : a31
                a11 ≤ a22 || continue

                a11*a22+a21*a31-a11*a31 ≠ 0 || continue
                a41, a41rem = divrem(a11*a22*a31, a11*a22+a21*a31-a11*a31)
                a41rem == 0 || continue

                a42, a42rem = divrem(a21*a41, a22)
                a42rem == 0 || continue
                a11 ≤ max(a41,a42)

                A = SMatrix{4,2,T}(a11,a21,a31,a41,a11,a22,a31,a42)
                all(a -> a > 0, A) || continue
                push!(res, A)

            end

            # Subsubcase a31 < a22. Then also a42 < a11
            for a42 = 1 : a11-1
                a22, a22rem = divrem(a31*(a11*a21+a11*a42-a21*a42), a11*a42)
                a22rem == 0 || continue
                a21 < a22 || continue
                a31 < a22 || continue
                a11 ≤ a22 || continue

                a41, a41rem = divrem(a22*a31*(a11-a42), a11*(a22-a31))
                a41rem == 0 || continue
                a42 < a41 || continue
                a11 ≤ a41 || continue

                A = SMatrix{4,2,T}(a11,a21,a31,a41,a11,a22,a31,a42)
                all(a -> a > 0, A) || continue
                push!(res, A)
            end
        end

    end

    return res

end


function get_degree_matrix_minors(ι :: T) where {T <: Integer}

    res = SVector{6,T}[]

    for A ∈ get_gorenstein_coefficients(ι)
        m12, m13, m14, m23, m24, m34 = gorenstein_coefficients_to_degree_matrix_minors(ι,A)

        # check necessary condition for almost freeness
        gcd(m34, m14, m24) == 1 || continue
        gcd(m34, m13, m23) == 1 || continue
        gcd(m24, m12, m23) == 1 || continue
        gcd(m14, m12, m13) == 1 || continue

        push!(res, SVector{6,T}(m12, m13, m14, m23, m24, m34))
    end

    return res

end

divisors(n :: Integer) = filter(k -> n % k == 0, 1 : abs(n))

# Returns all toric surface with given prescribed gorenstein index
# and free part of the degree matrix.
function classify_quadrilaterals_by_gorenstein_index(ι :: T) where {T <: Integer}

    result = Set{RationalPolygon{T,4,8}}()

    for (m12, m13, m14, m23, m24, m34) in get_degree_matrix_minors(ι)
        
        for ι4 in divisors(ι)
            ι4*(m12 - m24 - m14) % m14 == 0 || continue
            q = ι4*(m12 - m24 - m14) ÷ m14

            for c1 in divisors(q)
                b = q ÷ c1 # b = ι1*d4 + ι4*d1
                for ι1 in divisors(ι)
                    b % gcd(ι1, ι4) == 0 || continue
                    y2 = ι1*c1
                    for d1 = ceil(T, -1//c1) : floor(T, ι1 - 1//c1)
                        # validate that ι1 is the correct local gorenstein index
                        gcd(d1,ι1) == 1 || continue

                        (b - ι4*d1) % ι1 == 0 || continue
                        d4 = (b - ι4*d1) ÷ ι1
                        # validate that ι4 is the correct local gorenstein index
                        gcd(d4,ι4) == 1 || continue

                        x2 = 1 + c1*d1
                        gcd(x2, y2) == 1 || continue

                        (m23 + x2*m13) % m12 == 0 || continue
                        x3 = -(m23 + x2*m13) ÷ m12
                        (y2*m13) % m12 == 0 || continue
                        y3 = -(y2*m13) ÷ m12
                        gcd(x3, y3) == 1 || continue

                        (m24 + x2*m14) % m12 == 0 || continue
                        x4 = (m24 + x2*m14) ÷ m12
                        (y2*m14) % m12 == 0 || continue
                        y4 = (y2*m14) ÷ m12
                        gcd(x4, y4) == 1 || continue
                        
                        P = RationalPolygon(SMatrix{2,4,T}(1,0,x2,y2,x3,y3,x4,y4),one(T))
                        gorenstein_index(P) == ι || continue

                        push!(result, unimodular_normal_form(P))
                    end
                end
            end
        end
    end

    return result
end

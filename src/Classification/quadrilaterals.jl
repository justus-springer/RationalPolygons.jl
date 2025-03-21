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
        a11 :: T, a12 :: T, a21 :: T, a22 :: T,
        a31 :: T, a32 :: T, a41 :: T, a42 :: T) where {T <: Integer}

    # The gorenstein matrix
    G = T[ι ι-a21 ι-a22 ι ; ι ι ι-a31 ι-a32 ; ι-a42 ι ι ι-a41 ; ι-a11 ι-a12 ι ι]

    # An ad hoc way of computing the kernel of G, given that we already
    # know that it is two-dimensional.
    U = T[1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1]
    zero_column_after!(G, U, 1, 1)
    if iszero(G[2,2])
        i = findfirst(i -> !iszero(G[i,2]), 3 : 4)
        swaprows!(G, 2, 2 + i)
        swaprows!(U, 2, 2 + i)
    end
    zero_column_after!(G, U, 2, 2)
    
    # The lower two rows of U are now a lattice basis of the kernel of G,
    # i.e. can be viewed as the free part of the degree matrix Q.
    
    # We safeguard this computation by using big integers, since this seems to
    # be the place where overflow happens first. The resulting minors however
    # generally not very big, so we can switch to a fixed-size integer type
    # again afterwards.
    m12 = T(big(U[3,3]) * big(U[4,4]) - big(U[4,3]) * big(U[3,4]))
    m13 = T(big(U[3,2]) * big(U[4,4]) - big(U[4,2]) * big(U[3,4]))
    m14 = T(big(U[3,2]) * big(U[4,3]) - big(U[4,2]) * big(U[3,3]))
    m23 = T(big(U[3,1]) * big(U[4,4]) - big(U[4,1]) * big(U[3,4]))
    m24 = T(big(U[3,1]) * big(U[4,3]) - big(U[4,1]) * big(U[3,3]))
    m34 = T(big(U[3,1]) * big(U[4,2]) - big(U[4,1]) * big(U[3,2]))

    return (m12, m13, m14, m23, m24, m34)

end


function get_degree_matrix_minors(ι :: T) where {T <: Integer}

    res = SVector{6,T}[]

    for a11 = ι+1 : 3ι, a12 = 1 : a11-1
        for (a32, a21) in modified_unit_fraction_solutions(1 // ι - 1 // a11, a12, a11),
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

            m12, m13, m14, m23, m24, m34 = gorenstein_coefficients_to_degree_matrix_minors(ι,a11,a12,a21,a22,a31,a32,a41,a42)

            # check necessary condition for almost freeness
            gcd(m34, m14, m24) == 1 || continue
            gcd(m34, m13, m23) == 1 || continue
            gcd(m24, m12, m23) == 1 || continue
            gcd(m14, m12, m13) == 1 || continue

            push!(res, SVector{6,T}(m12, m13, m14, m23, m24, m34))
        end
    end

    # Treat a11 = a12 seperately
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

            m12, m13, m14, m23, m24, m34 = gorenstein_coefficients_to_degree_matrix_minors(ι,a11,a12,a21,a22,a31,a32,a41,a42)

            # check necessary condition for almost freeness
            gcd(m34, m14, m24) == 1 || continue
            gcd(m34, m13, m23) == 1 || continue
            gcd(m24, m12, m23) == 1 || continue
            gcd(m14, m12, m13) == 1 || continue

            push!(res, SVector{6,T}(m12, m13, m14, m23, m24, m34))

        end
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

            for c in divisors(q)
                b = q ÷ c
                for ι1 in divisors(ι)
                    b % gcd(ι1, ι4) == 0 || continue
                    y2 = ι1*c
                    for d = ceil(T, -1//c) : floor(T, ι1 - 1//c)
                        # checks if ι1 is really the local gorenstein
                        # index of v1 and v2
                        gcd(d,ι1) == 1 || continue

                        x2 = 1 + c*d
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
                        
                        # check if ι4 is really the local gorenstein
                        # index of v4 and v1
                        ι4*(1-x4) % y4 == 0 || continue
                        gcd(ι4, ι4*(1-x4) ÷ y4) == 1 || continue

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

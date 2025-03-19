export is_almost_free_grading_matrix
export get_gorenstein_coefficient_solutions
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

function get_gorenstein_coefficient_solutions(ι :: T) where {T <: Integer}

    res = SMatrix{4,2,T}[]

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

        push!(res, SMatrix{4,2,T}(a11,a21,a31,a41,a12,a22,a32,a42))

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

            push!(res, SMatrix{4,2,T}(a11,a21,a31,a41,a12,a22,a32,a42))

        end
    end

    return res

end

function get_degree_matrix_solutions(ι :: T) where {T <: Integer}

    res = SMatrix{2,4,T}[]

    for A ∈ get_gorenstein_coefficient_solutions(ι)
        G = SMatrix{4,4,T}([ι-A[1,1] ι-A[1,2] ι ι ; ι ι-A[2,1] ι-A[2,2] ι ; ι ι ι-A[3,1] ι-A[3,2] ; ι-A[4,2] ι ι ι-A[4,1]])
        U = hnfr(G).U
        Q = SMatrix{2,4,T}(U[3:4,:])
        is_almost_free_grading_matrix(Q) || continue
        push!(res, Q)
    end

    return res

end

divisors(n :: Integer) = filter(k -> n % k == 0, 1 : abs(n))

# Returns all toric surface with given prescribed gorenstein index
# and free part of the degree matrix.
function classify_quadrilaterals_by_gorenstein_index(ι :: T) where {T <: Integer}

    result = Set{RationalPolygon{T,4,8}}()

    for Q in get_degree_matrix_solutions(ι)
        w11,w12,w21,w22,w31,w32,w41,w42 = Q

        mu12 = w31 * w42 - w32 * w41
        mu23 = w11 * w42 - w12 * w41
        mu34 = w11 * w22 - w12 * w21
        mu14 = w21 * w32 - w22 * w31
        mu13 = w21 * w42 - w22 * w41
        mu24 = w11 * w32 - w12 * w31
        
        for ι4 in divisors(ι)
            ι4*(mu12 - mu24 - mu14) % mu14 == 0 || continue
            q = ι4*(mu12 - mu24 - mu14) ÷ mu14

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

                        (mu23 + x2*mu13) % mu12 == 0 || continue
                        x3 = -(mu23 + x2*mu13) ÷ mu12
                        (y2*mu13) % mu12 == 0 || continue
                        y3 = -(y2*mu13) ÷ mu12
                        gcd(x3, y3) == 1 || continue

                        (mu24 + x2*mu14) % mu12 == 0 || continue
                        x4 = (mu24 + x2*mu14) ÷ mu12
                        (y2*mu14) % mu12 == 0 || continue
                        y4 = (y2*mu14) ÷ mu12
                        gcd(x4, y4) == 1 || continue
                        
                        # check if ι4 is really the local gorenstein
                        # index of v4 and v1
                        ι4*(1-x4) % y4 == 0 || continue
                        gcd(ι4, ι4*(1-x4) ÷ y4) == 1 || continue

                        P = RationalPolygon(SMatrix{2,4,T}(1,0,x2,y2,x3,y3,x4,y4),1)
                        gorenstein_index(P) == ι || continue

                        push!(result, unimodular_normal_form(P))
                    end
                end
            end
        end
    end

    return result
end

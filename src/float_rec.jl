
using LinearAlgebra
using SparseArrays
using Nemo
using Distributions


# "Rational numbers with bound denominators are discrete!"

#------------------------------------------------------------------------------

# for the given float `y` and natural `M`
# computes rational fraction `p/q` approximating `y` such that
#   |p/q - y| < 1/(2M²)
#       and
#   q ≤ M
# i.e the fraction p/q with the simplest denominator
function float_reconstruction(y::T, M) where {T<:AbstractFloat}
    # the function implements algorithm derived from Corollary 6.3a
    # in Theory of Linear and Integer Programming, Schrijver, 1986
    # the notation is taken from there

    ε = eps(T)

    αβ = [abs(y), 1]
    A = [BigInt[1, 0], BigInt[0, 1]]

    k = 1
    # standard Extended Euclidean algorithm until q < M
    while true
        # determine the order of columns
        i, j = iseven(k) ? (1, 2) : (2, 1)

        abs(A[j][1]) > M && break

        k += 1

        # in case the float `y` is rational exactly (w.r.t to ε)
        isapprox(αβ[j], T(0), atol=ε, rtol=0) && break

        # EEA step
        t = floor(BigInt, αβ[i] / αβ[j])
        A[i] -= t * A[j]
        αβ[i] -= t * αβ[j]
    end

    i, j = iseven(k) ? (1, 2) : (2, 1)
    return -Int(sign(y)) * A[i][2] // A[i][1]
end


function solve_exact_wan(A, b)
    # the function implements Algorithm 1 from
    #   An algorithm to solve integer linear systems exactly
    #   using numerical methods, Zhendong Wan
    # the notation is taken from there

    # note that arithmetic here is short except for the last reconstruction step

    m = size(A, 1)

    # compute the inverse of A numerically
    F = lu(A)
    Ainv = inv(UpperTriangular(F.U)) * inv(UnitLowerTriangular(F.L)) * F.P

    # the common denominator of the answer vector
    d = BigInt(1)
    # the residual
    r = b
    i = 0

    # the hadamard bound for A #
    B = (BigInt(m)^(m/2) * BigInt(norm(A, Inf))^m)

    # the accumulated numerical answers
    xs = []
    αs = BigInt[]

    # @info "" det(A) B

    # while the necessary accuracy is not obtained
    while d <= 2m * B^2 * (2.0^(-i)*norm(b, Inf) + norm(A, Inf))
        i += 1
        @info "" i=i d=d

        # compute numerical solution for the residual
        # in standard precision w.r.t to machine epsilon eps(Float64)
        x̄ = Ainv * r

        # α, "the scale" of the error |r| / |(r - Ax̄)|
        α = floor(Int, min(
                    2^30,
                    2^floor(log2(norm(r, Inf) / norm(r - A*x̄, Inf)))
                )
        )

        if α < 2
            error("insufficient numerical accuracy")
        end

        x = round.(Int, x̄ * α)

        # update the denominator
        d *= α
        # and the residual
        r = α * r - A * x

        push!(xs, x)
        push!(αs, α)
    end

    n = zeros(BigInt, m)
    if i == 1
        ds = [1]
    else
        ds = [ prod(αs[j+1:end]) for j in 1:i-1 ]
        push!(ds, one(last(ds)))
    end
    # @info "" ds
    for j in 1:i
        n += ds[j] * xs[j]
    end

    float_reconstruction.(n / d, B)
end




function check()
    C = 50

    n = 10
    ok = 0
    for i in 1:C
        A = floor.(Int, 2*randn(n, n))
        B = floor.(Int, 2*randn(n, 1))
        Y = solve_exact(A, B)
        if A * Y == B
            ok += 1
        end
    end
    oks = (100 * ok / C)
    @info "n=$n, integers max=$(maximum(abs.(A)))" oks

    #############

    n = 10
    ok = 0
    for i in 1:C
        A = floor.(Int, 1000*randn(n, n))
        B = floor.(Int, 1000*randn(n, 1))
        Y = solve_exact(A, B)
        if A * Y == B
            ok += 1
        end
    end

    oks = (100 * ok / C)
    @info "n=$n, integers max=$(maximum(abs.(A)))" oks

    #############

    n = 10
    ok = 0
    for i in 1:C
        A = floor.(Int, 1000000*randn(n, n))
        B = floor.(Int, 1000000*randn(n, 1))
        Y = solve_exact(A, B)
        if A * Y == B
            ok += 1
        end
    end

    oks = (100 * ok / C)
    @info "n=$n, integers max=$(maximum(abs.(A)))" oks

end

# check()



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



#=
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
=#
# check()

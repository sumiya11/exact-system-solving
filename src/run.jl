
include("rational_rec.jl")
include("parser.jl")


REALS = []
NAIVE = []
MODULAR = []
DIXON = []

ALL = (REALS, NAIVE, MODULAR, DIXON)

function uwu()

    for (name, sz, count, coos) in load_COO_if(from_dim=40, to_dim=50)
        @info "loaded " name
        for coo in coos
            A = coo2matrix(coo...)
            A *= convert(BigInt, lcm(BigInt.(denominator.(A))))
            A *= rand(1:2^31-1)
            A = convert.(BigInt, A)
            b = convert.(BigInt, rand(Int, sz, 1))

            for i in 1:size(A, 1)
                A[i, i] = 1
            end

            start = time_ns()
            x1 = solve_in_reals(A, b)
            push!(REALS, (sz, (time_ns() - start) * 1e-9))

            start = time_ns()
            x2 = solve_exact_naive(A, b)
            push!(NAIVE, (sz, (time_ns() - start) * 1e-9))

            start = time_ns()
            x3 = solve_exact_modular(A, b)
            push!(MODULAR, (sz, (time_ns() - start) * 1e-9))

            start = time_ns()
            x4 = solve_exact_dixon(A, b)
            push!(DIXON, (sz, (time_ns() - start) * 1e-9))
        end
    end

end


#uwu()

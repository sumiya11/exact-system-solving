

include("float_rec.jl")


function myroot(f)
    x = 0
    ∇f = derivative(f)

    i = 0
    while i < 10
        x = x - f(x) / ∇f(x)
        i += 1
    end

    x
end

function myroots(f)

    n = degree(f)
    roots = []

    x = gen(parent(f))

    for i in 1:n
        f_over_real = change_base_ring(AbstractAlgebra.RealField, f)

        root = myroot(f_over_real)
        root_rec = 1

        for M in [100, 1000, 10000, 100000]
            root_rec = float_reconstruction(root, M)

            try
                f = divexact(f, x - AbstractAlgebra.QQ(root_rec))
            catch ArgumentError
                continue
            end
            break

        end

        push!(roots, root_rec)
    end

    roots
end


function test2()

    R, x = AbstractAlgebra.QQ["x"]

    # f = (x - 3//2)*(x - 11//15)*(x - 8)*(x - 2343//31)*(x + 23)*(x + 232323)

    n = 100
    f = R(1)

    for _ in 1:n
        f *= x - rand(1:5) // rand(1:5)
    end

    f
end

f = test2()

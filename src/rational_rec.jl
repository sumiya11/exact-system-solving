

using LinearAlgebra
using AbstractAlgebra
using SparseArrays
using Nemo
using Distributions
using Primes

include("utils.jl")

#------------------------------------------------------------------------------

# rational number reconstruction implementation borrowed from CLUE
# and modified a bit to suit the 'Modern Computer Algebra' definitions
# returns a rational r // h of QQ field in a canonical form such that
#   r // h ≡ a (mod m)
#
# let n = max( λ(a), λ(m) ) , where λ(x) is a number of bits for x
# O(n^2)
function rational_reconstruction(a::I, m::J) where {I<:Union{Int, BigInt}, J<:Union{Int, BigInt}}
    a = mod(a, m)
    if a == 0 || m == 0
        return Nemo.QQ(0, 1)
    end
    if m < 0
        m = -m
    end
    if a < 0
        a = m - a
    end
    if a == 1
        return Nemo.QQ(1, 1)
    end
    bnd = sqrt(float(m) / 2)

    U = I[1, 0, m]
    V = I[0, 1, a]
    while abs(V[3]) >= bnd
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end

    t = abs(V[2])
    r = V[3] * sign(V[2])
    # changed from `<= bnd` to `<= m / bnd`
    if t <= m / bnd && gcd(r, t) == 1
        return Nemo.QQ(r, t)
    end

    throw(DomainError(
        :($a//$m), "rational reconstruction of $a (mod $m) does not exist"
    ))
end

#------------------------------------------------------------------------------

function modular_reduction(x, field)
    n, d = field(numerator(x)), field(denominator(x))
    if iszero(d)
        throw(DomainError(
            :($x), "modular reduction of $x (to $field) does not exist"
        ))
    end
    n // d
end

function modular_reduction(x::Int, field)
    field(x)
end

#------------------------------------------------------------------------------

function solve_exact_dixon(A, b)
    n = size(A, 1)

    # prime p, such that it does not divide det(A)
    p = 2^31 - 1

    @assert mod(det(A), p) != 0

    F = GF(p)

    MSpace = MatrixSpace(F, n, n)
    VSpace = MatrixSpace(F, n, 1)
    Ar = MSpace(A)
    Ainv = inv(Ar)

    # upper bound for denominators of x
    B = 2*BigFloat(norm(A, 2))^(2n - 1) * norm(b, 2)

    x̂ = zeros(BigInt, n)
    d = b

    i = 0

    while BigInt(p)^i < B
        y = Ainv * VSpace(d)                # solve mod p
        ybig = (x->BigInt(x.data)).(Array(y))
        x̂ += ybig * BigInt(p)^i                    # add new digit
        d = Int.((d - A*ybig) / p)         # update residue
        i += 1
    end

    rational_reconstruction.(x̂, BigInt(p)^i)
end

function solve_in_reals(A, b)
    convert.(Float64, A) \ convert.(Float64, b)
end

function solve_exact_naive(A, b)
    A \ b
end

function solve_exact_reduce(A, b)
    n = size(A, 1)

    # upper bound for denominators of x
    B = 2*BigFloat(norm(A, 2))^(2n - 1) * norm(b, 2)
    p = Nemo.ZZ(nextprime(ceil(BigInt, B)))

    F = GF(p)
    MSpace = MatrixSpace(F, n, n)
    VSpace = MatrixSpace(F, n, 1)
    Ar = MSpace(A)
    Ainv = inv(Ar)

    x̂ = Ainv * VSpace(b)

    rational_reconstruction.(BigInt.(map(x -> x.data, x̂)), BigInt(p))
end

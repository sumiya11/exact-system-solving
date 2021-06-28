
using LinearAlgebra
using SparseArrays
using Nemo
using Distributions
using AbstractAlgebra


# "Rational numbers with bound denominators are discrete!"

#=
    The pipeline:
        1. solve the problem over reals to obtain real solution
        2. reconstruct real solution to rationals

    Notes:
        ◐ the accuracy of numerical solution matters a lot
        ◕ estimating bitsize of rational solution in advance is useful
        ◔ we want a numerical method that converges fast, e.g power iteration
=#


#------------------------------------------------------------------------------

# Descriprion for the function below.
# Suppose we have a float number `y` we want to convert to rational.
# It is well known there exists no rational number closest to `y`.
#
# However, one may consider the class of rationals with denominators bound by integer `M`.
# This class *contains a rational, closest to `y`* w.r.t each other rational in the class.
# Moreover, if there exist a suitable rational with bound denominator approximating `y`,
# it is *unique*, and there exists a logarithmic algorithm for finding the one.
#
# function `float_reconstruction` implements this algorithm.
# Intrinsically, the algorithm consists of performing some iterations of
# Extended Euclidean Algorithm for numbers `y` and `1`

# for the given float `y` and natural `M`
# computes rational fraction `p/q` approximating `y` such that
#   |p/q - y| < 1/(2M²)
#       and
#   q ≤ M
# i.e the fraction `p/q` with the simplest denominator closest to `y`
function float_reconstruction(y::T, M) where {T<:AbstractFloat}
    # the function implements algorithm derived from Corollary 6.3a of
    #   Theory of Linear and Integer Programming, Schrijver, 1986,
    # the notation is taken from there
    #
    # note that no formal correctness proof of this
    # particular implementation was carried out

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

# Example:

a = 4 / 3  # this sets a to 1.33..
b = float_reconstruction(a, 10)
@assert b == 4 // 3  # b is a rational exactly equal to 4 // 3

#------------------------------------------------------------------------------

function inner(x, y)
    sum(x .* y)
end

function rayleigh(A, x)
    inner(A * x, x) / inner(x, x)
end

function power_iteration(A; err=200)
    # vanilla version, just for illustration

    n = size(A, 1)
    M = BigFloat(err)
    x = ones(BigFloat, n)
    r = x

    i = 0
    while i < 50 && maximum(abs.(r)) > 1 / (2err^2)
        @info "$i th: $( maximum(abs.(r)) )"

        Ax = A * x
        x = Ax / norm(Ax)
        r = rayleigh(A, x)*x - A * x

        i += 1
    end

    rayleigh(A, x), x
end

# Example
A = [
    15//7 0 0;
    0 1 0;
    0 0 1//4;
]

# obtain numerical solution
# λ = 2.1428571427228059...
λ, x = power_iteration(A)

λ_true = float_reconstruction(λ, 100)
@assert λ_true == 15//7

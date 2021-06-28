

function bitsize(x::Union{Int, BigInt})
    ceil(BigInt, log2(abs(x)))
end

function bitsize(x::Rational)
    max(bitsize(numerator(x)), bitsize(denominator(x)))
end

function bitsize(x::fmpq)
    bitsize(Rational(x))
end

function bitsize(x::AbstractArray)
    maximum(bitsize, x)
end


function coo2matrix(n, m, COO)
    A = zeros(Rational, n, m)
    for (i, j, x) in COO
        A[i, j] = x
    end
    A
end

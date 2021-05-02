

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

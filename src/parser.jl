
import JSON
import Nemo: gen, terms, monomials
import Base.showerror

include("myeval.jl")

#------------------------------------------------------------------------------

struct ParseException <: Exception
    msg::String
end

Base.showerror(io::IO, e::ParseException) = print(io, e.msg)


#------------------------------------------------------------------------------

# returns a model parsed from the given `filename`
# (with matrices in dense format)
# each returned matrix is represented with a dense array
function load_dense_set(filename::String)
    ans = []

    filepath = "$(normpath(joinpath(@__FILE__, "..", "..")))/src/data/matrices/$(filename)"

    open(filepath, "r") do input_stream
        for arr in JSON.Parser.parse(input_stream)
            push!(ans, arr)
        end
    end

    return ans
end

#------------------------------------------------------------------------------

# returns a model parsed from the given `filename`
# (with matrices in coordinate format)
# each returned matrix is represented with an array of
# items (i, j, x) where typeof(x) is Rational
function load_COO_set(filename::String)
    model = []

    filepath = "$(normpath(joinpath(@__FILE__, "..", "..")))/src/data/matrices/$(filename)"

    str_to_nom = (x -> parse(BigInt, first(split(x, '/'))))
    str_to_den = (x -> (if occursin('/', x) parse(BigInt, last(split(x, '/'))) else 1 end))

    open(filepath, "r") do input_stream
        for arr in JSON.Parser.parse(input_stream)
            normal_repr = [
                arr[1], arr[2],
                map(x -> (x[1], x[2], str_to_nom(x[3]) // str_to_den(x[3])), arr[3])
            ]

            push!(model, normal_repr)
        end
    end

    return model
end

#------------------------------------------------------------------------------

# returns an array of all models from data/matrices that satisfy
# the given criterias:
#   `from_dim`  - lower bound for the model's matrices one-side dimension
#   `to_dim`    - upper bound -//-
#   `from_size` - lower bound for the model's total number of matrices
#   `to_size`   - upper bound -//-
#
# each returned model is represented with an array:
#   [filename, dim, size, matrices-data]
# each returned matrix is represented with an array of
# items (i, j, x) of nonzeroes where typeof(x) is Rational
function load_COO_if(;from_dim, to_dim, from_size=0, to_size=Inf)
    models = []

    for filename in readdir("src/data/matrices/")
        model_data = load_COO_set(filename)
        model_dim = first(first(model_data))
        model_size = length(model_data)

        if from_dim <= model_dim <= to_dim &&
                from_size <= model_size <= to_size
            push!(models, [filename, model_dim, model_size, model_data])
        end
    end

    @info "loaded models: $(length(models))"

    return models
end

#------------------------------------------------------------------------------

# returns a model parsed from the given `group\name` file
# (with matrices in MTX format)
# each returned matrix is represented with an array of
#   [filename, dim, size, data]
# where data is stored as items (i, j, x) where typeof(x) is Rational
#
# In case of any parse error occurred, e.g inexact rationals encountered,
# returs an empty array []
function load_MTX(group, name)
    model = []

    # seems to be a default separator
    sep = r"%-+\r?\n"

    filepath = replace(
        "$(normpath(joinpath(@__FILE__, "..", "..")))src/data/sscollection/$group/$name",
        "\\" => "/"
    )

    open(filepath, "r") do input_stream
        # some models have >3 sections, taking the last
        data = last(split(read(input_stream, String), sep))

        # should be ok
        data = split(data, r"\r?\n")
        meta = map(x -> parse(Int, x), split(first(data), ' '))
        # strange matrices, discard them
        # NEED TO BE HANDLED IN ANOTHER WAY
        if length(meta) != 3
            return []
        end

        m, n, nnz = meta
        # we want only square ones
        if m != n
            return []
        end

        try
            model = [
                group * '/' * name, # filename
                m,                  # dim
                1,                  # size
                collect(            # data
                    map(
                        # `Any` disallows auto-promoting Int to Rational
                        X -> Any[
                            parse(Int, X[1]),
                            parse(Int, X[2]),
                            length(X) > 2 ? Rational(parse(Float64, X[3])) : Rational(1)
                        ] ,
                        filter(!isempty, map(split, view(data, 1:length(data))))
                    )
                )
            ]
        catch e
            if isa(e, InexactError)
                # occurrs if floats are impossibly small/big,
                # discard the model then
                return []
            else
                rethrow(e)
            end
        end
    end

    return model
end

#------------------------------------------------------------------------------

# returns an array of all models from data/sscollection that satisfy
# the given criterias:
#   `from_dim`  - lower bound for the model's matrices one-side dimension
#   `to_dim`    - upper bound -//-
#
# each returned model is represented with an array:
#   [filename, dim, size, matrices-data]
# each returned matrix is represented with an array of
# items (i, j, x) of nonzeroes where typeof(x) is Rational
function load_MTX_if(;from_dim=0, to_dim=Inf)
    models = []

    extension = ".mtx"
    dir = "src/data/sscollection/"

    for group in readdir(dir)
        for name in readdir("$dir/$group/")
            if endswith(name, extension)
                model = load_MTX(group, name)
                if ! isempty(model)
                    fn, dim, sz, data = model
                    if from_dim <= model[2] <= to_dim
                        push!(models, [ fn, dim, sz, [ dim, dim, data ] ] )
                    end
                end
            end
        end
    end

    @info "loaded models: $(length(models))"

    return models
end

#------------------------------------------------------------------------------

function load_ODEs_if(;from_size=-Inf, to_size=Inf)
    testspath = "src/data/ODEs/testing"
    ext = ".ode"

    models = []

    for modeltype in readdir(testspath)
        for filename in readdir("$testspath/$modeltype")
            if endswith(filename, ext)
                try
                    model = load_ODEs("$testspath/$modeltype/$filename")
                    if from_size <= length(model) <= to_size
                        push!(models, [filename, model])
                    end
                catch ex
                    !isa(ex, ParseException) && rethrow(ex)
                    @info ex.msg
                end
            end
        end
    end

    @info "loaded models: $(length(models))"

    models
end


# filepath - absolute model file path
#
function load_ODEs(filepath)
    # *.ode --> polynomial form of ode system

    lines = []
    open(filepath, "r") do inputs
        lines = map(strip, readlines(inputs))
    end

    lines = filter(!isempty, map(strip, lines[
        (findfirst(startswith("begin"), lines) + 1
        :
        findlast(startswith("end"), lines) - 1)
    ]))
    # by default we consider the last begin/end block to yield reactions.
    # If this conjecture is false, we do not work
    reactions = filter(!isempty, map(strip, lines[
        (findlast(startswith("begin"), lines) + 1
        :
        findlast(startswith("end"), lines) - 1)
    ]))

    #=
        Parsing for the two ODE formats is implemented:
            - ERODE "arrows" format for biological processes
            - Ordinary "equations" format

        These are mainly distinct by the presence of the '->' token,
        and we will try to dispatch on it
    =#

    sep = occursin("->", prod(lines)) ? "->" : "="

    # scan reactions to discover all variables

    sep1 = Regex("($sep)|(,)|([+-/ \t\\*])")
    # assuming there is at least one utf letter in any variable
    strings = map(String, unique(filter!(
        x -> !isempty(findall(isletter, x)),
        split(join(reactions, ' '), sep1)
    )))

    if length(strings) > 1000
        throw(ParseException("too long system encountered at $filepath, skipping it"))
    end

    # we want to agree with CLUE in variable order
    strings = sort(strings)
    S, xs = QQ[strings...]

    # symbol :x to x from QQ[x]
    mapping = Dict{Symbol, fmpq_mpoly}(
        # Gleb: Meta.parse is a bit overshoot here, you can do Symbol(x)
        Symbol(x) => gen(S, i)
        for (i, x) in enumerate(strings)
    )

    ODEs = Dict(x => S(0) for x in xs)

    if sep === "->"
        sep2 = r"(->)|(,)"
        for reaction in reactions
            # divide each line into parts,
            # expected to be three of them
            lhs, rhs, speed = map(
                Meta.parse ∘ String ∘ strip,
                split(reaction, sep2)
            )

            reagents, products, speed = map(
                x -> myeval(x, mapping),
                (lhs, rhs, speed)
            )

            concentration = prod(terms(reagents))

            for reagent in monomials(reagents)
                ODEs[reagent] -= concentration * speed
            end
            for product in monomials(products)
                ODEs[product] += concentration * speed
            end
        end
    else
        sep2 = r"(=)|(,)"
        for reaction in reactions
            # divide each line into lhs, rhs
            lhs, rhs = map(
                Meta.parse ∘ String ∘ strip,
                split(reaction, sep2)
            )

            ẋ, rhs = map(
                x -> myeval(x, mapping),
                (lhs, rhs)
            )

            ODEs[ẋ] = rhs
        end
    end

    @info "loaded a system of $(length(ODEs)) ODEs from $filepath"

    return ODEs
end

#formulas for the design (loaded with StatsModels): X and the names of its columns from a formula and a table of
#samples, contrasts given by column names
module PLNRUVStatsModelsExt

using PLNRUV, StatsModels

function PLNRUV.design_matrix(f::FormulaTerm, data)
    rhs = apply_schema(f, schema(f, data), StatisticalModel).rhs #with the intercept, unless the formula removes it
    X = Float64.(modelcols(rhs, StatsModels.Tables.columntable(data)))
    names = coefnames(rhs)
    return (X = X isa AbstractVector ? reshape(X, :, 1) : X, names = names isa AbstractString ? [names] : collect(names))
end

#a contrast given by numbers, by a column name, by name => weight pairs, or a list of these (several contrasts)
contrast_spec(names, c::AbstractVector{<:Real}) = float.(c)
contrast_spec(names, c::AbstractMatrix{<:Real}) = float.(c)
contrast_spec(names, c::AbstractString) = contrast_vector(names, c)
contrast_spec(names, c::AbstractVector{<:Pair}) = contrast_vector(names, c)
contrast_spec(names, c::AbstractVector) = reduce(hcat, [contrast_spec(names, x) for x in c])

function PLNRUV.pln_mouthwash(Y::AbstractMatrix{<:Integer}, f::FormulaTerm, data, q::Integer; contrast = nothing, kwargs...)
    X, names = PLNRUV.design_matrix(f, data)
    c = contrast === nothing ? nothing : contrast_spec(names, contrast)
    return merge(pln_mouthwash(Y, X, q; contrast = c, kwargs...), (X = X, coefnames = names))
end

PLNRUV.choose_q(Y::AbstractMatrix{<:Integer}, f::FormulaTerm, data; kwargs...) = choose_q(Y, PLNRUV.design_matrix(f, data).X; kwargs...)

end

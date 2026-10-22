"""
    log_ratio(counts, reference = all genes; pseudocount = 0.5)

log-ratio of counts (samples in rows, genes in columns) to the geometric mean of a set of reference genes
(the multi-additive log-ratio, MALR); with all genes as reference it is the centred log-ratio
zeros are handled by a pseudocount
"""
function log_ratio(counts::AbstractMatrix{<:Real}, reference::AbstractVector{<:Integer} = axes(counts, 2); pseudocount::Real = 0.5)
    Z = log.(counts .+ pseudocount)
    Z .-= mean(view(Z, :, reference); dims = 2)
    return Z
end

#QR rotation of Wang J. et al., Ann. Stat. (2017): https://doi.org/10.1214/16-AOS1511
#with X = [X₁ X₂] (X₂ the covariates of interest) and X = Q R, the rows of Q'Y split into
#k₁ nuisance rows (dropped), k₂ rows carrying R₂₂ β₂ + confounding, and n - k residual rows free of X
function rotate(Y::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}, cov_of_interest::AbstractVector{<:Integer})
    n, k = size(X)
    @assert size(Y, 1) == n
    k₂ = length(cov_of_interest)
    k₁ = k - k₂
    nuisance = setdiff(1:k, cov_of_interest)
    F = qr([X[:, nuisance] X[:, cov_of_interest]])
    Yrot = F.Q' * Y
    return (
        Yint = Yrot[(k₁ + 1):k, :], Yres = Yrot[(k + 1):n, :],
        R22 = UpperTriangular(F.R[(k₁ + 1):k, (k₁ + 1):k]),
    )
end

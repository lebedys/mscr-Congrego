
export best_cheapest_generator

# metrics

#=
xlogx(x) = x > 0 ? -x * log(x) : zero(x)

entropy(p::AbstractVecOrMat{F}) where {F <: AbstractFloat} = sum(xlogx, p)

function entropy(c::AbstractVecOrMat{I}) where {I <: Integer}
    s = sum(c)
    return sum(cᵢ->xlogx(cᵢ/s), c)
end
=#

"""
    best_cheapest_generator(archive::Dict, pc::PartCatalogue)

Return the `(fitness, generator, components)` entry from a MAP-Elites archive with the highest fitness.
Ties are broken by lowest synthesis cost (via `get_price`).
"""
function best_cheapest_generator(archive::Dict, pc::PartCatalogue)
    argmax(((f, g, _),) -> (f, -get_price(g, pc)), values(archive))
end
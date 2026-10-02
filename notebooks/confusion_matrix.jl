# Confusion matrices, shared by ms_experimental_kl.ipynb and ms_experimental_cello.ipynb.

using Plots

const CONFUSION_AGREE    = colorant"#5E7CA0"
const CONFUSION_DISAGREE = colorant"#B8836A"
const CONFUSION_INVALID  = colorant"#B8B8B8"
const CONFUSION_TEXT     = colorant"#3A3A3A"
const CONFUSION_RULE     = colorant"#DCDCDC"

"""
    confusion_labels(observed, predictions)

Canonical phenotype order, restricted to the labels either side uses. One shared set
across every prediction keeps the diagonal true and the panels comparable.
"""
function confusion_labels(observed, predictions)
    used = Set{Symbol}()
    for call in skipmissing(observed)
        push!(used, Symbol(call))
    end
    for prediction in predictions, call in prediction
        push!(used, call)
    end
    return [label for label in [pheno_labels(2); :invalid] if label in used]
end

"""
    confusion_counts(observed, predicted, labels)

`[observed, predicted]` counts. Circuits without an observed call are skipped.
"""
function confusion_counts(observed, predicted, labels)
    index = Dict(label => i for (i, label) in enumerate(labels))
    matrix = zeros(Int, length(labels), length(labels))
    for (call, prediction) in zip(observed, predicted)
        ismissing(call) && continue
        matrix[index[Symbol(call)], index[prediction]] += 1
    end
    return matrix
end

"""
    confusion_scores(matrix, labels)

Agreement, Cohen's kappa, and the fraction of circuits the model settled on.
"""
function confusion_scores(matrix, labels)
    total = sum(matrix)
    agree = sum(matrix[i, i] for i in axes(matrix, 1))
    rows, cols = vec(sum(matrix, dims=2)), vec(sum(matrix, dims=1))
    chance = sum(rows .* cols) / total^2
    kappa = (agree / total - chance) / (1 - chance)
    invalid = findfirst(==(:invalid), labels)
    coverage = isnothing(invalid) ? 1.0 : 1 - sum(matrix[:, invalid]) / total
    return (; agree, total, kappa, coverage)
end

"""
    confusion_panel(matrix, labels; xlabel, ylabel, kwargs...)

Blue on the diagonal, brown off it, grey where the model did not settle.
`kwargs` are forwarded to `Plots.plot`.
"""
function confusion_panel(matrix, labels; xlabel="Predicted", ylabel="Observed", kwargs...)
    n = length(labels)
    hi = maximum(matrix)
    ypos(i) = n - i + 1

    p = plot(; xlims=(0.5, n + 0.5), ylims=(0.5, n + 0.5),
             xticks=(1:n, string.(labels)), yticks=(1:n, string.(reverse(labels))),
             xrotation=45, aspect_ratio=1, legend=false, framestyle=:box,
             tick_direction=:none, xlabel=xlabel, ylabel=ylabel,
             foreground_color_axis=CONFUSION_TEXT, foreground_color_text=CONFUSION_TEXT,
             foreground_color_border=CONFUSION_RULE, kwargs...)

    for k in 0:n
        plot!(p, [0.5, n + 0.5], [k + 0.5, k + 0.5];
              color=CONFUSION_RULE, linewidth=0.4, label="")
        plot!(p, [k + 0.5, k + 0.5], [0.5, n + 0.5];
              color=CONFUSION_RULE, linewidth=0.4, label="")
    end

    for i in 1:n, j in 1:n
        count = matrix[i, j]
        count == 0 && continue
        fill = labels[j] === :invalid ? CONFUSION_INVALID :
               (i == j ? CONFUSION_AGREE : CONFUSION_DISAGREE)
        y = ypos(i)
        plot!(p, Shape([j - 0.5, j + 0.5, j + 0.5, j - 0.5],
                       [y - 0.5, y - 0.5, y + 0.5, y + 0.5]);
              fillcolor=fill, fillalpha=0.30 + 0.70count / hi,
              linecolor=CONFUSION_RULE, linewidth=0.5, label="")
        # drawn twice a hair apart: GR offers Plots no bold weight
        for dx in (-0.022, 0.022)
            annotate!(p, j + dx, y, text(string(count), 9, :center, :vcenter, CONFUSION_TEXT))
        end
    end
    return p
end

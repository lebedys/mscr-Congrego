# Phenotype-frequency CSV output, shared by ms.ipynb and ms_misassembly.ipynb.

using CSV, DataFrames

"""
    write_phenotype_csv(path, entries)

Phenotype frequencies as one row per `(part_set, generator, labels, counts)` entry.
All entries must carry the same labels, which become the columns after
`part_set` and `generator`.
"""
function write_phenotype_csv(path, entries)
    labels = string.(entries[1][3])
    df = DataFrame(part_set=String[], generator=String[])
    for label in labels
        df[!, label] = Int[]
    end
    for (part_set, generator, entry_labels, counts) in entries
        string.(entry_labels) == labels ||
            error("$part_set/$generator carries different phenotype labels")
        push!(df, (part_set, generator, counts...))
    end
    CSV.write(path, df)
    return df
end

# Congrego

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Congrego

It is authored by Michiel Stock, Sim Castle, Thomas Gorochowski, and Adomas Lebedys.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "Congrego"
```
which auto-activate the project and enable local path handling from DrWatson.

## Combinatorial Synthesis of Diverse Genetic Circuits

Four notebooks under `notebooks/` were used for the thesis investigation.

### `ms.ipynb`
Phenotype distribution of the library as designed and as assembled.
- input: `scripts/parts/ms_original.jl`, `scripts/parts/ms_roadblock_pTac.jl` (part catalogues)
- input: the unique assembled genotypes, pasted inline in the notebook (see [Data](#data))
- output: `data/ms_phenotypes.csv` (phenotype frequencies per part set and generator)
- output: `plots/ms_phenotypes_original_vs_roadblock.{pdf,png}` (original against roadblocked)

### `ms_misassembly.ipynb`
Phenotype distribution of the misassembled library under cryptic-site removal, pTra disabling and pTac roadblocking.
- input: `scripts/parts/ms_{original,roadblock_pTac,disable_Tra,disable_Tra_roadblock_pTac}.jl` (part catalogues)
- output: `data/ms_misassembly_phenotypes.csv` (phenotype frequencies per part set and generator)
- output: `plots/ms_misassembly_phenotypes.{pdf,png}` (designed, assembled and both interventions, four panels)
- output: `plots/ms_misassembly_original.{pdf,png}` (the assembled library on its own)
- output: `plots/ms_misassembly_cryptic_disableTra.{pdf,png}` (cryptic sites removed, pTra disabled)
- output: `plots/ms_misassembly_cryptic_disableTra_roadblockpTac.{pdf,png}` (above, with pTac roadblocked)

### `ms_experimental_kl.ipynb`
Measured phenotypes, called by KL divergence, against the simulated prediction for each assembled genotype.
- input: `seq/out/michiel_lib_genotype_filtered.csv` (from `michiel_lib_genotype.ipynb`, see [Data](#data))
- input: `nb/out/analyse_pr_classify.csv` (from `analyse_pr_classify.ipynb`, see [Data](#data))
- input: `scripts/parts/ms_disable_Tra.jl`, `scripts/parts/ms_disable_Tra_roadblock_pTac.jl` (part catalogues)
- output: `data/ms_experimental_predictions_kl.csv` (observed and predicted phenotype per genotype)
- output: `plots/ms_experimental_confusion_kl.png` (observed against predicted, without and with roadblocking)
- output: `plots/ms_experimental_distributions_kl.png` (predicted vs observed distributions)

### `ms_experimental_cello.ipynb`
The same comparison, with the measured phenotypes called by Cello score instead.
- input: `seq/out/michiel_lib_genotype_filtered.csv` (from `michiel_lib_genotype.ipynb`, see [Data](#data))
- input: `nb/out/analyse_pr_classify.csv` (from `analyse_pr_classify.ipynb`, see [Data](#data))
- input: `scripts/parts/ms_disable_Tra.jl`, `scripts/parts/ms_disable_Tra_roadblock_pTac.jl` (part catalogues)
- output: `data/ms_experimental_predictions_cello.csv` (observed and predicted phenotype per genotype)
- output: `plots/ms_experimental_confusion_cello.png` (observed against predicted, without and with roadblocking)
- output: `plots/ms_experimental_distributions_cello.png` (predicted vs observed distributions)

### Data

- Part catalogues are included in this repository, under `scripts/parts/`.
- The genotype and classification CSVs come from
  [mscr-data-analysis](https://github.com/lebedys/mscr-data-analysis). Both `ms_experimental`
  notebooks resolve them through `analysis_dir`, which defaults to `../dev` and should be pointed
  at the cloned repository.
- The assembled genotypes pasted into `ms.ipynb` are from
  [Zenodo record 18378852](https://zenodo.org/records/18378852), under
  `masters_thesis_supplementary/data/genotype/mseg_output_michiel_cleaned.tsv`.
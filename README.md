# cd_predictor
Pipeline to generate SnoBIRD, a C/D box snoRNA predictor that is usable across eukaryotes. This also reproduces all figures present in Fafard-Couture et al. 2025

## Required dependencies

- **Snakemake** version >= 7.18.2<br>
- **Conda** version >= 4.12.0<br>
- **Mamba** version >= 0.15.3<br>

## Other dependencies
Most packages and tools used in this pipeline are automatically installed via conda/mamba. If one wants to recreate the model comparison figures, they must install manually the following tool:
- **snoReport2**: It must be installed manually following exactly the 
installation guidelines found [here](https://joaovicers.github.io/snoreport2/index.html#principal).
- **snoscan** via conda
- **infernal (cmscan)** with Rfam's CM models via conda
<br>

## 1. Downloading the required data
You **must** necessarily download the following data and datasets for this pipeline to work entirely, using the following command:
```bash
snakemake all_downloads --profile ../profile_local/
```
## 2. Generating the required datasets
To generate the tuning, training and test sets used for SnoBIRD training and comparison with other tools (i.e. a mandatory step to be able to run the following steps), run the following command:
```bash
snakemake --profile ../profile_local/
```
## 3. Training and testing SnoBIRD
To optimize, train and test SnoBIRD's models, you must run the following command:
```bash
snakemake all_snoBIRD_training --profile ../profile_local/
```
## 4. Recreating the main figures
To recreate automatically most of the main figures in the manuscript, you should run the following command (after having run the previous three commands): 
```bash
snakemake all_figures --profile ../profile_local/
```
## 5. Recreating other downstream figures that need specific complex input
To recreate the remaining figures, you cannot run the following command directly, because some of the dependent dataframes needed to recreate these remaining figures necessitate specific installation (e.g. installing and running the other tools, running SnoBIRD on a HPC on specific genome sequences, etc.). Please refer to rule `all_figures_independent` in the Snakefile to see which other jobs should be run before running this one. Once all these dependent jobs are done, you should then run the following command:
```bash
snakemake all_figures_independent --profile ../profile_local/
```




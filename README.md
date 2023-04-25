# cd_predictor
Pipeline to produce a C/D box snoRNA predictor that is usable across eukaryotes. 

## Tested with the following dependencies

- **Snakemake** version >= 7.18.2<br>
- **Conda** version >= 4.12.0<br>
- **Mamba** version >= 0.15.3<br>

## Other dependencies
Most packages and tools used in this pipeline are automatically installed via conda/mamba. One can also skip the models comparison 
part to avoid installing the following tools.
- **snoReport2**: Unfortunately, the tool snoReport2 is not available via conda and must be installed manually following exactly the 
installation guidelines found [here](https://joaovicers.github.io/snoreport2/index.html#principal).<br>








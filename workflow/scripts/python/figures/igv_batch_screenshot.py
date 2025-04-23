#!/usr/bin/python3
from pybedtools import BedTool
import subprocess as sp

snoBIRD_bed = BedTool(snakemake.input.snoBIRD_bed)

snoBIRD_bed.igv(name=True, path="./results/figures/screenshots", img="svg")

sp.call(f'mv {snoBIRD_bed.igv_script} {snakemake.output.batch_script}', shell=True)


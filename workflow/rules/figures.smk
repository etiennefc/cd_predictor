rule length_selection_positives_density:
    """ Create a density plot showing the length 
        distribution of positive examples. Also 
        show different thresholds based on 
        percentile/quartiles to choose a cutoff. 
        The chosen percentile is 95% (< 183 nt)"""
    input:
        positives = rules.get_sno_sequences.output.df
    output:
        density = 'results/figures/density/length_selection_positives.svg'
    params:
        percent_colors = ['#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#b10026']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/length_selection_positives_density.py"



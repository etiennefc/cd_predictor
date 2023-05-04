#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

def density(df, xlabel, ylabel, title, path, **kwargs):
    """
    Creates a simple density plot.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.kdeplot(df, fill=True, ax=ax, **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 35})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 35})
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    fig.suptitle(title, fontsize=65, weight='bold', x=0.36, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')


def density_percentile(df, xlabel, ylabel, title, path, percentile_list, percent_colors, percent_label, **kwargs):
    """
    Creates a density plot with vertical lines 
    correpsonding to percentiles in percentile_list.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.kdeplot(df, fill=True, ax=ax, **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 35})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 35})
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    legend_list = []
    for i, perc in enumerate(percentile_list):
        ax.vlines(x=perc, ymin=0, ymax=0.02, linestyles='dashed', colors=percent_colors[i])
        legend_element = mpatches.Patch(color=percent_colors[i], label=percent_label[i])
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(1,1),
                fontsize=20)
    fig.suptitle(title, fontsize=65, weight='bold', x=0.36, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

def lineplot(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path, **kwargs):
    """ 
    Create a vertical connected dot plot or lineplot.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    rc = {'ytick.labelsize': 20, 'xtick.labelsize': 20}
    plt.rcParams.update(**rc)    
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.lineplot(df, x=x_col, y=y_col, hue=hue_col, palette=color_dict, 
                    marker='o', markeredgecolor='grey', **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 25})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 25})
    ax.set_ylim(0, 1.05)
    fig.suptitle(title, fontsize=65, weight='bold', x=0.36, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

def pie_multiple(y, x, count_list, labels, colors, ax_title, title, legend_title, path, **kwargs):
    """
    Creates 8 pie charts from a list of list (count_list) where each global
    element corresponds to a list of local elements (ex: percent per rank across
    8 model intersection).
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(y, x, figsize=(25, 14))
    plt.subplots_adjust(hspace=0.5)
    ax = axes.flatten()
    for i, element in enumerate(count_list):
        count_per_element = count_list[i][:]
        ax[i].set_title(ax_title[i], fontdict={'fontsize': 25}, x=0.5, y=0.8)
        ax[i].pie(count_per_element, colors=colors, textprops={'fontsize': 19},
                    **kwargs)
        ax[i].axis('equal')
        white_circle = plt.Circle((0, 0), 0.4, color='white') #to create a donut chart
        ax[i].add_artist(white_circle) #to create a donut chart

    fig.suptitle(title, x=0.5, y=0.9, fontsize=25)
    fig.legend(labels=labels, loc='upper right', bbox_to_anchor=(1.1, 0.5),
                prop={'size': 20}, title=legend_title)
    plt.savefig(path, dpi=600)

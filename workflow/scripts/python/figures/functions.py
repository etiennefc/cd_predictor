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
    ax.set_xlabel(xlabel, fontdict={'fontsize': 30})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 30})
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    fig.suptitle(title, fontsize=40, weight='bold', x=0.5, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

def density_x(df_list, xlabel, ylabel, xscale, title, colors, crit_list, path, xvline=0, yminvline=0, ymaxvline=0, **kwargs):
    """
    Creates a density plot with a number x of dfs to represent on the same ax.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    for i, df in enumerate(df_list):
        sns.kdeplot(data=df, fill=True, ax=ax, color=colors[i], **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 35})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 35})
    ax.set_xscale(xscale)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.vlines(x=xvline, ymin=yminvline, ymax=ymaxvline, 
            linestyles='dashed', colors='black')

    legend_list = []
    for i, crit in enumerate(crit_list):
        legend_element = mpatches.Patch(color=colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0,1.1),
                fontsize=20)

    fig.suptitle(title, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=500)

def density_x_mod(df_list, xlabel, ylabel, xscale, title, colors, crit_list, path, xvline=0, yminvline=0, ymaxvline=0, **kwargs):
    """
    Creates a density plot with a number x of dfs to represent on the same ax.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    #fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    for i, df in enumerate(df_list):
        sns.kdeplot(df, fill=True, color=colors[i], **kwargs)
    #plt.xlabel(xlabel, fontdict={'fontsize': 35})
    #plt.ylabel(ylabel, fontdict={'fontsize': 35})
    #ax.set_xscale(xscale)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    #plt.spines['right'].set_linewidth(0)
    #plt.spines['top'].set_linewidth(0)
    #ax.vlines(x=xvline, ymin=yminvline, ymax=ymaxvline, 
    #        linestyles='dashed', colors='black')

    legend_list = []
    for i, crit in enumerate(crit_list):
        legend_element = mpatches.Patch(color=colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0,1.1),
                fontsize=20)

    plt.title(title, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=500)

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

def scatter_multiple(y, x, df, hue_col, hue_list, col_list_y, col_x, colors, title, path, **kwargs):
    """
    Creates n scatter plots from a list columns.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(y, x, figsize=(50, 14))
    for i, element in enumerate(col_list_y):
        for j, hue in enumerate(hue_list): 
            sns.scatterplot(data=df[df[hue_col] == hue], x=col_x, y=col_list_y[i], ax=axes[j, i], color=colors[j])
    

    fig.suptitle(title, x=0.5, y=0.95, fontsize=25)
    legend_list = []
    for i, crit in enumerate(hue_list):
        legend_element = mpatches.Patch(color=colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(2.25, 1.25),
                fontsize=18)
    plt.savefig(path, dpi=600)

def density_multiple(y, x, df, hue_col, hue_list, second_hue_col, col_list_x, colors, row_titles, title, path, **kwargs):
    """
    Creates n density plots from a list columns.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(y, x, figsize=(50, 14))
    for i, element in enumerate(col_list_x):
        for j, hue in enumerate(hue_list): 
            if isinstance(hue, list):
                temp_df = df[df[hue_col].isin(hue)]
            else:
                temp_df = df[df[hue_col].isin([hue])]
            sns.kdeplot(data=temp_df, hue=second_hue_col, x=col_list_x[i], ax=axes[j, i], color=colors, **kwargs)
    
    for row_idx, row_title in enumerate(row_titles):
        fig.text(0.04, 1 - (row_idx + 0.5) / y, row_title, va='center', ha='right', fontsize=15, weight='bold')

    fig.suptitle(title, x=0.5, y=0.99, fontsize=25)
    plt.tight_layout(rect=[0.05, 0, 1, 1])
    plt.savefig(path, dpi=600)

def shap_lineplot(x_val, y_val, xtick_labels, xlabel, ylabel, title, path, **kwargs):
    """ 
    Create a vertical connected dot plot or lineplot.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    rc = {'ytick.labelsize': 25, 'xtick.labelsize': 25}
    plt.rcParams.update(**rc)    
    plt.subplots(1, 1, figsize=(30, 8))
    plt.plot(x_val, y_val)
    plt.xlabel(xlabel, fontdict={'fontsize': 30})
    plt.ylabel(ylabel, fontdict={'fontsize': 30})
    plt.xticks(range(len(xtick_labels)), xtick_labels, fontsize=10, rotation=90)
    plt.margins(x=0)
    #ax.set_ylim(0, 1.05)
    plt.title(title, fontsize=30, x=0.5, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

def lineplot(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path, **kwargs):
    """ 
    Create a vertical connected dot plot or lineplot.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    rc = {'ytick.labelsize': 15, 'xtick.labelsize': 15, 'xtick.alignment': 'right'}
    plt.rcParams.update(**rc)    
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.lineplot(df, x=x_col, y=y_col, hue=hue_col, palette=color_dict, 
                    marker='o', markeredgecolor='grey', **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 20})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 20})
    xtick_labels = [tick.get_text() for tick in ax.get_xticklabels()]
    ax.set_xticks(list(range(0, len(xtick_labels))))
    ax.set_xticklabels(xtick_labels, rotation=25)
    ax.set_ylim(0, 1.05)
    plt.legend(fontsize=35)
    for handle in plt.gca().get_legend().legend_handles:
        handle.set_linewidth(16)
    fig.suptitle(title, fontsize=30, x=0.5, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

def lineplot_errorbars(df, x_col, y_col, std_col, hue_col, xlabel, ylabel, title, color_dict, path, **kwargs):
    """ 
    Create a vertical connected dot plot or lineplot with errorbars.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    rc = {'ytick.labelsize': 40, 'xtick.labelsize': 40}
    plt.rcParams.update(**rc)    
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.lineplot(df, x=x_col, y=y_col, hue=hue_col, color=color_dict, 
                    marker='o', markeredgecolor='grey', **kwargs)
    lower = df[y_col] - df[std_col]
    upper = df[y_col] + df[std_col]                
    ax.plot(df[x_col], lower, color='tab:grey', alpha=0.1)
    ax.plot(df[x_col], upper, color='tab:grey', alpha=0.1)
    ax.fill_between(df[x_col], lower, upper, color='lightgrey', alpha=0.2)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 40})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 40})
    ax.set_ylim(0, 1.05)
    fig.suptitle(title, fontsize=30, x=0.5, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')


def lollipop(df, color_dict, title, path):
    plt.rcParams['svg.fonttype'] = 'none'
    rc = {'ytick.labelsize': 20, 'xtick.labelsize': 20}
    plt.rcParams.update(**rc)  
    score_order = ['accuracy', 'precision', 'recall', 'f1_score']
    predictors = df['predictor'].unique()
    x = np.arange(len(score_order))

    # Set up the plot
    plt.figure(figsize=(12, 10))
    width = 0.15  # distance between side-by-side sticks

    # Define colors
    palette = color_dict

    # Draw lollipops for each predictor
    for i, predictor in enumerate(predictors):
        subset = df[df['predictor'] == predictor]
        values = subset.set_index('score_name').loc[score_order]['score_value']
        offsets = x + (i - len(predictors)/2) * width + width/2
        plt.vlines(offsets, 0, values, color=palette[predictor], linewidth=10)
        plt.plot(offsets, values, 'o', label=predictor, color=palette[predictor], markersize=20)

    # Formatting
    plt.xticks(x, score_order)
    plt.ylabel("Metrics value", fontsize=25)
    plt.xlabel("Metrics", fontsize=25)
    plt.ylim(0, 1.05)
    plt.legend()
    plt.tight_layout()
    plt.title(title, fontsize=25)
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
    fig.legend(labels=labels, loc='upper right', bbox_to_anchor=(0.5, 0.25),
                prop={'size': 10}, title=legend_title)
    plt.savefig(path, dpi=600)

def donut(count_list, labels, colors, ax_title, title, legend_title, path, **kwargs):
    """
    Creates a donut chart from a list.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(25, 14))
    plt.subplots_adjust(hspace=0.5)
    ax.set_title(ax_title, fontdict={'fontsize': 25}, x=0.5, y=0.8)
    ax.pie(count_list, colors=colors, textprops={'fontsize': 19},
                **kwargs)
    ax.axis('equal')
    white_circle = plt.Circle((0, 0), 0.4, color='white') #to create a donut chart
    ax.add_artist(white_circle) #to create a donut chart

    fig.suptitle(title, x=0.5, y=0.9, fontsize=25)
    fig.legend(labels=labels, loc='upper right', bbox_to_anchor=(1.1, 0.5),
                prop={'size': 20}, title=legend_title)
    plt.savefig(path, dpi=600)

def stacked_bar(lists, x_tick_labels, labels, title, xlabel, ylabel, colors, min_y, max_y, optional_annot, path, **kwargs):
    """
    Create a stacked bar chart from a list of lists ('lists').
    """
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    df = pd.DataFrame(lists, index=x_tick_labels, columns=labels)
    print(df)
    ax = df.plot.bar(stacked=True, figsize=(30, 12), color=colors, **kwargs)
    ax.set_xticklabels(x_tick_labels, rotation=0)
    
    # Add optional annotation above bars (ex: number of sno in each bar)
    n = len(labels)
    
    def chunker(list_, size):
        # Get chunks of n elements from list_ where n=size
        return (list_[pos:pos + size] for pos in range(0, len(list_), size))
    
    # Get the cumulative height of each bar until the before last stacked bar
    previous_heigths = [0] * len(x_tick_labels)
    for i, chunks_index in enumerate(chunker(list(range(0, n*len(x_tick_labels))[:-len(x_tick_labels)]), len(x_tick_labels))):  # from 0 to the before last bars
        for index_ in chunks_index:
            bar = list(ax.patches)[index_]
            previous_heigths[index_ - len(x_tick_labels) * i] += bar.get_height()
    
    # Add the cumulative height of previous bars to the last stacked bar of each stack
    last_bars = [bar_ for j, bar_ in enumerate(ax.patches) if j in list(range(0, n*len(x_tick_labels))[-len(x_tick_labels):])]  # last (i.e. topmost) bars
    for i, bar in enumerate(last_bars):
        ax.text((bar.get_x() + bar.get_width()/(len(x_tick_labels)/2) - 0.1), 
                (bar.get_height() + previous_heigths[i] + 10), optional_annot[i], fontsize=35)
    plt.legend(fontsize=30, loc='center right', bbox_to_anchor=(1.3, 0.5))


    plt.title(title, fontsize=38)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.autoscale()
    plt.margins(0.02)
    plt.ylim(min_y, max_y)
    plt.savefig(path, bbox_inches='tight', dpi=600)

def stacked_bar2(lists, x_tick_labels, labels, title, xlabel, ylabel, colors, min_y, max_y, optional_annot, path, **kwargs):
    """
    Create a stacked bar chart from a list of lists ('lists').
    """
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 30}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    df = pd.DataFrame(lists, index=x_tick_labels, columns=labels)
    print(df)
    ax = df.plot.bar(stacked=True, figsize=(35,8), color=colors, **kwargs)
    ax.set_xticklabels(x_tick_labels, rotation=0)
    
    # Add optional annotation above bars (ex: number of sno in each bar)
    for rect, annotation in zip(ax.patches, [f'({i})' for i in optional_annot] * len(labels)):
        x = rect.get_x() + rect.get_width() / 2  # Center above the bar
        ax.text(x, 106, annotation, ha='center', va='bottom', fontsize=15)

    plt.legend(fontsize=30, loc=5, bbox_to_anchor=(1, -1.25))
    plt.title(title, fontsize=40, y=1.1, x=0.5)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.autoscale()
    plt.margins(0.02)
    plt.ylim(min_y, max_y)
    plt.savefig(path, bbox_inches='tight', dpi=600)

def stacked_bar3(lists, x_tick_labels, labels, title, xlabel, ylabel, colors, min_y, max_y, optional_annot, path, **kwargs):
    """
    Create a stacked bar chart from a list of lists ('lists').
    """
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    df = pd.DataFrame(lists, index=x_tick_labels, columns=labels)
    print(df)
    print(df['FN'].sum())
    ax = df.plot.bar(stacked=True, figsize=(40, 10), color=colors, **kwargs)
    #ax.set_xticklabels(x_tick_labels, rotation=45)
    
    # Add optional annotation above bars (ex: number of sno in each bar)
    n = len(labels)
    
    def chunker(list_, size):
        # Get chunks of n elements from list_ where n=size
        return (list_[pos:pos + size] for pos in range(0, len(list_), size))
    
    # Get the cumulative height of each bar until the before last stacked bar
    previous_heigths = [0] * len(x_tick_labels)
    for i, chunks_index in enumerate(chunker(list(range(0, n*len(x_tick_labels))[:-len(x_tick_labels)]), len(x_tick_labels))):  # from 0 to the before last bars
        for index_ in chunks_index:
            bar = list(ax.patches)[index_]
            previous_heigths[index_ - len(x_tick_labels) * i] += bar.get_height()
    
    # Add the cumulative height of previous bars to the last stacked bar of each stack
    last_bars = [bar_ for j, bar_ in enumerate(ax.patches) if j in list(range(0, n*len(x_tick_labels))[-len(x_tick_labels):])]  # last (i.e. topmost) bars
    for i, bar in enumerate(last_bars):
        ax.text((bar.get_x() + bar.get_width()/(len(x_tick_labels)/2) - 0.1), 
                551, optional_annot[i], fontsize=35)
    plt.legend(fontsize=50, loc='center right', bbox_to_anchor=(1, 0.9))


    plt.title(title, fontsize=38)
    plt.xlabel(xlabel, fontsize=50)
    plt.ylabel(ylabel, fontsize=50)
    plt.autoscale()
    plt.margins(0.02)
    plt.ylim(min_y, max_y)
    plt.savefig(path, bbox_inches='tight', dpi=600)

def grouped_stacked_bar2(lists, x_tick_labels, labels, title, xlabel,
                        ylabel, colors, legend_title, optional_annot, path, **kwargs):
    """
    Create a grouped stacked bar chart. Lists is a list of two lists of lists;
    labels is the labels of the stacked variables.
    """
    plt.rcParams['svg.fonttype'] = 'none'

    # Create dfs from the lists in 'lists'
    df1 = pd.DataFrame(lists[0], index=x_tick_labels, columns=labels)
    df2 = pd.DataFrame(lists[1], index=x_tick_labels, columns=labels)

    # Create the bar plot
    fig, ax = plt.subplots(1, 1, figsize=(35, 8))
    df2.plot.bar(ax=ax, position=0, width=.3, color=colors, stacked=True,
                edgecolor='black', **kwargs)
    df1.plot.bar(ax=ax, position=1, width=.3, color=colors, stacked=True,
                edgecolor='black', **kwargs)
    plt.autoscale()
    plt.title(title, fontsize=35, x=0.5, y=1.15)
    plt.xlabel(xlabel, fontsize=35)
    plt.ylabel(ylabel, fontsize=35)
    plt.margins(0.02)
    plt.ylim(0, 110)

    for i, tick_label in enumerate(x_tick_labels):
        # Calculate the position of the first group (grouped stacked bar)
        x_left = i - 0.3  # Adjust based on the position and width
        x_right = i   # Adjust based on the width of the bar
        y_top = df1.loc[tick_label].sum()  # Total height of the stacked bar
        rectangle = mpatches.Rectangle(
            (x_left, 0), x_right - x_left, y_top,
            fill=False, edgecolor='black', linewidth=4.5
        )
        ax.add_patch(rectangle)

    # Add optional annotation above bars (ex: number of sno in each bar)
    for rect, annotation in zip(ax.patches, [f'({i})' for i in optional_annot] * len(labels)):
        x = rect.get_x() + rect.get_width() / 2 - 0.3  # Center above the bar
        ax.text(x, 104, annotation, ha='center', va='bottom', fontsize=15)

    # Add legend
    legend_list = []
    for i, crit in enumerate(list(reversed(labels))):
        legend_element = mpatches.Patch(color=colors[crit], label=crit)
        legend_list.append(legend_element)
    legend_list.append(
        mpatches.Patch(facecolor='white', edgecolor='black', linewidth=4.5, label="SnoBIRD predictions"))
    legend_list.append( 
        mpatches.Patch(facecolor='white', edgecolor='black', label="Actual genome proportion"))
    legend = ax.legend(handles=legend_list, bbox_to_anchor=(1,-0.7), fontsize=25)
    legend.set_title(legend_title,prop={'size':18})
    ax.add_artist(legend)

    plt.savefig(path, bbox_inches='tight', dpi=600)


def percent_count(count_list):
    """
    Create from a list of lists a percentage list of lists.
    """
    percent_list = []

    for i, cat in enumerate(count_list):
        temp = []
        for j, crit in enumerate(cat):
            total = sum(cat)
            if total != 0:
                percent = crit/total * 100
                percent = round(percent, 2)
            else:
                percent = 0
            temp.append(percent)
        percent_list.append(temp)

    return percent_list

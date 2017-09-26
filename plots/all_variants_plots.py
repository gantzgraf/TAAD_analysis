''' Contains functions utilised for producing plots from an all variants df
'''
import matplotlib as mpl
mpl.use('Agg')   # allows one to run matplotlib and seaborn headless
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import data_cleaning.simple_filters as sf
from PIL import Image
from io import BytesIO


def all_variants_barplot(df, outfile):
    ''' Produces a split barplot showing the percentage of
        pathogenic or likley pathogenic variants amongst 
        all variants disovered per gene within the given
        df

    NOTE: 
        THIS IS FOR MOST DAMAGING VARS (SKI EXON1 & AB < 0.3 REMOVED) ONLY
    '''
    clean_all_variants = df
    variant_counts = variant_counts_df(clean_all_variants, 'Symbol', 'All Genes')
    split_barplot_variants(variant_counts, outfile)


def variant_counts_df(df, gene_column, column_to_sort_by='All Genes'):
    ''' Get pathogenic/likely pathogenic and all variant counts for all genes
        within a df

    Args:
        gene_column: name of the column contaning the gene name
        column_to_sort_by: sort by All Gene or by Pathogenic Genes 

    Returns:
        df containing the number of total variants and pathogenic/likely
        pathogenic variants
    '''
    counts_new = df[df['AB'] > 0.3][gene_column].value_counts().to_dict()
    cond = ((df['New Category'] == 'Pathogenic/Likely Pathogenic') & (df['validation'] == 1))
    path_new = df[cond][gene_column].value_counts().to_dict()
    
    compare_table = pd.DataFrame([counts_new, path_new]).transpose().fillna("-")


    # produce a table comparing all gene counts againts pathogenic gene counts
    compare_table.columns = ['All Genes', 'Pathogenic Genes']
    compare_table = compare_table.apply(pd.to_numeric, errors='coerce')
    compare_table = compare_table.sort_values(column_to_sort_by, ascending=False)
    
    return compare_table


def split_barplot_variants(df, outfile, colour="b", left_extend=0.15):
    ''' Produces a split barplot showing the percentage of
        pathogenic or likley pathogenic variants amongst 
        all variants disovered per gene within the given
        df

    Args:
        outfile: desired path and name of the outputted PNG file
        colour: refers to the bar colours, defaulted to blue
        left_extend: amount of etrax space to give the y-labels
    '''

    all_column = df.columns[0] 
    pathogenic_column  = df.columns[1]
    
    #convert entire df to numeric type
    df = df.apply(pd.to_numeric, errors='coerce')

    # set up colors etc.
    f, ax = plt.subplots(figsize=(9,9))
    sns.set(font_scale=1.05)
    sns.set_style("whitegrid")
    
    # ALL
    sns.set_color_codes("muted")
    sns.barplot(x=all_column, y='index', data=df.reset_index(),
                label="All Variants", color=colour)   
    
    # Pathogenic
    sns.set_color_codes("pastel")
    sns.barplot(x=pathogenic_column, y='index', data=df.reset_index(),
                label="Pathogenic or Likely\nPathogenic Variants", color=colour)
    
    # Add a legend and informative axis label, swap labels
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1],loc="lower right", frameon=False)
    ax.set_title("")

    # set x limits and labels
    ax.set_xlabel("Variants Called", fontsize=15)
    ax.set_ylabel("")
    ax.tick_params(labelsize=12)
    
    # calculate percent of each genes that are pathogenic
    percent_pathogenic = [p/h*100 for p, h in zip(df[pathogenic_column].tolist(), df[all_column].tolist())]

    for p, total in zip(ax.patches, percent_pathogenic):
        h = p.get_width()
        if np.isnan(total) == True:
            annotation = "0.0%"
        else:
            annotation = str(round(total, 1))+"%"
        ax.text(h+3, p.get_y()+0.5, annotation, ha='left')
    
    
    # adjust y-axis if needed
    plt.gcf().subplots_adjust(left=left_extend)
    fig = ax.get_figure()
    #fig.savefig(outfile+".png")
     # adjust y-axis if needed
    
    # Workaround to save in TIFF format
    png1 = BytesIO()
    fig.savefig(png1, format='png', dpi=300)
    # fig.savefig(outfile+'.png', dpi=300)

    png2 = Image.open(png1)
    png2.save(outfile+".tiff")
    png1.close()
    


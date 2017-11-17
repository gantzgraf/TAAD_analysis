''' A collection of plotting functions that display phenotype features against variant classifications '''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats
from data_cleaning.conversion import convert2category
import plots.plot_manipulations as pm
import data_cleaning.simple_filters as sf

def age_v_family_history(df, column, outfile=None):
    ''' Boxplot showing median age at diagnosis between
        pateints with and without a family history.
    '''
    df = df.dropna(subset = [column, 'family_history']) 
    df = df[df['family_history'] != 'unknown']
    df['family_history'] = convert2category(df['family_history'], ['yes', 'no'])

    fig, ax = plt.subplots(figsize=(15, 12))
    sns.set(font_scale=1.5)
    sns.set_style("whitegrid")
    sns.set_palette('Greys')

    ax = sns.boxplot(data=df, x='family_history', y=column)
    ax.set(title="", ylabel='Age at Diagnosis', xlabel='', ylim=(0,100))
    name_list = pm.rename_xtick(df, 'family_history')
    ax.set_xticklabels(name_list)
    
    family_history = df[df['family_history'] == 'yes'][column]
    no_family_history = df[df['family_history'] == 'no'][column]
    fam_p_val = str(pm.RoundToSigFigs(stats.ranksums(family_history, no_family_history)[1], 2 ))
    pm.line_between_plots(ax, x1=0, x2=1, height=90, string='p = '+fam_p_val, fontsize=18)

    if outfile:
        ax.figure.savefig(outfile)

def variant_class_violin(df, column, title='', outfile=None):
    ''' Produces a violin plot of the age at surgery vs the variant class.'''
    df = df.dropna(subset = [column, 'New Category'])   
    y_label = column.replace("age at diagnosis", "Age at Diagnosis")
    name_list = pm.rename_xtick(df, 'New Category', 
                                NA=True, validated_path=True)

    fig, ax = plt.subplots(figsize=(15, 12))
    sns.set(font_scale=1.5)
    sns.set_style("whitegrid")
    sns.set_palette('Greys')
    ax = sns.violinplot(data=df, x='New Category', y=column, cut=0)
    ax.set(title=title, ylabel=y_label, xlabel='', xticklabels=name_list,
           ylim=(0, 105))

    # Seperate and clean each variant class Series
    pathogenic = df[df['New Category'] == "Pathogenic/Likely Pathogenic"]
    pathogenic = sf.validated_only(pathogenic)[column].dropna()
    damaging = df[df['New Category'] == "VUS"][column].dropna()
    benign = df[df['New Category'] == "Likely Benign / No Variant"][column].dropna()

    # Get unpaired ranksum wilcoxon p-values between each new category and place them upon the plot
    #path_ben_p_val = str(pm.RoundToSigFigs(stats.ranksums(pathogenic, benign)[1] ,3))
    path_ben_p_val = str.format("{:.3g}", stats.ranksums(pathogenic, benign)[1])
    pm.line_between_plots(ax, x1=0, x2=2, height=benign.max()+15, 
                        string="p = "+path_ben_p_val, fontsize=15)

    #path_dam_p_val = str(pm.RoundToSigFigs(stats.ranksums(pathogenic, damaging)[1], 3))
    path_dam_p_val = str.format("{:.3g}", stats.ranksums(pathogenic, damaging)[1])
    pm.line_between_plots(ax, x1=0, x2=1,height=benign.max()+10, 
                        string="p = "+path_dam_p_val, fontsize=15)

    #dam_ben_p_val = str(pm.RoundToSigFigs(stats.ranksums(damaging, benign)[1], 3))
    dam_ben_p_val = str.format("{:.3g}", stats.ranksums(damaging, benign)[1])
    pm.line_between_plots(ax, x1=1, x2=2, height=benign.max()+5, 
                        string="p = "+dam_ben_p_val, fontsize=15)

    if outfile:
        ax.figure.savefig(outfile)

    return ax

def age_group_v_pathogenic_piechart(df, outfile=None):
    '''  Two side-by-side piecharts displaying patients under and over 50
         and the percentage of each variant class is present in each group.
    '''
    # remove non validated PLP variants and NaN values in Age Group
    df = df[~((df['validation'] == 0) & (df['New Category'] == 'Pathogenic/Likely Pathogenic'))]
    df = df.dropna(subset=['Age Group', 'New Category']) 
    
    # create percentage lists for piechart
    table = df.groupby(['Age Group', 'New Category']).size()
    percentage_table = table.groupby(level=0).apply(lambda x: x / float(x.sum())).unstack().T
    young_percent = percentage_table['Under 50'].tolist()[::-1] # [PLP, VUS, Benign]
    old_percent = percentage_table['Over 50'].tolist()[::-1]    
    
    # create values for xtick labels
    total_counts = df['Age Group'].value_counts().to_dict()
    young_xlabel = "Under 50\nn = {}".format(total_counts.get('Under 50'))
    old_xlabel = "Over 50\nn = {}".format(total_counts.get('Over 50'))
    
    # sharing an axis allows the use of pm.line_between_plots
    fig, ax = plt.subplots(figsize=(15, 12))
    mpl.rcParams['font.size'] = 12 # cannot find an ax method for this
    ax.pie(young_percent, labels=['Likely Benign / No Variant', 'VUS', 'Pathogenic/Likely Pathogenic'],
           colors=['grey', 'silver', 'white'], shadow=False,  autopct='%1.1f%%', 
           center=(0,0), startangle=90, labeldistance=1.2, pctdistance=0.7) 
    ax.pie(old_percent, labels=['Likely Benign / No Variant', 'VUS', 'Pathogenic/Likely Pathogenic'],
           colors=['grey', 'silver', 'white'], shadow=False,  autopct='%1.1f%%', 
           center=(2.5,0), startangle=90, labeldistance=1.2, pctdistance=0.7)

    # axis modifications
    ax.axis('equal')
    ax.set_xticks([0, 2.5])
    ax.set_xticklabels([young_xlabel, old_xlabel])
    ax.tick_params(labelsize=15)

    # STATS
    # As I am an analysing a large number of samples, the chi-square test is approriate
    # and will yield a simliar result as the Fisher-Freeman-Halton test
    contingency_table = table.unstack().values.tolist()
    chi2, pvalue, dof, expected = stats.chi2_contingency(contingency_table)
    pm.line_between_plots(axs=ax, x1=0, x2=2.5, height=1.5, fontsize=15, extend=0.2,
                          string="p = {:.2g}".format(pvalue))
    if outfile:
        ax.figure.savefig(outfile)
    
    return ax

def fh_vs_genetic_diagnosis(df, outfile=None):
    ''' Countplot displaying counts for each variant classification
        split by family history.
    '''
    df = df.dropna(subset=['New Category'])
    df = df[df['family_history'] != 'unknown']
    df['family_history'] = convert2category(df['family_history'], ['yes', 'no'])
    
    sns.set_style("whitegrid")
    sns.set_palette('Greys')
    fig, ax = plt.subplots(figsize=(14, 12))
    name_list = pm.rename_xtick(df, 'family_history', counts=False)
    ax = sns.countplot(data=df, x='family_history', hue='New Category')
    ax.set(title="  .  ", xlabel="", ylabel="Samples", xticklabels=name_list,ylim=(0,500))    
    ax.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=8, ncol=3, mode="expand", borderaxespad=0.)
    
    # label percentages above bars
    fh_total = len(df[df['family_history'] == 'yes'])
    no_fh_total = len(df[df['family_history'] == 'no'])
    for p, total in zip(ax.patches, [fh_total, no_fh_total, fh_total, no_fh_total, fh_total, no_fh_total]):
        h = p.get_height()
        ax.text(p.get_x()+0.15, h+3, str('{:.1f}'.format((h/total)*100)+"%"), ha='center') 
    
    # STATS: As I am an anlasying a large number of samples, the chi-square test is approriate
    # and will yield a simliar result as the Fisher-Freeman-Halton test
    table = df.groupby(['family_history', 'New Category']).size().unstack(1).fillna(0.0)
    array = table.values.tolist()
    chi2, pvalue, dof, expected = stats.chi2_contingency(array)
    pm.line_between_plots(axs=ax, x1=0, x2=1, height=array[1][2]+20, extend=10, 
                          string="p = "+str(pm.RoundToSigFigs(pvalue,3)),
                          fontsize=14)

    if outfile:
        ax.figure.savefig(outfile)

def gender_vs_genetic_diagnosis(df, outfile=None):
    ''' Countplot displaying counts for each variant classification
        split by gender
    '''
    df = df.dropna(subset=['Gender', 'New Category'])

    sns.set(font_scale=1.2)
    sns.set_style("whitegrid")
    sns.set_palette("Greys")
    
    fig, ax = plt.subplots(figsize=(14, 12))
    ax = sns.countplot(data=df, x='Gender', hue='New Category')
    name_list = pm.rename_xtick(df, 'Gender', counts=False)
    ax.set_xticklabels(name_list)
    ax.set(title='  .  ', xlabel="", ylabel="Samples", ylim=(0,700))
    ax.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=8, ncol=3, mode="expand", borderaxespad=0.)

    # STATS
    gen_total = len(df[df['Gender'] == 'Male'])
    no_gen_total = len(df[df['Gender'] == 'Female'])
    table = df.groupby(['Gender', 'New Category']).size().unstack(1)
    array = table.values.tolist()
    chi2, pvalue, dof, expected = stats.chi2_contingency(array)
    pm.line_between_plots(axs=ax, x1=0, x2=1, height=array[0][2]+45, extend=10, 
                          string="p = "+str(round(pvalue,3)), fontsize=14)

    # label percentages above bars
    for p, total in zip(ax.patches, [gen_total, no_gen_total, gen_total, no_gen_total,
                        gen_total, no_gen_total]):
        h = p.get_height()
        ax.text(p.get_x()+0.15, h+3, str(round((h/total), 3)*100)[:4]+"%", ha='center')

    if outfile:
        fig = ax.get_figure()
        fig.savefig(outfile, bbox_inches='tight')

    return ax

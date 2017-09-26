''' A collection of plotting functions which compare phenotype features against validated PLP genes'''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import data_cleaning.simple_filters as sf
from data_cleaning.conversion import convert2category

def fh_v_genes_facetgrid(df, outfile=None):
    ''' MultiAxis countplot of validated PLP
        variant gene counts between patients with 
        and without a family history.
    '''
    df = sf.truly_pathogenic(df)
    df = df[df['family_history'] != "Unknown Family History"]

    sns.set(font_scale=1.5, style="whitegrid")
    g = sns.FacetGrid(df, row='family_history', hue='Symbol',
                      aspect=1.9, sharey=False, size=6.5, palette='Greys',
                      row_order=['yes', 'no'])
    symbol_order = df['Symbol'].value_counts().index.values.tolist()
    g = g.map(sns.countplot, "Symbol", order=symbol_order)
    g.facet_axis(0,0).set(ylim=(0,16), title='Family History')
    g.facet_axis(1,0).set(ylim=(0,16), title='No Family History')
    g.set_axis_labels("", "Samples")

    if outfile:
        g.savefig(outfile+'.png')

def age_diagnosis_v_genes(df, outfile=None):
    ''' Stripplot displaying age at diagnosis against 
        genes (PLP most damaging)
    '''
    df = df.dropna(subset = ['Symbol', 'age at diagnosis'])
    pathogenic_df = sf.truly_pathogenic(df)
    genes = list(pathogenic_df['Symbol'].dropna().unique())
    pathogenic_df['Symbol'] = convert2category(pathogenic_df['Symbol'], genes)

    sns.set(style='whitegrid')
    g = sns.factorplot(data=pathogenic_df, x='Symbol',y='age at diagnosis',
                          kind='strip', size=10, palette='Greys', linewidth=0.3)
    g.set(ylabel='Age at Diagnosis', xlabel='')
    # add median dashes to the plot - adding a Line2D object would perhaps be better
    median_age = [pathogenic_df[pathogenic_df['Symbol'] == x]['age at diagnosis'].median() for x in genes]
    for pos, median_val in zip(range(0, len(genes)), median_age):
        plt.text(pos - 0.1, median_val - 1, "—", color="black", fontsize=14)

    if outfile:
        g.savefig(outfile)

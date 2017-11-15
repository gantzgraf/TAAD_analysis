import all_variant_dataframe as av
import most_damaging_dataframe as md
import tables.variant_table as vt
import tables.variant_summary as vs
import tables.risk_ratio as rr
import tables.demographics as demo
# from plots.all_variants_plots import all_variants_barplot
import plots.phenotype_gene_plots as pgp
import plots.phenotype_variant_plots as pvp
import plots.all_variants_plots as avp
import pandas as pd
import os

pd.set_option('display.max_columns', 500)

def main(yale_phenotype, yale_all_variants, yale_most_damaging,
         uk_phenotype, uk_all_variants, uk_most_damaging, yale_survival,
         FILE_PATH):
    ''' Create, merge and clean variant CSV files and utilise the resulting
        DataFrames to produce cleaned data, tables and plots within the output
        directory.
    '''
    # create output dirs
    if not os.path.exists(FILE_PATH+'/output/'):
        for sub_dir in ['cleaned_data', 'plots', 'tables']:
            os.makedirs(FILE_PATH+'/output/'+sub_dir)
    # create and clean all variants DataFrame
    all_tuple = av.create_all_variants(yale_phenotype, yale_all_variants, 
                                       uk_phenotype, uk_all_variants)
    UK_all_variants, Yale_all_variants, all_variants = all_tuple
    all_variants.reset_index(inplace=True)
    # create and clean most damaging DataFrame
    most_damaging = md.create_most_damaging(UK_all_variants, uk_most_damaging, uk_phenotype,
                                            Yale_all_variants, yale_most_damaging, yale_phenotype,
                                            yale_survival, FILE_PATH)
    # need to filter on depth so that we only calculate risks etc. on samples 
    # we have sequenced successfully
    most_damaging = most_damaging[most_damaging['Depth'] != 'LOW']
    # output both DataFrames as CSV file
    all_variants.to_csv(FILE_PATH+"/output/cleaned_data/All_Variants.csv")
    most_damaging.to_csv(FILE_PATH+"/output/cleaned_data/Most_Damaging.csv")
    # used most damaging DataFrame to produce plots and tables
    tables(most_damaging, all_variants)
    plots(most_damaging)

def tables(most_damaging, all_variants):
    ''' Create all the tables for the manuscript '''
    demo.demographics_table(df=most_damaging, 
                            outfile=FILE_PATH+"output/tables/Patient_Demographics.csv")
    vt.variant_table(df=all_variants, 
                     outfile=FILE_PATH+"output/tables/"
                     "Pathogenic & Likely Pathogenic Variants Detected by NGS Panel.csv")
    vt.variant_table(df=all_variants, 
                     pathogenic=False, 
                     outfile=FILE_PATH+"output/tables/VUS Variants Detected by NGS Panel.csv")
    vs.variant_summary_table(df=all_variants, 
                             outfile=FILE_PATH+"output/tables/Summary_of_Variants.csv")
    rr.risk_ratio_table(df=most_damaging, 
                        outfile=FILE_PATH+"output/tables/RR_table.csv") 

def plots(most_damaging):
    ''' Generate all plots associated with the manuscript'''
    plot_path = FILE_PATH+'output/plots/'
    no_mfs = most_damaging[most_damaging['Known Syndrome'] != 'Marfan']

    avp.all_variants_barplot(df=most_damaging,
                             outfile=plot_path+'All Most Damaging Variant Counts.png')
    pvp.age_v_family_history(df=most_damaging, 
                             column='age at diagnosis', 
                             outfile=plot_path+'Age at Diagnosis Vs Family History.png')
    pvp.variant_class_violin(df=most_damaging, 
                             column='age at diagnosis',
                             title='Age at Diagnosis Vs Variant Class',
                             outfile=plot_path+\
                             'Age at Diagnosis Vs Variant Class.png')
    pvp.variant_class_violin(df=no_mfs,
                             column='age at diagnosis',
                             title='Age at Diagnosis Vs Variant Class - No MFS',
                             outfile=plot_path+\
                             'Age at Diagnosis Vs Variant Class - No MFS.png')
    pvp.age_group_v_pathogenic_piechart(df=most_damaging,
                                        outfile=plot_path+\
                                        'Age Group Vs Variant Class.png')
    pvp.age_group_v_pathogenic_piechart(df=no_mfs,
                                        outfile=plot_path+\
                                        'Age Group Vs Variant Class'
                                        '- No MFS.png')
    pvp.fh_vs_genetic_diagnosis(df=most_damaging,
                                outfile=plot_path+\
                                'Family History Vs Variant Class.png')
    pvp.gender_vs_genetic_diagnosis(df=most_damaging,
                                    outfile=plot_path+\
                                    'Gender Vs Variant Class.png')
    pgp.fh_v_genes_facetgrid(df=most_damaging,
                             outfile=plot_path+\
                             'Family Vs PLP Genes.png')
    pgp.age_diagnosis_v_genes(df=most_damaging,
                              outfile=plot_path+\
                              'Age at Diagnosis Vs PLP Genes.png')


if __name__ == '__main__':
    FILE_PATH = os.path.dirname(os.path.abspath("__file__"))+"/TAAD_analysis/"
    ipath = FILE_PATH+'input_files/'

    yale_phenotype = ipath+'Yale_Phenotype_Data.csv'
    yale_all_variants = ipath+'Yale_All_Variants_Data.csv'
    yale_most_damaging = ipath+'Yale_Most_Damaging_Data.csv'
    yale_survival = FILE_PATH+"/input_files/Yale_Survival_Data_Clean.csv"
   
    uk_phenotype = ipath+'UK_Phenotype_Data.csv'
    uk_all_variants = ipath+'UK_All_Variants_Data.csv'
    uk_most_damaging = ipath+'UK_Most_Damaging_Data.csv'

    main(yale_phenotype, yale_all_variants, yale_most_damaging, 
         uk_phenotype, uk_all_variants, uk_most_damaging, yale_survival,
         FILE_PATH)

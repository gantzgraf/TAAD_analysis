import data_cleaning.genotype_phenotype as gp
import data_cleaning.simple_filters as sf
import data_cleaning.new_columns as nc
import pandas as pd
from data_cleaning import conversion
from data_cleaning import rename

def create_all_variants(yale_phenotype, yale_genotype, uk_phenotype, uk_genotype):
    ''' Clean and concatenate the Yale and UK all 
        variants data.
    
    Args:
        yale_phenotype: path to yale phenotype data
        yale_genotype: path to yale genotype data
        uk_phenotype: path to uk phenotype data
        uk_genotype: path to uk genotype data
    
    Returns:
        a tuple of cleaned and merged uk all variants, 
        yale all variants and combined all variants data
    '''    
    UK_all_variants_clean = cohort_all_variants(uk_phenotype, uk_genotype, 'UK')
    Yale_all_variants_clean = cohort_all_variants(yale_phenotype, yale_genotype, 'Yale')
    all_variants = pd.concat([UK_all_variants_clean,
                              Yale_all_variants_clean])
    return (UK_all_variants_clean, Yale_all_variants_clean, all_variants)

def cohort_all_variants(phenotype, genotype, cohort):
    ''' Merge and clean a cohorts phenotype and genotype data'''
    variants = gp.merge_genotype_phenotype(phenotype, genotype)
    variants['cohort'] = cohort
    clean_variants = clean_all_var_df(variants)
    return clean_variants

def clean_all_var_df(df, three_categories=True):
    ''' Clean up of the all variants data. '''
    df = df[df['Symbol'] != 'SMAD4'] #SMAD4 should be ignored
    df['Dup'] = df.apply(lambda x: mark_duplicate_samples(x, df), axis=1)
    df = df[~df['Dup'].str.contains("Duplicate")]
    df = conversion.convert2numeric(df, ['age at diagnosis'])
    # rename_columns performed in merge_genotype_phenotype
    df = rename.rename_columns(df)
    df = rename.rename_entries(df)
    df = nc.create_new_columns(df, three_categories)
    df = df[~df['Sample'].str.contains("Blank|blank|ddH20|dH2O|H2O|BLANK|ddh2o")]
    df = sf.no_SKI_exon1(df)
    return df

def mark_duplicate_samples(x, df):
    ''' Mark all duplicate samples so they can be later remove'''
    if x['Sample'][:-2] in df['Sample'].tolist() or x['Sample'].endswith('_pool7A'):
        return "Duplicate"
    else:
        return "-"

''' merge_genotype_phenotype() and it's helper functions merge the genotype and phenotype datasheets into one DataFrame'''
import pandas as pd
from data_cleaning import rename

def merge_genotype_phenotype(phenotype, genotype):
    ''' Clean genotype and phenotype data and merge them
        on sample
    Args:
        phenotype: path to phenotype data
        genotype: path to genotype data
    '''
    phenotype_clean = clean_phenotype_data(phenotype)
    genotype_clean = clean_genotype_data(genotype)
    merged = pd.merge(genotype_clean, phenotype_clean, on=['Sample'])
    return merged

def clean_phenotype_data(phenotype):
    ''' Clean the phenotype data.
    Args: 
        phenotype: path to data
    '''
    unwanted_char = {' ':'', '-':'', '\'':''}
    name_corrections = {'24SA1565': '21SA1565', '24SS1575': '21SS1575', 
                        '24GC1574': '21GC1574', '24DR1571': '21DR1571',
                        '24FP1566': '21FP1566', '24GC1574': '21GC1574', 
                        '24AS1570': '21AS1570', '24KS0915': '24ZS0915',
                        '1328': '21RL1328', '24DW932': '24DW0932', 
                        '926': '24GN0926', '931': '24SG0931', 
                        '937': '24CB0937', '1327': '21SN1327', 
                        '1374': '24AS1374', 'MK3598': 'MK_35_98'}
    pdf = pd.read_csv(phenotype, encoding='iso-8859-1')
    pdf_clean = rename.rename_columns(pdf)
    pdf_clean = rename.replace_series_strings(pdf_clean, 'Sample', 
                                              unwanted_char, 
                                              substring=True)
    pdf_clean = rename.replace_series_strings(pdf_clean, 'Sample', 
                                              name_corrections, 
                                              substring=False)
    return pdf_clean

def clean_genotype_data(genotype):
    ''' Clean the genotype data and filter for the
        genotype columns of interest
    Args:
        genotype: path to genotype
    '''
    genotype_columns = ['Sample', 'AD', 'AB', 'UID', 'validation', 
                        'Category', 'Score','Symbol', 'HGVS', 'Chrom', 
                        'Pos', 'Ref', 'Alt', 'Consequence', 'HGVSc', 
                        'HGVSp', 'Exon', 'Intron']
    gdf = pd.read_csv(genotype, encoding='iso-8859-1')
    gdf_clean = rename.rename_columns(gdf)
    mask = PLP2VUS(gdf_clean)
    gdf_clean.ix[mask, 'Category'] = 'Uncertain Significance'
    gdf_clean.ix[mask, 'New Category'] = 'VUS'
    gdf_filtered = gdf_clean[genotype_columns]
    return gdf_filtered

def PLP2VUS(df):
    ''' Mark variants to be converted to VUS classification.
        These variants had insufficient evidence in HGMD to 
        be truly classified as VUS.
    '''
    mask = (
           ((df['Chrom'] == 2) & (df['Pos'] == 189851842)) |
           ((df['Chrom'] == 3) & (df['Pos'] == 123401086)) |
           ((df['Chrom'] == 9) & (df['Pos'] == 101908876)) |
           ((df['Chrom'] == 15) & (df['Pos'] == 48760242)) |
           ((df['Chrom'] == 15) & (df['Pos'] == 48776056)) |
           ((df['Chrom'] == 15) & (df['Pos'] == 48644711)) |
           ((df['Chrom'] == 16) & (df['Pos'] == 15844048)) |
           ((df['Chrom'] == 15) & (df['Pos'] == 48782270)) |
           ((df['Chrom'] == 15) & (df['Pos'] == 48800841)) |
           ((df['Chrom'] == 15) & (df['Pos'] == 48829865))
    )
    return mask

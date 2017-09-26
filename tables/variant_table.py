''' variant_table() and its helper functions returns a table describing all validated pathogenic variants or VUS identified within the study along with selected phenotype and genotype features'''
import data_cleaning.simple_filters as sf

def variant_table(df, pathogenic=True, outfile=None):
    ''' Returns table of all P/LP or VUS variants
    Args:
        df: DataFrame
        pathogenic: filter for PLP (True) or VUS (False)
    '''
    table = filter_table(df, pathogenic)
    table = hgvs_variant(table)
    table = further_anonymise(table)
    table = clean_table(table)
    table = rename_columns(table)
    table = sort_reorder(table)
    table.reset_index(inplace=True, drop=True)
    if outfile:
        table.to_csv(outfile, index=False)
    return table

def filter_table(df, pathogenic=True):
    ''' Filter DataFrame for columns of interest'''
    coi = ['Sample', 'primary diagnosis', 'Known Syndrome', 
           'family_history', 'Symbol', 'HGVSc', 'HGVSp', 
           'Consequence', 'Category', 'cohort']
    if pathogenic:
        df = sf.truly_pathogenic(df)
    else:
        df = df[df['New Category'] == 'VUS']
    table = df[coi]
    return table

def hgvs_variant(table):
    '''Combine HGVSc and HGVSp nomenclature'''
    HGVSc = table.HGVSc.str.split(':').str[1]
    HGVSp = table.HGVSp.str.split(':').str[1]
    table['Variant'] = HGVSc + ":" + HGVSp.fillna('')    
    return table.drop(['HGVSc', 'HGVSp'], axis=1)

def further_anonymise(table):
    ''' Get rid of the sample/patient intials and prepend cohort'''
    new_sample = table.Sample.str.replace("[^0-9]", "_")
    yale = table['cohort'] == 'Yale'
    uk = table['cohort'] == 'UK'
    table.ix[uk, 'Sample'] = 'UK_' + new_sample[uk]
    table.ix[yale, 'Sample'] = 'Y' + new_sample[yale]
    table['Sample'] = table.Sample.str.replace('[_]+', '_')
    return table

def clean_table(table):
    table['Consequence'] = table['Consequence'].str.replace("_", " ")
    table['Known Syndrome'] = table['Known Syndrome'].fillna('N/A')
    table = table.drop('cohort', axis=1)
    return table

def rename_columns(table):
    table.rename(columns={'primary diagnosis':'Primary Diagnosis', 'Sample': 'ID',
                          'Symbol': 'Gene Affected', 'Consequence': 'Functional Category',
                          'Category': 'Classification', 'family_history': 'Family History',
                          'Known Syndrome': 'Clinical Diagnosis'}, inplace=True)
    return table

def sort_reorder(table):
    order = ['ID', 'Primary Diagnosis', 'Clinical Diagnosis', 
             'Family History', 'Gene Affected', 'Variant',
             'Functional Category', 'Classification']
    table = table[order]
    return table.sort_values('ID', ascending=True)

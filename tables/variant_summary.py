''' Generates a table of variant counts by gene for PLP variants and VUS.'''

def variant_summary_table(df, outfile=None):
    ''' Gene variant counts of all VUS and validated pathogenic
        and likely pathogenic variants return as a table.
    '''
    df = plp_vus(df)
    gp = df.groupby(['Symbol', 'Category'])['Category'].size()
    table = gp.unstack().fillna(0)
    table = rename_sort_columns(table)
    if outfile:
        table.to_csv(outfile)
    return table

def plp_vus(df):
    ''' Filter for VUS and validated PLP variants'''
    plp = ((df.Category.str.contains('Pathogenic')) & (df['validation'] == 1))
    vus = (df['Category'] == 'Uncertain Significance')
    return df[plp | vus]

def rename_sort_columns(df):
    df = df.rename(columns={'Uncertain Significance': 'VUS'})
    return df[['Pathogenic', 'Likely Pathogenic', 'VUS']]

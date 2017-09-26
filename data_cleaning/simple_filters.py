''' A collection of filtering functions. '''

def no_SKI_exon1(df):
    ''' Remove SKI exon 1 variants from the dataframe'''
    no_SKI = df[df.apply(lambda x: x['Symbol'] != "SKI" or x['Exon'] != "1/7", axis=1)]
    no_SKI = no_SKI[no_SKI.apply(lambda x: x['Symbol'] != "SKI" or x['Exon'] != "01-Jul", axis=1)]
    return no_SKI

def truly_pathogenic(df):
    ''' filter for validated pathogenic and likely pathogenic variants'''
    return validated_only(pathogenic_only(df))

def pathogenic_only(df):
    ''' Filter for pathogenic and likely pathogenic variants only '''
    pathogenic_new = df[(df['New Category'] == "Pathogenic") | 
                        (df['New Category'] == "Likely Pathogenic") | 
                        (df['New Category'] == "Pathogenic/Likely Pathogenic")]
    return pathogenic_new

def validated_only(df):
    ''' Filter out non validated rows from the df'''
    valid = df[df.apply(lambda x: x['validation'] == 1, axis=1)]
    return valid

def check_for_unwanted(df):
    ''' Print the number of samples containing SKI exon 1 and low AB within a given df
    '''
    num_ski = df[((df['Symbol'] == 'SKI') & (df['Exon'] == '1/7')) |
                 ((df['Symbol'] == 'SKI') & (df['Exon'] == "01-Jul"))].shape[0]
    num_ab = df[df['AB'] < 0.3].shape[0]
    print("{} of SKI exon 1 identified and {} of variants with a low AB".format(num_ski, num_ab))

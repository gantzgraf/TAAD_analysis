''' create_new_most_damaging() and its helper functions are designed to replace a samples/patients most damaging variant that is considered to be a false positive with the next most highly ranked damaging variant. '''
import pandas as pd
from data_cleaning import rename

def create_new_most_damaging(old_most_dam, all_vars, AB=0.3, Gene="SKI", Exon="1/7", Date="01-Jul"):
    ''' Replace the most damaging variant for each patients variant whom
        does not pass the allele balance threshold or whoms variant is
        within a known false positive gene and exon. If the existing most 
        damaging variant is the only variant for that patient, then said 
        variant will remain as the most damaging.
        
    Args:
        old_most_dam: existing dataframe which details the most damaging variant for each patient
        all_vars: a dataframe which contains all variants associated with the patients detailed in old_most_dam
        AB: allele balance minimum threshold
        Gene: gene in which a known false positive lies within 
        Exon: exon of said gene in which a known false positive lies within
        Date: converting a xsxl to csv results in exon nums turning into dates i.e. 1/7 becomes 01-Jul. This ensures said exons are filtered if this is the case.

    Returns:
        The old_most_dam df where the next most damaging variant has been 
        put in place of the old most damaging variant that did not pass 
        the allele balance threshold or was within a known false positive
    '''
    # drop duplicates (dropping is fine as they have been sorted by score)
    all_alt_vars = get_other_variants(old_most_dam, all_vars, AB, Gene, Exon, Date)
    alt_most_dam = all_alt_vars.drop_duplicates(['Sample'])
    alt_most_dam['new_md'] = "Y"   # mark sample/variants 

    # rename columns so they match with alt_most_dam 
    old_most_dam = rename.rename_columns(old_most_dam)

    # append the newly selected most damaging, sort so these variants appear above old most damaging vars and drop duplicates so only new most damaging remain.
    append_most_dam = old_most_dam.append(alt_most_dam)
    append_most_dam = append_most_dam.sort_values(['Sample', 'new_md'])
    
    new_most_dam = append_most_dam.drop_duplicates('Sample')
    
    return new_most_dam


def get_other_variants(most_damaging, all_var, AB, Gene, Exon, Date):
    ''' Get a list of all sample names that contain a given false positive variant 
        or a variant which does not pass the threshold of the allele balance and 
        use it to get all other variants asociated with said sample.

    Args:   
        most_damaging: csv containg most damaging variants per sample
        all_var: all called variants assocaited with each sample referred to in most_damaging patients

    Returns:
        A modified all_vars df which has all the alternative variants for each
        sample within the most_damaging df except variants below the AB threshold
        and within false positives. Returned df is sorted by samle name and variant
        score.
        
    '''
    # copy and rename columns
    df = most_damaging.copy()
    df = df.fillna("-")
    df = rename.rename_columns(df)

    # remove samples that have all NaN entries in the fields of interest
    df['all_nan'] = df.apply(lambda x: "Y" if (x['Symbol'] == "-" and x['Exon'] == "-" 
                             and x['AB'] == "-") else "N", axis=1)
    df = df[~df['all_nan'].str.contains("Y")]

    # convert AB to numeric
    df['AB'] = pd.to_numeric(df['AB'], errors='coerce')

    # filter for variants with SKI exon1 and AB < 0.3
    df = identify_unwanted(df, AB, Gene, Exon, Date)   
    df = df[df['TEST'].str.contains("LOW", na=False)]
                                               
    # list of all sample names 
    l = df['Sample'].tolist()
    
    # copy and rename columns
    all_vars = all_var.copy()
    all_vars = rename.rename_columns(all_vars)
    
    # filter for only rows that contain sample name in the given list
    all_vars['cross'] = all_vars['Sample'].isin(l)
    all_vars = all_vars[all_vars['cross'] == True]
    
    # filter for variants with AB > 0.3 or aren't SKI exon 1
    all_vars = identify_unwanted(all_vars, AB, Gene, Exon, Date)
    all_vars = all_vars[~all_vars['TEST'].str.contains("LOW", na=False)]

    # sort in order of sample name and score
    all_vars = all_vars.sort_values(['Sample','Score'], ascending=False)

    return all_vars 
 
def identify_unwanted(df, AB, Gene, Exon, Date):
    ''' Identify and mark rows with variants which have a allele 
        balance less than the given threshold or a variant within 
        a known false positive
    '''
    mask = (((df.Symbol.str.contains(Gene)) & (df.Exon.str.contains(Exon))) 
            | ((df.Symbol.str.contains(Gene)) & (df.Exon.str.contains(Date)))
            | (df.AB < AB))
    df.ix[mask, 'TEST'] = "LOW"
    
    return df

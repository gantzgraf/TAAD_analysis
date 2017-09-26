''' A collection of functions to rename columns and cell entries'''
import numpy as np

def rename_columns(df):
    ''' Clean up the columns names.'''
    df = df.rename(columns={'validated?(1=yes,2=no,0=not_done)': 'validation',
                            
                            'proven family history of aortic disease '
                            '(yes/no)': 'proven family_history',
                            'Proven Family History of Aortic Disease '
                            '(Yes/No)': 'proven family_history', 
                            
                            'probable family history of aortic '
                            'disease (yes/no)': 'probable family_history',
                            'Probable Family History of Aortic Disease '
                            '(Yes/No)': 'probable family_history',


                            'gender (m/f)': 'Gender', 'gender (male/female)': 'Gender', 
                            'gender (male/ female)': 'Gender', 'Gender (Male/Female)': 'Gender', 'Gender (Male/ Female)': 'Gender',

                            'primary diagnosis – indication for surgery: aneurysm / '
                            'dissection / rupture / transection / imh / pau': 'primary diagnosis',
                            'Primary Diagnosis  Presenting Indication (for surgery): Aneurysm / Dissection /  Transection / IMH / PAU': 'primary diagnosis',
                            'Primary Diagnosis  indication for surgery: Aneurysm / Dissection / Rupture / Transection / IMH / PAU': 'primary diagnosis',
                            'primary diagnosis – presenting indication (for surgery):'
                            ' aneurysm / dissection /  transection / imh / pau': 'primary diagnosis',
                            'Primary Diagnosis – indication for surgery: Aneurysm / '
                            'Dissection / Rupture / Transection / IMH / PAU': 'primary diagnosis',

                            'age at diagnosis (any aortic disease)': 'age at diagnosis', 
                            'Age at Diagnosis (any aortic disease)': 'age at diagnosis',

                            'location of primary diagnosis – ascending/arch/descending'
                            '/thoracoabdominal': 'location of primary diagnosis',
                            'Location of Primary Diagnosis – Ascending, Arch, Descending,'
                            ' Thoracoabdominal, Infrarenal': 'location of primary diagnosis',
                            'location of primary diagnosis – ascending, arch, descending, '
                            'thoracoabdominal, infrarenal': 'location of primary diagnosis',
                            'Location of Primary Diagnosis  Ascending, Arch, Descending, Thoracoabdominal, Infrarenal': 'location of primary diagnosis',
                            
                            'aortic size at diagnosis (primary location/earliest'
                            'measurement) (cm)': 'aortic size at diagnosis',
                            
                            'aortic size at diagnosis (primary location, earliest measureme'
                            'nt) (cm)': 'aortic size at diagnosis (cm)',
                            
                            '"aortic size at diagnosis (primary'
                            'location': 'aortic size at diagnosis (cm)', 
                            
                            'age at time of surgery': 'age_at_surgery',
                            'Age at Time of Surgery': 'age_at_surgery',
                            
                            'maximal aortic size  (cm) ': 'maximal aortic size (cm)',
                            'Maximal Aortic Size  (cm) ': 'maximal aortic size (cm)', 
                            'maximal aortic size (cm) ': 'maximal aortic size (cm)',
                            'Maximal Aortic Size (cm) ': 'maximal aortic size (cm)', 

                            '"aortic size at diagnosis (primary location': 'aortic size at diagnosis (cm)',
                            'aortic size at diagnosis (primary location/earliest measurement) (cm)': 'aortic size at diagnosis (cm)',
                            'Aortic Size at Diagnosis (primary location, earliest measurement) (cm)': 'aortic size at diagnosis (cm)',

                            'known syndrome - marfan / lds / eds': 'Known Syndrome',
                            'known syndrome - Marfan / LDS / EDS': 'Known Syndrome',

                            'extra-aortic  aneurysmal disease': 'extra-aortic aneurysmal disease',

                            'Mendelian ID': 'Sample',
                            'YALE Coding': 'Sample'
                            })
    return df

def rename_entries(df):
    ''' Rename many of the cell entries to increase consistancy'''
    df['Gender'] = df['Gender'].replace("male", "Male").replace('female', 'Female').fillna('-').replace("M","Male").replace("F", "Female")
    df['proven family_history'] = df['proven family_history'].fillna('unknown').replace("-", 'unknown').replace("N", "no"). replace("Y", "yes").replace('y', 'yes')
    df['probable family_history'] = df['probable family_history'].fillna('unknown').replace("-", 'unknown').replace("N", "no").replace("Y", "yes").replace('y', 'yes')
    df = replace_series_strings(df, col='primary diagnosis', dic={'aneurysm': 'Aneurysm', 'dissection': 'Dissection'}, substring=False)
    df = replace_series_strings(df, col='location of primary diagnosis', dic={'hemi': 'Arch', 'Hemi': 'Arch', 'thoracoabdominal': 'Thoracoabdominal'}, substring=True) 
    df['location of primary diagnosis'] = df['location of primary diagnosis'].replace("?", np.nan)
    return df                                    

def replace_series_strings(df, col, dic, substring):
    ''' Replace the keys with the items of the given 
        dictionary for all strings or substrings in a
        given column

    Args:
        col: column name to replace strings
        dic: dictionary where the key is the string to replace with the item
        substrings: search and replace for either substrings (True) or exact strings (False)

    Returns:
        dataframe with the given column having all the
        entries identified as the key in the given dict
        replaced with the item in said dict
    '''
    if not isinstance(substring, bool):
        raise TypeError("substring argument must equal True or False")

    for string, correction in dic.items():
        if substring is True:
            df[col] = df[col].str.replace(string, correction)
        elif substring is False:
            df[col] = df[col].replace(string, correction, regex=True)

    return df

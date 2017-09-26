''' create_new_columns() and its helper functions are used to produce new columns.'''

def create_new_columns(df, three_categories=True):
    ''' Utilise the below functions to alter existing and create 
        new columns required for the TAAD analysis.
    '''
    df['family_history'] = df.apply(determine_family_history, axis=1)
    df['Age Group'] = df.apply(determine_young_old, axis=1)
    df['location of primary diagnosis'] = df.apply(plus2slash,axis=1)
    df['location of primary diagnosis'] = df['location of primary diagnosis'].astype(str)
    df['simple location of primary diagnosis'] = df.apply(simple_location_primary_diagnosis, axis=1)
    if three_categories:
        df['New Category'] = df.apply(lambda x: determine_new_category(x, three_categories=True), axis=1)
        df['New Category code'] = df['New Category'].replace({'Likely Benign / No Variant': 1, 
                                                              'VUS': 2, 
                                                              'Pathogenic/Likely Pathogenic': 3})
    else:
        df['New Category'] = df.apply(determine_new_category, axis=1)
        df['New Category code'] = df['New Category'].replace({'Likely Benign / No Variant': 1, 
                                                              'Pathogenic/Likely Pathogenic': 2})
    return df

def determine_family_history(x):
    ''' Decide what the family_history column will contain for 
        each row by evaluating both probable and proven history.
    '''
    both = (x['probable family_history'].lower(), x['proven family_history'].lower())
    if 'yes' in both or 'y' in both:
        return 'yes'
    elif 'Marfan' in both or 'marfan' in both:
        return 'yes'
    elif 'unknown' in both and both[0] != 'no' and both[1] != 'no':
        return 'unknown'
    elif both[0] == 'no' or both[1] =='no':
        return 'no'
    else:
        print("Unable to determine family_history for{}\n{}".format(x['Sample'], both))
     
def determine_young_old(x):
    ''' Designate a sample young or old depending upon the age
        of the patient when diagnosed
    '''
    x = x['age at diagnosis']
    if x < 50:
        return 'Under 50'
    elif x >= 50:
        return 'Over 50'
    else:
        return "-"
    
def determine_new_category(x, three_categories=False):
    ''' Depending upon what the Category value, 
        decide which New Category value to return
    Args:
        three_categories: if True then P/LP, VUS & Likely Benign / No Variant
    '''    
    if three_categories is True:
        if (x['Symbol'] == 'SKI' and x['Exon'] == "1/7") or (x['Symbol'] == 'SKI' and x['Exon'] == "01-Jul"):
            return "Likely Benign / No Variant"
        elif x['Category'] == "Pathogenic" or x['Category'] == "Likely Pathogenic":
            return "Pathogenic/Likely Pathogenic"
        elif x['Category'] == 'Not Classified':
            return "Likely Benign / No Variant"
        elif x['Category'] == 'Uncertain Significance':
            return 'VUS'
        elif str(x['Category']) == 'nan' or x['Category'] == 'No Variant':
            return "Likely Benign / No Variant"
        else:
            return x['Category']
    else:
        if x['Symbol'] == 'SKI' and x['Exon'] == "01-Jul":
            return "Likely Benign / No Variant"
        elif x['Category'] in ('Uncertain Significance', 'Not Classified'):
            return 'Likely Benign / No Variant'
        elif x['Category'] == "Pathogenic" or x['Category'] == "Likely Pathogenic":
            return "Pathogenic/Likely Pathogenic"
        else:
            return x['Category']

def simple_location_primary_diagnosis(x):
    ''' Simplify the location of the primary diagnosis
        to two options: Ascending and Not Ascending
    '''
    x = x['location of primary diagnosis']
    if "Ascending" in x:
        return "Ascending"
    elif "Ascending" not in x:
        return "Not Ascending"
    else:
        return x
        
def plus2slash(x):
    ''' Rename a string so that any plus characters or ambigouse names are converted to an approriate
        charcter/string in the location of primary diagnosis fields.
    '''
    location = x['location of primary diagnosis']
    if isinstance(location, str):
        # This variables contents is cancer for the eyes - change it perhaps using a dict
        convert = location.replace("+", "/").replace(" ", "").replace("ascending", "Ascending").replace("descending", "Descending").replace("arch", "Arch").replace("hemi", "Hemi").replace("Thoracoabdomen","Thoracoabdominal").replace("Descending/PAU", "Descending").replace("Arch-stageI", "Arch").replace("aneurysm", "-").replace("Descendingthoracic", "Descending").replace("AorticValveReplacement", "-")
        return convert
    else:
        return location

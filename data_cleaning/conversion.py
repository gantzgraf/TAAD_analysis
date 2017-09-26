''' Functions which convert pandas DataFrame columns to various dtypes.'''
import pandas as pd

def convert2numeric(df, cols):
    ''' convert the parsed columns to numeric type'''
    df[cols] = df[cols].apply(lambda x: pd.to_numeric(x, errors='coerce'))
    return df

def convert_these_category(df, three_categories):
    ''' Convert multiple columns to categories. This is required to maintain
        the ordering of the xtick labels in downstream visualisations.
    Args:
        df: DataFrame
        three_categories: boolean
    '''
    df['family_history'] = convert2category(df['family_history'], ['yes', 'no', 'unknown'])
    df['Gender'] = convert2category(df['Gender'], ['Male', 'Female'])
    df['Age Group'] = convert2category(df['Age Group'], ['Under 50', 'Over 50'])
    labels = get_pathogenicity_labels(three_categories)
    df['New Category'] = convert2category(df['New Category'], labels) 
    return df

def get_pathogenicity_labels(three_categories):
    if three_categories:
        labels = ['Pathogenic/Likely Pathogenic', 
                  'VUS', 
                  'Likely Benign / No Variant']
    else:
        labels = ['Pathogenic/Likely Pathogenic', 
                  'Likely Benign / No Variant']
    return labels

def convert2category(column, label_order):
    ''' convert column to a category and set the order of the labels'''
    convert = column.astype("category")
    structured_cat = convert.cat.set_categories(label_order)
    return structured_cat

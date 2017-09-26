''' demogrphics_table() and its helper function retrieve all demographics data for each cohort and outputs as a csv table'''
import numpy as np
import pandas as pd
import collections
import operator

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_colwidth', -1)

OPS = { 
    ">": operator.gt,
    "<": operator.lt,
    ">=": operator.ge,
    "<=": operator.le,
    "==": operator.eq,
    "=": operator.eq,
    "!=": operator.ne,
    "&": operator.and_,
    "|": operator.or_ 
}

def demographics_table(df, outfile=None):
    ''' Get all demographic data for each cohort as a
        Series and combine them into a DataFrame.
    Args:
        df: cleaned most damaging dataframe
    '''
    series_storage = []
    for cohort_name in ['Yale', 'UK', None]:
        cohort_passed, df_passed = passed_depth_threshold(df, cohort_name)
        cohort_dict = get_demographics_dict(cohort_passed, df_passed)
        cohort_series = pd.Series(cohort_dict)
        series_storage.append(cohort_series)
    demographics = pd.concat(series_storage, axis=1)
    demographics.columns = ['Yale Cohort', 'UK Cohort', 'Whole Cohort']
    if outfile:
        demographics.to_csv(outfile)
    return  demographics

def passed_depth_threshold(df, cohort_name=None):
    ''' Returns most damagining dataframe that passed the
        sequencing depth threshold for the whole cohort
        and for a given sub-cohort (UK or Yale)
    '''
    passed = df[df['Depth'] != 'LOW']
    if cohort_name:
        cohort = passed[passed['cohort'] == cohort_name]
    else:
        cohort = passed
    return (cohort, passed)

def get_demographics_dict(cohort, passed):
    ''' Get all demographic data for the given cohort DataFrame
        and store inside a dict.
    Args:
        cohort: most damaging dataframe 
    '''
    d = collections.OrderedDict()
    d['Demographics'] = ''
    d['Numbers (%)'] = get_patient_numbers(cohort, passed)
    d['Age at Diagnosis, Median'] = cohort['age at diagnosis'].median()
    d['Age at Diagnosis, Min'] = cohort['age at diagnosis'].min()
    d['Age at Diagnosis, Max'] = cohort['age at diagnosis'].max()
    d['Male (%)'] = get_numbers(cohort, 'Gender', 'Male')
    d['Female (%)'] = get_numbers(cohort, 'Gender', 'Female')
    # ETHNICITY
    d['Probable/Proven Family History(%)'] = get_numbers(cohort, 'family_history', 'yes')
    d['No Family History (%)'] = get_numbers(cohort, 'family_history', 'no')
    d['Undergone Aortic Surgery'] = percent_operated(cohort)
    d['Primary Aortic Pathology'] = ''
    d['Aneurysm (%)'] = get_numbers(cohort, 'primary diagnosis', 'Aneurysm')
    d['Dissection (%)'] = get_numbers(cohort, 'primary diagnosis', 'Dissection')
    d['IMH/PAU (%)'] = get_numbers(cohort, 'primary diagnosis', 'IMH|PAU')
    d['Rupture (%)'] = get_numbers(cohort, 'Rupture (Y/N)', 'Y')
    d['Primary Anatomical Presentation'] = ''
    d['Ascending/Arch (%)'] = get_numbers(cohort, 
                                          'location of primary diagnosis',
                                          ".*Ascending.*|.*Arch.*")
    d['Descending/Thoracoabdominal (%)'] = get_numbers(cohort,
                                                       'location of primary diagnosis',
                                                       'Descending|Thoracoabdominal'
                                                       '|Infrarenal')
    d['Aortic Size'] = ''
    d['Maximum Aortic Diameter (cm), Median'] = cohort['maximal aortic size (cm)'].median()
    d['Maximum Aortic Diameter (cm), Min'] = '{0:1f}'.format(cohort['maximal aortic size (cm)'].min())
    d['Maximum Aortic Diameter (cm), Max'] = cohort['maximal aortic size (cm)'].max()
    d['Maximum Aortic Diameter < 5.5cm (%)'] = get_numbers(cohort,
                                                           'maximal aortic size (cm)',
                                                           5.5, op="<")
    d['Known Syndrome'] = ''
    d['MFS (%)'] = get_numbers(cohort, 'Known Syndrome', 'Marfan')
    d['LDS (%)'] = get_numbers(cohort, 'Known Syndrome', 'LDS')
    d['EDS (%)'] = get_numbers(cohort, 'Known Syndrome', 'EDS')
    d['Other (%)'] = get_other_syndromic_numbers(cohort)
    return d

def get_patient_numbers(cohort, passed):
    '''Get the number and percentage of patients which 
       passed the depth thershold from all cohorts
    '''
    return '{0:d}({1:.1f})'.format(len(cohort), len(cohort)/len(passed)*100)

def get_numbers(df, col, val, op='='):
    ''' Get the number and percentage of patints with/without
        a given phenotype feature

    Args: 
        df: dataframe of a given cohort
        col: column to query
        val: value to filter for in column
        op: operator symbol

    Return:
        str containing number of vals in col and 
        percentage of vals against total patient number

    Notes:
        The below example searches for samples/patients who have 
        either Marfan, LDS or EDS in their Known Syndrome column
        and returns the corresponding number and percentage.
        
            get_numbers(df, 'Known Syndrome', 'Marfan|LDS|EDS')
    '''
    calc = OPS[op]
    if isinstance(val, float) or (isinstance(val, str) and '|' not in val):
        df_num = len(df[calc(df[col], val)])
    elif '|' in val:
        if calc != operator.eq:
            raise TypeError('Operators other than eq will not'
                            'work when "|" is present in val')
        df_num = np.sum(df[col].str.match(val))
    return '{0:d} ({1:.1f})'.format(df_num, df_num/len(df)*100)

def percent_operated(df):
    ''' Get the number and percentage of patients that have had some 
        form of aortic operation.
    '''
    operated = df[['No.of Aortic Operations - Endovascular',
                   'No.of Aortic Operations - Open',
                   'No.of Aortic Operations - Hybrid']]
    operated_patients = len(df[(operated > 0).any(axis=1)])
    return '{0:d} ({1:.1f})'.format(operated_patients, operated_patients/len(df)*100)

def get_other_syndromic_numbers(df):
    ''' Get number and percentage of patients/samples
        without MFS, EDS or LDS.
    '''
    normal = df['Known Syndrome'].fillna('-')
    other = df[~normal.str.contains('Marfan|LDS|EDS|-')]
    return '{0:d} ({1:.1f})'.format(len(other), len(other)/len(df)*100)


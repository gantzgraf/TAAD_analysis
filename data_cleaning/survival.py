import pandas as pd
import os

FILE_PATH = os.path.dirname(os.path.abspath("__file__"))+"/"

def merge_survival_data(df, yale_survival):
    ''' Merge most_damaging data with Yale survival data'''
    survival = pd.read_csv(yale_survival, encoding="ISO-8859-1")
    # only use Sample data presenet within parsed df
    shared_survival = survival[survival['Sample'].isin(df['Sample'].unique())]
    shared_survival = shared_survival[['Sample', 'Long-term mortality (0=no, 1=yes)',
                                       'Type of surgery (0=elective, 1=urgent/emergent)',
                                       'Peri-operative morality (0=no, 1=yes)']]
    merged = pd.merge(df, shared_survival, how='outer')
    return merged

''' Output a list of individuals with three or more risk factors. '''

import data_cleaning.simple_filters as sf
import pandas as pd
import operator

OPS = {
    ">": operator.gt,
    "<": operator.lt,
    ">=": operator.ge,
    "<=": operator.le,
    "==": operator.eq,
    "=": operator.eq,
    "!=": operator.ne,
    "&": operator.and_,
    "|": operator.or_,
    "notna": pd.notna,
}
def at_risk_indvidiuals(df, outfile=None):
    ''' 
        Output CSV of patients with three or more risk factors 
        suggesting a genetic cause.
    '''
    risk_groups = (
                    ('Age Group', '=', 'Under 50'),
                    ('family_history', '=', 'yes'),
                    ('location of primary diagnosis', '=', 'Ascending'),
                    ('Known Syndrome', 'notna', None)
    )
    scores = [0] * len(df) 
    for i in range(len(df)): #can we do this without looping over dataframe?
        for rg in risk_groups:
            pheno, op, val = rg
            if val is not None:
                if OPS[op](df[pheno][i], val):
                    scores[i] += 1
            else:
                if OPS[op](df[pheno][i]):
                    scores[i] += 1
    df['risk_score'] = scores
    cols_to_keep = ['New Category', 'Category', 'Age Group', 'family_history', 
                    'location of primary diagnosis', 'Known Syndrome', 
                    'risk_score']
    at_risk = df[df['risk_score'] >= 3][cols_to_keep]
    if outfile:
        at_risk.to_csv(outfile)

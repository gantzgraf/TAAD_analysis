''' risk_ratio_table() and its helper functions construct a table detailing total patients, '''
import data_cleaning.simple_filters as sf
from scipy import stats
import pandas as pd
import numpy as np
import operator
import math

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

def risk_ratio_table(df, outfile=None):
    ''' Construct a risk ratio table for phenotypes
        associated with increased mortality in TAAD
        patients.
    '''
    # should include LDS with Marfan
    test_groups = [
        ('Age Group', '=', 'Under 50', 'Over 50'), 
        ('family_history', '=', 'yes', 'no'),
        ('Gender', '=', 'Male', 'Female'),
        ('location of primary diagnosis', '=', 'Ascending', None),
        ('primary diagnosis', '=', 'Dissection', None),
        ('maximal aortic size (cm)', '<=', 5, None),
        ('Long-term mortality (0=no, 1=yes)', '=', 0, 1),
        ('Known Syndrome', 'notna', None, None)
    ]

    all_table_rows = get_table_rows(df, test_groups)
    rr_table = pd.DataFrame(all_table_rows)
    rr_table.columns = ['','Total', 
                        'Number/Percentage with a Pathogenic or Likely '
                        'Pathogenic Variant Validated by Sanger', 
                        'RR(95% CI)', 'P-Value']
    rr_table = rr_table.set_index('')
    rr_table = rename_index(rr_table)
    print('\nINFO: The risk ratio row for syndromic is calculated for MFS only')
    if outfile:
        rr_table.to_csv(outfile)
    return rr_table

def get_table_rows(df, test_groups):
    ''' Get the total patient number, pathogenic variant number/percentage,
        risk ratio and p-value for each phenotype feature and conditions
        given and return as a nested list. 

    Args:
        df: DataFrame
        test_groups: list of tuples containing:
            (phenotype, str operator, exposed phenotype group, 
             non-exposed phenotype group)
    '''
    table_rows = []
    for group in test_groups:
        phenotype, operation, exposed, not_exposed = group
        op = OPS[operation]
        rr, pval = pathogenic_risk_ratio(df, phenotype, op, exposed, not_exposed)
        plp_num = pathogenic_number(df, phenotype, op, exposed)
        if exposed is None:
            total = len(df[op(df[phenotype])])
        else:
            total = len(df[op(df[phenotype], exposed)])
        table_rows.append([phenotype, total, plp_num, rr, pval])
    return table_rows

def pathogenic_number(df, col, op, value):
    ''' Returns the number and percentage of pathogenic 
        variants that meet the given phenotype conditions.
    Args:
        df: DataFrame
        col: phenotype column
        op: operator string
        value: value col will be filtered for
    '''
    if value is None:
        group = df[op(df[col])]
    else:
        group = df[op(df[col], value)]
    plp_group = sf.validated_only(sf.pathogenic_only(group))
    percent_plp = len(plp_group) / len(group) * 100
    return '{0:d} ({1:.1f})'.format(len(plp_group), percent_plp)

def pathogenic_risk_ratio(df, col, op, exposed_val=None, not_exposed_val=None):
    ''' Determine whether a phenotypic trait increases the 
        liklihood of carrying a causitive mutation for the
        disorder by performing riask ratio assessment
        and fishers-exact test.
    Args:
        df: DataFrame
        col: phenotype column
        op: operator function
        exposed_val: value of exposed group in column (optional - will 
                     test for True/False if None)
        not_exposed_val: value of non-exposed group (optional)
    Notes:
        The below example input will return the risk ratio
        of whether having a family history increase the likelihood
        of carrying a causitive TAAD mutation compared to this
        without a family history of the disorder.
          get_risk_ratio(df, 'family_history', '=', 'yes, 'no')
    Returns:
        Risk ratio and p-value in a tuple
    '''    
    # exposed and non-exposed groups
    if exposed_val is None:
        exposed = df[op(df[col])]
        not_exposed = df[~op(df[col])]
    else:
        exposed = df[op(df[col], exposed_val)]
        if not_exposed_val:
            not_exposed = df[op(df[col], not_exposed_val)]
        else:
            not_exposed = df[~op(df[col], exposed_val)]
        
    # Disease exposed and non-exposed groups
    d_exposed = len(sf.validated_only(sf.pathogenic_only(exposed)))
    d_not_exposed = len(sf.validated_only(sf.pathogenic_only(not_exposed)))
    
    # Non-Disease exposed and non-exposed groups
    nd_exposed = len(exposed) - d_exposed
    nd_not_exposed = len(not_exposed) - d_not_exposed
    
    rr = risk_ratio(a=d_exposed,
                    b=nd_exposed,
                    c=d_not_exposed,
                    d=nd_not_exposed)
    pval = stats.fisher_exact([[d_exposed, nd_exposed],
                              [d_not_exposed, nd_not_exposed]])[1]
    return (rr, pval)

def risk_ratio(a, b, c, d):
    ''' Calculate the risk ratio and confidence intervals.

      Algorithm extracted from the following sources:
          http://tinyurl.com/jsjvyra
          http://tinyurl.com/hr3lapx

    Args:
                         Disease     No Disease
      Exposed Group         a             b
      Not Exposed Group     c             d

    Returns:
      String containing the risk ratio value along with
      the lower and upper confidence intervals in brackets

    Meaning:
      RR > 1.0: the exposed group is at a higher risk of
      developing the disorder.

      CI range does not overlap 1.0, the result is statistically
      significant.

    Example:
      We have an exposed group (family history) composed of 100
      patients and non-exposed group (no known family history)
      composed of 200 patients. 20 patients in the exposed group
      have the disease and 25 patients in the non-exposed group
      have the disease. Thus: a=20, b=80, c=25, d=175.
    '''
    # calculate the risk ratio
    RR = (a/(a+b)) / (c/(c+d))

    # natural log of the RR
    LN = np.log(RR)

    # get the standard error
    SE = math.sqrt((1/a)+(1/c)-(1/(a+b))-(1/(c+d)))

    # calculate the lower and upper confidence intervals
    CI_low = np.exp(LN-1.96*SE)
    CI_up = np.exp(LN+1.96*SE)

    return "{:.2f} ({:.2f}-{:.2f})".format(RR, CI_low, CI_up)


def rename_index(df):
    df = df.rename(index={'Age Group': 'Young Age, < 50',
                          'family_history': 'Known or Probable Family History',
                          'Gender': 'Male',
                          'location of primary diagnosis': 'Ascending Aorta',
                          'primary diagnosis': 'Presence of Dissection',
                          'maximal aortic size (cm)': 'Large Aortic Diameter (>5cm)',
                          'Long-term mortality (0=no, 1=yes)': 'Short-Term Mortality',
                          'Known Syndrome': 'Syndromic'})
    return df

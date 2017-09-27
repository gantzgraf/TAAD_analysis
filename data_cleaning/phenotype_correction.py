''' phenotype_resolver() and it's helper functions checks and resolves duplicate samples that have different values in their phenotype fields '''
import numpy as np

def phenotype_resolver(df, phenotype_columns):
    ''' Differences in phenotype data between duplicate samples are corrected/reported
        and any samples that subsequently have no phenotype data are removed.

    Args:
        phenotype_columns: a list of phenotype columns 
    '''
    # duplicate samples with different phenotype information are dealt with here
    df = duplicate_column_checker(df, phenotype_columns)

    # get index numbers for all phenotype columns and filter out rows that have no phenotype data
    pheno_col_ix = [df.columns.get_loc(x) for x in phenotype_columns]
    df['Any Data'] = df.apply(lambda x: phenotype_data(x, pheno_col_ix), axis=1)
    print("\nINFO: Number of samples without phenotype data: {}\n".format(len(df[df['Any Data'] == "no"])))
    df = df[~df['Any Data'].str.contains("no")]

    return df

def duplicate_column_checker(df, columns_names, order=False, recurs=2, 
                             dup_ends=['_2','_3', '_pool7A', '_pool10A'],
                             column="Sample"):
    ''' Identify duplicate samples and verify whether they have the same data stored in the given columns.
        If one duplicate has NaN in its phenotype fields then copy the phenotype data from the
        other duplicate sample. If there are still differences between them, then report to user.

    Args:
      column_names: a list of column names to be investigated for differences between duplicates
      order: ascending - True or False
      recurs: utilised to stop the recursive function
      dup_end: characters which seperate the original and duplicate sample
      column: column in which to search & identify whether samples are duplicates

    Returns:
      a dataframe in which the duplicates differences in the given columns have been resolved
      and/or communicated to the user
    '''
    if recurs < 1:
        return df

    # replace all dashes with NaN
    df = df.replace('-', np.nan)

    # Identify which samples are duplicates and fill in a new column with the original samples name. This way all duplicates 
    # have the orginal sample name in its row.
    for dup in dup_ends:
        cond = ((df[column].str.endswith(dup)) & (df[column].str.len() > 4))
        df.ix[cond, 'same'] = df.ix[cond, column].str[:-(len(dup))]

    # replace NaN in same by entries in column
    df.same = df.same.combine_first(df[column])

    # sort columns first on same then column. column is reversed upon the funcs recursion
    df = df.sort_values(by=['same', column], ascending=[False, order]).reset_index()

    # get the clumn indexes for the given column names
    col_ix = [df.columns.get_loc(x) for x in columns_names]

    # assess whether the sample names in a row and the following row are duplicates and, if so, forward fill the fields 
    # described in column names
    for num in range(0, df.shape[0]-1):

        current_sam = df['same'].shift(-(num))[0]
        next_sam = df['same'].shift(-(num+1))[0]

        if current_sam == next_sam:

            df.iloc[num:num+2, col_ix] = df.iloc[num:num+2, col_ix].ffill(limit=1)

    # report duplicates that still have different phenotype data
    if recurs == 1:
        duplicates_column_diff(df, col_ix, column)

    # remove column index
    df.drop(['index'], axis=1, inplace=True)

    # do the same but reverse the order of secondary order
    return duplicate_column_checker(df, columns_names, order=True, recurs=recurs-1, dup_ends=dup_ends)

def duplicates_column_diff(df, col_ix, column):
    ''' Find duplicate samples and communicate whether they have different in values in the given column indexes.

    Args:
      col_ix: a list of column indexes to be investigated for differences between duplicate samples
      column: column name which will used to search and identify whether samples are duplicates

    Returns:
      a CSV file containing differnces in the given column indexes along with the column names in which there are 
      differences present
    '''
    # store names of duplicates that have differences in phenotypes
    diff_entries = []

    # unfortunately this would not work within duplicate_column_checker()
    for num in range(0, df.shape[0]-1):

        current_sam = df[column].shift(-(num))[0]
        next_sam = df[column].shift(-(num+1))[0]

        current_entries = df.iloc[num, col_ix].tolist()
        next_entries = df.iloc[num+1, col_ix].tolist()

        header_entries = df.columns[col_ix]

        if current_sam == next_sam:

            # something here to report phenotype differences between duplicates
            if num+2 <= df.shape[0]:

                # required as NaN does not equal NaN
                df = df.fillna("-")

                if current_entries != next_entries:
                    differences = [c for x, z, c in zip(current_entries, next_entries, header_entries)
                                   if str(x) != str(z)]
                    diff_entries.append([current_sam, next_sam] + differences)

    # save a file which describes the duplicate samples that have different info in the given column indexes
    # pd.DataFrame(data=diff_entries).to_csv(file_path+"error_diff_duplicate_entries_data.csv")

def phenotype_data(x, col_ix):
    ''' Check whether we have any data in the fields in the given column indexes.
    Args:
      col_ix: a list of column indexes to be investigated for differences between duplicate samples
    '''
    # start from column (index) and convert to a list and filter out duplicate values and NaNs
    pheno_cells = x[col_ix].tolist()
    all_row_cell_values = set([z for z in pheno_cells if str(z) != 'nan'])

    if str(x['validation']) == "1" and len(all_row_cell_values) == 0:
        print('{} has {} phenotype information but will not be included '
              'even though it has a validated variant: {}\n'.format(x['Sample'], len(all_row_cell_values),
                                                                    x['validation']))

    if len(all_row_cell_values) == 1 and "-" in all_row_cell_values:
        return "no"
    elif len(all_row_cell_values) == 0:
        return "no"
    else:
        return "yes"  

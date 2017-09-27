''' filter_by_depth() and it's helper functions allows one to filter genotype data by a given sequencing depth threshold'''
import numpy as np
import pandas as pd
import os

def filter_by_depth(df, depth_path, sample_column, depth_column, threshold, excluded_columns):
    ''' Alter the genotype columns to NaN for the samples that do not meet 
        the minimum depth threshold and those that still contain false positives 
        as their most damaging variants.
      
    Args:
      depth_path: dir containing depth files
      sample_column: sample column name in depth_df
      depth_column: depth column name in depth_df
      threshold: cut-off for depth number, anything below this number will have
                 its fields filled with NaN
      excluded_columns: fields in which one doesn't want to be filled with NaN
                        i.e. phenotype fields

    Returns:
      altered df where samples that are not meeting depth threshold are given np.nan within their fields
    '''
    depth_df = prepare_depth_df(depth_path)
    filtered_df = genotype_by_depth(df, depth_df, sample_column, 
                                    depth_column, threshold, 
                                    excluded_columns)
    return filtered_df

def prepare_depth_df(file_path):
    ''' Construct a dataframe which combines all of the CSVs containing
        sequencing depth information for the TAAD cohorts and recalulates 
        the depth data of interest.

    Args:
      file_path: absolute path to the input directory
    '''
    concat_depth_uk = merge_depth_data(file_path)
    concat_depth_yale = merge_depth_data(file_path, UK=False)
    depth_df_uk = recalculate_depth(concat_depth_uk)
    depth_df_yale = recalculate_depth(concat_depth_yale)

    depth_df_uk['Cohort'] = depth_df_uk.apply(lambda x: "UK", axis=1)
    depth_df_yale['Cohort'] = depth_df_yale.apply(lambda x: "Yale", axis=1)
    depth_df = pd.concat([depth_df_uk, depth_df_yale])

    # For all duplicate values in the sample_id column, get the mean values of the subsequent columns.
    # This combines duplicate samples depth info into mean values.
    depth_df = depth_df.reset_index().groupby('sample_id').mean().reset_index()

    return depth_df

def merge_depth_data(file_path, UK=True):
    ''' Find, open and concatenate the depth data for each sample in a
        given cohort.

    Args:
      UK: if True then merge the UK data, else Yale data

    Returns:
      concatenated DataFrame where each sample is reffered to multiple
      times in the index.
    '''
    if UK:
        direct = 'UK_Depth'
    else:
        direct = 'Yale_Depth'

    # open all data frames of interest and store in a list
    taadx = file_path + direct+"/depth_vs_taadx/"
    taadz = file_path + direct+"/depth_vs_taadz/"
    all_taad = [taadx+x for x in os.listdir(taadx) if x.endswith("sample_summary")] + \
               [taadz+x for x in os.listdir(taadz) if x.endswith("sample_summary")]
    all_pd = [pd.read_csv(x, sep="\t").set_index('sample_id') for x in all_taad]
    all_pd = [x.apply(pd.to_numeric) for x in all_pd]

    # concatenate the dataframes and get the mean value for columns between duplicate sample_ids
    df = pd.concat(all_pd).sort_index()

    return df

def recalculate_depth(df):
    ''' Recalulate the total reads, and %_bases_above_x for each sample
        in the depth df. This function assumes the given df contains the
        index refers to the sampl id multiple times (at least twice; one
        for the X and Z assay).
    '''
    # calculate the total number of reads for both X & Z assays for each sample; ABS Total
    total = df['total'].groupby(df.index).sum().to_frame().rename(columns={'total': 'ABS_Total'})
    df = df.join(total)
    df = df[df['ABS_Total'] != 0]

    # Use ABS_Total to calculate the percentage of reads the each assay contributed to the total 
    # num of reads for each sample
    df['Read_%'] = df.apply(lambda x: int(x['total'])/int(x['ABS_Total']), axis=1)

    # use the Read_% to recalulate the %_bases_above_x and join all the desired elements together
    df = df[['total', '%_bases_above_49', '%_bases_above_99', 'ABS_Total', 'Read_%']]
    above_49 = df.apply(lambda x: x['%_bases_above_49']*x['Read_%'], axis=1).groupby(df.index).sum().to_frame().rename(columns={0: '%_bases_above_49'})
    above_99 = df.apply(lambda x: x['%_bases_above_99']*x['Read_%'], axis=1).groupby(df.index).sum().to_frame().rename(columns={0: '%_bases_above_99'})
    total = total[total['ABS_Total'] != 0]
    df = total.join(above_49).join(above_99)

    return df

def genotype_by_depth(df, depth_df, sample_column, depth_column, threshold, excluded_columns):
    ''' Alter the columns in a row to NaN if the sample does not meet the minimum
      depth threshold or still have false positive variants as their most damaging.
    '''
    # create a depth column that details whether the depth is above or below the threshold
    df['Depth'] = df.apply(lambda x: depth_status(x, depth_df, threshold, sample_column, depth_column), axis=1)

    # get all non-phenotype column names in a list
    genotype_columns = [x for x in df.columns if x not in excluded_columns]

    # Place np.nan in the genotype fields for all samples which did not pass the x%_above_49 reads.
    # This way they will be included in the phenotype plots and demographics but not in the plots
    # concerning variant info
    df.loc[df.Depth == 'LOW', genotype_columns] = np.nan
    print("\nINFO: {} have not passed the % above 49 reads".format(df[df['Depth'] == 'LOW'].shape[0]))

    # SKI EXON 1 & AB < 0.3
    # change genotype to np.nan for these samples (this is after a next most damaging variant has been sought)
    cond = (df.AB < 0.3) | ((df.Symbol == "SKI") & (df.Exon == "1/7")) | ((df.Symbol == "SKI") & (df.Exon == "01-Jul"))
    df.loc[cond, genotype_columns] = np.nan
    df.loc[cond, 'New Category'] = "Likely Benign / No Variant"

    return df

def depth_status(most_damaging_df, depth_df, threshold, sample_column, depth_column):
    ''' Return a value depending upon whether the depth for the given sample
        meets the given threshold
    '''
    # create a list of samples that have a depth lower than the given threshold
    low_depth_list = list(depth_df[depth_df[depth_column] <= threshold][sample_column])

    # return LOW or HIGH if the sample from most damaging df is in the low depth list
    if most_damaging_df['Sample'] in low_depth_list:
        return 'LOW'
    else:
        return 'HIGH'

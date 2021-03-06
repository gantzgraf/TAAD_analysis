import data_cleaning.new_columns as nc
import data_cleaning.get_next_most_damaging as nmd
import data_cleaning.genotype_phenotype as gp
import data_cleaning.phenotype_correction as pc
import data_cleaning.filter_by_depth as fd
from data_cleaning import conversion
from data_cleaning import survival
from data_cleaning import rename
import pandas as pd

def create_most_damaging(uk_all, uk_most_damaging, uk_phenotype,
                         yale_all, yale_most_damaging, yale_phenotype,
                         yale_survival, file_path, depth_threshold=80):
    ''' Merge the most damaging genotype, phenotype and survival
        data from both cohorts and clean the merged data.

    Args:
        uk_all: All UK variants DataFrame
        yale_all: All Yale variants DataFrame
        uk_most_damaging: Most Damaging UK CSV file
        yale_most_damaging: Most Damaging Yale CSV file
        uk_phenotype: phenotype data for the UK cohort
        yale_phenotype: phenotype data for the Yale cohort
        yale_survival: survival data for Yale patients
        file_path: path to the TAAD_analysis directory
        depth_threshold: determines read depth threshold for filtering

    Returns:
        a cleaned most damaging variants dataframe which includes
        genotype, phenotype and survival data for all patients in
        all cohorts.
    '''
    # Merge Genotype-Phenotype
    uk_md, yale_md = most_damaging_dataframes(uk_most_damaging, uk_phenotype,
                                              yale_most_damaging, yale_phenotype)
    # Next Most Damaging Variant
    df = next_most_damaging_combine(uk_md, yale_md, uk_all, yale_all)
    # Merge Survival Data
    df = survival.merge_survival_data(df, yale_survival)
    survival_columns = ['Sample', 'Long-term mortality (0=no, 1=yes)', 
                        'Type of surgery (0=elective, 1=urgent/emergent)',
                        'Peri-operative morality (0=no, 1=yes)']
    phenotype_columns = get_phenotype_columns(yale_phenotype)
    md_phenotype_columns = phenotype_columns + survival_columns
    # Dtype Conversion
    numeric = ['validation', 'AB', 'age_at_surgery', 'age at diagnosis',
               'maximal aortic size (cm)', 'aortic size at diagnosis (cm)', 
               'No.of Aortic Operations - Endovascular', 
               'No.of Aortic Operations - Open', 
               'No.of Aortic Operations - Hybrid']
    df = conversion.convert2numeric(df, numeric)
    # Resolve Phenotype Difference
    df = pc.phenotype_resolver(df, phenotype_columns)
    # Correct Typos
    df = rename.rename_entries(df)
    # Resolve Data Ambiguity
    df = nc.create_new_columns(df, three_categories=True)
    # Remove Duplicates & Negative Controls
    df = df.drop_duplicates(['Symbol', 'Exon', 'Category', 'same',
                             'age at diagnosis', 'primary diagnosis',
                             'Gender', 'location of primary diagnosis',
                             'proven family_history', 'maximal aortic size (cm)',
                             'probable family_history', 'validation'])
    ### drop the remaining duplicates that have the least pathogenic variant
    df = df.sort_values('New Category code').groupby('same').last()
    df = df[~df['Sample'].str.contains("Blank|blank|ddH20|dH2O|H2O|BLANK|ddh2o")]
    # Filter by Sequencing Depth
    exclude = phenotype_columns + ['Sample', 'Depth', 'cohort', 
                                   'simple location of primary diagnosis', 
                                   'Age Group', 'family_history']
    df = fd.filter_by_depth(df=df, 
                            depth_path=file_path+"input_files/", 
                            sample_column='sample_id', 
                            depth_column='%_bases_above_49', 
                            threshold=depth_threshold, 
                            excluded_columns=exclude)
    # Recategorise Pathogenicity
    df['New Category'] = df['New Category'].fillna("Likely Benign / No Variant")
    df['Category'] = df['Category'].fillna('No Variant')
    # Convert to Categories
    df = conversion.convert_these_category(df, three_categories=True)
    return df

def most_damaging_dataframes(uk_most_damaging, uk_phenotype, yale_most_damaging, yale_phenotype):
    ''' Concatenate the Yale and UK Most Damaging variants
        data with their respective phenotype data.
    '''
    Yale_most_damaging = gp.merge_genotype_phenotype(
        yale_phenotype, yale_most_damaging)
    UK_most_damaging = gp.merge_genotype_phenotype(
        uk_phenotype, uk_most_damaging)
    Yale_most_damaging['cohort'] = "Yale"
    UK_most_damaging['cohort'] = "UK"
    return (UK_most_damaging, Yale_most_damaging)

def next_most_damaging_combine(uk_md, yale_md, uk_all, yale_all):
    ''' Replace known false positive variants in applicable
        samples in both cohorts and combine them.

    Args:
        uk_md: UK most damaging data in a DataFrame format
        yale_md: Yale most damaging data in a DataFrame format
        uk_all: UK all variants data in a DataFrame format
        yale_all: Yale all variants data in a DataFrame format
    
    Returns:
        a combined DataFrame of both cohorts most damaging variants
        with false positive variants replaced with the next most damaging
        variants
    '''
    uk_real_md = nmd.create_new_most_damaging(uk_md, uk_all)
    yale_real_md = nmd.create_new_most_damaging(yale_md, yale_all)
    combined_md = pd.concat([uk_real_md, yale_real_md])
    return combined_md

def get_phenotype_columns(p):
    ''' Open a phenotype file and return the
        column names as a list.
    Args:
        p - phenotype file
    '''
    p = rename.rename_columns(pd.read_csv(p, encoding='iso-8859-1'))
    phenotype_columns = list(p.columns.values)
    phenotype_columns.remove('Sample')
    return phenotype_columns

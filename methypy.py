import pandas as pd
import numpy as np

# D:\Python\ZBiB\Renia\PZRP 5K Ex.xlsx

input_file = input("Please select your file: ")

EVENTS = [
    '0000', '0001', '0010', '0011',
    '0100', '0101', '0110', '0111',
    '1000', '1001', '1010', '1011',
    '1100', '1101', '1110', '1111'
]


def generate_sequences(df):
    """
    Generate methylation sequences by combining values from specific columns.
    
    Parameters:
        df (pd.DataFrame): Input dataframe containing the data.
        
    Returns:
        pd.DataFrame: Dataframe with additional columns containing generated sequences.
    """
    # Find columns that start with 'A_R'
    ar_columns = [col for col in df.columns[1:-1] if col.startswith('A_R')]
    
    # Create an empty dataframe for storing sequences
    sequences_df = pd.DataFrame(index=df.index, columns=[f'R{i}' for i in range(1, len(ar_columns) + 1)])
    
    # Iterate over each column of 'A_R#' to generate sequences
    for idx, col in enumerate(ar_columns):
        # Combine values from 'Da', 'A_R#', 'Dk', and corresponding 'K_R#' columns
        sequence = df['Da'].astype(str) + df[col].astype(str) + df['Dk'].astype(str) + df[col.replace('A_R', 'K_R')].astype(str)
        
        # Update the corresponding column in the sequences dataframe
        sequences_df[f'R{idx+1}'] = sequence
    
    # Attach the sequences to the original dataframe
    df = pd.concat([df, sequences_df], axis=1)
    
    return df


def update_methylation_df(methylation_df, df, methyl_type, half_cols):
    """
    Update the methylation dataframe with counts of sequences for a given methylation type.
    
    Parameters:
        methylation_df (pd.DataFrame): Dataframe to update with counts.
        df (pd.DataFrame): Dataframe containing the sequences.
        methyl_type (str): The methylation type to filter on (e.g., 'CG', 'CXG', 'CXX').
        half_cols (int): Number of half-columns (half the total columns excluding the first column).
        
    Returns:
        pd.DataFrame: Updated methylation dataframe with sequence counts.
    """
    # Select rows that match the given methylation type
    methyl_df = df[df['MethylType'] == methyl_type]  

    for index, row in methylation_df.iterrows():
        sequence = row['Sequence']

        for i in range(1, half_cols):
            # Count occurrences of the sequence in the methylation type-specific dataframe
            count = methyl_df[methyl_df[f'R{i}'] == sequence].shape[0]
            methylation_df.at[index, f'R{i}'] = count

    return methylation_df


def add_summary_rows(df):
    """
    Add summary rows to the dataframe based on specific categories.

    Parameters:
    - df (DataFrame): The input dataframe containing methylation sequence counts.

    Returns:
    - DataFrame: The dataframe with added summary rows.
    """
    # Calculate the sum of all rows
    total_all_events = df.sum()

    # Calculate sums for specific categories
    sum_specific_events = df.loc[['0001', '0010', '0101', '0110', '1001', '1010', '1101', '1110']].sum()
    sum_dme = df.loc[['0110', '0111']].sum()
    sum_dnme = df.loc[['1001', '1011']].sum()
    sum_ce = df.loc[['0100', '1000']].sum()
    sum_snms = df.loc[['1100', '1101', '1111']].sum()
    sum_sms = df.loc[['0011', '1101']].sum()

    # Add the summary rows to the dataframe
    df.loc['CXX_E'] = total_all_events
    df.loc['CXX_SE'] = sum_specific_events
    df.loc['CXX_DME'] = sum_dme
    df.loc['CXX_DNME'] = sum_dnme
    df.loc['CXX_CE'] = sum_ce
    df.loc['CXX_SNMSs'] = sum_snms
    df.loc['CXX_SMSs'] = sum_sms

    # Calculate the total of the total counts in each category
    total_total_count_in_each_category = df.loc[['CXX_SE', 'CXX_DME', 'CXX_DNME', 'CXX_CE', 'CXX_SNMSs', 'CXX_SMSs']].sum()
    df.loc['CXX_TTCIE'] = total_total_count_in_each_category

    return df


def methypy(input_file):
    """
    Main function to process the input file and generate methylation dataframes.
    
    Parameters:
        input_file (str): Path to the input Excel file.
        
    Returns:
        tuple: Updated dataframes for each methylation type and the original dataframe with sequences.
    """
    # Read the input Excel file
    df = pd.read_excel(input_file)
    half_cols = len(df.columns) // 2

    # Prepare column names for the dataframe
    df.columns = ['MethylType'] + ['Da'] + [f'A_R{i}'for i in range(1, half_cols)] + ['Dk'] + [f'K_R{i}' for i in range(1, half_cols)]

    # Generate methyl sequences based on the prepared dataframe
    df_with_sequences = generate_sequences(df)

    # Prepare dataframes to keep track of each methylation type's sequences
    cxx_df = pd.DataFrame({'Sequence': EVENTS})
    cg_df = pd.DataFrame({'Sequence': EVENTS})
    cxg_df = pd.DataFrame({'Sequence': EVENTS})

    # Initialize the sequence count columns with 0
    for df in [cxx_df, cg_df, cxg_df]:
        for i in range(1, half_cols):
            df[f'R{i}'] = 0

    # Update dataframes for each methylation type with sequence counts
    cg_df_updated = update_methylation_df(cg_df, df_with_sequences, 'CG', half_cols)
    cxg_df_updated = update_methylation_df(cxg_df, df_with_sequences, 'CXG', half_cols)
    cxx_df_updated = update_methylation_df(cxx_df, df_with_sequences, 'CXX', half_cols)

    # Set the 'Sequence' column as the index for each updated dataframe
    cg_df_updated.set_index('Sequence', inplace=True)
    cxg_df_updated.set_index('Sequence', inplace=True)
    cxx_df_updated.set_index('Sequence', inplace=True)

    # Add summary rows to each updated dataframe
    cg_df_updated_with_summary = add_summary_rows(cg_df_updated)
    cxg_df_updated_with_summary = add_summary_rows(cxg_df_updated)
    cxx_df_updated_with_summary = add_summary_rows(cxx_df_updated)

    return cg_df_updated_with_summary, cxg_df_updated_with_summary, cxx_df_updated_with_summary, df_with_sequences


# Execute the main function and get the updated dataframes
cg_df, cxg_df, cxx_df, df_with_sequences = methypy(input_file)

# Save the updated dataframes to CSV files
df_with_sequences.to_csv('sekwencje.csv', index=False)
cg_df.to_csv('cg_df.csv', index=False)
cxg_df.to_csv('cxg_df.csv', index=False)
cxx_df.to_csv('cxx_df.csv', index=False)
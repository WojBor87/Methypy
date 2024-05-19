import pandas as pd
import os

# D:\Python\ZBiB\Renia\PZRP 5K Ex.xlsx

# Define the valid file extensions for input Excel files
VALID_FILE_EXTENSION = ['.xlsx']

EVENTS = [
    '0000', '0001', '0010', '0011',
    '0100', '0101', '0110', '0111',
    '1000', '1001', '1010', '1011',
    '1100', '1101', '1110', '1111'
]


def is_valid_file(file_path):
    """
    Check if the file path has a valid extension.
    
    Parameters:
        file_path (str): The file path to check.
        
    Returns:
        bool: True if the file path has a valid extension, False otherwise.
    """
    _, file_ext = os.path.splitext(file_path)
    return file_ext.lower() in VALID_FILE_EXTENSION


def read_input_file(input_file):
    """
    Read the input Excel file if it exists and has a valid format.
    
    Parameters:
        input_file (str): The path to the input Excel file.
        
    Returns:
        pd.DataFrame or None: The DataFrame read from the Excel file, or None if reading failed.
    """
    if not os.path.exists(input_file):
        print("Error: File not found.")
        return None
    
    if not is_valid_file(input_file):
        print("Error: Invalid file format. Supported formats: .xlsx, .xls")
        return None
    
    try:
        df = pd.read_excel(input_file)
        return df
    except Exception as e:
        print(f"Error: Failed to read the Excel file. {e}")
        return None


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
    if methyl_type != "T":
        methyl_df = df[df['MethylType'] == methyl_type]  
    else:
        methyl_df = df

    for index, row in methylation_df.iterrows():
        sequence = row['Sequence']

        for i in range(1, half_cols):
            # Count occurrences of the sequence in the methylation type-specific dataframe
            count = methyl_df[methyl_df[f'R{i}'] == sequence].shape[0]
            methylation_df.at[index, f'R{i}'] = count

    return methylation_df


def add_summary_rows(df, methyl_type):
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
    df.loc[f'{methyl_type}_E'] = total_all_events
    df.loc[f'{methyl_type}_SE'] = sum_specific_events
    df.loc[f'{methyl_type}_DME'] = sum_dme
    df.loc[f'{methyl_type}_DNME'] = sum_dnme
    df.loc[f'{methyl_type}_CE'] = sum_ce
    df.loc[f'{methyl_type}_SNMSs'] = sum_snms
    df.loc[f'{methyl_type}_SMSs'] = sum_sms

    # Calculate the total of the total counts in each category
    total_total_count_in_each_category = df.loc[[
        f'{methyl_type}_SE', 
        f'{methyl_type}_DME', 
        f'{methyl_type}_DNME', 
        f'{methyl_type}_CE', 
        f'{methyl_type}_SNMSs', 
        f'{methyl_type}_SMSs'
        ]].sum()
    df.loc[f'{methyl_type}_TTCIE'] = total_total_count_in_each_category

    # Calculate D1
    d1 = df.loc[[f'{methyl_type}_TTCIE', '0000']].sum()
    df.loc[f'{methyl_type}_D1'] = d1 

    return df


def calculations(summary_df, methyl_type, dataframe_d1):
    summary_df.loc['CXX_D1'] = dataframe_d1['CXX']
    summary_df.loc['CG_D1'] = dataframe_d1['CG']
    summary_df.loc['CXG_D1'] = dataframe_d1['CXG']
    summary_df.loc['D1'] = dataframe_d1['T']
    summary_df.loc[f'{methyl_type}_DMV'] = 100 * summary_df.loc[f'{methyl_type}_DME'] / summary_df.loc['D1']
    summary_df.loc[f'{methyl_type}_DNMV'] = 100 * summary_df.loc[f'{methyl_type}_DNME'] / summary_df.loc['D1']
    summary_df.loc[f'{methyl_type}_SV'] = 100 * summary_df.loc[f'{methyl_type}_SE'] / summary_df.loc['D1']
    summary_df.loc[f'{methyl_type}_CV'] = 100 * summary_df.loc[f'{methyl_type}_CE'] / summary_df.loc['D1']
    summary_df.loc[f'{methyl_type}_TTCIV'] = 100 * summary_df.loc[f'{methyl_type}_TTCIE'] / summary_df.loc['D1']

    # SV+DMV+DNMV
    SUM = summary_df.loc[f'{methyl_type}_SV'] + summary_df.loc[f'{methyl_type}_DMV'] + summary_df.loc[f'{methyl_type}_DNMV']
    # CV / SUM
    CV_multiplier = summary_df.loc[f'{methyl_type}_CV'] / SUM

    summary_df.loc[f'{methyl_type}_SV_CN'] = summary_df.loc[f'{methyl_type}_SV'] + (summary_df.loc[f'{methyl_type}_SV'] * CV_multiplier) / summary_df.loc[f'{methyl_type}_E']
    summary_df.loc[f'{methyl_type}_DMV_CN'] = summary_df.loc[f'{methyl_type}_DMV'] + (summary_df.loc[f'{methyl_type}_DMV'] * CV_multiplier) / summary_df.loc[f'{methyl_type}_E']
    summary_df.loc[f'{methyl_type}_DNMV_CN'] = summary_df.loc[f'{methyl_type}_DNMV'] + (summary_df.loc[f'{methyl_type}_DNMV'] * CV_multiplier) / summary_df.loc[f'{methyl_type}_E']
    summary_df.loc[f'{methyl_type}_dMET_CN'] = summary_df.loc[f'{methyl_type}_DNMV_CN'] - summary_df.loc[f'{methyl_type}_DMV_CN']

    return summary_df


def methypy(input_file):
    """
    Main function to process the input file and generate methylation dataframes.
    
    Parameters:
        input_file (str): Path to the input Excel file.
        
    Returns:
        tuple: Updated dataframes for each methylation type and the original dataframe with sequences.
    """
    # Read the input Excel file
    df = read_input_file(input_file)
    if df is None:
        return None, None, None, None

    half_cols = len(df.columns) // 2

    # Prepare column names for the dataframe
    df.columns = ['MethylType'] + ['Da'] + [f'A_R{i}'for i in range(1, half_cols)] + ['Dk'] + [f'K_R{i}' for i in range(1, half_cols)]

    # Generate methyl sequences based on the prepared dataframe
    df_with_sequences = generate_sequences(df)

    # Prepare dataframes to keep track of each methylation type's sequences
    cxx_df = pd.DataFrame({'Sequence': EVENTS})
    cg_df = pd.DataFrame({'Sequence': EVENTS})
    cxg_df = pd.DataFrame({'Sequence': EVENTS})
    total_df = pd.DataFrame({'Sequence': EVENTS})

    # Initialize the sequence count columns with 0
    for df in [cxx_df, cg_df, cxg_df, total_df]:
        for i in range(1, half_cols):
            df[f'R{i}'] = 0

    # Define a dictonary to store DataFrames
    dataframes = {
        "CXX": cxx_df,
        "CG": cg_df,
        "CXG": cxg_df,
        "T": total_df
    }

    dataframe_d1 = {
    }

    for methyl_type, df in dataframes.items():
        # Update dataframe with sequence counts
        updated_df = update_methylation_df(df, df_with_sequences, methyl_type, half_cols)

        # Set the 'Sequence' column as the index
        updated_df.set_index('Sequence', inplace=True)

        # Add summary rows
        updated_df_with_summary = add_summary_rows(updated_df, methyl_type)

        # Update the dataframe in the dictionary
        dataframes[methyl_type] = updated_df_with_summary
        dataframe_d1[methyl_type] = updated_df_with_summary.loc[f'{methyl_type}_D1']
    
    for methyl_type, df in dataframes.items():
        # Update dataframe with calculations
        updated_df_with_calculations = calculations(df, methyl_type, dataframe_d1)

        # Update the dataframe in the dictionary
        dataframes[methyl_type] = updated_df_with_calculations

    return dataframes['CG'], dataframes['CXG'], dataframes['CXX'], dataframes['T']


def main():
    # Prompt user to select the input file
    input_file = input("Please enter the path to your Excel file: ")
    
    # Process the input file
    cg_df, cxg_df, cxx_df, total_df = methypy(input_file)
    
    if cg_df is not None and cxg_df is not None and cxx_df is not None and total_df is not None:
        # Save the updated dataframes to CSV files
        cg_df.to_csv('results\cg_df.csv', index=True)
        cxg_df.to_csv('results\cxg_df.csv', index=True)
        cxx_df.to_csv('results\cxx_df.csv', index=True)
        total_df.to_csv('results\\total_df.csv', index=True)
        print("Processing completed successfully.")
    else:
        print("Processing failed. Please check the input file path and format.")


if __name__ == "__main__":
    main()
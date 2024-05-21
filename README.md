# Methypy

## Description
Methypy is a Python script designed for analyzing methylation sequences from an Excel dataset. It processes the input data to generate specific methylation sequences and updates dataframes for different methylation types (CG, CXG, CXX) with counts of these sequences. The results are saved as CSV files for further analysis.

## Features
- **Data Input:** Reads methylation data from an Excel file.
- **Sequence Generation:** Combines values from specific columns to generate methylation sequences.
- **Sequence Counting:** Updates dataframes with counts of sequences for different methylation types.
- **Output:** Saves the updated dataframes to CSV files.

## Files
- **methypy.py:** Main script for processing the methylation data.
- **results/cg_df.csv:** CSV file with counts of CG methylation sequences.
- **results/cxg_df.csv:** CSV file with counts of CXG methylation sequences.
- **results/cxx_df.csv:** CSV file with counts of CXX methylation sequences.
- **results/total_df.csv:** CSV file with counts of all methylation sequences.

## Usage
1. **Install Requirements: Install the required Python packages using the requirements.txt file.
   ```sh
   pip install -r requirements.txt
2. **Run the Script: Execute the script with the provided Excel file path.
   ```sh
   python methypy.py
3. **Provide File Path: When prompted, enter the path to your Excel file.
4. Analyze Results: Review the generated CSV files in the results directory for methylation sequence analysis.

## Author
- GitHub: [WojBor87](https://github.com/WojBor87)

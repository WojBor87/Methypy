# Methypy

## Description
Methypy is a Python script designed for analyzing methylation sequences from an Excel dataset. It processes the input data to generate specific methylation sequences and updates dataframes for different methylation types (CG, CXG, CXX) with counts of these sequences. The results are saved as CSV files for further analysis.

## Features
- **Data Input:** Reads methylation data from an Excel file.
- **Sequence Generation:** Combines values from specific columns to generate methylation sequences.
- **Sequence Counting:** Updates dataframes with counts of sequences for different methylation types.
- **Output:** Saves the updated dataframes to CSV files.

## Files
- **methylation_analysis.py:** Main script for processing the methylation data.
- **sekwencje.csv:** CSV file containing the original data with generated sequences.
- **cg_df.csv:** CSV file with counts of CG methylation sequences.
- **cxg_df.csv:** CSV file with counts of CXG methylation sequences.
- **cxx_df.csv:** CSV file with counts of CXX methylation sequences.

## Usage
1. **Install Dependencies:** Ensure you have pandas and numpy installed in your Python environment.
2. **Run the Script:** Execute the script with the provided Excel file path.
3. **Analyze Results:** Review the generated CSV files for methylation sequence analysis.

## Author
- GitHub: [WojBor87](https://github.com/WojBor87)
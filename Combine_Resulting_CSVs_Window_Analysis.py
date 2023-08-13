import pandas as pd
import os

# Folder path containing CSV files
folder_path = "Z:\Desktop\PbPbC_Analysis\Exposure_Analysis_N2Z_MediaSeperated\Chem_347_M2\R_Results\Pb_Median_Pd_Accum_Results"
def preprocess_dataframe(df, filename):
    df['File_Name'] = filename  # Add a new column with the file name
    df = df[['File_Name'] + [col for col in df.columns if col != 'File_Name']]  # Move the File_Name column to the front

    # Make other changes to the DataFrame here
    # Example: df['column_name'] = df['column_name'].apply(your_function)
    return df


combined_data = []

# Iterate through all CSV files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith('.csv'):
        csv_path = os.path.join(folder_path, filename)
        # Read the CSV file skipping the first two rows
        df = pd.read_csv(csv_path, skiprows=2)
        df['Index_Row'] = df.index  # Add a new column with index row names
        df = df[['Index_Row'] + [col for col in df.columns if col != 'Index_Row']]  # Move the Index_Row column to the front
        df = preprocess_dataframe(df, filename)
        combined_data.append(df)

if combined_data:
    combined_df = pd.concat(combined_data, ignore_index=True)

    # Save the combined DataFrame to a new CSV file
    combined_df.to_csv('Z:\Desktop\PbPbC_Analysis\Exposure_Analysis_N2Z_MediaSeperated\Chem_347_M2\R_Results\Pd_Median_Pd_Accum_CompiledResults.csv', index=False)
    print("Combined CSV files successfully.")
else:
    print("No CSV files were processed.")
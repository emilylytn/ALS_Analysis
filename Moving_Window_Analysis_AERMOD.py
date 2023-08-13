"""
Code created by Emily Luy Tan
Feb 29, 2023

Purpose: all code in one place
    - new way to deal with moving window based on population density readings of the cases
"""
import os
import pandas as pd
from pathlib import Path
import numpy as np
import arcpy


# Adding a column and filling it with the correct year's exposure reading
# HC
# In this line, specify the input .shp filepath (point shp appended with all year chemical exposures) 
# <-----------------USER INPUT
fc = r"C:\Users\f002w4d\Desktop\Oct_24_Analysis\Exposure_readings\Analysis_2_23\AERMOD_pb\Mortality\Mortality_OH_Expanded.shp"
# region ArcPy
# add empty field to be filled correct exposure reading and name it "expsr_val"
arcpy.AddField_management(fc, "expsr_val", "FLOAT")
# HC
# Creating a temporary sub-table that will allow us to match the "year" column with the correct year of exposure data
fields = ['study_ID', 'year', 'expsr_val']
# appending all the chemical concentration fields to this sub-table without having to list them all out manually
# HC
for year in range (2001, 2021):  # <--------USER INPUT
    fields.append('pb_' + str(year))  # <------------- USER INPUT: range of years(?), specify chemical abbreviation as found in input file attribute field names. -------------->

# update expsr_val field with correct exposure value for the corresponding year in each row if year >= 2000
with arcpy.da.UpdateCursor(fc, fields) as cursor:   # cursor reads through shapefile attribute table row by row
    for row in cursor:
        year = row[1]   # isolate column of sub-table with year information
        if year >= 2001:
            exposure_field = 'pb_' + str(year)    # <------------- USER INPUT: specify chemical abbreviation as found in input file attribute field names. Same as above. ------------------->
            exposure_field_index = fields.index(exposure_field)
            exposure = row[exposure_field_index]    # index into column with matching year in column name
            row[2] = row[exposure_field_index]
            cursor.updateRow(row)       # Add data to exmpty "expsr_val" column
        elif year < 2001:
            row[2] = -9999    # If the year of data is beyond our scope, fill the exposure column in with N/A
            cursor.updateRow(row)
del cursor
#endregion

#region Clean chem data (get one chem reading for each row)
# create a Pandas dataframe (df) that isolates each case's diagnosis year
df_diagnosis = pd.DataFrame(arcpy.da.FeatureClassToNumPyArray(fc, ['study_ID', 'index_year', 'expsr_val']))
df_diag = df_diagnosis.replace(-9999, 0)   # Replacing years that have exposure values of -9999 with 0 so they don't interfere with calculations. \n This was deemed appropriate since the majority of area within the AERMOD Chemical layers is not filled within chemical data
# arbitraily use median to make sure there is only one record of the index year in this dataframe
diagnosis_grouped = df_diag.groupby('study_ID')['index_year'].median()
diagnosis_grouped = diagnosis_grouped.reset_index()

# create df of exposure values for each year per patient
df = pd.DataFrame(arcpy.da.FeatureClassToNumPyArray(fc, ['study_ID', 'year', 'expsr_val'], skip_nulls=True))
df_yrs_out_of_scope_removed =  df.drop(df[df['year'] <= 2000].index) # remove study_ids that read -9999 because there were not those raster years. 
df_yrs_out_of_scope_removed = df_yrs_out_of_scope_removed.reset_index(drop=True)
# df of patient exposure per year set as the median of all location exposures for that year (getting rid of redundant years for each patient)
df_2 = df_yrs_out_of_scope_removed.replace(-9999, 0)
grouped = df_2.groupby(['study_ID', 'year'])['expsr_val'].mean() # df with multilevel indecies... lvl0=patientID, lvl1=year, data=exposure (for year)
# if there are mulitple locations for a single year within a patient's record (i.e. the patient moved once or more wihtin the year), their chemical exposure is averaged between these locations within the year
patient_entire_df = grouped.reset_index()
patient_entire_df 
#endregion

#region Clean popdens data (get one popdens reading for each row)
# repeat this same process with the Population Density exposure readings
# HC
popd = r"C:\Users\f002w4d\Desktop\Oct_24_Analysis\Exposure_readings\Analysis_2_23\Averaged_Xun_PopDens\Mortality\Mortality_OH_ExpandedAll.shp"


# create a df of study_IDs, Year, and the population density exposure
df_pd = pd.DataFrame(arcpy.da.FeatureClassToNumPyArray(popd, ['study_ID', 'year', 'pd_7x7'], skip_nulls=True))
df_yrs_out_of_scope_removed_pd =  df_pd.drop(df_pd[df_pd['year'] <= 2000].index) # remove rows that read -9999 because there were not those raster years. 
df_yrs_out_of_scope_removed_pd = df_yrs_out_of_scope_removed_pd.reset_index(drop=True)
# df of patient exposure per year set as the mean of all location exposures for that year (getting rid of redundant years for each patient)
df_2_pd = df_yrs_out_of_scope_removed_pd.replace(-9999, 0)
grouped_pd = df_2_pd.groupby(['study_ID', 'year'])['pd_7x7'].mean() # df with multilevel indecies... lvl0=patientID, lvl1=year, data=exposure (for year)
patient_entire_df_pd = grouped_pd.reset_index()
patient_entire_df_pd 


# merge the chemical exposure df and the population density dfs using both the study_ID and year as keys
df_merged = pd.merge(patient_entire_df, patient_entire_df_pd, how='outer', left_on=['study_ID','year'], right_on = ['study_ID','year'])
print(df_merged)
#endregion

#region calculate the median chemical/accumulated popdens exposure over various time frames
# functions to loop through the merged df and calculate the median chemical exposure over various time frames
def find_yrs_prior_median(num_years, patient_all_expsr_df, patient_ID, diagnosis_yr):
    earliest_yr = diagnosis_yr - num_years
    patient_expsr_subframe= patient_all_expsr_df[patient_all_expsr_df['study_ID']==patient_ID] # create a subframe of a single study_ID's full data
    filteredbyyear= patient_expsr_subframe[patient_expsr_subframe['year'] >= earliest_yr] # filter this study_ID's subframe so that is only includes the years within a time frame
    med_expsr=filteredbyyear['expsr_val'].median() # take the median of these year's exposure values
    return med_expsr # output of fumtion is this median

# functions to loop through the merged df and calculate the accumulated population density exposure over various time frames
def find_yrs_prior_accum_pd(num_years, patient_all_expsr_df, patient_ID, diagnosis_yr):
    earliest_yr = diagnosis_yr - num_years
    patient_expsr_subframe= patient_all_expsr_df[patient_all_expsr_df['study_ID']==patient_ID]
    filteredbyyear= patient_expsr_subframe[patient_expsr_subframe['year'] >= earliest_yr]
    pd_expsr=filteredbyyear['pd_7x7'].sum()
    return pd_expsr

# set the number of years you want in your time frames
# Run the functions to create a df where each patient has a study ID, index year, and calculated exposures for 0-i years
for i in range(13):
    new_col_name = f'med_expsr_{i}_yrs'
    diagnosis_grouped[new_col_name] = diagnosis_grouped.apply(lambda x: find_yrs_prior_median(i, df_merged, x['study_ID'], x['index_year']), axis=1)   

for k in range(13):
    new_col_name_pd = f'accum_pd_expsr_{k}_yrs'
    diagnosis_grouped[new_col_name_pd] = diagnosis_grouped.apply(lambda x: find_yrs_prior_accum_pd(k, df_merged, x['study_ID'], x['index_year']), axis=1)   

patient_median_exprs_df = diagnosis_grouped
#endregion

#export the final df as a .csv file
patient_median_exprs_df.to_csv (r"C:\Users\f002w4d\Desktop\Oct_24_Analysis\Exposure_readings\Analysis_2_23\Summary_exposures\Mortality\Pb_Median_Pop_Accum_15.csv", index = False, header = True)
print(patient_median_exprs_df)


#___Do it for the controls too__________________________________________________________


controls = "/Users/Shashwat/Desktop/Pb_Median_Pop_Accum_15 2.csv"
df_control = pd.read_csv(controls)

cases = "/Users/Shashwat/Desktop/Pb_Median_Pop_Accum_15.csv"
patient_median_exprs_df = pd.read_csv(cases)


#region Function: sliding window
def sliding_window(iterable, size, overlap=0):
    """
        >>> list(sliding_window([1, 2, 3, 4], size=2))
        [(1, 2), (3, 4)]
        >>> list(sliding_window([1, 2, 3], size=2, overlap=1))
        [(1, 2), (2, 3)]
        >>> list(sliding_window([1, 2, 3, 4, 5], size=3, overlap=1))
        [(1, 2, 3), (3, 4, 5)]
        >>> list(sliding_window([1, 2, 3, 4], size=3, overlap=1))
        [(1, 2, 3), (3, 4)]
        >>> list(sliding_window([1, 2, 3, 4], size=10, overlap=8))
        [(1, 2, 3, 4)]
    """
    start = 0
    end = size
    step = size - overlap
    if step <= 0:
        ValueError("overlap must be smaller then size")
    length = len(iterable)
    windows95=[]
    while end < length:
        output = iterable.iloc[start:end]
        windows95.append(output)
        start += step
        end += step
    return windows95
#endregion

#region function: rejoin attribute data
"""
-------------------------------------------------------------
Function: Rejoin CDCP and Case data with age and sex data
-------------------------------------------------------------
"""

def add_attribute_data(attribute_data_file, join_to_data_df):

    # Read csv files as pandas df
    df_attributes=pd.read_csv(attribute_data_file)
    
    joined_df=pd.merge(join_to_data_df, df_attributes, on="study_ID", how="inner")
    joined_df.insert(0, 'ID', joined_df.index + 1)
    joined_df=joined_df.drop(columns=['study_ID'])
    final_window_df=joined_df[['ID','disease', new_chem_col_name, new_pop_col_name, 'AGE', 'SEX']]
    return final_window_df
#endregion

#region Loop: For each year, for each case pop dens window, match appropriate control pop dens
for i in range(13):
    new_df_name = f'Med_pb_Accum_pd_{i}_yrs'
    new_chem_col_name = f'med_expsr_{i}_yrs'
    new_pop_col_name = f'accum_pd_expsr_{i}_yrs'
    new_df = patient_median_exprs_df.filter(['study_ID', new_chem_col_name, new_pop_col_name], axis=1)
    sorted_df=new_df.sort_values(by=[f'accum_pd_expsr_{i}_yrs'],  ascending=True)
    sorted_df['rank'] =sorted_df[f'accum_pd_expsr_{i}_yrs'].rank(axis=0, method='first') # rank population density from smallest to biggest
    df_window=sliding_window(sorted_df,50,overlap=40) 
    
    window_count = 50
    for d in df_window:
        low_pop=d.iat[0,2] # save population density associated with the lowest rank
        high_pop=d.iat[49,2] #save population density associated with the highest rank (together, these are the bound for the range od population densities from the control file)
        filtered_control_df=df_control.loc[df_control[f'accum_pd_expsr_{i}_yrs']>=low_pop].loc[df_control[f'accum_pd_expsr_{i}_yrs']<=high_pop] #make a df of the control data that falls within the same population density range as that of the current window
        # filter out correct year of controls here (study_id, chem year, pop year.... then concat)
        matched_controls_df=filtered_control_df[['study_ID', new_chem_col_name, new_pop_col_name]]
        d=d.drop(columns=['rank'])
        matched_df=pd.concat([d,matched_controls_df]) # combine the control df and case df one on top of another
        final_wind_df=add_attribute_data(r"/Users/Shashwat/Desktop/Emily GIS work/total_attribute_data.csv", matched_df)
        output_file_name = f'{new_df_name}' + f'_Window{window_count}' 
        final_wind_df.to_csv("/Users/Shashwat/Desktop/Emily GIS work/R_Inputs/Raw_Value_Outputs/Mean_Pb_Accum_Pop_12" + "/" + output_file_name + ".csv", index = False, header = True)
        
        # log transformation...
        df_orig=final_wind_df
        # isolate colums you want to log transform and log to transform
        df_log = df_orig.iloc[:,[2]].applymap(lambda x: (x * 1000000))
        df_log.columns = 'mil_' + df_log.columns # rename these columns to show they have been transformed
        # replace non-transformed columns with the transformed ones
        df_orig.drop(df_orig.iloc[:,[2]], axis=1, inplace=True)
        df_orig.insert(2, df_log.columns[0], df_log.iloc[:,[0]])
        df_log = df_orig.iloc[:,[2]].applymap(lambda x: np.log10(x+1))
        df_log.columns = 'log_' + df_log.columns # rename these columns to show they have been transformed
        # replace non-transformed columns with the transformed ones
        df_orig.drop(df_orig.iloc[:,[2]], axis=1, inplace=True)
        df_orig.insert(2, df_log.columns[0], df_log.iloc[:,[0]])
        df_orig.to_csv("/Users/Shashwat/Desktop/Emily GIS work/R_Inputs/Raw_Value_Outputs/Mean_Pb_Accum_Pop_12_LogMill" + "/" + output_file_name + "_LogM.csv", index = False, header = True)
        
        window_count += 10

        #within this loop, deal with each individual window and export as individual file.

#endregion
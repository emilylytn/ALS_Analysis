import os
import arcpy
from arcpy import env
from arcpy.sa import *
import pandas as pd
from pathlib import Path
import numpy as np


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = "toolbox"

        # List of tool classes associated with this toolbox
        self.tools = [Tool1, Tool2, Tool3, Tool5, Tool6, Tool7]


class Tool1(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Exposures"
        self.description = "The Extract Exposures tool reads the exposure for shapefile points across multiple rasters. Exposure values will be added to the original input .shp attribute table." 
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Environmental Rasters Location (.gdb)",
                            name="environ",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="Exposure Points Folder",
                            name="caseFolder",
                            datatype="DEFolder",
                            parameterType="Required",
                            direction="Input"),
                                        
            arcpy.Parameter(displayName="Raster Name string segment that contains the year",
                            name="seg_num",
                            datatype="GPString",
                            parameterType="Required",
                            direction="Input"),   

            arcpy.Parameter(displayName="Raster Name string segment that PRECEEDS year",
                            name="pre_year",
                            datatype="GPString",
                            parameterType="Required",
                            direction="Input"),    

            arcpy.Parameter(displayName="Raster Name string segment that FOLLOWS year",
                            name="post_year",
                            datatype="GPString",
                            parameterType="Required",
                            direction="Input"),                   
        ]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        environ = parameters[0].valueAsText
        caseFolder = parameters[1].valueAsText
        seg_num = parameters[2].valueAsText
        pre_year = parameters[3].valueAsText
        post_year = parameters[4].valueAsText

        # Set the current workspace (your gdb of environmental rasters)
        arcpy.env.workspace = environ  

        # Set the case folder including all the case files (.shp)
        # This file path name should be a FOLDER that contains .shp files.
        print(caseFolder)

        # Read all raster layers in current workplace
        rasters = arcpy.ListRasters("*", "all")

        inRasterList = []

        for raster in rasters:
            year = raster.split("_")[int(seg_num)]  # Assuming the file name format is "name_year.shp"
            rasterName = pre_year + year + post_year  # Change the string to match the names of your raster files' names
            fieldName = "pb_" + year  # This is the name of the columns of exposures that will be created
            inRasterList.append([rasterName, fieldName])
            arcpy.AddMessage("Current raster name: " + str(rasterName))

        for file in os.listdir(caseFolder):
            if file.endswith(".shp"):
                CasesSet = caseFolder + "/" + file
                print(CasesSet)
                
                ExtractMultiValuesToPoints(CasesSet, inRasterList, "NONE")

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return


class Tool2(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Clean Exposure Data"
        self.description = "This tool adds an exposure value column for each case/control, extracting the correct year of exposure data for each row in data table. The data also gets cleaned, eliminating duplicates and any years of data that fall outside of a temporal limit. WARNING: This code permenantly changes input .shp"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Point .shp with all raw exposure readings",
                            name="fc",
                            datatype="DEShapeFile",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="First year of Environmental Data",
                            name="startyr",
                            datatype="GPString",
                            parameterType="Required",
                            direction="Input"),
                                        
            arcpy.Parameter(displayName="Last year of Environmental Data",
                            name="endyr",
                            datatype="GPString",
                            parameterType="Required",
                            direction="Input")
        ]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        fc = parameters[0].valueAsText
        startyr = parameters[1].valueAsText
        endyr = parameters[2].valueAsText

        # In this line, specify the input .shp filepath (point shp appended with all year chemical exposures)
        # <-----------------USER INPUT
        # ArcPy
        # add empty field to be filled correct exposure reading and name it "expsr_val"
        arcpy.AddField_management(fc, "expsr_val", "FLOAT")
        # HC
        # Creating a temporary sub-table that will allow us to match the "year" column with the correct year of exposure data
        fields = ['study_ID', 'year', 'expsr_val']
        # appending all the chemical concentration fields to this sub-table without having to list them all out manually
        # HC

        for year in range(int(startyr), (int(endyr) + 1)):  # <--------USER INPUT
            fields.append('pb_' + str(year))  # <------------- USER INPUT: range of years(?), specify chemical abbreviation as found in input file attribute field names. -------------->

        # update expsr_val field with correct exposure value for the corresponding year in each row if year >= 2000
        with arcpy.da.UpdateCursor(fc, fields) as cursor:   # cursor reads through shapefile attribute table row by row
            for row in cursor:
                year = row[1]   # isolate column of sub-table with year information
                if year >= int(startyr):
                    exposure_field = 'pb_' + str(year)    # <------------- USER INPUT: specify chemical abbreviation as found in input file attribute field names. Same as above. ------------------->
                    exposure_field_index = fields.index(exposure_field)
                    exposure = row[exposure_field_index]    # index into column with matching year in column name
                    row[2] = row[exposure_field_index]
                    cursor.updateRow(row)       # Add data to exmpty "expsr_val" column
                elif year < int(startyr):
                    row[2] = -9999    # If the year of data is beyond our scope, fill the exposure column in with N/A
                    cursor.updateRow(row)
        del cursor

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class Tool3(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Summary Exposure Calculation"
        self.description = "This tool adds population density 'exposure' data for points. It also calculates Median, Sum,  Max, and Min summary statistics across a set number of timeframes. "
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Point .shp with all exposure readings",
                            name="fc",
                            datatype="DEShapeFile",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="First year of Environmental Data",
                            name="startyr",
                            datatype="GPString",
                            parameterType="Required",
                            direction="Input"),
                                        
            arcpy.Parameter(displayName="Last year of Environmental Data",
                            name="endyr",
                            datatype="GPString",
                            parameterType="Required",
                            direction="Input"),

            arcpy.Parameter(displayName="Population Density exposures .shp",
                            name="popd",
                            datatype="DEShapeFile",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="Number of timeframe years",
                            name="timeframes_num",
                            datatype="GPLong",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="Output Folder",
                            name="output_folder",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),

            arcpy.Parameter(displayName="Output path (for summary statistic .csv)",
                            name="outpath_base",
                            datatype="GPString",
                            parameterType="Required",
                            direction="Output"),
            
            arcpy.Parameter(displayName="Output path (for TEST.csv)",
                            name="outpath_test",
                            datatype="DEFile",
                            parameterType="Required",
                            direction="Output"),

                            
        ]
        return params


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        fc=parameters[0].valueAsText
        startyr=parameters[1].valueAsText
        endyr=parameters[2].valueAsText
        popd=parameters[3].valueAsText
        timeframes_num=parameters[4].valueAsText
        output_folder=parameters[5].valueAsText
        outpath_base=parameters[6].valueAsText
        outpath_test=parameters[7].valueAsText

        summary_types = ["median", "sum", "max", "min"]

        #region Clean chem data (get one chem reading for each row)
        # create a Pandas dataframe (df) that isolates each case's diagnosis year
        df_diagnosis = pd.DataFrame(arcpy.da.FeatureClassToNumPyArray(fc, ['study_ID', 'index_year', 'expsr_val']))
        df_diag = df_diagnosis.replace(-9999, 0)   # Replacing years that have exposure values of -9999 with 0 so they don't interfere with calculations. \n This was deemed appropriate since the majority of area within the AERMOD Chemical layers is not filled within chemical data
        # arbitraily use median to make sure there is only one record of the index year in this dataframe
        diagnosis_grouped = df_diag.groupby('study_ID')['index_year'].median()
        diagnosis_grouped = diagnosis_grouped.reset_index()

        # create df of exposure values for each year per patient
        df = pd.DataFrame(arcpy.da.FeatureClassToNumPyArray(fc, ['study_ID', 'year', 'expsr_val'], skip_nulls=True))
        df_yrs_out_of_scope_removed =  df.drop(df[df['year'] <= (int(startyr)-1)].index) # remove study_ids that read -9999 because there were not those raster years. 
        df_yrs_out_of_scope_removed = df_yrs_out_of_scope_removed.reset_index(drop=True)
        # df of patient exposure per year set as the median of all location exposures for that year (getting rid of redundant years for each patient)
        df_2 = df_yrs_out_of_scope_removed.replace(-9999, 0)
        grouped = df_2.groupby(['study_ID', 'year'])['expsr_val'].mean() # df with multilevel indecies... lvl0=patientID, lvl1=year, data=exposure (for year)
        # if there are mulitple locations for a single year within a patient's record (i.e. the patient moved once or more wihtin the year), their chemical exposure is averaged between these locations within the year
        patient_entire_df = grouped.reset_index()
        patient_entire_df.to_csv (outpath_test, index = False, header = True)
        #endregion

        #region Clean popdens data (get one popdens reading for each row)
        # repeat this same process with the Population Density exposure readings
        # HC

        # create a df of study_IDs, Year, and the population density exposure
        df_pd = pd.DataFrame(arcpy.da.FeatureClassToNumPyArray(popd, ['study_ID', 'year', 'pd_7x7'], skip_nulls=True))
        df_yrs_out_of_scope_removed_pd =  df_pd.drop(df_pd[df_pd['year'] <= (int(startyr)-1)].index) # remove rows that read -9999 because there were not those raster years. 
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
        
        for calc_type in summary_types:
            outpath = arcpy.os.path.join(output_folder, f"{outpath_base}_{calc_type}.csv")

            if calc_type == "median":
                arcpy.AddMessage("Currently summarizing: Median")
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
                for i in range(int(timeframes_num)+1):
                    new_col_name = f'med_expsr_{i}_yrs'
                    diagnosis_grouped[new_col_name] = diagnosis_grouped.apply(lambda x: find_yrs_prior_median(i, df_merged, x['study_ID'], x['index_year']), axis=1)   

                for k in range(int(timeframes_num)+1):
                    new_col_name_pd = f'accum_pd_expsr_{k}_yrs'
                    diagnosis_grouped[new_col_name_pd] = diagnosis_grouped.apply(lambda x: find_yrs_prior_accum_pd(k, df_merged, x['study_ID'], x['index_year']), axis=1)   

                patient_median_exprs_df = diagnosis_grouped
                #endregion

                #export the final df as a .csv file (summary statistic)
                patient_median_exprs_df.to_csv (outpath, index = False, header = True)
                arcpy.AddMessage(f"Median: {outpath} csv created")

            elif calc_type =="sum":
                arcpy.AddMessage("Currently summarizing: Accum")
                #region calculate the median chemical/accumulated popdens exposure over various time frames
                # functions to loop through the merged df and calculate the median chemical exposure over various time frames
                def find_yrs_prior_median(num_years, patient_all_expsr_df, patient_ID, diagnosis_yr):
                    earliest_yr = diagnosis_yr - num_years
                    patient_expsr_subframe= patient_all_expsr_df[patient_all_expsr_df['study_ID']==patient_ID] # create a subframe of a single study_ID's full data
                    filteredbyyear= patient_expsr_subframe[patient_expsr_subframe['year'] >= earliest_yr] # filter this study_ID's subframe so that is only includes the years within a time frame
                    med_expsr=filteredbyyear['expsr_val'].sum() # take the median of these year's exposure values
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
                for i in range(int(timeframes_num)+1):
                    new_col_name = f'accum_expsr_{i}_yrs'
                    diagnosis_grouped[new_col_name] = diagnosis_grouped.apply(lambda x: find_yrs_prior_median(i, df_merged, x['study_ID'], x['index_year']), axis=1)   

                for k in range(int(timeframes_num)+1):
                    new_col_name_pd = f'accum_pd_expsr_{k}_yrs'
                    diagnosis_grouped[new_col_name_pd] = diagnosis_grouped.apply(lambda x: find_yrs_prior_accum_pd(k, df_merged, x['study_ID'], x['index_year']), axis=1)   

                patient_median_exprs_df = diagnosis_grouped
                #endregion

                #export the final df as a .csv file (summary statistic)
                patient_median_exprs_df.to_csv (outpath, index = False, header = True)
                arcpy.AddMessage(f"Accum: {outpath} csv created")
                

            elif calc_type == "max":
                arcpy.AddMessage("Currently summarizing: Max")
                #region calculate the median chemical/accumulated popdens exposure over various time frames
                # functions to loop through the merged df and calculate the median chemical exposure over various time frames
                def find_yrs_prior_median(num_years, patient_all_expsr_df, patient_ID, diagnosis_yr):
                    earliest_yr = diagnosis_yr - num_years
                    patient_expsr_subframe= patient_all_expsr_df[patient_all_expsr_df['study_ID']==patient_ID] # create a subframe of a single study_ID's full data
                    filteredbyyear= patient_expsr_subframe[patient_expsr_subframe['year'] >= earliest_yr] # filter this study_ID's subframe so that is only includes the years within a time frame
                    med_expsr=filteredbyyear['expsr_val'].max() # take the median of these year's exposure values
                    return med_expsr # output of fumtion is this median

                # functions to loop through the merged df and calculate the accumulated population density exposure over various time frames
                def find_yrs_prior_accum_pd(num_years, patient_all_expsr_df, patient_ID, diagnosis_yr):
                    earliest_yr = diagnosis_yr - num_years
                    patient_expsr_subframe= patient_all_expsr_df[patient_all_expsr_df['study_ID']==patient_ID]
                    filteredbyyear= patient_expsr_subframe[patient_expsr_subframe['year'] >= earliest_yr]
                    pd_expsr=filteredbyyear['pd_7x7'].sum()
                    return pd_expsr

                for i in range(int(timeframes_num)+1):
                    new_col_name = f'max_expsr_{i}_yrs'
                    diagnosis_grouped[new_col_name] = diagnosis_grouped.apply(lambda x: find_yrs_prior_median(i, df_merged, x['study_ID'], x['index_year']), axis=1)   

                for k in range(int(timeframes_num)+1):
                    new_col_name_pd = f'accum_pd_expsr_{k}_yrs'
                    diagnosis_grouped[new_col_name_pd] = diagnosis_grouped.apply(lambda x: find_yrs_prior_accum_pd(k, df_merged, x['study_ID'], x['index_year']), axis=1)   

                patient_median_exprs_df = diagnosis_grouped
                #endregion

                #export the final df as a .csv file (summary statistic)
                patient_median_exprs_df.to_csv (outpath, index = False, header = True)
                arcpy.AddMessage(f"Max: {outpath} csv created")

            elif calc_type == "min":
                arcpy.AddMessage("Currently summarizing: Min")
                #region calculate the median chemical/accumulated popdens exposure over various time frames
                # functions to loop through the merged df and calculate the median chemical exposure over various time frames
                def find_yrs_prior_median(num_years, patient_all_expsr_df, patient_ID, diagnosis_yr):
                    earliest_yr = diagnosis_yr - num_years
                    patient_expsr_subframe= patient_all_expsr_df[patient_all_expsr_df['study_ID']==patient_ID] # create a subframe of a single study_ID's full data
                    filteredbyyear= patient_expsr_subframe[patient_expsr_subframe['year'] >= earliest_yr] # filter this study_ID's subframe so that is only includes the years within a time frame
                    med_expsr=filteredbyyear['expsr_val'].min() # take the median of these year's exposure values
                    return med_expsr # output of fumtion is this median

                # functions to loop through the merged df and calculate the accumulated population density exposure over various time frames
                def find_yrs_prior_accum_pd(num_years, patient_all_expsr_df, patient_ID, diagnosis_yr):
                    earliest_yr = diagnosis_yr - num_years
                    patient_expsr_subframe= patient_all_expsr_df[patient_all_expsr_df['study_ID']==patient_ID]
                    filteredbyyear= patient_expsr_subframe[patient_expsr_subframe['year'] >= earliest_yr]
                    pd_expsr=filteredbyyear['pd_7x7'].sum()
                    return pd_expsr

                for i in range(int(timeframes_num)+1):
                    new_col_name = f'min_expsr_{i}_yrs'
                    diagnosis_grouped[new_col_name] = diagnosis_grouped.apply(lambda x: find_yrs_prior_median(i, df_merged, x['study_ID'], x['index_year']), axis=1)   

                for k in range(int(timeframes_num)+1):
                    new_col_name_pd = f'accum_pd_expsr_{k}_yrs'
                    diagnosis_grouped[new_col_name_pd] = diagnosis_grouped.apply(lambda x: find_yrs_prior_accum_pd(k, df_merged, x['study_ID'], x['index_year']), axis=1)   

                patient_median_exprs_df = diagnosis_grouped
                #endregion

                #export the final df as a .csv file (summary statistic)
                patient_median_exprs_df.to_csv (outpath, index = False, header = True)
                arcpy.AddMessage(f"Min: {outpath} csv created")

                # Run the functions to create a df where each patient



class Tool5(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Window Analysis and Export"
        self.description = "This tool takes case and control summary exposure csvs and runs a rolling window analysis on them. Creates both log and non-log transformed export files. Output should be a folder as there will be many csv files."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Folder that contains case and control summary folders",
                            name="aaa_folder",
                            datatype="DEFolder",
                            parameterType="Required",
                            direction="Input"),

            arcpy.Parameter(displayName="Number of timeframe years",
                            name="timeframes_num2",
                            datatype="GPLong",
                            parameterType="Required",
                            direction="Input"),

            arcpy.Parameter(displayName="Attribute Data File (containing controls and cases)",
                            name="attrib_file",
                            datatype="DEFile",
                            parameterType="Required",
                            direction="Input"),

            arcpy.Parameter(displayName="Output folder for Window Analysis files",
                            name="output_folder",
                            datatype="DEFolder",
                            parameterType="Required",
                            direction="Input"),

        ]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed. This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        aaa_folder = parameters[0].valueAsText
        timeframes_num2 = parameters[1].valueAsText
        attrib_file = parameters[2].valueAsText
        output_folder = parameters[3].valueAsText

        import os
        import pandas as pd
        import numpy as np
        import glob

        control_folder = os.path.join(aaa_folder, "Controls")
        mortality_folder = os.path.join(aaa_folder, "Mortality")

        calc_types = ["median", "sum", "max", "min"]

        abbreviation = {
            "median": "med",
            "sum": "accum",
            "max": "max",
            "min": "min"
        }

        for calc_type in calc_types:
            arcpy.AddMessage(f"Processing {calc_type} files...")
            control_files = glob.glob(os.path.join(control_folder, f"*{calc_type}*.csv"))
            mortality_files = glob.glob(os.path.join(mortality_folder, f"*{calc_type}*.csv"))

            calc_abb = abbreviation.get(calc_type, calc_type)

            arcpy.AddMessage(f"Expected control files: {control_files}")
            arcpy.AddMessage(f"Expected mortality files: {mortality_files}")

            existing_control_files = [file_name for file_name in control_files if os.path.exists(
                os.path.join(control_folder, file_name.replace("\\", "/")))]
            existing_mortality_files = [file_name for file_name in mortality_files if os.path.exists(
                os.path.join(mortality_folder, file_name.replace("\\", "/")))]

            if existing_control_files:
                control_file2 = str(existing_control_files[0])
            else:
                arcpy.AddMessage(f"No {calc_type} control files found or there is an error.")

            if existing_mortality_files:
                mortality_file2 = str(existing_mortality_files[0])
            else:
                arcpy.AddMessage(f"No {calc_type} mortality files found or there is an error.")

            if os.path.join(control_folder, control_file2) in control_files and os.path.join(mortality_folder,
                                                                                              mortality_file2) in mortality_files:
                arcpy.AddMessage(f"There's a match! {control_file2} and {mortality_file2}")
                control_path = os.path.join(control_folder, control_file2)
                arcpy.AddMessage(f"Importing {control_path} into dataframe...")
                mortality_path = os.path.join(mortality_folder, mortality_file2)
                arcpy.AddMessage(f"Importing {control_path} into dataframe...")

                df_control = pd.read_csv(control_path)
                patient_median_exprs_df = pd.read_csv(mortality_path)

                # The rest of your existing code for window analysis goes here

                # region Function: sliding window
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
                    windows95 = []
                    while end < length:
                        output = iterable.iloc[start:end]
                        windows95.append(output)
                        start += step
                        end += step
                    return windows95
                # endregion

                # region function: rejoin attribute data
                """
                -------------------------------------------------------------
                Function: Rejoin CDCP and Case data with age and sex data
                -------------------------------------------------------------
                """

                def add_attribute_data(attribute_data_file, join_to_data_df):

                    # Read csv files as pandas df
                    df_attributes = pd.read_csv(attribute_data_file)

                    joined_df = pd.merge(join_to_data_df, df_attributes, on="study_ID", how="inner")
                    joined_df.insert(0, 'ID', joined_df.index + 1)
                    joined_df = joined_df.drop(columns=['study_ID'])
                    final_window_df = joined_df[['ID', 'disease', new_chem_col_name, new_pop_col_name, 'AGE', 'SEX']]
                    return final_window_df
                # endregion

                # Loop: For each year, for each case pop dens window, match appropriate control pop dens
                for i in range(int(timeframes_num2) + 1):
                    new_df_name = f'{calc_abb}_pb_Accum_pd_{i}_yrs'
                    new_chem_col_name = f'{calc_abb}_expsr_{i}_yrs'
                    new_pop_col_name = f'accum_pd_expsr_{i}_yrs'
                    new_df = patient_median_exprs_df.filter(['study_ID', new_chem_col_name, new_pop_col_name],
                                                             axis=1)
                    sorted_df = new_df.sort_values(by=[f'accum_pd_expsr_{i}_yrs'], ascending=True)
                    sorted_df['rank'] = sorted_df[f'accum_pd_expsr_{i}_yrs'].rank(axis=0, method='first')
                    df_window = sliding_window(sorted_df, 50, overlap=40)

                    window_count = 50

                    for d in df_window:
                        low_pop = d.iat[0, 2]
                        high_pop = d.iat[49, 2]
                        filtered_control_df = df_control.loc[df_control[f'accum_pd_expsr_{i}_yrs'] >= low_pop].loc[
                            df_control[f'accum_pd_expsr_{i}_yrs'] <= high_pop]
                        matched_controls_df = filtered_control_df[['study_ID', new_chem_col_name, new_pop_col_name]]
                        d = d.drop(columns=['rank'])
                        matched_df = pd.concat([d, matched_controls_df])
                        final_wind_df = add_attribute_data(str(attrib_file), matched_df)

                        arcpy.AddMessage(f"Exporting {calc_type} results...")

                        # Create a folder for each calc_type within the specified parent output folder
                        calc_type_folder = os.path.join(output_folder, f"{calc_type}_Output")
                        os.makedirs(calc_type_folder, exist_ok=True)
                        output_file_name = f'{new_df_name}' + f'_Window{window_count}'
                        output_file_path = os.path.join(calc_type_folder, output_file_name + ".csv")
                        final_wind_df.to_csv(output_file_path, index=False, header=True)

                        arcpy.AddMessage(f"Exporting log results for {calc_type}...")

                        # log transformation...
                        df_orig = final_wind_df
                        df_log = df_orig.iloc[:, [2]].applymap(lambda x: (x * 1000000))
                        df_log.columns = 'mil_' + df_log.columns
                        df_orig.drop(df_orig.iloc[:, [2]], axis=1, inplace=True)
                        df_orig.insert(2, df_log.columns[0], df_log.iloc[:, [0]])
                        df_log = df_orig.iloc[:, [2]].applymap(lambda x: np.log10(x + 1))
                        df_log.columns = 'log_' + df_log.columns
                        df_orig.drop(df_orig.iloc[:, [2]], axis=1, inplace=True)
                        df_orig.insert(2, df_log.columns[0], df_log.iloc[:, [0]])

                        log_folder_path = os.path.join(output_folder, "Log")
                        os.makedirs(log_folder_path, exist_ok=True)
                        calc_type_log_folder = os.path.join(log_folder_path, f"{calc_type}_Log_Output")
                        os.makedirs(calc_type_log_folder, exist_ok=True)

                        log_output_file_path = os.path.join(calc_type_log_folder, output_file_name + "_LogM.csv")
                        df_orig.to_csv(log_output_file_path, index=False, header=True)

                        window_count += 10
        return


class Tool6(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "MUST FIX R Logistic Regression"
        self.description = "Log Regression"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Folder containing all window analysis files",
                            name="input_WA_folder",
                            datatype="DEFolder",
                            parameterType="Required",
                            direction="Input"), 
            
            arcpy.Parameter(displayName="Output folder for all logistic regression files",
                            name="output_folder",
                            datatype="DEFolder",
                            parameterType="Required",
                            direction="Input"), 
        ]
        return params


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        input_WA_folder = parameters[0].valueAsText
        output_folder = parameters[1].valueAsText


        import os
        import pandas as pd
        import statsmodels.api as sm
        import csv

        # Folder name variations
        folder_variations = ["Accum", "Median", "Max", "Min"]

        # Iterate over each folder variation
        for variation in folder_variations:
            # Set up the working directory
            input_folder = os.path.join(input_WA_folder, variation)
            output_folder = os.path.join(output_folder, variation)

            os.chdir(input_folder)
        
            # Specify the output folder for results with variation in the folder path
            output_folder = output_folder
            
            # Get a list of all CSV files in the specified path
            file_list = [file for file in os.listdir() if file.endswith("LogM.csv")]
            
            for file_name in file_list:
                # Read CSV file
                temp_data = pd.read_csv(file_name)
                temp_data_new = temp_data.iloc[:, 1:]  # Exclude the first column
                temp_data_new.columns = ["disease"] + list(temp_data_new.columns[1:])
                
                # Perform logistic regression
                temp_regression = sm.Logit(temp_data_new["disease"], temp_data_new.iloc[:, 1:]).fit()
                
                # Construct the output file path with "_Results.csv" at the end
                output_file_path = os.path.join(output_folder, os.path.splitext(file_name)[0] + "_results.csv")
                
                # Write basic results to the file
                df_basic = pd.DataFrame({
                    "variable": temp_data.columns[1],
                    "num_sample": len(temp_data),
                    "call": str(temp_regression.model),
                    "aic": temp_regression.aic,
                    "iteration": temp_regression.nit
                }, index=[0])

                df_basic.to_csv(output_file_path, sep=",", index=False, header=True, mode='a', line_terminator='\n')
                
                # Write coefficient results to the same file
                temp_regression_coefficient = temp_regression.summary2().tables[1]
                temp_regression_coefficient["oddsRatio"] = temp_regression_coefficient["Coef."].apply(lambda x: round(pow(2, x), 3))
                temp_regression_coefficient.to_csv(output_file_path, sep=",", index=True, header=True, mode='a', line_terminator='\n')


                return

        def postExecute(self, parameters):
            """This method takes place after outputs are processed and
            added to the display."""
            return
    

class Tool7(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Subject Location Preparation"
        self.description = "This tool adds a row of data for each location of each subject for each year"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="CSV file containing untouched, unduplicated subject locations",
                            name="input",
                            datatype="DEFile",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="Output CSV file name containing duplicated subject locations",
                            name="output",
                            datatype="DEFile",
                            parameterType="Required",
                            direction="Input")
        ]
        return params


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        input_file = parameters[0].valueAsText
        output_file = parameters[1].valueAsText
        
        import pandas as pd
        import os
        import csv

        # Read CSV file
        kidca_hmaprdngs_big = pd.read_csv(input_file)

        # Data wrangling
        kidca_hmaprdngs_xlong_pat = kidca_hmaprdngs_big.apply(lambda row: pd.DataFrame({
            'ID': row['ID'],
            'study_ID': row['study_ID'],
            'ALS_status': row['ALS_status'],
            'SEX': row['SEX'],
            'AGE': row['AGE'],
            'source_type': row['source_type'],
            'source': row['source'],
            'seq': row['seq'],
            'latitude': row['latitude'],
            'longitude': row['longitude'],
            'yr_addr_end': row['yr_addr_end'],
            'yr_addr_start': row['yr_addr_start'],
            'index_year': row['index_year'],
            'num_yrs_ad': row['num_yrs_ad'],
            'index_yr_minus15': row['index_minus15'],
            'year': list(range(row['yr_addr_start'], row['yr_addr_end'] + 1))
        }), axis=1)

        # Concatenate the list of DataFrames into a single DataFrame
        kidca_hmaprdngs_xlong_pat = pd.concat(kidca_hmaprdngs_xlong_pat.values, ignore_index=True)

        # Write to CSV
        kidca_hmaprdngs_xlong_pat.to_csv(output_file, index=False)
    
        return

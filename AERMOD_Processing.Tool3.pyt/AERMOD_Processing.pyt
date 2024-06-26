# -*- coding: utf-8 -*-

import arcpy
import os
import pandas as pd
import numpy as np
import glob
from glob import iglob


class Toolbox(object):
    def __init__(self):
        """This toolbox creates raster files from raw EPA AERMOD data."""
        self.label = "Toolbox"
        self.alias = "toolbox"

        # List of tool classes associated with this toolbox
        self.tools = [Tool1, Tool2, Tool3, Tool4, Tool5, Tool6, Tool7, Tool8]


class Tool1(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Aggregate Media"
        self.description = "Aggregate media files: for each unique Id (cell), take the sum of the conc values so that there is only one conc value for each cell per chemical"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Year folders directory location",
                            name="Year_dir",
                            datatype="GPString",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="Chemical Number",
                            name="chem_num",
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
        
        Year_dir=parameters[0].valueAsText
        chem_num=parameters[1].valueAsText

        Year_list= [f.path for f in os.scandir(Year_dir) if f.is_dir()]
        for file in Year_list:
            chem_folder=glob.glob(os.path.join(file, f'Chem_{chem_num}'))
            media_csv_files=glob.glob(os.path.join(file, f'Chem_{chem_num}', '*.csv'))
            if chem_folder:  # Check if chem_folder is not empty
                # Extract the first element from the list
                chem_folder = chem_folder[0]
                agg_dir= chem_folder + '\\' + 'Aggregated_Media'
                os.makedirs(agg_dir)
                arcpy.AddMessage(f"Created folder: {agg_dir}")
                for media in media_csv_files:
                    current_file=os.path.basename(os.path.normpath(media))
                    df= pd.read_csv(media)
                    ID_group=df.groupby([f'uID_{chem_num}'], as_index=False).agg({f'rels_{chem_num}': 'first', f'chem_{chem_num}': 'first', f'facl_{chem_num}': 'first', f'media_{chem_num}': 'first', f'conc_{chem_num}': 'sum', f'toxc_{chem_num}': 'sum', f'score_{chem_num}': 'sum', f'canc_{chem_num}': 'sum', f'nocan_{chem_num}': 'sum', f'pop_{chem_num}': 'first'})
                    ID_group.to_csv(agg_dir + '\\' + current_file, index=False, header=True)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class Tool2(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Copy Aggregate Media into Individual GeoDatabses"
        self.description = "batch wise copy rows of aggregated media into appropriate year gdbs (example with Diazinon)"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Aggregated Media Directory Path",
                            name="agg_dir",
                            datatype="DEFolder",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="Folder where new geodatabses for each year will be created",
                            name="output_gdb",
                            datatype="DEFolder",
                            parameterType="Required",
                            direction="Input"),

            arcpy.Parameter(displayName="First Year of Data",
                            name="startyr",
                            datatype="GPLong",
                            parameterType="Required",
                            direction="Input"),

            arcpy.Parameter(displayName="Last Year of Data",
                            name="endyr",
                            datatype="GPLong",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(
                displayName="Chemical number",
                name="chem_num",
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
        
        agg_dir=parameters[0].valueAsText
        output_gdb=parameters[1].valueAsText
        startyr=parameters[2].valueAsText
        endyr=parameters[3].valueAsText
        chem_num=parameters[4].valueAsText


        # Define the base path
        base_path = agg_dir

        # Iterate over the years 2001-2020
        for year in range(int(startyr), (int(endyr)+1)):
            year_str = str(year)
            agg_med_dir = os.path.join(base_path, year_str, f"Chem_{chem_num}", "Aggregated_Media")

            # Create a year-specific geodatabase
    
            year_geodatabase = os.path.join(output_gdb, year_str + ".gdb")
            out_year_GDB= year_str + ".gdb"
            if not arcpy.Exists(year_geodatabase):
                arcpy.management.CreateFileGDB(output_gdb, out_year_GDB)

                arcpy.AddMessage(f"Database {year_geodatabase} created.")

            agg_file_list = [f for f in iglob(agg_med_dir, recursive=True) if os.path.isdir(f)]

            for file in agg_file_list:
                individual_med_file = glob.glob(os.path.join(file, '*.csv'))
                for med_file in individual_med_file:
                    in_table = str(med_file)
                    med_file_norm = os.path.normpath(in_table)
                    med_file_piece = med_file_norm.split(os.sep)
                    current_media = str(med_file_piece[-1])[:-4]
                    specific_out = os.path.join(year_geodatabase, "Aggregated" + '_' + current_media)

                    # Check if the specific_out table already exists
                    if not arcpy.Exists(specific_out):
                        arcpy.AddMessage(f"{in_table}, {specific_out} will be created.")
                        arcpy.CopyRows_management(in_table, specific_out)
                        arcpy.AddMessage(f"Table {specific_out} created.")
                    else:
                        arcpy.AddMessage(f"Table {specific_out} already exists. Skipping.")

        arcpy.AddMessage("Processing completed.")


    ###    gdb_dir = r"E:\AERMOD\Joined_shps_Diazinon\USA_*"
    ###    folder_list = [f for f in iglob(gdb_dir, recursive=True) if os.path.isdir(f)]
    ###    for file in folder_list:
    ###        arcpy.FeatureClassToGeodatabase_conversion(["D:\AERMOD\poly_gc14_conus_810m_Orig.shp"], file)


        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return
    
class Tool3(object): 

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "QueryTableCreation"
        self.description = "Joining empty shp with media chem data and exporting to the correct gdb that holds all years, all media for one chemcial, LOOP THROUGH  join via makiing query table (M2 only... must do again for M1)"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(
                displayName="Input Geodatabase parent folder",
                name="input_gdb_root",
                datatype="DEWorkspace",
                parameterType="Required",
                direction="Input"
                ),
                # Output geodatabase (required, workspace)
            arcpy.Parameter(
                displayName="Output Geodatabase Root path",
                name="output_gdb",
                datatype="DEWorkspace",
                parameterType="Required",
                direction="Input"
            ),
            
            # Table names (required, string)
            arcpy.Parameter(
                displayName="M1 or M2",
                name="media",
                datatype="GPString",
                parameterType="Required",
                direction="Input"
            ),
            
            arcpy.Parameter(
                displayName="Chemical number",
                name="chem_num",
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
        
        # Gather user inputs for paths and parameters
        input_gdb_root = parameters[0].valueAsText  # Get root directory of input geodatabases
        output_gdb = parameters[1].valueAsText  # Get output geodatabase path
        media = parameters[2].valueAsText  # Get comma-separated table names
        chem_num = parameters[3].valueAsText # Get common field name for joining

        # Iterate through input geodatabases
        for gdb_folder in iglob(os.path.join(input_gdb_root, "*"), recursive=True):
            if os.path.isdir(gdb_folder):
                arcpy.env.workspace = gdb_folder

                # Extract year from folder name
                year = os.path.basename(gdb_folder).split("_")[-1][:4]

                # Construct query table name and output feature class path
                qt_name = f"QueryTable_{year}_{media}"
                output_fc = os.path.join(output_gdb, f"USA_{year}_{media}")

                if media == "M1":
                    media_full = "Media_1"
                else:
                    media_full = "Media_2"

                join_shp = f"Aggregated_{media_full}"

                # Print debug information
               

                arcpy.AddMessage(f"Processing: {gdb_folder}")
                arcpy.AddMessage(f"Year: {year}")
                arcpy.AddMessage(f"Query table name: {qt_name}")
                arcpy.AddMessage(f"Output feature class path: {output_fc}")
                arcpy.AddMessage(f"Join shapefile: {join_shp}")
                arcpy.AddMessage(f"Fields in poly_gc14_conus_810m_Orig: {arcpy.ListFields('poly_gc14_conus_810m_Orig')}")
                arcpy.AddMessage(f"Fields in {join_shp}: {arcpy.ListFields(join_shp)}")

                
                if arcpy.Exists(output_fc):
                    arcpy.AddMessage(f"Query table {output_fc} already exists. Skipping to the next iteration.")
                else:
                    # Create query table with user-specified table names and join field
                    arcpy.management.MakeQueryTable(
                        f"poly_gc14_conus_810m_Orig;{join_shp}",
                        qt_name,
                        "USE_KEY_FIELDS",
                        None,
                        None,
                        f"poly_gc14_conus_810m_Orig.UniqueID = {join_shp}.uID_{chem_num}"
                    )

                    arcpy.AddMessage(f"Query table {qt_name} created successfully.")
                    # Copy features if the query table was created
                    arcpy.management.CopyFeatures(qt_name, output_fc)
                    arcpy.AddMessage(f"Query table: {qt_name} saved to {output_fc}")
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class Tool4(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Empty Shape Grid Copying"
        self.description = "Copies empty shp grid into each individual yearly gdb"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Base Folder",
                            name="base_folder",
                            datatype="DEFolder",
                            parameterType="Required",
                            direction="Input"),

            arcpy.Parameter(displayName="Input Shapefile",
                            name="input_shapefile",
                            datatype="DEShapefile",
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
        
        base_folder = parameters[0].valueAsText
        input_shapefile = parameters[1].valueAsText

        # Iterate over each year folder (GDB) and perform the FeatureClassToGeodatabase_conversion
        for year_folder in os.scandir(base_folder):
            if year_folder.is_dir():
                gdb_path = year_folder.path
                arcpy.FeatureClassToGeodatabase_conversion([input_shapefile], gdb_path)
                arcpy.AddMessage(f"Conversion completed for GDB in year: {os.path.basename(gdb_path)}")
        arcpy.AddMessage("Tool execution completed.")
        
        return


    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class Tool5(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Set Null Data to Zero"
        self.description = "This tool iterates over a geodatabase of rasters, identifies null data cells in each raster, and sets these Null values to Zero."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Folder of Rasters (with null values)",
                            name="input_gdb",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),
                                        
            arcpy.Parameter(displayName="Output Location for Zero Rasters",
                            name="Con_Output",
                            datatype="DEWorkspace",
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
        input_gdb = parameters[0].valueAsText
        Con_Output = parameters[1].valueAsText

        arcpy.env.workspace = input_gdb

        rasters = arcpy.ListRasters()

        for raster in rasters:
            # Process: Is Null
            is_null_raster = arcpy.sa.IsNull(raster)

            # Process: Con
            output_raster = arcpy.sa.Con(is_null_raster, 0, raster, "Value = 1")
            output_path = f"{Con_Output}" + "/" + f"{raster}" + "_Zero"
            output_raster.save(output_path)
            arcpy.AddMessage(f"{raster} set to zero.")

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class Tool6(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add Media 1 and Media 2 Rasters"
        self.description = "***"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Folder of Rasters (with M1)",
                            name="input_gdb",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),
                                        
            arcpy.Parameter(displayName="Folder of Rasters (with M2)",
                            name="input_gdb_2",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),
                                        
            arcpy.Parameter(displayName="Output Location for Added Rasters",
                            name="Add_Output",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="First Year of Data",
                            name="startyr",
                            datatype="GPLong",
                            parameterType="Required",
                            direction="Input"),

            arcpy.Parameter(displayName="Last Year of Data",
                            name="endyr",
                            datatype="GPLong",
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
        input_gdb = parameters[0].valueAsText
        input_gdb_2 = parameters[1].valueAsText
        Add_Output = parameters[2].valueAsText
        startyr = parameters[3].valueAsText
        endyr = parameters[4].valueAsText

        import arcpy
        import os
        import arcpy.sa

        # Set the workspaces for M1 and M2 rasters
        workspace_m1 = input_gdb
        workspace_m2 = input_gdb_2

        # List all M1 and M2 rasters in their respective workspaces
        arcpy.env.workspace = workspace_m1
        rasters_m1 = arcpy.ListRasters("*", "ALL")
        arcpy.env.workspace = workspace_m2
        rasters_m2 = arcpy.ListRasters("*", "ALL")

        # Loop through each year
        for year in range(int(startyr), int(endyr) + 1):
            # Filter M1 and M2 rasters for the current year
            m1_rasters_for_year = [raster for raster in rasters_m1 if f"_{year}_" in raster]
            m2_rasters_for_year = [raster for raster in rasters_m2 if f"_{year}_" in raster]

            # Check if there are rasters for both M1 and M2 for the current year
            if len(m1_rasters_for_year) > 0 and len(m2_rasters_for_year) > 0:
                # Construct the full paths to the input rasters
                full_path_m1 = os.path.join(workspace_m1, m1_rasters_for_year[0])
                arcpy.AddMessage(f"{full_path_m1}")
                full_path_m2 = os.path.join(workspace_m2, m2_rasters_for_year[0])
                arcpy.AddMessage(f"{full_path_m2}")

                # Perform raster addition
                outrast = arcpy.sa.Plus(arcpy.sa.Raster(full_path_m1), arcpy.sa.Raster(full_path_m2))

                # Save the output raster
                output_name = f"USA_{year}_M1M2"
                output_path = os.path.join(Add_Output, output_name)
                arcpy.management.CopyRaster(outrast, output_path)  # Save the output in M1M2 workspace (you can change this if desired)
                arcpy.AddMessage(f"{output_path} has been created.")
            else:
                arcpy.AddMessage(f"Not enough rasters found for year {year} in both M1 and M2 workspaces. Skipping.")


        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return
    

class Tool7(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multiply all rasters by 1 Million"
        self.description = "Multiplies the rasters in a gdb by 1 million, Creates a output gdb based on the media number, and sorts the outputs into the appropriate media gdb."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Folder of Rasters (with M1 and M2 to be multiplied)",
                            name="input_workspace",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),
                                        
            arcpy.Parameter(displayName="Output workspace where M1_Million and M2_Million geodatabases will be created",
                            name="output_workspace",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="Chemical Number",
                            name="chem_num",
                            datatype="GPLong",
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
        
        input_workspace = parameters[0].valueAsText
        output_workspace = parameters[1].valueAsText
        chem_num = parameters[2].valueAsText


        
        # Set input workspace
        arcpy.env.workspace = input_workspace
        
        # List all rasters in the geodatabase
        rasters = arcpy.ListRasters()

        # Create M1_Million and M2_Million geodatabases if they don't exist
        M1_gdb = os.path.join(output_workspace, f"Orig_{chem_num}_M1_Mill.gdb")
        M2_gdb = os.path.join(output_workspace, f"Orig_{chem_num}_M2_Mill.gdb")

        if not arcpy.Exists(M1_gdb):
            arcpy.CreateFileGDB_management(M1_gdb)
            arcpy.AddMessage(f"{M1_gdb} created.")

        if not arcpy.Exists(M2_gdb):
            arcpy.CreateFileGDB_management(M2_gdb)
            arcpy.AddMessage(f"{M2_gdb} created.")

        # Loop through each raster
        for raster in rasters:
            # Define the output geodatabase based on raster name
            year = raster.split("_")[1]
            media = raster.split("_")[2]
            output_gdb = M1_gdb if media == "M1" else M2_gdb

            # Construct the output raster path
            output_raster_name = f"Chem_{chem_num}_{year}_{media}_Mill"  
            output_raster = os.path.join(output_gdb, output_raster_name)
            
            # Multiply raster values by 1,000,000
            out_raster = arcpy.sa.Times(raster, 1000000)
            
            # Save the multiplied raster
            out_raster.save(output_raster)

            arcpy.AddMessage(f"{raster} processed and saved to {output_raster}")

        return


    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return


class Tool8(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert datatype to INT"
        self.description = "This script creates new geodatabases within a folder named INT at the same level as the input M1_Million and M2_Million geodatabases. It then copies the contents of M1_Million and M2_Million geodatabases into these new geodatabases. Finally, it converts float fields to integer fields within the new geodatabases."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = [
            arcpy.Parameter(displayName="Folder of Rasters (with M1 or M2 to be converted to INT)",
                            name="input_workspace",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),
                                        
            arcpy.Parameter(displayName="Output workspace where M1_INT and M2_INT geodatabases will be created",
                            name="output_workspace",
                            datatype="DEWorkspace",
                            parameterType="Required",
                            direction="Input"),
            
            arcpy.Parameter(displayName="Chemical Number",
                            name="chem_num",
                            datatype="GPLong",
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
        
        input_workspace = parameters[0].valueAsText
        output_workspace = parameters[1].valueAsText
        chem_num = parameters[2].valueAsText
        

        # Create INT folder within the input workspace
        int_folder = os.path.join(output_workspace, "INT")
        os.makedirs(int_folder, exist_ok=True)

        media = os.path.basename(input_workspace).split("_")[2]
        

        # Create new geodatabase Orig_chemnum_INT
        int_gdb = os.path.join(int_folder, f"Orig_{chem_num}_{media}_INT.gdb")

        if not arcpy.Exists(int_gdb):
            arcpy.CreateFileGDB_management(int_folder, f"Orig_{chem_num}_{media}_INT.gdb")
            arcpy.AddMessage(f"Orig_{chem_num}_{media}_INT.gdb created.")

        # Define input and output geodatabases
        arcpy.env.workspace = input_workspace
        rasters = arcpy.ListRasters()

        for raster in rasters:
            # Convert raster to integer using raster calculator
            year = raster.split("_")[2]
            output_raster_path = os.path.join(int_gdb, f"Chem_{chem_num}_{year}_{media}_INT")
            out_raster = arcpy.sa.Int(raster)
            out_raster.save(output_raster_path)
            arcpy.AddMessage(f"{out_raster} successfully converted to INT and saved to {output_raster_path}")

        return


    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return
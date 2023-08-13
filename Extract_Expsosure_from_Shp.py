import os
import arcpy
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")

# Set the current workspace (your gdb of environmental rasters)
arcpy.env.workspace = "D:\AERMOD\AllYears_346_347_Sum.gdb"   

# set the case folder including all the case files (.shp)
# This file path name should be a FOLDER that contains .shp files. 

caseFolder = "D:\Exposure_Analyses\ExposureAnalysis_07072023_PbPbC\ExposureReadings\Mortality"

# Read all raster layers in current workplace
rasters = arcpy.ListRasters("*", "all")

inRasterList = []

for raster in rasters:
    year = raster.split("_")[1]  #Assuming the file name format is "name_year.shp"
    rasterName = "Sum_" + year  #Change the string to match the names of your raster files' names
    fieldName = "pb_" + year  #This is the name of the columns of exposures that will be created
    inRasterList.append([rasterName, fieldName])


for file in os.listdir(caseFolder):
    if file.endswith(".shp"):
        CasesSet = caseFolder + "/" + file
        print(CasesSet)

        ExtractMultiValuesToPoints(CasesSet,inRasterList,"NONE")

	
	

 
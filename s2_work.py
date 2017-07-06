"""
Description of this script here.
Rather than making a complicated script which does a lot of things
going to try and keep it simple, limit it to as few tasks as possible.

So:
1. Open a S2 data product from a zip or .SAFE format. 
Uses a slightly modified tiff-generator.py script for this task.
Has a dependency on the Scripts folder which contains gdal_merge.py

2. Retrieve the metadata which is commonly used in atmospheric correction such as the sun zenith angle. 

3. Handle the image arrays from the merged data files. 
Adjust the DN within the array so that the newly created tiff files contain pixels representing
TOA radiance instead of reflectance.

4. Define AC (Atmospheric Correction) models here.

5. A main function for whenever this script is called independently.
"""

import tiffgenerator as tg #Necessary for the running of: 1,2
import os #Necessary for the running of: 2
import numpy as np #Necessary for the running of: 2, 5
import rastercr as rc #Necessary for the running of: 3
import gdal #Necessary for the running of 5
from gdalconst import * #Necessary for the running of: 3, 5
import time #Necessary for the running of: 5

#==============================================================================================
#1: Opening of a zipped S2 L1C Product.
def open_s2():
	"""
	Requests an input location for some S2L1C product, a location to store the product and then
	generates a merged .tif file with all the bands.
	Modified the original tiff-generator code to work when provided with an extracted .SAFE file
	Returns the requested paths.
	"""
	inputPath = input("Please enter the location of the S2A L1C product to be extracted: \n")
	outputPath = input("Please choose a location to store the extracted and processed data: \n")
	if not outputPath[len(outputPath)-1] == "/":
		#log that an improper output path name was initally chosen
		outputPath = outputPath + "/" #Makes sure that the data is stored in a directory.
	tg.generate_geotiffs(inputPath,outputPath) #Having added an if __main__ check, the generate_geotiffs now works as expected, for a zip input.
	return inputPath, outputPath

#==============================================================================================

#==============================================================================================
#2: Retrieve the metadata that is contained within the MTD_TL and MTD_MSI xml files.
def metadata_get(filename, metadata):
	"""	
	Opens a file which contains some metadata that the user wishes to extract from an .xml file.
	For example, sun azimuth and zenith angles.
	Within the structure of the xml file, the relevant information will be contained within a
	section of code that will be written to a file when found.
	"""
	mtd_filname = metadata + "_meta_1.txt"
	i=1
	while os.path.exists(os.getcwd()+"/" + mtd_filname):
		#While a metadata file exists with a matching name, loop through some integers until a 
		#unique identifier is found to work as a suffix for the filename.
		#Prevents overwriting old data or bloating the metadata files.	
		x = int(mtd_filname[len(metadata)+6:-4])
		i+=1
		mtd_filname = metadata + "_meta_" + str(i) + ".txt"
	with open(filename, "r") as meta_file:
		pmt = 1 #Sets a parameter that can be turned on/off when the relevant data is found.
		for lines in meta_file: #Reads through the metadata file provided.
			if metadata in lines:
				pmt=pmt*-1 #Switches the parameter when it is found.
			if pmt<0:
				#print("Writing " + metadata +" to file.")
				with open(mtd_filname,"a") as out_file:
					out_file.write(lines)
	#if not args.quiet:	
	print("Finished creating metadata relating to %s file." %(metadata))

def selected_metadata(inputPath, outputPath):
	"""
	This function uses the metadata_get function to create the desired meta-files.
	"""
	#Get the L1C granule subdirectory path necessary to access that metadata file.
	GranuleSubDir = tg.get_immediate_subdirectories(outputPath+"/"+inputPath[:60] +".SAFE/GRANULE")
	#The following are the desired metadata files that need to be created
	metadata_get(outputPath+"/" + inputPath[:60] + ".SAFE/MTD_MSIL1C.xml","Reflectance_Conversion") #Gets the solar reflectance coefficient and the solar irradiances.
	metadata_get(outputPath+"/" + inputPath[:60] + ".SAFE/GRANULE/"+GranuleSubDir[0]+"/MTD_TL.xml","Mean_Sun") #Gets the solar zenith and azimuth angles.
	metadata_get(outputPath+"/" + inputPath[:60] + ".SAFE/GRANULE/"+GranuleSubDir[0]+"/MTD_TL.xml","Mean_Viewing_Incidence_Angle_List") #Obtains the average viewing angles.


def sol_irr(solar_file):
	"""
	This function opens a 'solar file', which contains the information produced for the "Reflectance Conversion" .txt file above.
	The function iterates through the file and produces corrected solar irradiances for each band.
	"""
	sol_irr = []
	with open(solar_file,"r") as sol_info:
		for lines in sol_info:
			if "<U>" in lines:
				conv_coef = float(lines.strip().replace("<U>","</U>").strip("</U>"))
			elif "SOLAR_IRR" in lines:
				x = lines.replace("W/m\xc2\xb2/\xc2\xb5m","x").strip("SOLAR_IRRADIANCE")
				x = lines.replace("W/m²/µm","x").strip("SOLAR_IRRADIANCE")
				sol_irr = np.append(sol_irr,(x.strip().replace('<SOLAR_IRRADIANCE bandId="',"Band").replace('" unit="x">',"|").replace("</SOLAR_IRRADIANCE>","")))
	return sol_irr, conv_coef

def sun_ang(sun_file):
	"""
	Short function that cycles through the small text files which contain the sun zenith and azimuth angles.
	"""
	angles = []
	with open(sun_file,"r") as sun_info:
		for line in sun_info:
			if "ANGLE" in line:
				x = line.split('"deg">')
				x = x[1].split('</')
				angles = np.append(angles,float(x[0]))
	return angles #Returns zenith and azimuth angles respectively

def view_ang(view_file):
	"""
	Produce a dictionary with keys that correspond to the S2 bands and items:
	band, zenith angle and azimuthal angle.
	Read in a dictionary with zenith and azimuth viewing angles.
	Returns the pair of values corresponding to the input bandId, in the format zenith, azimuth.
	"""
	view_arr = []
	with open(view_file, "r") as view_info:
		for lines in view_info:
			view_arr = np.append(view_arr, lines.strip()) #Remove the whitespace for each line.
	view_arr = view_arr[1:-1] #Removes the first and last lines which did not contain desired data.
	means = []
	means_dict = {}
	for term in view_arr:
		if "bandId" in term:
			#means[term[30:-1]]=term[37:-1] #assigns a key with the band id
			#This is only for S2 at the moment, might not apply to other metadata
			y=int(term[38:-2]) #stores the band id
			if y < 8:
				y+=1
			elif y == 8: 
				y=str(y)+"A"
			means = np.append(means,("B"+str(y)).replace('"',''))
		elif "ZENITH" in term:
			#adds zenith viewing angle to that band key
			means = np.append(means,term[25:-15])
		else:
			#adds azimuth viewing angle to that band key
			means= np.append(means,term[26:-16])
		for i in range(0,len(means),4):
			means_dict[means[i].replace('"','')] = means[i:i+3]
	#Set i to 1 to iterate through the bands and return a sorted list of zeniths and azimuths
	i = 1
	view_zen = [] #A list to store the view zenith angles.
	view_azi = [] #A list to store the view azimuth angles.
	while str(i) != "13":
		zen, azi = means_dict["B"+str(i)][1], means_dict["B"+str(i)][2]
		view_zen = np.append(view_zen,zen)
		view_azi = np.append(view_azi, azi)
		if i == 8:
			i = "8A"
		elif i == "8A":
			i = 9
		else:
			i+=1
	return view_zen, view_azi
#==============================================================================================

#==============================================================================================
#3: Create an array of top-of-the-atmosphere pixel values.
def create_arr(tif_data, band_no,cols, rows, ):
	"""
	Reads in a data file such as the merged tif created in Part 1.
	Uses gdal to convert this to an array for some band within the tif.
	"""
	band = tif_data.GetRasterBand(band_no) #Obtains the band information
	data = band.ReadAsArray(0,0, cols,rows) #Sets it as an array.
	return data

def S2_adjust(band_id, solar_corrected, sun_zenith):
	"""
	Reads in the array for some S2 band and adjusts the pixel data to return TOA radiance.
	rho = pi*L/E_s * cos(Theta_s)
	"""
	scale_factor = 10000 #Scale factr used for the reflectance values.
	band_id = np.array(band_id)
	band_id = band_id * float(solar_corrected) *np.cos(sun_zenith*(np.pi/180)) / (np.pi * scale_factor)
	return band_id

def create_tiff(data,name,rasterOrigin,pixelWidth,pixelHeight):
	"""
	Re-normalize to GS range (0-255).
	Then creates a tiff for some input data.
	"""	
	max_data = 0
	for i in range(len(data)):
		if max(data[i]) > max_data:	
			max_data = max(data[i])
	data = data/(max_data/255)
	newRasterfn = 'input_' + name + '.tif'
	rc.main(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,data)
	print("Created tiff")
	return data

#Should add function to handle the creation of bands that don't come from the .zip or .SAFE files, e.g. re-merging band images after correction. Would be useful, but not too sure about how to do it just yet.
#==============================================================================================

#==============================================================================================
#4.A: The SMAC (Simplified Method for Atmospheric Correction) function:
def smac_func(band_no,phi_v,phi_s,theta_v,theta_s):
	"""
	Declare source of the programme.
	Requires access to the folder "COEFS" which contains the data files for the 
	smac coefficients.
	"""
	nom_smac = 'COEFS/Coef_S2A_CONT_B' + band_no + '.dat'
	coefs=smac.coeff(nom_smac)
	r_surf = smac.smac_inv(r_toa , theta_s, phi_s, theta_v[band_no], phi_v[band_no],1013,0.1,0.3,0.3, coefs)
	r_toa2 = smac.smac_dir(r_surf, theta_s, phi_s, theta_v[band_no], phi_v[band_no],1013,0.1,0.3,0.3, coefs)
	return r_surf, r_toa2

#==============================================================================================

#==============================================================================================
#5: Combine all the functions to create a sample output when this script is run by itself.
def main():
	inputPath, outputPath = open_s2() #Opens a Sentinel-2 product
	selected_metadata(inputPath, outputPath) #Create metadata files.
	#Corrects the solar irradiance for the reflectance based on Earth's position at measurement:
	sol_irr_at_sat, correction_coefficient = sol_irr("Reflectance_Conversion_meta_1.txt") #still using first metadata files created
	sol_cor = []
	for i in sol_irr_at_sat:
		sol_cor = np.append(sol_cor,i.split("|")[0])
		sol_cor = np.append(sol_cor,(float(i.split("|")[1]))*correction_coefficient) #sol_cor is stored in the form of a list [BAND_ID, VALUE, BAND_ID, VALUE,...]
	#Produce a tuple for the average solar zenith and azimuth angles.
	sun_zen, sun_azi = sun_ang("Mean_Sun_meta_1.txt") #still using first metadata files created
	#Produce a tuple for the average viewing zenith and azimuth angles.
	view_zen, view_azi = view_ang("Mean_Viewing_Incidence_Angle_List_meta_1.txt") #still using first metadata files created
	print(view_zen)
	#Correct each band from reflectance to radiance.
	data_file = (outputPath+"/" + inputPath[:60] + "_PROCESSED/merged.tif")
	dataset = gdal.Open(data_file, GA_ReadOnly) #Creates the gdal osgeo class object.
	cols = dataset.RasterXSize #Retrieve the number of columns
	rows = dataset.RasterYSize #Retrieve the number of rows
	bands = dataset.RasterCount #Retrieve the number of bands; 13 for S2 product.
	geotransform = dataset.GetGeoTransform() #Tuple with geoinformation for product.
	originX = geotransform[0]
	originY = geotransform[3]
	rasterOrigin = (originX, originY)
	pixelWidth = geotransform[1]
	pixelHeight = geotransform[5]
	for i in range(bands,0,-1):
		band = create_arr(dataset, i ,cols, rows)
		band_rad = S2_adjust(band, sol_cor[i*2 - 1], sun_zen)
		band_gs = create_tiff(band_rad, "band_" +str(i),rasterOrigin,pixelWidth,pixelHeight)
#==============================================================================================

#==============================================================================================
if __name__ == "__main__":
	start_time = time.time()
	main()
	print("\n === Completed after %s seconds ===" % (time.time() - start_time))
#==============================================================================================

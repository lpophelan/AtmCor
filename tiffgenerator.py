import sys
import subprocess
sys.path.append('Scripts/')
import gdal_merge
import zipfile  
import os
import time
import readline, glob
from pathlib import Path

def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]


def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]


def generate_geotiffs(inputProductPath, outputPath):

	#Trying to handle the case where tiffgenerator is fed a .SAFE file
	if inputProductPath[-4:] == "SAFE":
		basename = os.path.basename(inputProductPath[:-1])
		iPP = inputProductPath + "/"
		subprocess.call(["mkdir",outputPath])
		subprocess.call(["cp","-R",iPP,outputPath])
	elif inputProductPath[-3:] == "zip":
		basename =  os.path.basename(inputProductPath)
		#print(basename[:-3])
		if os.path.isdir(outputPath + basename[:-3] + "SAFE") :
			print('Already extracted')
		else:
			zip = zipfile.ZipFile(inputProductPath) 
			zip.extractall(outputPath) 
			print("Extracting Done") 
	else:
		print("File format supplied not currently supported.")
	

	directoryName = outputPath + basename[:-3] + "SAFE/GRANULE"
	#productName = os.path.basename(inputProductPath)[:-4]
	productName = basename[:-4]
	outputPathSubdirectory = outputPath + productName + "_PROCESSED"

	if not os.path.exists(outputPathSubdirectory):
		os.makedirs(outputPathSubdirectory)

	subDirectorys = get_immediate_subdirectories(directoryName)

	results = []

	for granule in subDirectorys:
		unprocessedBandPath = outputPath + productName + ".SAFE/GRANULE/" + granule + "/" + "IMG_DATA/"
		results.append(generate_all_bands(basename,unprocessedBandPath, granule, outputPathSubdirectory))
	
	#gdal_merge.py -n 0 -a_nodata 0 -of GTiff -o /home/daire/Desktop/merged.tif /home/daire/Desktop/aa.tif /home/daire/Desktop/rgbTiff-16Bit-AllBands.tif
	merged = outputPathSubdirectory + "/merged.tif"
	params = ['',"-of", "GTiff", "-o", merged]

	for granule in results:
		params.append(granule)

	gdal_merge.main(params)


def generate_all_bands(basename,unprocessedBandPath, granule, outputPathSubdirectory):

	granuleBandTemplate = basename[38:-20]+basename[10:27]
	
	outputPathSubdirectory = outputPathSubdirectory 
	if not os.path.exists(outputPathSubdirectory+ "/IMAGE_DATA"):
		os.makedirs(outputPathSubdirectory+ "/IMAGE_DATA")
	
	outPutTiff = '/'+granule[:-6]+'16Bit-AllBands.tif'
	outPutVRT = '/'+granule[:-6]+'16Bit-AllBands.vrt'

	outPutFullPath = outputPathSubdirectory + "/IMAGE_DATA/" + outPutTiff
	outPutFullVrt = outputPathSubdirectory + "/IMAGE_DATA/" + outPutVRT
	inputPath = unprocessedBandPath + granuleBandTemplate
	#print(inputPath)
	print("In progress.... \n")

	bands = {"B1" :  inputPath + "B01.jp2",
	"B2" :  inputPath + "B02.jp2",
	"B3" :  inputPath + "B03.jp2",
	"B4" :  inputPath + "B04.jp2",
	"B5" :  inputPath + "B05.jp2",
	"B6" :  inputPath + "B06.jp2",
	"B7" :  inputPath + "B07.jp2",
	"B8" :  inputPath + "B08.jp2",
	"B8A" :  inputPath + "B8A.jp2",
	"B9" :  inputPath + "B09.jp2",
	"B10" :  inputPath + "B10.jp2",
	"B11" :  inputPath + "B11.jp2",
	"B12" :  inputPath + "B12.jp2"}

	#print(bands)

	cmd = ['gdalbuildvrt', '-resolution', 'user', '-tr' ,'20', '20', '-separate' ,outPutFullVrt]


	for band in sorted(bands.values()):
		cmd.append(band)
           
	my_file = Path(outPutFullVrt)
	if not my_file.is_file():
		# file exists
		subprocess.call(cmd)

	#, '-a_srs', 'EPSG:3857'
	cmd = ['gdal_translate', '-of' ,'GTiff', outPutFullVrt, outPutFullPath]

	my_file = Path(outPutTiff)
	if not my_file.is_file():
		# file exists
		subprocess.call(cmd)



	#params = ['', '-o', outPutFullPath, '-separate', band_01, band_02, band_03, band_04, band_05, band_06, band_07, band_08, band_8A, band_09, band_10, band_11, band_12]

	#gdal_merge.main(params)
	
	return(outPutFullPath)

if __name__ == "__main__":
	start_time = time.time()
	outputPath = '../Output/'
	readline.set_completer_delims(' \t\n;')
	readline.parse_and_bind("tab: complete")
	readline.set_completer(complete)
	inputPath = input("Input Path? ")
	generate_geotiffs(inputPath, outputPath)
	print("--- %s seconds ---" % (time.time() - start_time))

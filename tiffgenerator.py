import sys
import subprocess
sys.path.append('Scripts/')
import gdal_merge
import zipfile  
import os
import time
import readline, glob
from pathlib import Path
import logging

logging.basicConfig(filename="S2A-AtmCor.log",format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# create logger
logger = logging.getLogger('Tiff Generator')
logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)

def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]


def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def generate_geotiffs(inputProductPath, outputPath):
    """
    Rewrite to match PEP8 guidelines, i.e. if inputProductPath ends with zip 
    to catch zip/SAFE format.
    Use .startswith() to catch L1 or L2 products.
    """

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
            logger.info('Already extracted')
        else:
            zip = zipfile.ZipFile(inputProductPath) 
            zip.extractall(outputPath) 
            logger.info("Extracting Done") 
    else:
        logger.warning("File format supplied not currently supported.")
        raise NotImplementedError

    directoryName = outputPath + basename[:-3] + "SAFE/GRANULE"
    #productName = os.path.basename(inputProductPath)[:-4]
    productName = basename[:-4]
    outputPathSubdirectory = outputPath + productName + "_PROCESSED"

    if not os.path.exists(outputPathSubdirectory):
        os.makedirs(outputPathSubdirectory)

    subDirectorys = get_immediate_subdirectories(directoryName)

    results = []

    for granule in subDirectorys:
        unprocessedBandPath = ''.join([
            outputPath,productName,".SAFE/GRANULE/",granule,"/IMG_DATA/"])
        results.append(generate_all_bands(
            basename,unprocessedBandPath, granule, outputPathSubdirectory))
    #gdal_merge.py -n 0 -a_nodata 0 -of GTiff -o /home/daire/Desktop/merged.tif /home/daire/Desktop/aa.tif /home/daire/Desktop/rgbTiff-16Bit-AllBands.tif
    merged = outputPathSubdirectory + "/merged.tif"
    params = ['',"-of", "GTiff", "-o", merged]
    logger.debug("Params: %s", params)
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

    outPutFullPath = ''.join([outputPathSubdirectory,"/IMAGE_DATA/",outPutTiff])
    outPutFullVrt = ''.join([outputPathSubdirectory,"/IMAGE_DATA/",outPutVRT])
    inputPath = ''.join([unprocessedBandPath,granuleBandTemplate])
    logger.debug("In progress.... \n")

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

    cmd = ['gdalbuildvrt', '-resolution', 'user', '-tr' ,'20', '20', '-separate' ,outPutFullVrt]
    logger.debug("Command: %s", cmd)

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
    logger.info("--- %s seconds ---", time.time() - start_time)
    print("--- %s seconds ---" % (time.time() - start_time))

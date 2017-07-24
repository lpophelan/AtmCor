#!/usr/bin/env python
# encoding=utf-8
"""
Description of this script here.
Rather than making a complicated script which does a lot of things
going to try and keep it simple, limit it to as few tasks as possible.

So:
1. Open a S2 data product from a zip or .SAFE format. 
Uses a slightly modified tiff-generator.py script for this task.
Has a dependency on the Scripts folder which contains gdal_merge.py

2. Retrieve the metadata which is commonly used in atmospheric correction such 
as the sun zenith angle. 

3. Handle the image arrays from the merged data files. 
Adjust the DN within the array so that the newly created tiff files contain 
pixels representing
TOA radiance instead of reflectance.

4. Define AC (Atmospheric Correction) models here. Currently it includes an 
implementation of the CNES SMAC tool and Py6S.

5. A main function for whenever this script is called independently.
"""

import logging
import os #Necessary for the running of: 2
import time #Necessary for the running of: 5
import gdal #Necessary for the running of: 5
import pickle #Necessary for the running of: 4
import numpy as np #Necessary for the running of: 2, 5
import subprocess as sp #Necessary for the running of: 5

from gdalconst import * #Necessary for the running of: 3, 5
from Py6S import * #Necessary for the running of: 4.B

import tiffgenerator as tg #Necessary for the running of: 1,2
import rastercr as rc #Necessary for the running of: 3
import smac #Necessary for the running of: 4.A

#Authorship information next:
__author__ = "Liam Phelan"
__version__ = "1.0"
__status__ = "Prototype"

logging.basicConfig(filename="S2A-2017_05_08.log",format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# create logger
logger = logging.getLogger('Main Script')
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


#=============================================================================
#1: Opening of a zipped S2 L1C Product.
def open_s2():
    """
    Requests an input location for some S2L1C product, a location to store the 
    product and then generates a merged .tif file with all the bands.
    Modified the original tiff-generator code to work when provided with an 
    extracted .SAFE file.
    Returns the requested paths.
    """
    inputPath = input(\
        "Enter the location of the S2A L1C product to be extracted: \n")
    outputPath = input(\
        "Choose a location to store the extracted and processed data: \n")
    if not outputPath.endswith("/"):
        logger.debug("Appended / to %s", outputPath)
        outputPath = ''.join([outputPath,"/"]) #Data stored in directory.
    tg.generate_geotiffs(inputPath,outputPath)
    logger.debug("The input file was: %s",inputPath)
    logger.debug("The processed data was stored in: %s", outputPath) 
    return inputPath, outputPath

#=============================================================================

#=============================================================================
#2: Retrieve metadata that's contained within the MTD_TL & MTD_MSI xml files.
def metadata_get(filename, metadata):
    """	
    Opens a file which contains some metadata that the user wishes to extract 
    from an .xml file.
    For example, sun azimuth and zenith angles.
    Within the xml file, the relevant information will be contained within a
    section of code that will be written to file when found.
    """
    mtd_filname = ''.join([metadata,"_meta_1.txt"])
    i = 1
    while os.path.exists('/'.join([os.getcwd(),mtd_filname])):
        #Prevents overwriting old data or bloating the metadata files.	
        x = int(mtd_filname[len(metadata)+6:-4])
        logger.debug("%s already exists, checking for v.%s",mtd_filname,i+1) 
        i += 1
        mtd_filname = ''.join([metadata,"_meta_",str(i),".txt"])
    with open(filename, "r",encoding='utf-8') as meta_file:
        pmt = 1 #Parameter to turn on/off when the relevant data is found.
        for lines in meta_file:
            if metadata in lines:
                pmt *= -1 #Switches the parameter when it is found.
            if pmt < 0:
                #print("Writing " + metadata +" to file.")
                with open(mtd_filname,"a",encoding='utf-8') as out_file:
                    out_file.write(lines)
    #if not args.quiet:	
    print("Finished creating metadata relating to %s file." %(metadata))
    logger.info("Created file containing %s metadata.", metadata)

def selected_metadata(inputPath, outputPath):
    """
    Function uses the metadata_get function to create the desired meta-files.
    """
    #Get L1C granule subdirectory path necessary to access that metadata file.
    GranuleSubDir = tg.get_immediate_subdirectories(''.join([
        outputPath,"/",inputPath[:60],".SAFE/GRANULE"]))
    #The following are the desired metadata files that need to be created:
    #Solar reflectance and irradiances.
    metadata_get(''.join(
        [outputPath,"/",inputPath[:60],".SAFE/MTD_MSIL1C.xml"]),
        "Reflectance_Conversion")
    #The global footprint
    metadata_get(''.join(
        [outputPath,"/",inputPath[:60],".SAFE/MTD_MSIL1C.xml"]),
        "Global_Footprint")
    #The wavelengths for the different bands.
    metadata_get(''.join(
        [outputPath,"/",inputPath[:60],".SAFE/MTD_MSIL1C.xml"]),
        "Wavelength")  
    #Solar zenith and azimuth angles.
    metadata_get(''.join([outputPath,"/",inputPath[:60],".SAFE/GRANULE/",
        GranuleSubDir[0],"/MTD_TL.xml"]),"Mean_Sun") 
    #Mean viewing angles.
    metadata_get(''.join([outputPath,"/",inputPath[:60],".SAFE/GRANULE/",
        GranuleSubDir[0],"/MTD_TL.xml"]),"Mean_Viewing_Incidence_Angle_List") 


def sol_irr(solar_file):
    """
    This function opens a 'solar file', which contains the information 
    produced for the "Reflectance Conversion" .txt file above.
    The function iterates through the file and produces corrected solar 
    irradiances for each band.
    """
    sol_irr = []
    with open(solar_file,"r",encoding='utf-8') as sol_info:
        for lines in sol_info:
            if "<U>" in lines:
                conv_coef = float(lines.strip().replace("<U>","</U>"). \
                   strip("</U>"))
            elif "SOLAR_IRR" in lines:
                x = lines.replace("W/m\xc2\xb2/\xc2\xb5m","x"). \
                    strip("SOLAR_IRRADIANCE")
                x = lines.replace("W/m²/µm","x").strip("SOLAR_IRRADIANCE")
                sol_irr = np.append(sol_irr,(x.strip(). \
                    replace('<SOLAR_IRRADIANCE bandId="',"Band").replace(
                    '" unit="x">',"|").replace("</SOLAR_IRRADIANCE>","")))
    logger.debug("Solar irradiance: %s", sol_irr)
    logger.debug("Conversion coefficient: %s", conv_coef)
    return sol_irr, conv_coef

def band_wave(wavelength_file):
    """
    Opens the provided wavelength file to convert the information from .xml 
    format to useable data.
    """
    with open(wavelength_file,"r",encoding='utf-8') as wave_info:
        min_wave = []
        max_wave = []
        mid_wave = []
        for line in wave_info:
            lam_wave = line.strip(" ")
            if lam_wave.endswith("MIN>\n"):
                mini = ''
                for i in lam_wave:
                    if i.isdigit() or i == ".":
                        mini = ''.join([mini,str(i)])
                min_wave = np.append(min_wave,mini)
            elif lam_wave.endswith("MAX>\n"):
                maxi = ''
                for i in lam_wave:
                    if i.isdigit() or i == ".":
                        maxi = ''.join([maxi,str(i)])
                max_wave = np.append(max_wave,maxi)
            elif lam_wave.endswith("CENTRAL>\n"):
                midi = ''
                for i in lam_wave:
                    if i.isdigit() or i == ".":
                        midi = ''.join([midi,str(i)])
                mid_wave = np.append(mid_wave,midi)
    logger.debug("Minimum wavelength for some band: %s", min_wave)
    logger.debug("Central wavelength for some band: %s", mid_wave)
    logger.debug("Maximum wavelength for some band: %s", max_wave)
    return min_wave, mid_wave, max_wave

def lat_and_long(global_file):
    """
    Opens the global footprint file to convert the metadata from .xml format
    into useable data.
    """
    with open(global_file, "r",encoding='utf-8') as foot_info:
        for line in foot_info:
            line = line.strip(" ") #First line just contains text.
            line = line.replace("<EXT_POS_LIST>","")
            foot = line.replace("</EXT_POS_LIST>\n","")
        foot = foot.split(" ")
        lat, lon = [], []
        for i in range(len(foot)-1): #Avoid adding last blank space.
            if i%2 == 0:
                lat = np.append(lat, foot[i])
            else:
                lon = np.append(lon, foot[i])
    logger.debug("Latitude: %s %s", lat, type(lat))
    logger.debug("Longitude: %s %s", lon, type(lon))
    return lat, lon

def sun_ang(sun_file):
    """
    Short function that cycles through the small text files which contain 
    the sun zenith and azimuth angles.
    """
    angles = []
    with open(sun_file,"r",encoding='utf-8') as sun_info:
        for line in sun_info:
            if "ANGLE" in line:
                x = line.split('"deg">')
                x = x[1].split('</')
                angles = np.append(angles,float(x[0]))
    logger.debug("Sun zenith and azimuth angles: %s", angles)
    return angles #Returns zenith and azimuth angles respectively

def view_ang(view_file):
    """
    Produce a dictionary with keys that correspond to the S2 bands and items:
    band, zenith angle and azimuthal angle.
    Read in a dictionary with zenith and azimuth viewing angles.
    Returns the pair of values corresponding to the input bandId, in the
    format zenith, azimuth.
    """
    view_arr = []
    with open(view_file, "r",encoding='utf-8') as view_info:
        for lines in view_info:
            view_arr = np.append(view_arr, lines.strip()) #Remove whitespace.
    view_arr = view_arr[1:-1] #Removes lines which did not have desired data.
    means = []
    means_dict = {}
    for term in view_arr:
        if "bandId" in term:
            #means[term[30:-1]]=term[37:-1] #assigns a key with the band id
            #This is only for S2 right now, does not apply to other metadata.
            y=int(term[38:-2]) #stores the band id
            if y < 8:
                y += 1
            elif y == 8: 
                y = ''.join([str(y),"A"])
            means = np.append(means,(''.join(["B",str(y)]).replace('"','')))
        elif "ZENITH" in term:
            #adds zenith viewing angle to that band key
            means = np.append(means,term[25:-15])
        else:
            #adds azimuth viewing angle to that band key
            means= np.append(means,term[26:-16])
        for i in range(0,len(means),4):
            means_dict[means[i].replace('"','')] = means[i:i+3]
    #Iterate through the bands and return sorted list of zeniths and azimuths
    i = 1
    view_zen = [] #Store the view zeniths.
    view_azi = [] #Store the view azimuths.
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
    logger.debug("Viewing zenith angles: %s", view_zen)
    logger.debug("Viewing azimuth angles: %s", view_azi)
    return view_zen, view_azi
#=============================================================================

#=============================================================================
#3: Create an array of top-of-the-atmosphere pixel values.
def create_arr(tif_data, band_no,cols, rows, ):
    """
    Reads in a data file such as the merged tif created in Part 1.
    Uses gdal to convert this to an array for some band within the tif.
    """
    band = tif_data.GetRasterBand(band_no) #Obtains the band information
    data = band.ReadAsArray(0,0, cols,rows) #Sets it as an array.
    return data

def s2_adjust(band_id, solar_corrected, sun_zenith):
    """
    Reads in the array for some S2 band and adjusts the pixel data to return 
    TOA radiance.
    rho = pi*L/E_s * cos(Theta_s)
    """
    scale_factor = 10000 #Scale factor used for the reflectance values.
    band_id = np.array(band_id)
    band_id = (band_id 
        * float(solar_corrected) 
        * np.cos(sun_zenith*(np.pi/180))
        / (np.pi*scale_factor))
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
    newRasterfn = ''.join([name,'.tif'])
    rc.main(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,data)
    logger.debug("Raster file created.")
    return data

#Should add function to handle the creation of bands that don't come from the 
#.zip or .SAFE files, e.g. re-merging band images after correction. Would be 
#useful, but not too sure about how to do it just yet.
#=============================================================================

#=============================================================================
#4.A: The SMAC (Simplified Method for Atmospheric Correction) function:
def smac_func(band_no,r_toa,azimuth_view,azimuth_sun,zenith_view,zenith_sun, 
              coefficients, s28a):
    """
    Declare source of the programme.
    Requires access to the folder "COEFS" which contains the data files for 
    the smac coefficients.
    S2-8A parameter takes 8A into account.
    r_surf=smac_inv(r_toa,theta_s,phi_s,theta_v,phi_v,pressure,AOT,UO3,UH2O,
    coefs)
    For SMAC, Thetas are the zeniths and Phis are the azimuths
    AOT: Aerosol Optical Thickness at 550nm
    UO3: Ozone content
    UH20: Water vapour content in kg/m**2
    """
    if s28a is True:
        nom_smac = ''.join([coefficients,str(band_no),'a.dat'])
    else:
        nom_smac = ''.join([coefficients,str(band_no),'.dat'])
    coefs=smac.coeff(nom_smac)
    r_surf = smac.smac_inv(
        r_toa , zenith_sun, azimuth_sun, 
        float(zenith_view[band_no]),
        float(azimuth_view[band_no]),
        1013,0.1,0.3,0.3, coefs)
    logger.debug("Created surface reflectance using SMAC.")
    return r_surf

#4.B: The 6S function:
#(Second Simulation of a Satellite Signal in the Solar Spectrum Vector Code) 
def sixs_func(
    satellite_latitude, observation_date, month, day,
    solar_zenith, solar_azimuth, view_zenith, view_azimuth,
    target_altitude, band_wavelength, iLUT):
    """
    Do analysis with 6S function.
    """
    s = SixS()
    try:
        if s.test() != 0:
            logger.error("Py6S test failed to return correct response.")
            raise ValueError("6S test not functioning correctly, returned 0.")
    except NameError:
        logger.error("Python 6S failed to initialise.")
        raise NameError()
    #Run 6S simulation defined by SixS object across the whole VNIR range.
    #wavelengths, results = SixSHelpers.Wavelengths.run_vnir(
    #    s, output_name="pixel_radiance")
    logger.debug("AtmosProf params: Lat = %s | %s; ObsDate = %s | %s",
        satellite_latitude, type(satellite_latitude),
        observation_date, type(observation_date))
    s.atmos_profile = AtmosProfile.FromLatitudeAndDate(
        float(satellite_latitude), observation_date)
    logger.debug("6S Atmospheric Profile: %s ",s.atmos_profile)
    #Maritime profile suitable for Ireland
    s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Maritime)
    logger.debug("6S Aerosol Profile: %s ",s.aero_profile)
    #Need to properly choose a GR model.
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.3)
    logger.debug("6S Ground Reflectance: %s ",s.ground_reflectance)
    s.geometry = Geometry.User()
    #Solar and viewing zeniths and azimuths must be floats, while month
    #and day must be ints.
    logger.debug("Solar and Viewing paramater types: %s | %s | %s | %s",
        type(solar_zenith), type(solar_azimuth),
        type(view_zenith), type(view_azimuth))
    logger.debug("Month and day parameters: %s %s | %s %s",
        month, type(month), day, type(day))
    s.geometry.solar_z = float(solar_zenith)
    print(s.geometry.solar_z)
    logger.info("SixS Geometry: Solar Z = %s", s.geometry.solar_z)
    s.geometry.solar_a = float(solar_azimuth)
    s.geometry.view_z = float(view_zenith)
    s.geometry.view_a = float(view_azimuth)
    s.geometry.month = int(month)
    s.geometry.day = int(day)
    logger.debug("6S Geometry: %s ",s.geometry)
    s.altitudes = Altitudes()
    s.altitudes.set_target_custom_altitude(target_altitude)
    s.altitudes.set_sensor_satellite_level()
    logger.debug("6S Altitude: %s ",s.altitudes)
    logger.debug("6S Wavelength: %s", band_wavelength)
    print(band_wavelength, type(band_wavelength))
    s.wavelength = Wavelength(getattr(PredefinedWavelengths,band_wavelength))
    #Similarly to ground reflectance, can improve the choice here.
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(0.23)
    try:
        logger.debug("Debug report")
        s.produce_debug_report()
        logger.debug(s.produce_debug_report())
    except TypeError:
        logger.warning("Failed to produce debug report")
            
    s.run()
    file_out = ''.join(["6S_outputs_",str(band_wavelength),".txt"])
    s.outputs.write_output_file(file_out)
    sp.call(["cat",file_out])
    """
    An interpolated look-up table requires the following input variables 
    (in order) to provide atmospheric correction coefficients:
        solar zentith [degrees] (0 - 75)
        water vapour [g/m2] (0 - 8.5)
        ozone [cm-atm] (0 - 0.8)
        aerosol optical thickness [unitless] (0 - 3)
        surface altitude [km] (0 - 7.75)
    """
    #Should replace the hard values with a function which calls the correct
    #value from the 
    water_vapour = 1.349 #g/cm2
    ozone = 0.343 #cm-atm
    aerosol_optical_thickness = 0.5
    surface_altitude = target_altitude
    a, b = iLUT(
        solar_zenith, water_vapour, 
        ozone, aerosol_optical_thickness, 
        surface_altitude)
    print(a,b)
    logger.info("Corrrection coefficients for Py6S: %s %s", a, b)
    return a, b

def sixs_surf(cor_a, cor_b, sen_rad):
    """
    Reads in the correction coefficients, the radiance at the sensor and
    returns a value for the surface reflectance.
    """
    ρ = (sen_rad - cor_a) / cor_b
    return ρ
#=============================================================================

#=============================================================================
#5: Combine the functions for a sample output when this script's called alone.
def main():
    inputPath, outputPath = open_s2() #Opens a Sentinel-2 product.
    obs_date = (inputPath.split("_"))[2].split("T")[0] #Get observation date.
    selected_metadata(inputPath, outputPath) #Create metadata files.
    #Correct solar irradiance for reflectance at measurement Earth's position:
    sol_irr_at_sat, correction_coefficient = sol_irr(
        "Reflectance_Conversion_meta_1.txt")
    sol_cor = []
    #sol_cor is stored in the form of [BAND_ID, VALUE, BAND_ID, VALUE,...]
    for i in sol_irr_at_sat:
        sol_cor = np.append(sol_cor,i.split("|")[0])
        sol_cor = np.append(
            sol_cor,(float(i.split("|")[1]))*correction_coefficient) 
    #Create a function to collect wavelengths
    #band_lambdas = wave_fun("Wavelength_meta_1.txt")
    sat_lat, sat_lon = lat_and_long("Global_Footprint_meta_1.txt")
    #Produce a tuple for the average solar zenith and azimuth angles.
    sun_zen, sun_azi = sun_ang("Mean_Sun_meta_1.txt")
    #Produce a tuple for the average viewing zenith and azimuth angles.
    view_zen, view_azi = view_ang(
        "Mean_Viewing_Incidence_Angle_List_meta_1.txt")
    #Correct each band from reflectance to radiance.
    data_file = ''.join(
        [outputPath,"/",inputPath[:60],"_PROCESSED/merged.tif"])
    dataset = gdal.Open(data_file, GA_ReadOnly) #Make gdal osgeo class object.
    cols = dataset.RasterXSize #The number of columns; 5400 for S2 at 20m.
    rows = dataset.RasterYSize #The number of rows; 5400 for S2 at 20m.
    bands = dataset.RasterCount #the number of bands; 13 for S2 product.
    geotransform = dataset.GetGeoTransform() #Geoinformation for product.
    originX = geotransform[0]
    originY = geotransform[3]
    rasterOrigin = (originX, originY)
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    #Create arrays to append bands to.
    radiance_bands = ["gdal_merge.py", "-o", "s2_radiance.tif","-separate"] 
    smac_bands = ["gdal_merge.py", "-o", "smac.tif","-separate"]
    sixs_bands = ["gdal_merge.py", "-o", "smac.tif","-separate"]
    #The different band wavelengths which are necessary for 6S:
    sixs_s2_wavelengths = """
    S2A_MSI_01
    S2A_MSI_02
    S2A_MSI_03
    S2A_MSI_04
    S2A_MSI_05
    S2A_MSI_06
    S2A_MSI_07
    S2A_MSI_08
    S2A_MSI_09
    S2A_MSI_10
    S2A_MSI_11
    S2A_MSI_12
    S2A_MSI_13
    """
    sixs_s2_wavelengths = sixs_s2_wavelengths.strip("\n").strip(" ") \
        .replace("\n","").split("    ")
    #Ground altitude, would be useful to tie in equivalent S1 data.
    gnd_alt = 0.2 #Surface height above sea level (km)
    #Base filename for where the iLUTs are stored.
    #Assumes that working in fixed directory, needs to be adapted.
    base = "6S_emulator/files/LUTs/S2A_MSI/Continental/view_zenith_0/files/\
        iLUTs/S2A_MSI/Continental/view_zenith_0/"
    base = base.replace("        ",'')
    logger.debug(base)
    #Iterate through the 13 Sentinel bands, creating tiffs.
    for i in range(bands):
        #Initialise 6S
        if i+1 < 10:                       
            fpath = ''.join([base,"S2A_MSI_0",str(i+1),".ilut"])
        else:
            fpath = ''.join([base,"S2A_MSI_",str(i+1),".ilut"])
        logger.debug("Filepath: %s", fpath)
        with open(fpath, "rb") as ilut_file:
            iLUT = pickle.load(ilut_file)
        a, b = sixs_func(sat_lat[i], obs_date, obs_date[4:-2], obs_date[-2:],
            sun_zen, sun_azi, view_zen[i], view_azi[i],
            gnd_alt, sixs_s2_wavelengths[i], iLUT)
        band = create_arr(dataset, i+1 ,cols, rows)
        logger.debug("Created array for band %s", str(i+1))
        band_rad = s2_adjust(band, sol_cor[i*2 + 1], sun_zen)
        logger.debug("Changed reflectance to radiance for band %s", str(i+1))
        sixs_boa = sixs_surf(a, b, band_rad)
        logger.debug("Determiend BOA reflectance using 6S parameters.")
        if i < 8:
            name = ''.join(["B",str(i+1)]) #Band names as in original MTD data.
            band_gs = create_tiff(
                band_rad, name,rasterOrigin,pixelWidth,pixelHeight)
            band_6s = create_tiff(
                sixs_boa, ''.join(["6S", name]), rasterOrigin, pixelWidth,
                pixelHeight)
            band_smac = smac_func(
                i+1,np.array(band),view_azi,sun_azi,
                view_zen, sun_zen,'COEFS/Coef_S2A_CONT_B', False)
        elif i == 8:
            name = ''.join(["B",str(i),"A"]) #Band name in original MTD data.
            band_gs = create_tiff(
                band_rad, name,rasterOrigin,pixelWidth,pixelHeight)
            band_6s = create_tiff(
                sixs_boa, ''.join(["6S", name]), rasterOrigin, pixelWidth,
                pixelHeight)
            band_smac = smac_func(
                i,np.array(band),view_azi,sun_azi,
                view_zen, sun_zen,'COEFS/Coef_S2A_CONT_B', True)
        else:
            name = ''.join(["B",str(i)]) #Band names as in original MTD data.
            band_gs = create_tiff(
                band_rad, name,rasterOrigin,pixelWidth,pixelHeight)
            band_6s = create_tiff(
                sixs_boa, ''.join(["6S", name]), rasterOrigin, pixelWidth,
                pixelHeight)
            band_smac = smac_func(
                i,np.array(band),view_azi,sun_azi,
                view_zen, sun_zen,'COEFS/Coef_S2A_CONT_B', False)
        radiance_bands.append(''.join([name,".tif"]))
        name_sm = ''.join(["SMAC_",name])
        band_smac_gs1 = create_tiff(
            band_smac, name_sm,rasterOrigin,pixelWidth,pixelHeight)
        print("%.f" % (i*100/13), end="..", flush=True)
        smac_bands.append(''.join([name_sm,".tif"]))
        sixs_bands.append(''.join(["6S",name,".tif"]))

    print("100 - done.")
    sp.call(radiance_bands)
    sp.call(smac_bands)
    #sp.call(["gdalwarp","-s_srs","'EPSG:32629'","-t_srs","'EPSG:32629'","s2_radiance.tif","out.tif"])
    #sp.call(["gdalwarp","-s_srs","'EPSG:32629'","-t_srs","'EPSG:32629'","smac","smac_out.tif"])
    #Errors when gdalwarp called as above, but the equivalent 
    #gdalwarp -s_srs 'EPSG:32629' -t_srs 'EPSG:32629' s2_radiance.tif out.tif
    #is fine when input directly to terminal. 

#=============================================================================

#=============================================================================
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    logger.info("Run time was %.3f seconds", end_time-start_time)
    print("\n == Completed after %.3f seconds ==" % (end_time - start_time))
#=============================================================================

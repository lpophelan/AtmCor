#!/usr/bin/env python
# encoding=utf-8
"""
This python script contains some functions which can read in a file that is 
created by the Py6S outputs class. It currently finds some 
"""
import logging
import time

logging.basicConfig(filename="S2A-AtmCor.log",format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# create logger
logger = logging.getLogger('6S Output Reader')
logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)

#=============================================================================
def isempty(string):
    """
    Reads in a string and determines whether it contains any alphanumeric 
    characters or whether it does not contain informative text, for example 
    purely whitespace.
    """
    logger.debug(string)
    for char in string:
        if char.isalpha() or char.isnumeric():
            return False
    return True
#=============================================================================

#=============================================================================
def sixsread(filename):
    """
    Opens a provided file and uses the isempty function to extract information 
    from it.
    """
    with open(filename, "r", encoding="utf-8") as f:
        content = []    
        for line in f:
            this_line = line.replace("*","").strip(" ") #Prepares each line.
            logger.debug(isempty(this_line))
            if not isempty(this_line):
                content.append(this_line)
        logger.info("Content list created from %s", filename)
        for i in range(len(content)):
            content[i] = content[i].replace("\n","")
            logger.debug(content[i])
    return content
#=============================================================================

#=============================================================================
def wv_and_ozone(data):
    """    
    Fetches the values for water vapour content and ozone which can be 
    retrieved from the content list created by the sixsread function.
    """    
    for i in range(len(data)):
        if data[i].startswith("uh2o="):
            logger.debug(data[i])
            #Water vapour and ozone should match below format
            wv = data[i].split("        ")[0].split("=")[1].split(" ")[1]
            logger.info("Water vapour content: %s", wv)
            o3 = data[i].split("        ")[1].split("=")[1].split(" ")[1]
            logger.info("Ozone content: %s", o3)
            return wv, o3
    logger.warning("Water vapour and ozone values not found")
    #raise NameError
    return None, None
#=============================================================================

#=============================================================================
def main():
    #An example use of the script and its functions.
    outputs = sixsread("6S_outputs_S2A_MSI_02.txt")
    water_vapour, ozone = wv_and_ozone(outputs)
    logger.info("Water vapour (g/cm2) and ozone (cm-atm) values retrieved.")
#=============================================================================

#=============================================================================
if __name__ == "__main__":
    start_time = time.time()
    st_gmt = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))
    logger.info("=====LOG START: %s ===== ", st_gmt)
    main()
    end_time = time.time()
    logger.info("Completed after %.3f seconds", end_time-start_time)
    print("\n == Completed after %.3f seconds ==" % (end_time - start_time))
#=============================================================================


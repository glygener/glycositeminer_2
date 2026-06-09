import sys,os
import json
import glob
import subprocess
from optparse import OptionParser


def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())
   
 
    source = "pubtator_downloads"
    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + source) == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + source
        x = subprocess.getoutput(cmd)


    base_url = "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/"

    file_name = "bioconcepts2pubtatorcentral.offset.gz"
    out_file = config_obj["data_dir"] + "/" + source + "/" + file_name
    #cmd = "wget %s/%s -O %s" % (base_url, file_name, out_file)
    cmd = "curl %s/%s -o %s" % (base_url, file_name, out_file)
    #print (cmd)
    x = subprocess.getoutput(cmd)

    file_name = "bioconcepts2pubtatorcentral.gz"
    out_file = config_obj["data_dir"] + "/" + source + "/" + file_name
    #cmd = "wget %s/%s -O %s" % (base_url, file_name, out_file)
    cmd = "curl %s/%s -o %s" % (base_url, file_name, out_file)
    #print (cmd)
    x = subprocess.getoutput(cmd)

    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/" + source
    x = subprocess.getoutput(cmd)

    return



if __name__ == '__main__':
    main()



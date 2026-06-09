import sys,os
import json
import glob
import subprocess
from optparse import OptionParser


def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())
    if config_obj["mode"] == "dev":
        config_obj["data_dir"], config_obj["misc_dir"] = "reldir/", "miscdir/"

   
 
    source = "pubtator_downloads"
    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + source) == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + source
        x = subprocess.getoutput(cmd)


    #bioconcepts2pubtatorcentral.offset is an archive with titles/abstracts and the annotations in full texts. 
    out_file = config_obj["data_dir"] + "/" + source + "/bioconcepts2pubtatorcentral.offset.gz"
    cmd = "wget https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/bioconcepts2pubtatorcentral.offset.gz "
    cmd += "-O %s" % (out_file)           
    x = subprocess.getoutput(cmd)
    #print (cmd)

    #bioconcepts2pubtatorcental is a combination of all entity annotations in PubTatorCentral [1]. 
    out_file = config_obj["data_dir"] + "/" + source + "/bioconcepts2pubtatorcentral.gz"
    cmd = "wget https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/bioconcepts2pubtatorcentral.gz "
    cmd += "-O %s" % (out_file)
    x = subprocess.getoutput(cmd)
    #print (cmd)


    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/" + source
    x = subprocess.getoutput(cmd)



if __name__ == '__main__':
    main()



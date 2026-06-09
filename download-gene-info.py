import sys,os
import json
import glob
import subprocess
from optparse import OptionParser


def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())

    source = "gene_info"
    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + source) == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + source
        x = subprocess.getoutput(cmd)

    out_file = config_obj["data_dir"] + "glygen/species_info.csv"
    tax_id_dict = {}
    with open (out_file, "r") as FR:
        for line in FR:
            row = line.strip().split(",")
            if row[-2] == "yes":
                tax_id_dict[row[0].strip()] = True

    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz"
    
    out_file = config_obj["data_dir"] + "gene_info/tmp.gene_info.gz"
    cmd = "curl %s -o %s " % (url, out_file)
    x = subprocess.getoutput(cmd)
    cmd = "gunzip " + out_file
    x = subprocess.getoutput(cmd)


    in_file = config_obj["data_dir"] + "gene_info/tmp.gene_info"
    out_file = config_obj["data_dir"] + "gene_info/All_Data.gene_info"
    
    FW = open(out_file, "w") 
    with open (in_file, "r") as FR:
        for line in FR:
            tax_id = line.split("\t")[0].strip()
            if tax_id in tax_id_dict:
                FW.write(line)
    FW.close()

    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/" + source
    x = subprocess.getoutput(cmd)
   




if __name__ == '__main__':
    main()



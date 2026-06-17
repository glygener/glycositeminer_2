import os,sys
import string
import glob
import json
import datetime
import gzip
import subprocess

__version__="1.0"
__status__ = "Dev"




###############################
def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
 
    gz_file = config_obj["data_dir"] + "pubtator_downloads/bioconcepts2pubtatorcentral.offset.gz"
    medline_abstracts_dir = config_obj["data_dir"] + "medline_abstracts/"
    pubtator_extracts_dir = config_obj["data_dir"] + "pubtator_extracts/"
    log_file = config_obj["data_dir"] + "logs/pubtator_extracts.log"

    if os.path.isdir(pubtator_extracts_dir) == False:
        cmd = "mkdir -p " + pubtator_extracts_dir
        x = subprocess.getoutput(cmd)
  
 

    pmid_dict = {}
    for in_file in glob.glob(medline_abstracts_dir + "pmid.*.txt"):
        doc_id = in_file.split("/")[-1].split(".")[-2]
        pmid_dict[doc_id] = True

    total_common_count = len(list(pmid_dict.keys()))
    with open(log_file, "w") as FL:
        FL.write("")



    common_count, pubtator_count = 0, 0
    part = 0
    prev_pmid = ""
    buf = ""
    with gzip.open(gz_file,'r') as fin:        
        for line in fin:        
            line = line.decode()
            if line.strip() == "":
                continue
            pmid = line.split("|")[0].split()[0]
            if prev_pmid != "" and pmid != prev_pmid:
                pubtator_count += 1 
                if pmid in pmid_dict:
                    out_file = pubtator_extracts_dir + "pmid.%s.txt" % (prev_pmid)
                    with open(out_file, "w") as FW:
                        FW.write("%s\n" % (buf))
                    buf = "" 
                    common_count += 1
                buf += line
            if pmid != prev_pmid and pubtator_count%10000 == 0:
                with open(log_file, "a") as FL:
                    FL.write("parsed %s putator PMIDs, found %s/%s of medline extracts\n" % (pubtator_count,common_count,total_common_count))
            prev_pmid = pmid

    with open(log_file, "a") as FL:
        FL.write("finished parsing %s putator PMIDs, found %s/%s of medline extracts\n" % (pubtator_count,common_count,total_common_count))  
    
    if prev_pmid != "":
        out_file = pubtator_extracts_dir + "pmid.%s.txt" % (prev_pmid)
        with open(out_file, "w") as FW: 
            FW.write("%s\n" % (buf))
    
    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/pubtator_extracts" 
    x = subprocess.getoutput(cmd)
    cmd = "chmod -R 777 logs/*"
    x = subprocess.getoutput(cmd)

    return


            

if __name__ == '__main__':
	main()


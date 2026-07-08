import os,sys
import string
import glob
import json
import datetime
import gzip
import subprocess
from util import get_otext2taxid_reduced_dict

__version__="1.0"
__status__ = "Dev"




###############################
def main():


    config_obj = json.loads(open("conf/config.json", "r").read())


    misc_dir = config_obj["data_dir"] + "misc/"
    ent_dir = config_obj["data_dir"] + "llm_entities/*/"
    otext2taxid_dict = get_otext2taxid_reduced_dict(misc_dir, ent_dir)

    seen_cmb = {}
    file_list = glob.glob(config_obj["data_dir"] + "llm_entities/*/*.json")
    for in_file in file_list:
        doc_id = in_file.split("/")[-1].split(".")[1]
        doc = json.load(open(in_file))
        for obj in doc["glycosylation_sites"]:
            if "organism" not in obj:
                continue
            if obj["organism"] in ["", "none", None]:
                continue
            otext = obj["organism"].lower()
            if otext not in otext2taxid_dict:
                continue
            for tax_id in otext2taxid_dict[otext]:
                sp = otext2taxid_dict[otext][tax_id]
                cmb = "%s|%s" % (doc_id, sp)
                seen_cmb[cmb] = True

    out_file = config_obj["data_dir"] + "predicted/predicted_corrected.csv"
    FW = open(out_file, "w")
    in_file = config_obj["data_dir"] + "predicted/predicted.csv"
    with open(in_file, "r") as FR:
        for line in FR:
            row = line[1:-2].strip().split("\",\"")
            doc_id, sp = row[0], row[4]
            if doc_id == "evidence":
                newrow = row + ["species_flag"]
                FW.write("\"%s\"\n" % ("\",\"".join(newrow)))
            else:
                cmb = "%s|%s" % (doc_id, sp)
                flag = "species_match" if cmb in seen_cmb else "species_mismatch"
                newrow = row + [flag]
                FW.write("\"%s\"\n" % ("\",\"".join(newrow)))
    FW.close()
 

    return


            

if __name__ == '__main__':
	main()


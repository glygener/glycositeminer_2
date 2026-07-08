import os
import json
from sklearn import svm
from sklearn import datasets
import pickle
import util
import subprocess
import glob

def load_verification_dict(data_dir):

    verification_dict = {}
    file_list = glob.glob(data_dir + "/confirmation/pmid.*.json")
    for in_file in file_list:
        doc_id = in_file.split(".")[-2]
        doc = json.load(open(in_file))
        for obj in doc:
            site = "%s|%s|%s" % (doc_id, obj["canon"], obj["pos"])
            if site not in verification_dict:
                verification_dict[site] = "no"
            if obj["answer"] == "yes":
                verification_dict[site] = "yes"

    return verification_dict




def main():

    config_obj = json.loads(open("conf/config.json", "r").read())

    verification_dict = load_verification_dict(config_obj["data_dir"])

    in_file = config_obj["data_dir"] + "predicted/predicted_tmp.csv"
    out_file = config_obj["data_dir"] + "predicted/predicted.csv"
    FW = open(out_file, "w")
    row = ["evidence", "uniprotkb_ac", "glycosylation_site", "amino_acid", "tax_name", "glygen_status", "algorithm", "llm_verification_flag", "curation_flag"]
    FW.write("\"%s\"\n" % ("\",\"".join(row)))
    
    line_list = open(in_file, "r").read().split("\n")[0:-1]
    idx = 0
    f_list = []
    for line in line_list:
        idx += 1
        row = line.replace("\"", "").split(",")
        if idx == 1:
            f_list = row
        else:
            site = "%s|%s|%s" % (row[0], row[1], row[2])
            flag = "llm_" + verification_dict[site] if site in verification_dict else "llm_unknown" 
            row[f_list.index("llm_verification_flag")] = flag
            FW.write("\"%s\"\n" % ("\",\"".join(row))) 
    FW.close()




    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/predicted/"
    x = subprocess.getoutput(cmd)



    return


if __name__ == '__main__':
    main()





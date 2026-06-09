import sys,os
import json
import glob
import subprocess
from optparse import OptionParser


def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())
    if config_obj["mode"] == "dev":
        config_obj["rel_dir"], config_obj["misc_dir"] = "reldir/", "miscdir/"

    source = "glygen"
    # create dir if doesn't exist    
    if os.path.isdir(config_obj["rel_dir"] + source) == False:
        cmd = "mkdir -p " + config_obj["rel_dir"] + source
        x = subprocess.getoutput(cmd)

    ds_list = [
        "canonicalsequences.fasta",
      "genelocus.csv",
      "genenames_refseq.csv",
      "genenames_uniprotkb.csv",
      "xref_geneid.csv",
    ]

    #tax_id,short_name,long_name,common_name,glygen_name,nt_file,is_reference,sort_order


    glygen_dir = "/data/projects/glygen/generated/"

    out_file = config_obj["rel_dir"] + "glygen/species_info.csv"
    cmd = "cp /data/projects/glygen/generated/misc/species_info.csv %s" % (out_file)
    x = subprocess.getoutput(cmd)
    sp_list = []
    with open (out_file, "r") as FR:
        for line in FR:
            row = line.strip().split(",")
            if row[-3] == "yes":
                sp_list.append(row[1])


    for ds in ds_list:
        for sp in sp_list:
            file_name = "%s_protein_%s" % (sp, ds)
            out_file = config_obj["rel_dir"] + "glygen/%s" % (file_name) 
            cmd = "cp /data/projects/glygen/generated/datasets/unreviewed/%s reldir/glygen/" % (file_name)
            x = subprocess.getoutput(cmd)
            #print (cmd)
    

    cmd = "chmod -R 777 " + config_obj["rel_dir"] + "/" + source
    x = subprocess.getoutput(cmd)


if __name__ == '__main__':
    main()



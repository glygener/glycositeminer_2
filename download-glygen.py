import sys,os
import json
import glob
import subprocess
from optparse import OptionParser


def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())

    source = "glygen"
    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + source) == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + source
        x = subprocess.getoutput(cmd)

    ds_list = [
        "canonicalsequences.fasta",
      "genelocus.csv",
      "genenames_refseq.csv",
      "genenames_uniprotkb.csv",
      "xref_geneid.csv",
    ]

    #tax_id,short_name,long_name,common_name,glygen_name,nt_file,is_reference,sort_order


    glygen_rel = config_obj["glygen_rel"] 
    out_file = config_obj["data_dir"] + "glygen/species_info.csv"
    ftp_url = "https://data.glygen.org/ln2releases/data/v-%s/misc/species_info.csv" % (glygen_rel)
    cmd = "curl %s -o %s" % (ftp_url, out_file)
    x = subprocess.getoutput(cmd)
    sp_list = []
    with open (out_file, "r") as FR:
        for line in FR:
            row = line.strip().split(",")
            if row[-3] == "yes":
                sp_list.append(row[1])


    ftp_url = "https://data.glygen.org/ln2releases/data/v-%s/reviewed/" % (glygen_rel)
    

    for src in config_obj["glygen_src_list"]:
        for sp in sp_list:
            file_name = "%s_proteoform_glycosylation_sites_%s.csv" % (sp, src)
            out_file = config_obj["data_dir"] + "glygen/%s" % (file_name)
            cmd = "curl %s%s -o %s" % (ftp_url, file_name, out_file)
            x = subprocess.getoutput(cmd)
            cmd = "grep uniprotkb_canonical_ac %s |wc " % (out_file)
            x = subprocess.getoutput(cmd)
            parts = x.strip().split(" ")
            if x.strip().split(" ")[0] == "0":
                cmd = "rm -f %s" % (out_file)
                x = subprocess.getoutput(cmd)
 
    for ds in ds_list:
        for sp in sp_list:
            file_name = "%s_protein_%s" % (sp, ds)
            out_file = config_obj["data_dir"] + "glygen/%s" % (file_name) 
            cmd = "curl %s%s -o %s" % (ftp_url, file_name, out_file)
            x = subprocess.getoutput(cmd)
            cmd = "curl %s%s -o %s" % (ftp_url, file_name, out_file)


 
    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/" + source
    x = subprocess.getoutput(cmd)


if __name__ == '__main__':
    main()



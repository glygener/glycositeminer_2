import sys,os
import json
import glob
import subprocess
from optparse import OptionParser


def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())
    if config_obj["mode"] == "dev":
        config_obj["rel_dir"] = "/data/shared/glyds/releases/ftp/glycositeminer/"

    # create dir if doesn't exist    
    downloads_dir = config_obj["rel_dir"] + "downloads/"
    if os.path.isdir(downloads_dir) == False:
        cmd = "mkdir -p " + downloads_dir
        x = subprocess.getoutput(cmd)

    tarball_list = [ 
        "medline_abstracts.tar.gz",
        "medline_extracts.tar.gz",
        "pubtator_extracts.tar.gz",
        "llm_entities.tar.gz",
        "gene_info.tar.gz",
        "glygen.tar.gz", 
        "curation.tar.gz",
        "confirmation.tar.gz",
        "misc.tar.gz"
    ]
    for tarball in tarball_list:
        out_file = downloads_dir + tarball
        cmd = "wget https://data.glygen.org/ftp/glycositeminer/tarballs/%s -O %s"  % (tarball, out_file)
        cmd += " --no-check-certificate"
        x = subprocess.getoutput(cmd)
        cmd = "tar xfz %s -C %s --touch" % (out_file, config_obj["rel_dir"])
        x = subprocess.getoutput(cmd)
        cmd = "chmod -R 777 " + config_obj["rel_dir"] + "/" + tarball.split(".")[0]
        x = subprocess.getoutput(cmd)


    cmd = "chmod -R 777 " + downloads_dir
    x = subprocess.getoutput(cmd)
   


if __name__ == '__main__':
    main()



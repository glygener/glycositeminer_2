import sys,os
import json
import glob
import subprocess
from optparse import OptionParser


def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())

    source = "confirmation"
    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + source) == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + source
        x = subprocess.getoutput(cmd)

    out_file = config_obj["data_dir"] + "confirmation/confirmation.tar.gz"
    ftp_url = "https://data.glygen.org/ln2downloads/tmp/confirmation.tar.gz"
    cmd = "curl %s -o %s" % (ftp_url, out_file)
    x = subprocess.getoutput(cmd)

    cmd = "tar xvfz %s -C %s" % (out_file, config_obj["data_dir"])
    x = subprocess.getoutput(cmd)
    
    
    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/" + source
    x = subprocess.getoutput(cmd)

    return



if __name__ == '__main__':
    main()



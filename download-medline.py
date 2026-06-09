import sys,os
import json
import glob
import subprocess


def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())

    xml_folder = config_obj["data_dir"] + "medline_xml"
    if os.path.isdir(xml_folder) == False:
        cmd = "mkdir -p " + xml_folder
        x = subprocess.getoutput(cmd)

    # Make sure to update "conf/medline.json"
    doc = json.loads(open("conf/medline.json", "r").read())

    zeros = "0000000"
    cat_list = ["baseline", "updatefiles"]
    for category in cat_list:
        start = doc[category]["start"]
        end = doc[category]["end"]
        year = doc[category]["year"]
        ftp_url = doc[category]["url"]
        for idx in range(start, end + 1):
            padlen = 4 - len(str(idx))
            idx = zeros[:padlen] + str(idx)
            file_name = "pubmed%sn%s.xml.gz" % (year, idx)
            out_file_one = xml_folder + "/" + file_name
            #cmd = "wget %s%s -O %s" % (ftp_url, file_name, out_file_one)
            cmd = "curl %s%s -o %s" % (ftp_url, file_name, out_file_one)
            x = subprocess.getoutput(cmd)

    cmd = "chmod -R 777 " + xml_folder
    x = subprocess.getoutput(cmd)

    return


if __name__ == '__main__':
    main()



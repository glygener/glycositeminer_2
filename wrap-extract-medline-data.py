import os
import sys
import json
import glob


def main():

    DEBUG = False
    #DEBUG = True


    config_obj = json.loads(open("conf/config.json", "r").read())
 
    pattern = "*"
    if DEBUG:
        pattern = "pubmed23n000*"
    xml_dir = config_obj["data_dir"] + "medline_xml/"
    xml_file_list = sorted(glob.glob(xml_dir + "/%s.xml.gz" % (pattern)))
    n = len(xml_file_list)
    batch_size = int(len(xml_file_list)/10) if len(xml_file_list) > 10 else len(xml_file_list)
    range_list = [] 
    for i in range(0, 12):
        s = i*batch_size + 1
        e = s + batch_size
        if e >= n:
            e = n
            range_list.append({"s":s, "e":e})
            break
        range_list.append({"s":s, "e":e}) 

    for o in range_list:
        cmd = "nohup python3 extract-medline-data.py -s %s -e %s  & " % (o["s"], o["e"])
        os.system(cmd)
        #print (cmd)

    return




if __name__ == '__main__':
    main()



import sys
import string
import json
import subprocess
import glob



def main():

   
    config_obj = json.loads(open("conf/config.json", "r").read())
    
    glygen_dir = config_obj["data_dir"] + "glygen/"
    file_list = glob.glob(glygen_dir+"*_proteoform_glycosylation_sites_*.csv")
    if file_list == []:
        print ("\n\tNo *_proteoform_glycosylation_sites_*.csv found under %s" % (glygen_dir))
        exit()


    out_file = glygen_dir + "known_sites.csv"
    FW = open(out_file, "w") 
    newrow = ["doc_id","canon","pos","amino_acid","source"]
    FW.write("\"%s\"\n" % ("\",\"".join(newrow)))
 
    row_list = []
    seen = {}        
    for in_file in file_list:
        file_name = in_file.split("/")[-1]
        lcount = 0
        f_list = []
        with open(in_file, "r") as FR:
            for line in FR:
                lcount += 1
                if lcount == 1:
                    f_list = line.strip().replace("\"", "").split(",")
                else:
                    row = line.strip().split("\",\"")
                    canon = row[f_list.index("uniprotkb_canonical_ac")].replace("\"", "")
                    xref_key = row[f_list.index("xref_key")].replace("\"", "")
                    xref_id = row[f_list.index("xref_id")].replace("\"", "")
                    pos = row[f_list.index("glycosylation_site_uniprotkb")].replace("\"", "")
                    start_pos = row[f_list.index("start_pos")].replace("\"", "")
                    end_pos = row[f_list.index("end_pos")].replace("\"", "")
                    aa = row[f_list.index("amino_acid")].replace("\"", "")
                    mining_tools = row[f_list.index("mining_tools")] if "mining_tools" in f_list else ""
                    if mining_tools.lower() == "glycositeminer":
                        continue
                    if xref_key.find("xref_pubmed") == -1:
                        continue
                    if pos.strip() == "":
                        continue
                    if int(start_pos) != int(end_pos):
                        continue
                    newrow = [xref_id,canon,pos,aa,file_name]
                    s = json.dumps(newrow)
                    if s not in seen:
                        FW.write("\"%s\"\n" % ("\",\"".join(newrow)))
                        seen[s] = True
    FW.close()


    
    cmd = "chmod -R 777 %s" % (glygen_dir)
    x = subprocess.getoutput(cmd)


    return



if __name__ == '__main__':
    main()



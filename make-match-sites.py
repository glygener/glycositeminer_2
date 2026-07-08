import sys,os
import json
import util
import glob
import subprocess





def main():


    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())

    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + "match_sites/") == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + "match_sites/"
        x = subprocess.getoutput(cmd)


    known_site_dict = util.get_known_site_dict(config_obj["data_dir"])    
    species_dict = util.get_species_dict(config_obj["misc_dir"])
    seq_dict, canon2taxid = util.load_seq_dict(config_obj["data_dir"], config_obj["misc_dir"] )

    mismatch_dict = {}
    cls_dict = {}         
    tmp_dict = {}
    file_list = glob.glob(config_obj["data_dir"] + "integrated/*.json")
    for in_file in file_list:
        doc = json.load(open(in_file))
        for site in doc:
            for obj in doc[site]:
                doc_id, canon = obj["docid"], obj["canon"]
                gene_text = obj["gene_text"]
                species_list = obj["specieslist"]
                status = obj["status"]
                pos = site.split("|")[1]
                cmb = "%s|%s|%s" % (doc_id, canon, pos)
                if status == "aa_mismatch":
                    mismatch_dict[cmb] = True
                elif status == "aa_match":
                    if cmb not in tmp_dict:
                        tmp_dict[cmb] = {}
                    tmp_dict[cmb][gene_text] = True
                    cls_dict[cmb] = "x"
                    seq_cmb = "%s|%s" % (doc_id, canon)
                    if cmb in known_site_dict:
                        cls_dict[cmb] = "1"
                    tax_id = canon2taxid[canon]
                    tax_name = species_dict[tax_id]


    taken_dict = {}
    for cmb in cls_dict:
        doc_id,canon, pos = cmb.split("|")
        pos_cmb = "%s|%s" % (doc_id, pos)
        if cls_dict[cmb] == "1":
            taken_dict[pos_cmb] = True
    for cmb in cls_dict:
        doc_id,canon, pos = cmb.split("|")
        pos_cmb = "%s|%s" % (doc_id, pos)
        if cls_dict[cmb] == "x" and pos_cmb in taken_dict:
            cls_dict[cmb] = "0"

    out_file = config_obj["data_dir"] + "/match_sites/mismatch_sites.csv"
    FW = open(out_file, "w")
    FW.write("%s\n" % ("site"))        
    for cmb in mismatch_dict:
        FW.write("%s\n" % (cmb))
    FW.close()

    out_file = config_obj["data_dir"] + "/match_sites/match_sites.csv"
    FW = open(out_file, "w")
    FW.write("%s,%s\n" % ("site","cls"))
    for cmb in tmp_dict:
        doc_id,canon, pos = cmb.split("|")
        pos_cmb = "%s|%s" % (doc_id, pos)
        cls = cls_dict[cmb]
        FW.write("%s,%s\n" % (cmb, cls))
    FW.close()


    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/match_sites"
    x = subprocess.getoutput(cmd)
   
    return
 



if __name__ == '__main__':
    main()



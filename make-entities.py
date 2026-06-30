import sys,os
import json
import util
import glob
import entities
import subprocess

def get_doc_id_list():
    tmp_list = []
    #file_list = glob.glob(config_obj["data_dir"] + "/medline_extracts/pmid.*.json") 
    file_list = glob.glob(config_obj["data_dir"] + "/llm_entities/1/llm.*.json") 
    for in_file in file_list:
        doc_id = in_file.split(".")[-2]
        tmp_list.append(doc_id)
    return tmp_list




def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())


    DEBUG = False
    #DEBUG = True

    doc_id_list = get_doc_id_list()
    if DEBUG:
        doc_id_list = ["9931318"]

    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + "entities/") == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + "entities/"
        x = subprocess.getoutput(cmd)


    aa_dict = util.get_aa_dict(config_obj["misc_dir"])
    genename_dict = util.load_genename_dict(config_obj["data_dir"])


    seen_ent_file = {}
    known_aa = util.get_known_aa(config_obj["data_dir"], config_obj["misc_dir"])
    for doc_id in doc_id_list:
        extract_file = config_obj["data_dir"] + "medline_extracts/pmid.%s.json" % (doc_id)
        ent_file = config_obj["data_dir"] + "entities/site.%s.json" % (doc_id)
        entities.dump_glycositeminer_entities(doc_id, config_obj["data_dir"],genename_dict, aa_dict, known_aa)
        if os.path.isfile(ent_file):
            entities.dump_pubtator_entities(doc_id, config_obj["data_dir"])
            seen_ent_file[ent_file] = True


    llm_site_dict = util.load_llm_sites(config_obj["data_dir"], config_obj["misc_dir"])
    for in_file in seen_ent_file:
        doc_id = in_file.split(".")[-2]
        doc = json.load(open(in_file, "r"))
        for sent_idx in doc:
            for obj in doc[sent_idx]:
                pos = obj["text"].split("-")[-1].strip()
                cmb = "%s|%s" % (doc_id, pos)
                obj["llm_flag"] = cmb in llm_site_dict
        with open(in_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))


    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/entities/"
    x = subprocess.getoutput(cmd)
   

    return
 

if __name__ == '__main__':
    main()



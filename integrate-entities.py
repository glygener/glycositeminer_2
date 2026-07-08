import sys,os
import json
import util
import glob
import entities
import subprocess






def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())

    DEBUG = False
    #DEBUG = True

    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + "integrated/") == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + "integrated/"
        x = subprocess.getoutput(cmd)



    doc_id_list = []
    file_list = glob.glob(config_obj["data_dir"] + "entities/site.*.json")
    for in_file in file_list:
        doc_id = in_file.split(".")[-2]
        doc_id_list.append(doc_id)


    species_dict = util.get_species_dict(config_obj["misc_dir"])
    seq_dict, canon2taxid = util.load_seq_dict(config_obj["data_dir"], config_obj["misc_dir"] )
    geneid2canon = util.load_geneid2canon(config_obj["data_dir"])
    
    site_dict = {}
    entities.load_site_entities(site_dict, doc_id_list, config_obj["data_dir"], config_obj["misc_dir"])
    entities.load_glyco_entities(site_dict, doc_id_list, config_obj["data_dir"])
    entities.load_species_entities(site_dict, doc_id_list, config_obj["data_dir"])
    entities.load_gene_entities(site_dict, doc_id_list, config_obj["data_dir"], geneid2canon)
    entities.load_extragene_entities(site_dict, doc_id_list, config_obj["data_dir"], canon2taxid)
    entities.load_manualgene_entities(site_dict, config_obj["data_dir"], config_obj["misc_dir"])


    
    doc_id_list = list(site_dict.keys())
    if DEBUG:
        #doc_id_list = ["19369259"] 
        doc_id_list = ["9931318"]
  

    obj_dict = {}
    for doc_id in doc_id_list:
        obj_dict[doc_id] = []

        site_sent_idx_list = list(site_dict[doc_id]["site_sent_idx"].keys())
        glyco_sent_idx_list = list(site_dict[doc_id]["glyco_sent_idx"].keys())
        gene_sent_idx_list = list(site_dict[doc_id]["gene_sent_idx"].keys())
        extragene_sent_idx_list = list(site_dict[doc_id]["extragene_sent_idx"].keys())

        glyco_sent_idx_list = ["no_glyco_ents"] if glyco_sent_idx_list == [] else glyco_sent_idx_list
        gene_sent_idx_list = ["no_gene_ents"] if gene_sent_idx_list == [] else gene_sent_idx_list
        extragene_sent_idx_list = ["no_extragene_ents"] if extragene_sent_idx_list == [] else extragene_sent_idx_list


        gene_cmb_list = list(site_dict[doc_id]["gene"].keys())
        str_a = "_".join(site_sent_idx_list)
        str_b = "_".join(glyco_sent_idx_list)
        str_c = "_".join(gene_sent_idx_list)
        str_d = "_".join(extragene_sent_idx_list)


        seen_species = {}
        for tax_id in site_dict[doc_id]["species"]:
            tax_name = species_dict[tax_id] if tax_id in species_dict else tax_id
            seen_species[tax_name] = True
        str_sp = ";".join(list(seen_species.keys()))

        for site in site_dict[doc_id]["site"]:
            if gene_cmb_list == []:
                gene_id, gene_text, canon = "no_gene_id", "no_gene_text", "no_canon"
                obj = {
                    "docid":doc_id, "site":site, "geneid":gene_id, 
                    "gene_text":gene_text, "canon":canon,
                    "str_a":str_a, "str_b":str_b,"str_c":str_c,"str_d":str_d,"specieslist":str_sp
                }
                obj_dict[doc_id].append(obj)
            else:
                for gene_cmb in gene_cmb_list:
                    if gene_cmb.find("direct|") != -1:
                        gene_id, gene_text, canon = gene_cmb.split("|") 
                        obj = {
                            "docid":doc_id, "site":site, "geneid":gene_id, "gene_text":gene_text, 
                            "canon":canon,
                            "str_a":str_a, "str_b":str_b,"str_c":str_c,"str_d":str_d,
                            "specieslist":str_sp
                        }
                        obj_dict[doc_id].append(obj)
                    else:
                        gene_text, gene_id = gene_cmb.split("|")
                        if gene_id in geneid2canon:
                            for canon in geneid2canon[gene_id]:
                                obj = {
                                    "docid":doc_id, "site":site, "geneid":gene_id, 
                                    "gene_text":gene_text, "canon":canon,
                                    "str_a":str_a, "str_b":str_b,"str_c":str_c,"str_d":str_d,
                                    "specieslist":str_sp
                                }
                                obj_dict[doc_id].append(obj)
                        else:
                            canon = "no_canon"
                            obj = {
                                "docid":doc_id, "site":site, "geneid":gene_id, "canon":canon,
                                "gene_text":gene_text,
                                "str_a":str_a, "str_b":str_b,"str_c":str_c,"str_d":str_d,"specieslist":str_sp
                            }
                            obj_dict[doc_id].append(obj)

    for doc_id in obj_dict:
        tmp_dict = {}
        for obj in obj_dict[doc_id]:
            aa_one, pos = obj["site"].split("|")
            obj["status"] = "aa_mismatch"
            canon = obj["canon"]

            # Using full sequence targets
            if canon in seq_dict and pos.isdigit():
                seq_idx = int(pos)
                seq = seq_dict[canon]
                if seq_idx < len(seq):
                    obj["status"] = "aa_match" if aa_one.upper() == seq[seq_idx-1] else obj["status"] 
            site = obj["site"]
            if site not in tmp_dict:
                tmp_dict[site] = []
            tmp_dict[site].append(obj)
            


        out_file = config_obj["data_dir"] + "integrated/site.%s.json" % (doc_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(tmp_dict, indent=4)))
        if DEBUG:
            print ("created file:%s\n" % (out_file))

    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/integrated/"
    x = subprocess.getoutput(cmd)
   
        

if __name__ == '__main__':
    main()



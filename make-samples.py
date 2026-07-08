import os
import json
import util
import glob
import numpy
import subprocess




def main():

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())
    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + "/samples/") == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + "/samples/"
        x = subprocess.getoutput(cmd)

    seq_dict, canon2taxid = util.load_seq_dict(config_obj["data_dir"], config_obj["misc_dir"] )
    docid2taxid = util.load_docid2taxid(config_obj["data_dir"])
    glygen_species_dict = util.get_species_dict(config_obj["misc_dir"])
   
    docid2taxname = {} 
    for doc_id in docid2taxid:
        for tax_id in docid2taxid[doc_id]:
            if tax_id not in glygen_species_dict:
                continue
            tax_name = glygen_species_dict[tax_id]
            if doc_id not in docid2taxname:
                docid2taxname[doc_id] = {}
            docid2taxname[doc_id][tax_name] = True
   
 
    for in_file in glob.glob(config_obj["data_dir"] + "/sites/sites.*.json"):
        doc = json.load(open(in_file))
        for obj in doc:
            tax_name = obj["tax_name"]
            if tax_name == "unknown":
                continue
            if doc_id not in docid2taxname:
                docid2taxname[doc_id] = {}
            docid2taxname[doc_id][tax_name] = True

   


    in_file = config_obj["data_dir"] + "/match_sites/match_sites.csv"
    df = {}
    util.load_sheet(df, in_file, [], ",")
    f_list = df["fields"]
    canon2pos = {}
    pos2canon = {}
    docid2pos = {}
    for row in df["data"]:
        site = row[f_list.index("site")]
        cls = row[f_list.index("cls")]
        doc_id, canon, pos = site.split("|")
        pos = int(pos)
        if doc_id not in canon2pos:
            canon2pos[doc_id] = {}
        if canon not in canon2pos[doc_id]:
            canon2pos[doc_id][canon] = {}
        canon2pos[doc_id][canon][pos] = cls
        if doc_id not in pos2canon:
            pos2canon[doc_id] = {}
        if pos not in pos2canon[doc_id]:
            pos2canon[doc_id][pos] = {}
        pos2canon[doc_id][pos][canon] = cls

        if doc_id not in docid2pos:
            docid2pos[doc_id] = {}
        docid2pos[doc_id][pos] = True




    S = 10
    # FO =  1 if the organism or species corresponding to canon has been mentioned in the abstract
    # F0 =  0 if the abstract does not mention any organism
    # F0 = -1 if one or more organisms have been mentioned in the abstract and none matches the species of canon 
    # Fs =  1  if number of sequences matching the site is one
    # Fs = -1  if number of sequences matching the site is many

    FW1 = open(config_obj["data_dir"] + "/samples/samples_labeled.csv", "w")
    FW2 = open(config_obj["data_dir"] + "/samples/samples_all.csv", "w")
    row = ["site"]
    for i in range(1, S + 1):
        row.append("T1S%s" % (i))
    row += ["tax_id_match", "target_count", "cls"]
    FW1.write("\"%s\"\n" % ("\",\"".join(row)))
    FW2.write("\"%s\"\n" % ("\",\"".join(row)))

    for doc_id in canon2pos:
        canon_list = list(canon2pos[doc_id].keys())
        pos_list = sorted(list(docid2pos[doc_id].keys()))
        for canon in canon2pos[doc_id]:
            canon_full = canon.split("_")[0]
            for pos in canon2pos[doc_id][canon]:
                row_one, row_two = [], []
                for j in range(0, S):
                    flag = 0
                    s_pos = 0 
                    if j < len(pos_list):
                        s_pos = pos_list[j]
                        flag = 1 if s_pos in canon2pos[doc_id][canon] else -1
                    if pos == s_pos:
                        row_one.append(str(flag))
                    else:
                        row_two.append(str(flag))
                site = "%s|%s|%s" % (doc_id, canon, pos)
                row = [site] + row_one + row_two
                tax_id = canon2taxid[canon_full]
                tax_name = glygen_species_dict[tax_id]
                sp_flag = "0" 
                if doc_id in docid2taxname:
                    sp_flag = "1" if tax_name in docid2taxname[doc_id] else "-1"
                nseq_flag = "1" if len(list(pos2canon[doc_id][pos].keys())) == 1 else "-1"
                cls = canon2pos[doc_id][canon][pos]
                row += [sp_flag, nseq_flag, cls]
                if cls in ["1", "0"]:
                    FW1.write("\"%s\"\n" % ("\",\"".join(row)))
                FW2.write("\"%s\"\n" % ("\",\"".join(row)))
    FW1.close()


    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/samples" 
    x = subprocess.getoutput(cmd)

    return


if __name__ == '__main__':
    main()

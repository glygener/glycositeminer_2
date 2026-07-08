import os
import json
from sklearn import svm
from sklearn import datasets
import pickle
import util
import subprocess
import glob



def main():

    config_obj = json.loads(open("conf/config.json", "r").read())

    # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + "predicted/") == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + "predicted/"
        x = subprocess.getoutput(cmd)

    cutoff_file = config_obj["data_dir"] + "models/cutoff.json"
    cutoff_dict = json.load(open(cutoff_file))
    
    #curation_dict = json.load(open(config_obj["data_dir"] + "models/curated.json"))
    curation_dict = {}


    sample_file = config_obj["data_dir"] + "samples/samples_all.csv"
    model_file_svm = config_obj["data_dir"] + "models/svm.pkl"
    model_file_mlp = config_obj["data_dir"] + "models/mlp.pkl"

    species_dict = util.get_species_dict(config_obj["misc_dir"])
    seq_dict, canon2taxid = util.load_seq_dict(config_obj["data_dir"], config_obj["misc_dir"] )
    aa_dict = util.get_aa_dict(config_obj["misc_dir"]) 
    known_site_dict = util.get_known_site_dict(config_obj["data_dir"])

    X_tst, site_list = util.load_samples_unlabeled(sample_file)


    seen_cmb = {}
    count_dict = {"svm":{"total":0}, "mlp":{"total":0}, "either":{"total":0},"llm_either":{"total":0}}
    seen_row = {}
    with open(model_file_svm, 'rb') as f:
        clf_svm = pickle.load(f)
        y_prd_svm = clf_svm.predict(X_tst)
        proba_svm = clf_svm.predict_proba(X_tst)
        y_pred_prob_svm = proba_svm[:, 1] 
        for i in range(0, len(y_prd_svm)):
            site = site_list[i]
            doc_id, canon, pos = site.split("|")
            canon_full, new_pos = canon, pos
            if canon.find("_PDB") != -1:
                s = int(canon.split("_")[-3])
                new_pos = s + int(pos) - 1
                canon_full = canon.split("_")[0]
            aa_one =  seq_dict[canon_full][int(new_pos) -1]
            aa_three = aa_dict["one2three"][aa_one] if aa_one in aa_dict["one2three"] else ""
            glygen_status = "in_glygen_yes" if site in known_site_dict else "in_glygen_no"
            canon = site.split("|")[1]
            tax_id = canon2taxid[canon_full]
            tax_name = species_dict[tax_id] if tax_id in species_dict else tax_id
            flag = "pass" if y_pred_prob_svm[i] > cutoff_dict["svm"] else "fail"
            #print (y_pred_prob_svm[i], flag, "svm")
            if flag == "pass":
                row = [doc_id,canon,pos,aa_three,tax_name, glygen_status]
                row_str = json.dumps(row)
                if row_str not in seen_row:
                    seen_row[row_str] = {}
                seen_row[row_str]["svm"] = True
                if tax_name not in count_dict["svm"]:
                    count_dict["svm"][tax_name] = 0
                count_dict["svm"][tax_name] += 1
                count_dict["svm"]["total"] += 1
                cmb = "%s|%s|%s" % (doc_id,canon,pos)
                if cmb not in seen_cmb:
                    seen_cmb[cmb] = True
                    if tax_name not in count_dict["either"]:
                        count_dict["either"][tax_name] = 0
                    count_dict["either"][tax_name] += 1
                    count_dict["either"]["total"] += 1
                    ver_flag = "unknown"
                    if ver_flag == "yes":
                        if tax_name not in count_dict["llm_either"]:
                            count_dict["llm_either"][tax_name] = 0
                        count_dict["llm_either"][tax_name] += 1
                        count_dict["llm_either"]["total"] += 1
 
    with open(model_file_mlp, 'rb') as f:
        clf_mlp = pickle.load(f)
        y_prd_mlp = clf_mlp.predict(X_tst)
        proba_mlp = clf_mlp.predict_proba(X_tst)
        y_pred_prob_mlp = proba_mlp[:, 1]
        for i in range(0, len(y_prd_mlp)):
            site = site_list[i]
            doc_id, canon, pos = site.split("|")
            canon_full, new_pos = canon, pos
            if canon.find("_PDB") != -1:
                s = int(canon.split("_")[-3])
                new_pos = s + int(pos) - 1
                canon_full = canon.split("_")[0]
            aa_one =  seq_dict[canon_full][int(new_pos) -1]
            aa_three = aa_dict["one2three"][aa_one] if aa_one in aa_dict["one2three"] else ""
            glygen_status = "in_glygen_yes" if site in known_site_dict else "in_glygen_no"
            canon = site.split("|")[1]
            tax_id = canon2taxid[canon_full]
            tax_name = species_dict[tax_id] if tax_id in species_dict else tax_id
            flag = "pass" if y_pred_prob_mlp[i] > cutoff_dict["mlp"] else "fail"
            #print (y_pred_prob_svm[i], flag, "mlp")
            if flag == "pass":
                row = [doc_id,canon,pos,aa_three,tax_name, glygen_status]
                row_str = json.dumps(row)
                if row_str not in seen_row:
                    seen_row[row_str] = {}
                seen_row[row_str]["mlp"] = True
                if tax_name not in count_dict["mlp"]:
                    count_dict["mlp"][tax_name] = 0
                count_dict["mlp"][tax_name] += 1
                count_dict["mlp"]["total"] += 1 
                cmb = "%s|%s|%s" % (doc_id,canon,pos)
                if cmb not in seen_cmb:
                    seen_cmb[cmb] = True
                    if tax_name not in count_dict["either"]:
                        count_dict["either"][tax_name] = 0
                    count_dict["either"][tax_name] += 1
                    count_dict["either"]["total"] += 1
                    ver_flag = "unknown"
                    if ver_flag == "yes":
                        if tax_name not in count_dict["llm_either"]:
                            count_dict["llm_either"][tax_name] = 0
                        count_dict["llm_either"][tax_name] += 1
                        count_dict["llm_either"]["total"] += 1
    out_file = config_obj["data_dir"] + "predicted/predicted_tmp.csv"
    FW = open(out_file, "w")
    row = ["evidence", "uniprotkb_ac", "glycosylation_site", "amino_acid", "tax_name", "glygen_status", "algorithm", "llm_verification_flag", "curation_flag"]
    FW.write("\"%s\"\n" % ("\",\"".join(row)))
    for row_str in seen_row:
        algo = ";".join(list(seen_row[row_str].keys()))
        row = json.loads(row_str)
        site = "%s|%s|%s" % (row[0], row[1], row[2])
        ver_flag = "llm_unknown"
        cur_flag = "curation_" + curation_dict[site] if site in curation_dict else "not_curated" 
        FW.write("\"%s\"\n" % ("\",\"".join(row +  [algo,ver_flag, cur_flag]))) 
    FW.close()





    sp_list = list(set(list(count_dict["svm"].keys()) + list(count_dict["mlp"].keys())))
    sp_list.remove("total")
    sp_list.append("total")

    sp_obj_list = [
        {"sp":"human","lbl":"Human"},
        {"sp":"mouse","lbl":"Mouse"},
        {"sp":"rat","lbl":"Rat"},
        {"sp":"bovine","lbl":"Bovine"},
        {"sp":"sarscov2","lbl":"SARS-CoV-2"},
        {"sp":"yeast","lbl":"Yeast"},
        {"sp":"arabidopsis","lbl":"Arabidopsis"},
        {"sp":"pig","lbl":"Pig"},
        {"sp":"fruitfly","lbl":"Fruit fly"},
        {"sp":"chicken","lbl":"Chicken"},
        {"sp":"sarscov1","lbl":"HCoV-SARS"},
        {"sp":"zebrafish","lbl":"Zebrafish"},
        {"sp":"hamster","lbl":"Hamster"},
        {"sp":"hcv1a","lbl":"HCV-H77"},
        {"sp":"total","lbl":"Total"}
    ]
    tmp_obj_list = []
    seen = {}
    for obj in sp_obj_list:
        tmp_obj_list.append({"sp":obj["sp"],"lbl":obj["lbl"]})
        seen[obj["sp"]] = True
    for sp in sp_list:
        if sp not in seen:
            tmp_obj_list.append({"sp":sp,"lbl":""})
 
    out_file = config_obj["data_dir"] + "/predicted/stats.txt"
    
    FW = open(out_file, "w")
    row = ["tax_name","tax_name_lbl","SVM", "MLP","SVM or MLP", "LLM verified"]
    FW.write("%s\n" % (",".join(row)))
    for obj in tmp_obj_list:
        sp, lbl = obj["sp"], obj["lbl"]
        row = [sp, lbl]
        total = 0
        for alg in ["svm", "mlp", "either", "llm_either"]:
            n = count_dict[alg][sp] if sp in count_dict[alg] else 0
            total += n
            row.append(str(n))
        if total > 0:
            FW.write("%s\n" % (",".join(row)))
    FW.close()
 

    cmd = "chmod -R 777 " + config_obj["data_dir"] + "/predicted/"
    x = subprocess.getoutput(cmd)



    return


if __name__ == '__main__':
    main()





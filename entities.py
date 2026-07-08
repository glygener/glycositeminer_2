import sys,os
import json
import util
import glob
import subprocess
import numpy
import pickle



def write_obj_list(doc_id, obj_list, out_file):

    combo_id = "%s|%s|%s" % (doc_id, -1, -1)
    doc = {combo_id:obj_list}
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(doc, indent=4)))



    
def dump_pubtator_entities(doc_id, base_dir):   
    
    gene_obj_list = []
    species_obj_list = []
    in_file = base_dir + "pubtator_extracts/pmid.%s.txt" % (doc_id)
    print (in_file, os.path.isfile(in_file))

    if os.path.isfile(in_file) == False:
        return
    with open(in_file, "r") as FR:
        for line in FR:
            row = line.split("\t")
            if len(row) != 6:
                continue
            d, s, e = row[0].strip(), row[1].strip(), row[2].strip()
            text, ent_type = row[3].strip(), row[4].strip()
            if d != doc_id:
                continue
            
            if ent_type in ["Gene", "Species"]  and text != "":
                ent_id = row[5].strip()
                source = "Entrez" if ent_type == "Gene" else "NCBI"
                o = {"text":text, "id":[{"source": source, "idString": ent_id}], "s": s, "e": e}
                if ent_type == "Gene":
                    gene_obj_list.append(o)
                if ent_type == "Species":
                    species_obj_list.append(o)

    gene_dict, species_dict = {}, {}
    combo_id = "%s|%s|%s" % (doc_id, -1, -1)
    if len(gene_obj_list) > 0:
        gene_dict = {combo_id:gene_obj_list}
        out_file = base_dir + "/entities/gene.%s.json" % (doc_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(gene_dict, indent=4)))

    if len(species_obj_list) > 0:
        species_dict = {combo_id:species_obj_list}
        out_file = base_dir + "/entities/species.%s.json" % (doc_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(species_dict, indent=4)))
   
    return







def dump_glycositeminer_entities(doc_id, base_dir, gene_dict, aa_dict, known_aa):
    in_file = base_dir + "/medline_extracts/pmid.%s.json" % (doc_id)
    ignore_terms = ["the", "to", "and", "in", "a", "was", "eg", "in", "be", "how", "for", "vs"]
    seen_aa = {}
    for aa in known_aa:
        seen_aa[aa.lower()] = True
    
    site_dict, glyco_dict, extragene_dict = {}, {}, {}
    doc = json.loads(open(in_file, "r").read())
    #s_list, seen_lbl = [], {}
    all_ent_list = []
    for obj in doc["sent_list"]:
        sent = obj["sentence"]
        sentence_index = obj["sentence_idx"]
        g_list, e_list = [], []
        s_list, seen_lbl = [], {}
        sent_combo_id = "%s|%s|%s" % (sentence_index, 1, len(sent))
        for o in obj["entities"]:
            lbl, text, start, end = o["label"], o["text"], o["start"], o["end"]
            seen_lbl[lbl] = True
            if lbl == "SITE":
                parts = text.split("-")
                in_flag = False
                if len(parts) > 1:
                    aa, pos = parts[0].strip().lower(), parts[1].strip()
                    if aa in seen_aa and pos.isdigit():
                        s_list.append({"text":text,"s":start, "e":end})
                        in_flag = True
                all_ent_list.append({"text":text,"s":start, "e":end})
            if lbl == "GLYCOSYLATION":
                g_list.append({"text":text,"s":start, "e":end})
        if "SITE" in seen_lbl and s_list != []:
            site_dict[sent_combo_id] = s_list
        if g_list != []:
            glyco_dict[sent_combo_id] = g_list
    
        for term in sent.split(" "):
            term = term.replace(",", "").replace(".", "")
            if term.lower() == term:
                continue
            if len(term) < 3:
                continue
            if term.lower() in ignore_terms:
                continue
            if term in aa_dict["three"] or term in aa_dict["full"]:
                continue
            if term in gene_dict:
                s, e = -1, -1
                ent = {"text": term, "id": [], "s": s, "e": e}
                for canon in gene_dict[term]:
                    ent["id"].append({"source": "UniProtKB", "idString": canon})
                e_list.append(ent)
        if e_list != []:
            extragene_dict[sent_combo_id] = e_list

    if site_dict != {}:
        out_file = base_dir + "/entities/site.%s.json" % (doc_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(site_dict, indent=4)))
        
        if glyco_dict != {}:
            out_file = base_dir + "/entities/glyco.%s.json" % (doc_id)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(glyco_dict, indent=4)))
        
        if extragene_dict != {}:
            out_file = base_dir + "/entities/extragene.%s.json" % (doc_id)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(extragene_dict, indent=4)))

    return 









def integrate_entities(doc_id_list, base_dir, misc_dir):


    species_dict = util.get_species_dict(misc_dir)
    seq_dict, canon2taxid = util.load_seq_dict(base_dir, misc_dir)
    geneid2canon = util.load_geneid2canon(base_dir)
 
    site_dict = {}
    load_site_entities(site_dict, doc_id_list, base_dir, misc_dir)
    load_glyco_entities(site_dict, doc_id_list, base_dir)
    load_species_entities(site_dict, doc_id_list, base_dir)
    load_gene_entities(site_dict, doc_id_list, base_dir, geneid2canon)
    load_extragene_entities(site_dict, doc_id_list, base_dir, canon2taxid)


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


        gene_id_list = list(site_dict[doc_id]["gene"].keys())
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
            if gene_id_list == []:
                gene_id, canon = "no_gene_id", "no_canon"
                obj = {
                    "docid":doc_id, "site":site, "geneid":gene_id, "canon":canon,
                    "str_a":str_a, "str_b":str_b,"str_c":str_c,"str_d":str_d,"specieslist":str_sp
                }
                obj_dict[doc_id].append(obj)
            else:
                for gene_id in gene_id_list:
                    if gene_id.find("direct|") != -1:
                        gene_id, canon = gene_id.split("|") 
                        obj = {
                            "docid":doc_id, "site":site, "geneid":gene_id, "canon":canon,
                            "str_a":str_a, "str_b":str_b,"str_c":str_c,"str_d":str_d,"specieslist":str_sp
                        }
                        obj_dict[doc_id].append(obj)
                    elif gene_id in geneid2canon:
                        for canon in geneid2canon[gene_id]:
                            obj = {
                                "docid":doc_id, "site":site, "geneid":gene_id, "canon":canon,
                                "str_a":str_a, "str_b":str_b,"str_c":str_c,"str_d":str_d,"specieslist":str_sp
                            }
                            obj_dict[doc_id].append(obj)
                    else:
                        canon = "no_canon"
                        obj = {
                            "docid":doc_id, "site":site, "geneid":gene_id, "canon":canon,
                            "str_a":str_a, "str_b":str_b,"str_c":str_c,"str_d":str_d,"specieslist":str_sp
                        }
                        obj_dict[doc_id].append(obj)

    tmp_dict = {}
    for doc_id in obj_dict:
        tmp_dict[doc_id] = {}
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
            if site not in tmp_dict[doc_id]:
                tmp_dict[doc_id][site] = []
            tmp_dict[doc_id][site].append(obj)
            

    return tmp_dict




def load_site_entities(site_dict, doc_id_list, base_dir, misc_dir):

    aa_dict = util.get_aa_dict(misc_dir)
    for doc_id in doc_id_list:
        in_file = base_dir + "/entities/site.%s.json" % (doc_id)
        if os.path.isfile(in_file) == False:
            continue
        doc = json.loads(open(in_file, "r").read())
        for c in doc:
            sent_idx = c.split("|")[0]
            for obj in doc[c]:
                if obj["llm_flag"] == False:
                    continue
                aa, pos = obj["text"].split("-")[0].strip(), obj["text"].split("-")[1].strip()
                aa = aa[0].upper() + aa[1:].lower()
                pos = util.check_pos(pos)
                if len(aa) == 3:
                    aa = aa[0].upper() + aa[1:].lower()
                    aa = aa_dict["three2one"][aa]
                elif len(aa) > 3:
                    aa = aa_dict["full2one"][aa]
                site = "%s|%s" % (aa, pos)
                if doc_id not in site_dict:
                    site_dict[doc_id] = {
                        "site":{}, 
                        "gene":{},
                        "canon":{},
                        "species":{},
                        "site_sent_idx":{}, 
                        "glyco_sent_idx":{}, 
                        "gene_sent_idx":{},
                        "extragene_sent_idx":{}
                    }
                site_dict[doc_id]["site_sent_idx"][sent_idx] = True
                site_dict[doc_id]["site"][site] = True
    return


def load_glyco_entities(site_dict,doc_id_list, base_dir):

    for doc_id in doc_id_list:
        in_file = base_dir + "/entities/glyco.%s.json" % (doc_id)
        if os.path.isfile(in_file) == False:
            continue
        doc = json.loads(open(in_file, "r").read())
        for c in doc:
            sent_idx = c.split("|")[0]
            if doc_id in site_dict:
                site_dict[doc_id]["glyco_sent_idx"][sent_idx] = True

    return


def load_gene_entities(site_dict, doc_id_list, base_dir, geneid2canon):
    
    for doc_id in doc_id_list:
        in_file = base_dir + "/entities/gene.%s.json" % (doc_id)
        if os.path.isfile(in_file) == False:
            continue
        if doc_id not in site_dict:
            continue
        doc = json.loads(open(in_file, "r").read())
        for c in doc:
            sent_idx = c.split("|")[0]
            site_dict[doc_id]["gene_sent_idx"][sent_idx] = True
            for obj in doc[c]:
                gene_text = obj["text"]
                for o in obj["id"]:
                    if o["source"] == "Entrez":
                        gene_id = o["idString"]
                        gene_cmb = gene_text + "|" + o["idString"]
                        site_dict[doc_id]["gene"][gene_cmb] = True
                        if gene_id in geneid2canon:
                            for canon in geneid2canon[gene_id]:
                                site_dict[doc_id]["canon"][canon] = True
    return


def load_extragene_entities(site_dict, doc_id_list, base_dir, canon2taxid):
   
    for doc_id in doc_id_list:
        in_file = base_dir + "/entities/extragene.%s.json" % (doc_id)
        if os.path.isfile(in_file) == False:
            continue
        if doc_id not in site_dict:
            continue
        tax_id_list = list(site_dict[doc_id]["species"].keys())
        doc = json.loads(open(in_file, "r").read())
        for c in doc:
            sent_idx = c.split("|")[0]
            site_dict[doc_id]["extragene_sent_idx"][sent_idx] = True
            for obj in doc[c]:
                gene_text = obj["text"] if obj["text"].strip() != "" else "no_gene_text"
                for o in obj["id"]:
                    if o["source"] == "UniProtKB":
                        canon = o["idString"] if o["idString"].strip() != "" else "no_canon"
                        tax_id = str(canon2taxid[canon]) if canon in canon2taxid else ""
                        gene_cmb = "direct|%s|%s" % (gene_text,  canon)
                        if len(site_dict[doc_id]["canon"].keys()) == 0:
                            site_dict[doc_id]["gene"][gene_cmb] = True


    return




def load_manualgene_entities(site_dict, base_dir, misc_dir):

    in_file = misc_dir + "canon_map.csv"
    df = {}
    util.load_sheet(df, in_file, [], ",")
    f_list = df["fields"]
    for row in df["data"]:
        doc_id = row[f_list.index("doc_id")] 
        canon = row[f_list.index("canon")]
        gene_text = row[f_list.index("text")]
        gene_id = "direct|%s|%s" % (gene_text, canon)
        if doc_id in site_dict:
            site_dict[doc_id]["gene_sent_idx"]["0"] = True
            site_dict[doc_id]["gene"][gene_id] = True
    
    return




def load_species_entities(site_dict, doc_id_list, base_dir):

    for doc_id in doc_id_list:
        in_file = base_dir + "/entities/species.%s.json" % (doc_id)
        if os.path.isfile(in_file) == False:
            continue
        if doc_id not in site_dict:
            continue
        doc = json.loads(open(in_file, "r").read())
        for c in doc:
            sent_idx = c.split("|")[0]
            site_dict[doc_id]["gene_sent_idx"][sent_idx] = True
            for obj in doc[c]:
                gene_txt = obj["text"]
                for o in obj["id"]:
                    if o["source"] == "NCBI":
                        tax_id = o["idString"]
                        site_dict[doc_id]["species"][tax_id] = True

    return






def get_match_dict(entegrated_dict):

    match_dict, docid2pos, docid2canon = {}, {}, {}
    for doc_id in entegrated_dict:
        doc = entegrated_dict[doc_id]
        docid2pos[doc_id] = {}
        docid2canon[doc_id] = {}
        glyco_ent, direct, entrez = False, False, False
        for site in doc:
            pos = site.split("|")[1]
            if pos.find("bad_pos") != -1:
                continue
            for obj in doc[site]:
                if obj["geneid"] == "direct":
                    direct = True
                if obj["geneid"] not in  ["direct", "no_gene_id"]:
                    entrez = True
                v = obj["canon"].split("|")[0]
                if v in ["no_canon"]:
                    continue
                if obj["str_b"] != "no_glyco_ents":
                    glyco_ent = True
                if obj["str_b"] == "no_glyco_ents":
                    continue
                r_canon = util.get_r_canon(obj["canon"])
                c = "%s|%s|%s" % (doc_id, r_canon, pos)
                if c not in match_dict:
                    match_dict[c] = "1" if obj["status"] == "aa_match" else "-1"
                elif match_dict[c] == "-1":
                    match_dict[c] = "1" if obj["status"] == "aa_match" else "-1"
                docid2canon[doc_id][r_canon] = True
                docid2pos[doc_id][pos] = True

    return match_dict, docid2pos, docid2canon









def get_pos_list(doc_id, orig_pos_list, canon_list, match_dict, S):
    
    score2pos = {}
    for pos in orig_pos_list:
        score = 0
        for canon in canon_list:
            cmb = "%s|%s|%s" % (doc_id, canon, pos)
            score += int(match_dict[cmb]) if cmb in match_dict else 0
        if score not in score2pos:
            score2pos[score] = {}
        score2pos[score][pos] = True

    pos_list = []
    for score in sorted(score2pos, reverse=True):
        for pos in score2pos[score]:
            if len(pos_list) < S:
                pos_list.append(int(pos))

    return pos_list




def make_predictions(sample_rows, cutoff_svm, cutoff_mlp, base_dir):

    X_tst, site_list = [], []
    for row in sample_rows:
        site_list.append(row[0])
        newrow = []
        for v in row[1:-1]:
            newrow.append(float(v))
        X_tst.append(numpy.array(newrow))
    X_tst = numpy.array(X_tst).astype(float)

    svm_obj_list, mlp_obj_list = [], []
    model_file_svm = base_dir + "/models/svm.pkl"
    model_file_mlp = base_dir + "/models/mlp.pkl"
    with open(model_file_svm, 'rb') as f:
        clf = pickle.load(f)
        y_prd = clf.predict(X_tst)
        proba = clf.predict_proba(X_tst)
        y_prob = proba[:, 1] 
        for i in range(0, len(y_prd)):
            site = site_list[i]
            canon = site.split("|")[1]
            flag = "pass" if y_prob[i] > cutoff_svm else "fail"
            svm_obj_list.append({"site":site,"status":flag,"y_pred_prob":y_prob[i],"y_prd":int(y_prd[i])}) 
    

    with open(model_file_mlp, 'rb') as f:
        clf = pickle.load(f)
        y_prd = clf.predict(X_tst)
        proba = clf.predict_proba(X_tst)
        y_prob = proba[:, 1]
        for i in range(0, len(y_prd)):
            site = site_list[i]
            canon = site.split("|")[1]
            flag = "pass" if y_prob[i] > cutoff_mlp else "fail"
            mlp_obj_list.append({"site":site,"status":flag,"y_pred_prob":y_prob[i],"y_prd":int(y_prd[i])}) 


    return svm_obj_list, mlp_obj_list





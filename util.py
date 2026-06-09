import os
import json
import csv
import glob
import subprocess
import numpy
import requests


def write_output(out_file, row_list):

    with open(out_file, "w") as FW: 
        for row in row_list:
            FW.write("\"%s\"\n" % ("\",\"".join(row)))
    return



def get_species_dict(misc_dir):

    species_map = json.load(open(misc_dir + "species_map.json"))
    species_dict = {}
    for tax_id in species_map:
        tax_name = species_map[tax_id][0]
        species_dict[tax_id] = tax_name
        species_dict[tax_name] = tax_id
    return species_dict





def get_aa_dict(misc_dir):

    aa_dict = {"full2one":{}, "one2three":{}, "three2one":{}, "three2full":{}, "one2full":{}, "full2three":{}, "one":[], "three":[], "full":[]}
     # create amino acid patterns
    in_file = misc_dir + "aa.csv"
    line_list = open(in_file, "r").read().split("\n")
    for line in line_list[1:]:
        if line.strip() == "":
            continue
        one,three,full = line.strip().split(",")
        three = three[0].upper() + three[1:]
        aa_dict["one"].append(one)
        aa_dict["three"].append(three)
        aa_dict["full"].append(full)
        aa_dict["one"].append(one.lower())
        aa_dict["three"].append(three.lower())
        aa_dict["full"].append(full.lower())
        aa_dict["three2one"][three] = one
        aa_dict["one2three"][one] = three
        aa_dict["full2one"][full] = one
        aa_dict["three2full"][three] = full
        aa_dict["one2full"][one] = full
        aa_dict["full2three"][full] = three
    return aa_dict


def get_aa_three(aa, aa_dict):
    aa_three = "Xxx"
    aa_three = aa if len(aa) == 3 else aa_three
    aa_three = aa_dict["one2three"][aa.upper()] if aa.upper() in aa_dict["one2three"] else aa_three
    aa_three = aa_dict["full2three"][aa] if aa in aa_dict["full2three"] else aa_three
    return aa_three
            

def get_known_aa(base_dir, misc_dir):

    aa_dict = get_aa_dict(misc_dir)
    known_aa = {}
    in_file = base_dir + "glygen/known_sites.csv"
    df = {}
    load_sheet(df, in_file, [], ",")
    f_list = df["fields"]
    for row in df["data"]:
        pos = row[f_list.index("pos")]
        if pos == "":
            continue
        aa_three = row[f_list.index("amino_acid")]
        if aa_three not in aa_dict["three2one"]:
            continue
        aa_one = aa_dict["three2one"][aa_three] 
        aa_full = aa_dict["three2full"][aa_three] 
        known_aa[aa_one] = True
        known_aa[aa_three] = True
        known_aa[aa_full] = True
    return known_aa

def get_glyco_terms(misc_dir):

    tmp_dict = {}
    in_file = misc_dir + "glyco.csv"
    df = {}
    load_sheet(df, in_file, [], ",")
    f_list = df["fields"]
    for row in df["data"]:
        term = row[f_list.index("term")]
        tmp_dict[term] = True
        tmp_dict[term.lower()] = True
    return tmp_dict


def no_digits(input_string):
    return not any(char.isdigit() for char in input_string)


def transform_sent(sent, aa_dict):

    char_set = "abcdefghigklmnopqrstuvwxyz"

    in_w_list = sent.split(" ")
    w_list = []
    for w_idx in range(0, len(in_w_list)):
        w = in_w_list[w_idx].strip()
        orig = w.strip()
        w = w.strip()
        if w == "":
            continue
        if w[0] == "(":
            w = w[1:]
        if w != "":
            if w[-1] in [")",",","."]:
                w = w[:-1]
        if w == "":
            continue
        flag = False
        for aa in aa_dict["full"]:
            if w.lower().find(aa.lower()) != -1:
                one, two = aa, w.lower().replace(aa.lower(), "")
                two = two.replace("-", "").replace(",", "").replace(".", "").strip()
                three = in_w_list[w_idx + 1] if len(in_w_list) > w_idx + 1 else ""
                three = three.replace(".","").replace(",", "")
                #print ("|%s|%s|%s|%s|%s" % (w, one, two, three, two.isdigit()))
                if two.isdigit():
                    w = one + " - " + two
                    flag = True
                    #print (orig, "-->", w)
                    break
                elif three.isdigit():
                    #w = one + " - " + three
                    w = one + " -"
                    flag = True
                    #print (orig, "-->", w)
                    break
        if flag == False and w[:3] in aa_dict["three"]:
            #Asn123, Asn(123), Asn-123
            one, two = w[:3], w[3:].strip()
            three = in_w_list[w_idx + 1] if len(in_w_list) > w_idx + 1 else ""
            three = three.replace(".","")
            if two == "" and three.isdigit():
                #w = one + " - " + three.strip()
                w = one + " -"
            elif two != "":
                if two[0] in ["1","2","3","4","5","6","7","8","9"]:
                    w = one + " - "
                    for jj in range(0, len(two)):
                        if two[jj].isdigit() == False:
                            break
                        w += two[jj]
                    w = w.replace(".", "").replace(",", "")
                    flag = True
                    #print (one, two, "-->", w)
                elif len(two) > 1:
                    #print (w, one, two)
                    if two[0] in ["-", "("] and two[1] in ["1","2","3","4","5","6","7","8","9"]:
                        w = one + " - " 
                        for jj in range(1, len(two)):
                            if two[jj].isdigit() == False:
                                break
                            w += two[jj]
                        w = w.replace(".", "").replace(",", "")
                        flag = True
        if flag == False and w[0] in aa_dict["one"]:
            #N123, N(123), N-123
            tmp_flag = False
            if len(w) > 1:
                if w[1] in ["1","2","3","4","5","6","7","8","9"]:
                    tmp_flag = True
            if len(w) > 2:
                if w[1] in ["-","("] and w[2] in ["1","2","3","4","5","6","7","8","9"]:
                    tmp_flag = True
            if tmp_flag:
                w = w[0] + " - " + w[1:].replace(" ", "").replace("(", "").replace(")", "").replace("-","")
                w = w.replace(".", "").replace(",", "")
                flag = True
                #print (orig, "-->", w)
        w_list.append(w)

    s = " ".join(w_list)
    s = s.replace("-  ", "- ") 
    return s

def get_obj_dict(in_file):

    obj_dict = {}
    line_list = open(in_file, "r").read().split("\n")
    for line in line_list:
        if line.strip() == "":
            continue
        obj = json.loads(line)
        doc_id = obj["docId"]
        obj_dict[doc_id] = obj
    return obj_dict



def load_data(file):
    data = json.loads(open(file, "r", encoding="utf-8").read())
    return (data)


def save_data(file, data):
    with open (file, "w", encoding="utf-8") as FW:
        FW.write("%s\n" % (json.dumps(data,indent=4)))

    return

def process_ne(file):
    
    data = load_data(file)
    name_list = []
    for item in data:
        name_list.append(item)
    
    for item in data:
        item = item.replace("The", "").replace("the", "").replace("and", "").replace("And", "")
        names = item.split(" ")
        for name in names:
            name = name.strip()
            name_list.append(name)
        if "(" in item:
            names = item.split("(")
            for name in names:
                name = name.replace(")", "").strip()
                name_list.append(name)
        if "," in item:
            names = item.split(",")
            for name in names:
                name = name.replace("and", "").strip()
                if " " in name:
                    new_names = name.split()
                    for x in new_names:
                        x = x.strip()
                        name_list.append(x)
                name_list.append(name)
    
    final_name_list = []
    #titles = ["Dr.", "Professor", "Mr.", "Mrs.", "Ms.", "Miss", "Aunt", "Uncle", "Mr. and Mrs."]
    titles = []
    for character in name_list:
        if "" != character:
            final_name_list.append(character)
            for title in titles:
                titled_char = f"{title} {character}"
                final_name_list.append(titled_char)
    
    final_name_list = list(set(final_name_list))
    final_name_list.sort()
    return (final_name_list)





def load_sheet(sheet_obj, in_file, field_list, separator):

    seen = {}
    sheet_obj["fields"] = []
    sheet_obj["data"] = []
    field_ind_list = []
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=separator, quotechar='\"')
        row_count = 0
        ncols = 0
        for row in csv_grid:
            if json.dumps(row) in seen:
                continue
            seen[json.dumps(row)] = True
            row_count += 1
            for j in range(0, len(row)):
                row[j] = row[j].replace("\"", "`")
            if row_count == 1:
                ncols = len(row)
                for j in range(0, len(row)):
                    f = row[j].strip()
                    if field_list != [] and f in field_list:
                        field_ind_list.append(j)
                        sheet_obj["fields"].append(f)
                    else:
                        field_ind_list.append(j)
                        sheet_obj["fields"].append(f)
            else:
                #make sure every row has ncols columns
                if len(row) != ncols:
                    print ("bad structure of CSV file")
                    print (sheet_obj["fields"])
                    print (row)
                    exit()
                new_row = []
                for j in field_ind_list:
                    new_row.append(row[j].strip())
                sheet_obj["data"].append(new_row)



    return







            
def load_seq_dict(base_dir, misc_dir):

    species_dict = get_species_dict(misc_dir)
    seq_hash = {} 
    sp_hash = {}
    file_list = glob.glob(base_dir + "glygen/*.fasta")
    for in_file in file_list:
        tax_name = in_file.split("/")[-1].split("_")[0]
        tax_id = species_dict[tax_name]
        with open(in_file, "r") as FR:
            for line in FR:
                if line[0] == ">":
                    seq_id = line.split("|")[1]
                    seq_hash[seq_id] = ""
                    sp_hash[seq_id] = tax_id
                else:
                    seq_hash[seq_id] += line.strip()

    return seq_hash, sp_hash


def check_pos(pos):

    try:
        pos = int(pos)
        return str(pos)
    except:
        return "bad_pos_%s" % (pos)


def load_pmid2canon(in_file):

    tmp_dict = {}
    df = {}
    load_sheet(df, in_file, [], ",")
    f_list = df["fields"]
    for row in df["data"]:
        doc_id = row[f_list.index("doc_id")]
        gene_id = row[f_list.index("gene_id")]
        canon = row[f_list.index("canon")]
        flag_one = row[f_list.index("pubtator_status")]
        flag_two = row[f_list.index("glygen_status")]    
        if doc_id not in tmp_dict:
            tmp_dict[doc_id] = {}
        c = "%s|%s" % (gene_id, canon)
        tmp_dict[doc_id][c] = "%s|%s" % (flag_one, flag_two)
        
    return tmp_dict



def get_site_dict(in_file, aa_dict, known_aa):

    doc = json.loads(open(in_file, "r").read())
    site_dict = {}
    for combo in doc:
        sent_idx = combo.split("|")[0]
        for obj in doc[combo]:
            aa = obj["text"].split("-")[0].strip()
            #if len(aa) < 3:
            #    continue
            pos = check_pos(obj["text"].split("-")[1].strip())
            if pos.isdigit() == False:
                continue
            pos = int(pos)
            aa_one = aa
            if len(aa) == 3:
                aa_one = aa_dict["three2one"][aa] if aa in aa_dict["three2one"] else ""
            elif len(aa) > 3:
                aa_one = aa_dict["full2one"][aa] if aa in aa_dict["full2one"] else ""
            if aa_one not in known_aa:
                continue
            site = "%s|%s" % (aa_one, pos)
            if site not in site_dict:
                site_dict[site] = {}
            site_dict[site][sent_idx] = True

    return site_dict


def get_gene_sent_idx(in_file):

    tmp_dict = {}
    if os.path.isfile(in_file):
        doc = json.loads(open(in_file, "r").read())
        for combo in doc:
            sent_idx = combo.split("|")[0]
            for obj in doc[combo]:
                for o in obj["id"]:
                    gene_id = o["idString"]
                    if gene_id not in tmp_dict:
                        tmp_dict[gene_id] = {}
                    tmp_dict[gene_id][sent_idx] = True
    return tmp_dict

def get_glyco_sent_idx(in_file):

    tmp_dict = {}
    if os.path.isfile(in_file):
        doc = json.loads(open(in_file, "r").read())
        for combo in doc:
            sent_idx = combo.split("|")[0]
            tmp_dict[sent_idx] = True
   
 
    return list(tmp_dict.keys())



def get_sent_dist(list_one, list_two):

    min_d = 1000
    for idx_one in list_one:
        for idx_two in list_two:
            d = abs(int(idx_one) - int(idx_two))
            min_d = d if d < min_d else min_d

    return min_d



def load_geneid2ensg(base_dir):

    geneid2ensg = {}
    file_list = [base_dir + "gene_info/All_Data.gene_info"]
    for in_file in file_list:
        df = {}
        load_sheet(df, in_file, [], "\t")
        f_list = df["fields"]
        for row in df["data"]:
            gene_id = row[f_list.index("GeneID")]
            xrefs = row[f_list.index("dbXrefs")].split("|")
            for x in xrefs:
                if x.split(":")[0] == "Ensembl":
                    ensg = x.split(":")[1]
                    if gene_id not in geneid2ensg:
                        geneid2ensg[gene_id] = {}
                    geneid2ensg[gene_id][ensg] = True

    return geneid2ensg


def load_geneid2canon(base_dir):

    canon2sp = {}
    geneid2ensg = load_geneid2ensg(base_dir)
    ensg2canon = {}
    file_list = glob.glob(base_dir + "glygen/*_protein_genelocus.csv")
    for in_file in file_list:
        sp = in_file.split("/")[-1].split("_")[0]
        df = {}
        load_sheet(df, in_file, [], ",")
        f_list = df["fields"]
        for row in df["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            ensg = row[f_list.index("ensembl_gene_id")]
            if ensg not in ensg2canon:
                ensg2canon[ensg] = {}
            ensg2canon[ensg][canon] = True
            canon2sp[canon] = sp  

    geneid2canon = {}
    file_list = glob.glob(base_dir + "glygen/*_protein_xref_geneid.csv")
    for in_file in file_list:
        sp = in_file.split("/")[-1].split("_")[0]
        df = {}
        load_sheet(df, in_file, [], ",")
        f_list = df["fields"]
        for row in df["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            gene_id = row[f_list.index("xref_id")]
            if gene_id not in geneid2canon:
                geneid2canon[gene_id] = {}
            geneid2canon[gene_id][canon] = sp
            canon2sp[canon] = sp

    for gene_id in geneid2ensg:
        for ensg in geneid2ensg[gene_id]:
            if ensg in ensg2canon:
                for canon in ensg2canon[ensg]:
                    if gene_id not in geneid2canon:
                        geneid2canon[gene_id] = {}
                    geneid2canon[gene_id][canon] = canon2sp[canon]

    return geneid2canon








def get_known_site_dict(base_dir):

    tmp_dict = {}
    in_file = base_dir + "glygen/known_sites.csv"
    df = {}
    load_sheet(df, in_file, [], ",")
    f_list = df["fields"]
    for row in df["data"]:
        doc_id = row[f_list.index("doc_id")]
        canon = row[f_list.index("canon")]
        ac = row[f_list.index("canon")].split("-")[0]
        pos = row[f_list.index("pos")]
        aa_three = row[f_list.index("amino_acid")]
        source = row[f_list.index("source")]
        combo_one = "%s|%s|%s" % (doc_id, canon, pos)
        combo_two = "%s|%s|%s" % (doc_id, ac, pos)
        if combo_one not in tmp_dict:
            tmp_dict[combo_one] = {}
            tmp_dict[combo_two] = {}
        tmp_dict[combo_one][source] = True
        tmp_dict[combo_two][source] = True
    return tmp_dict






def write_output(row_list, out_file):    
    with open(out_file, "w") as FW:
        for row in row_list:
            FW.write("%s\n" % (",".join(row)))
    return


def load_docid2taxid(base_dir):

    tmp_dict = {}
    for in_file in glob.glob(base_dir + "/entities/species.*.json"):
        doc_id = in_file.split(".")[-2]
        sp_doc = json.loads(open(in_file, "r").read())
        for cmb in sp_doc:
            for obj in sp_doc[cmb]:
                for o in obj["id"]:
                    tax_id = o["idString"]
                    if doc_id not in tmp_dict:
                        tmp_dict[doc_id] = {}
                    tmp_dict[doc_id][tax_id] = True
    return tmp_dict
            


def load_rank_dict(base_dir):

    tmp_dict = {}
    for in_file in glob.glob(base_dir + "/ranked/sites.*.csv"):
        typ = in_file.split(".")[-2]
        df = {}
        load_sheet(df, in_file, [], ",")
        f_list = df["fields"]
        for row in df["data"]:
            doc_id = row[f_list.index("doc_id")]
            site = row[f_list.index("site")]
            canon = row[f_list.index("canon")].split("|")[1].split("_")[0]
            rank = row[f_list.index("rank")]
            c = "%s|%s|%s" % (doc_id, canon, site.split("|")[1])
            tmp_dict[c] =  rank
    return tmp_dict

def load_glygenstatus_dict(base_dir):

    tmp_dict = {}
    for in_file in glob.glob(base_dir + "/ranked/sites.*.csv"):
        typ = in_file.split(".")[-2]
        df = {}
        load_sheet(df, in_file, [], ",")
        f_list = df["fields"]
        for row in df["data"]:
            doc_id = row[f_list.index("doc_id")]
            site = row[f_list.index("site")]
            canon = row[f_list.index("canon")].split("|")[1].split("_")[0]
            rank = row[f_list.index("rank")]
            status = row[f_list.index("glygen_status")]
            c = "%s|%s|%s" % (doc_id, canon, site.split("|")[1])
            tmp_dict[c] =  status
    return tmp_dict



def load_class_dict(base_dir):
    tmp_dict = {}
    in_file = base_dir + "/samples/samples_labeled.csv"
    df = {}
    load_sheet(df, in_file, [], ",")
    f_list = df["fields"]
    for row in df["data"]:
        site = row[f_list.index("site")]
        tmp_dict[site] = row[f_list.index("cls")]
    return tmp_dict




def get_r_canon(canon):
    r_canon = canon.replace("PDB|", "").replace("MET|", "").replace("HOMOLOG|", "").split("|")[0]
    
    if canon.find("PDB|") != -1:
        r_canon = canon.split("|")[0] + "." + r_canon
    
    return r_canon



def get_match_dict(enteg_dir):

    site_count = 0
    docid2canon = {}
    docid2pos = {}
    match_dict = {}
    file_list = glob.glob(enteg_dir +  "/site.*.json")
    for in_file in file_list:
        doc_id = in_file.split(".")[-2]
        doc = json.loads(open(in_file, "r").read())
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
                r_canon = get_r_canon(obj["canon"])
                c = "%s|%s|%s" % (doc_id, r_canon, pos)
                if c not in match_dict:
                    match_dict[c] = "1" if obj["status"] == "aa_match" else "-1"
                elif match_dict[c] == "-1":
                    match_dict[c] = "1" if obj["status"] == "aa_match" else "-1"
                docid2canon[doc_id][r_canon] = True
                docid2pos[doc_id][pos] = True
        if glyco_ent:
            site_count += len(doc.keys())
            #print ("yyy", doc_id, entrez, direct, site_count)
    return match_dict, docid2pos, docid2canon

def load_site_dict(enteg_dir):

    docid2site = {}
    file_list = glob.glob(enteg_dir +  "/site.*.json")
    for in_file in file_list:
        doc_id = in_file.split(".")[-2]
        doc = json.loads(open(in_file, "r").read())
        docid2site[doc_id] = {}
        for site in doc:
            pos = site.split("|")[1]
            for obj in doc[site]:
                r_canon = get_r_canon(obj["canon"])
                s = "%s|%s" % (r_canon, pos)
                docid2site[doc_id][s] = True

    return docid2site


def load_samples_unlabeled(in_file):

    X, site_list  = [], []

    df = {}
    load_sheet(df, in_file, [], ",")
    f_list = df["fields"]
    row_list = []
    seen = {}
    for row in df["data"]:
        site = row[0]
        if site in seen:
            continue
        seen[site] = True
        row_list.append(row)
        site_list.append(site)
        newrow = []
        for v in row[1:-1]:
            newrow.append(float(v))
        X.append(numpy.array(newrow))
    X = numpy.array(X).astype(float)
    return X, site_list


def load_samples_labeled(in_file, exclude_site_list, cls):

    X, y  = [], []
    site_list = []
    
    df = {}
    load_sheet(df, in_file, [], ",")
    f_list = df["fields"]
    row_list = []
    for row in df["data"]:
        if row[-1] == cls:
            site = row[0]
            if site in exclude_site_list:
                continue
            row_list.append(row)
            site_list.append(site)
            newrow = []
            for v in row[1:-1]:
                newrow.append(float(v))
            X.append(numpy.array(newrow))
            y.append(int(row[-1]))

    X = numpy.array(X).astype(float)
    y = numpy.array(y).astype(int)
    return X, y, site_list







def load_curated_neg():

    seen = {}
    file_list = glob.glob("/data/shared/glycositeminer/userdata/curation/*")
    for in_file in file_list:
        doc = json.loads(open(in_file, "r").read())
        doc_id = doc["docid"]
        for site in doc["sitedict"]:
            if doc["sitedict"][site] == "not_glycosylated":
                pos = site.split("|")[-1].split("-")[1].strip()
                cmb = "%s|%s" % (doc_id, pos)
                seen[cmb] = True
    return seen






def load_target_count_dict(match_dict):

    tmp_dict = {}
    for cmb in match_dict:
        if match_dict[cmb] == "-1":
            continue
        doc_id, canon, pos = cmb.split("|")
        site = "%s|%s" % (doc_id, pos)
        if site not in tmp_dict:
            tmp_dict[site] = {}
        tmp_dict[site][canon] = True

    target_count_dict = {}
    for site in tmp_dict:
        n = len(tmp_dict[site].keys())
        target_count_dict[site] = n

    return target_count_dict


def get_otext2taxid_reduced_dict(misc_dir, ent_dir):

    out_dict = {}
    species_map = json.load(open(misc_dir + "species_map.json"))
    file_list = glob.glob(ent_dir + "/llm.*.json")
    seen = {}
    for in_file in file_list:
        doc = json.load(open(in_file))
        for obj in doc["glycosylation_sites"]:
            otext = "N/A"
            if "organism" not in obj:
                if "organism" in doc:
                    otext = doc["organism"]
            else:
                otext = obj["organism"]
            otext = "N/A" if otext in [None, ""] else otext
            otext = otext.lower()
            if otext in seen:
                continue
            tax_dict = {}
            for tax_id in species_map:
                tax_name_list = [name.lower() for name in species_map[tax_id]]
                if tax_id == "10090":
                    tax_name_list.append("mice")
                for name in tax_name_list:
                    if otext.find(name) != -1:
                        tax_dict[tax_id] = tax_name_list[0]
            if tax_dict != {}:
                out_dict[otext] = tax_dict
            seen[otext] = True

    return out_dict
 

    

def get_otext2taxid_dict(misc_dir, ent_dir):

    file_list = glob.glob(ent_dir + "/llm.*.json")
    doc_id2otext = {}
    idx = 0
    for in_file in file_list:
        idx += 1
        doc_id = in_file.split(".")[-2]
        doc = json.load(open(in_file))
        for obj in doc["glycosylation_sites"]:
            o_text = "N/A"
            if "organism" not in obj:
                if "organism" in doc:
                    o_text = doc["organism"]
            else:
                o_text = obj["organism"]
            o_text = "N/A" if o_text in [None, ""] else o_text
            o_text = o_text.lower()
            if doc_id not in doc_id2otext:
                doc_id2otext[doc_id] = {}
            doc_id2otext[doc_id][o_text] = True
        


    seen_otext_one = {}
    seen_otext_two = {}
    for doc_id in doc_id2otext:
        for o_text in doc_id2otext[doc_id]:
            seen_otext_one[o_text] = True
            w_list = o_text.split(" ")
            if len(w_list) > 1: 
                txt = w_list[0] + " " + w_list[1]
                seen_otext_two[txt] = True


    otext2taxid_three = json.load(open(misc_dir  + "curated_taxid_map.json"))
         
    in_file = misc_dir + "names.dmp"
    otext2taxid_one, otext2taxid_two = {}, {} 
    with open(in_file, "r") as FR:
        lcount = 0
        for line in FR:
            lcount += 1
            row = line.strip().split("|")
            tax_id, o_text = row[0].strip(), row[1].strip().lower()
            w_list = o_text.split(" ")
            if o_text in seen_otext_one:
                if o_text not in otext2taxid_one:
                    otext2taxid_one[o_text] = {}
                otext2taxid_one[o_text][tax_id] = True
            elif len(w_list) > 1:
                txt = w_list[0] + " " + w_list[1]
                if txt in seen_otext_two:
                    if txt not in otext2taxid_two:
                        otext2taxid_two[txt] = {}
                    otext2taxid_two[txt][tax_id] = True

    out_dict = {}
    for doc_id in doc_id2otext:
        for o_text in doc_id2otext[doc_id]:
            w_list = o_text.split(" ") 
            out_dict[o_text] = {"A":[], "B":[], "C":[], "D":[]}
            tax_id_list = []
            if o_text in otext2taxid_one:
                tax_id_list = list(otext2taxid_one[o_text].keys())
                out_dict[o_text]["A"] = tax_id_list
            elif len(w_list) > 1:
                txt = w_list[0] + " " + w_list[1]
                if txt in otext2taxid_two:
                    tax_id_list = list(otext2taxid_two[txt].keys())
                    if len(tax_id_list) > 10:
                        tax_id_list = tax_id_list[:10]
                    out_dict[o_text]["B"] = tax_id_list
            if tax_id_list == []:
                if o_text in otext2taxid_three:
                    tax_id_list = list(otext2taxid_three[o_text].keys())
                    out_dict[o_text]["C"] = tax_id_list


    return out_dict



def load_pdb_range_dict():

    range_dict = {}
    ptrn = " <http://purl.uniprot.org/core/chain> "
    file_list = glob.glob("/data/shared/glygen/downloads/ebi/current/uniprot-proteome-*.nt")
    for in_file in file_list:
        cmd = "grep \"%s\" %s" % (ptrn, in_file)
        line_list = subprocess.getoutput(cmd).split("\n")
        for line in line_list:
            row = line.strip().split(" ")
            #<http://purl.uniprot.org/isoforms/J3QS39-1#PDB_3PRP_tt1tt75>
            canon = row[0].split("/")[-1].split("#")[0]
            chain = row[0].split("/")[-1].split("#")[1].split("_")[1]
            r = row[-2].split("=")[-1].replace("\"","")
            s, e = r.split("-")
            if s.isdigit() and e.isdigit():
                # we care only about those that are shifted
                if int(s) > 1:
                    c = chain + "-" + r
                    if canon not in range_dict:
                        range_dict[canon] = {}
                    range_dict[canon][c] = True

    return range_dict





def parse_one_word(w, aa_dict):


    obj_list = []
    new_site = "Xxx"
    flag = "1x"
    if len(w) > 1:
        if w[0].lower() in aa_dict["one"] and w[1].isdigit():
            d = ""
            for j in range(1, len(w)):
                if w[j].isdigit() == False:
                    break
                d += w[j]
            new_site = aa_dict["one2three"][w[0].upper()] + "-" + d
            flag = "1a"
            obj_list.append({"newsite":new_site, "flag":flag})
    if len(w) > 3:
        if w[0:3].lower() in aa_dict["three"] and w[3].isdigit():
            d = ""
            for j in range(3, len(w)):
                if w[j].isdigit() == False:
                    break
                d += w[j]
            new_site = w[0].upper() + w[1:3].lower() + "-" + d
            flag = "1b"
            obj_list.append({"newsite":new_site, "flag":flag})

    for aa in aa_dict["full"]:
        aa = aa[0].upper() + aa[1:].lower()
        if w.lower().find(aa.lower()) != -1:
            d = w.lower().replace(aa.lower(), "")
            if d.isdigit():
                new_site = aa_dict["full2three"][aa] + "-" + d
                flag = "1c"
                obj_list.append({"newsite":new_site, "flag":flag})
    if obj_list == []:
        obj_list.append({"newsite":new_site, "flag":flag})

    seen = []
    for obj in obj_list:
        seen.append(json.dumps(obj))

    res = []
    for s in list(set(seen)):
        res.append(json.loads(s))

    return res



def parse_two_words(w, w_next, aa_dict):

    w, w_next = w.lower(), w_next.lower()
    
    
    obj_list = []
    new_site = "Xxx"
    rflag = False
    f = w in aa_dict["one"] or w in aa_dict["three"] or w in aa_dict["full"]
    if f:
        aa = ""
        aa = w[0].upper() if len(w) == 1 else aa
        aa = w[0].upper() + w[1:].lower() if len(w) >= 3 else aa
        aa = aa_dict["one2three"][aa] if len(aa) == 1 else aa
        aa = aa_dict["full2three"][aa] if len(aa) > 3 else aa
        if w_next.isdigit():
            new_site = aa + "-" + w_next
            obj_list.append({"flag":"2a", "newsite":new_site})
            rflag = True
        else:
            for p in w_next.split("/"):
                if p.isdigit():
                    new_site = aa + "-" + p
                    obj_list.append({"flag":"2b", "newsite":new_site})
                    rflag = True
    if rflag == False:
        o_list = parse_one_word(w, aa_dict)
        for o in o_list:
            if o["flag"] in ["1a","1b", "1c"]:
                obj_list.append({"flag":"2c_%s" % (o["flag"]), "newsite":o["newsite"]})
                rflag = True
        o_list = parse_one_word(w_next, aa_dict)
        for o in o_list:
            if o["flag"] in ["1a","1b", "1c"]:
                obj_list.append({"flag":"2d_%s" % (o["flag"]), "newsite":o["newsite"]})
                rflag = True
    if rflag == False:
        obj_list.append({"flag":"2x", "newsite":new_site})


    seen = []
    for obj in obj_list:
        seen.append(json.dumps(obj))

    res = []
    for s in list(set(seen)):
        res.append(json.loads(s))

    return res


def transform_llm_site(orig_site, aa_dict):

    site_list = []
    if orig_site in [None, "", "N/A"]:
        return site_list

    if no_digits(orig_site):
        return  []

    site = orig_site.strip()
    site = site.replace("-", " ").replace("(", " ").replace(")", " ").replace("'", "").replace(",", "")
    site = site.strip()
    w_list = site.split(" ")
    if w_list == [""]:
        return []


    obj_list = []
    if len(w_list) == 1:    
        w = w_list[0]
        o_list = parse_one_word(w, aa_dict)
        for o in o_list:
            obj_list.append(o)
    elif len(w_list) == 2:
        w, w_next = w_list[0].lower(), w_list[1].lower()
        o_list = parse_two_words(w, w_next, aa_dict)
        for o in o_list:
            obj_list.append(o)
    else:
        w_list.append("*****")
        for i in range(0, len(w_list)-1):
            w, w_next = w_list[i].lower(), w_list[i+1].lower()
            o_list = parse_two_words(w, w_next, aa_dict)
            for o in o_list:
                obj_list.append(o)
   
    site_list = [] 
    for obj in obj_list:
        if obj["newsite"] == "Xxx":
            continue
        site_list.append(obj["newsite"])

    return list(set(site_list))





def load_genename_dict(base_dir):

    seen = {}
    file_list = glob.glob(base_dir + "glygen/*_protein_genenames_uniprotkb.csv")
    for in_file in file_list:
        sp = in_file.split("/")[-1].split("_")[0]
        data_frame = {}
        load_sheet(data_frame, in_file, [], ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            for f in ["gene_symbol_recommended","gene_symbol_alternative","orf_name"]:
                gene_text = row[f_list.index(f)].strip()
                gene_text_lc = gene_text.lower()
                for txt in [gene_text, gene_text_lc]:
                    if txt == "":
                        continue
                    if txt not in seen:
                        seen[txt] = {}
                    seen[txt][canon] = sp

    file_list = glob.glob(base_dir + "glygen/*_protein_genenames_refseq.csv")
    for in_file in file_list:
        sp = in_file.split("/")[-1].split("_")[0]
        data_frame = {}
        load_sheet(data_frame, in_file, [], ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            for f in ["refseq_gene_name"]:
                for gene_text in row[f_list.index(f)].split(";"):
                    gene_text = gene_text.strip()
                    gene_text_lc = gene_text.lower()
                    for txt in [gene_text, gene_text_lc]:
                        if txt == "":
                            continue
                        if txt not in seen:
                            seen[txt] = {}
                        seen[txt][canon] = sp

    return seen







def load_proteinname_dict(base_dir):

    field_map = {
        "*_protein_recnames.csv":[
            "uniprotkb_canonical_ac","recommended_name_full","recommended_name_short","ec_name"
        ],
        "*_protein_altnames.csv":[
            "uniprotkb_canonical_ac","alternative_name_full","alternative_name_short","ec_name"
        ],
        "*_protein_submittednames.csv":[
            "uniprotkb_canonical_ac","submitted_name_full","submitted_name_short","ec_name"
        ],
        "*_protein_proteinnames_refseq.csv":[
            "uniprotkb_canonical_ac","refseq_ac","refseq_protein_name"
        ],
    }

    seen = {}
    for ptrn in field_map:
        file_list = glob.glob(base_dir + "glygen/%s" % (ptrn))
        for in_file in file_list:
            sp = in_file.split("/")[-1].split("_")[0]
            data_frame = {}
            load_sheet(data_frame, in_file, [], ",")
            f_list = data_frame["fields"]
            for row in data_frame["data"]:
                canon = row[f_list.index("uniprotkb_canonical_ac")]
                for f in field_map[ptrn]:
                    for protein_text in row[f_list.index(f)].strip().split(";"):
                        protein_text = protein_text.strip()
                        protein_text_lc = protein_text.lower()
                        for txt in [protein_text, protein_text_lc]:
                            if txt == "":
                                continue
                            if txt not in seen:
                                seen[txt] = {}
                            seen[txt][canon] = sp


    return seen





def load_proteinname_dict_by_sp(base_dir):

    field_map = {
        "*_protein_recnames.csv":[
            "uniprotkb_canonical_ac","recommended_name_full","recommended_name_short","ec_name"
        ],
        "*_protein_altnames.csv":[
            "uniprotkb_canonical_ac","alternative_name_full","alternative_name_short","ec_name"
        ],
        "*_protein_submittednames.csv":[
            "uniprotkb_canonical_ac","submitted_name_full","submitted_name_short","ec_name"
        ],
        "*_protein_proteinnames_refseq.csv":[
            "uniprotkb_canonical_ac","refseq_ac","refseq_protein_name"
        ],
    }

    seen = {}
    for ptrn in field_map:
        file_list = glob.glob(base_dir + "glygen/%s" % (ptrn))
        for in_file in file_list:
            sp = in_file.split("/")[-1].split("_")[0]
            if sp not in seen:
                seen[sp] = {}
            data_frame = {}
            load_sheet(data_frame, in_file, [], ",")
            f_list = data_frame["fields"]
            for row in data_frame["data"]:
                canon = row[f_list.index("uniprotkb_canonical_ac")]
                for f in field_map[ptrn]:
                    for val in row[f_list.index(f)].strip().split(";"):
                        val = val.strip()
                        if val == "":
                            continue
                        if val not in seen[sp]:
                            seen[sp][val] = {}
                        seen[sp][val][canon] = sp

    return seen






def dump_pubtator_entities(doc_id, base_dir):

    gene_obj_list = []
    species_obj_list = []
    in_file = base_dir + "pubtator_extracts/pmid.%s.txt" % (doc_id)
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
        out_file = base_dir + "pubtator/entities/gene.%s.json" % (doc_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(gene_dict, indent=4)))

    if len(species_obj_list) > 0:
        species_dict = {combo_id:species_obj_list}
        out_file = base_dir + "pubtator/entities/species.%s.json" % (doc_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(species_dict, indent=4)))

    return


def load_llm_sites(rel_dir, misc_dir):


    ent_dir = rel_dir + "llm_entities/*/"
    otext2taxid_dict = get_otext2taxid_reduced_dict(misc_dir, ent_dir)
    glygen_species_dict = get_species_dict(misc_dir)
    pubtator_docid2taxid = load_docid2taxid(rel_dir)
    
    aa_dict = get_aa_dict(misc_dir)
    tmp_dict = {}
    file_list = glob.glob(rel_dir + "llm_entities/*/*.json")
    for in_file in file_list:
        doc_id = in_file.split(".")[-2]
        doc = json.load(open(in_file))
        seen_site_text = {}
        for obj in doc["glycosylation_sites"]:
            site_text = obj["glycosylation_site"] if "glycosylation_site" in obj else ""
            otext = "N/A"
            if "organism" not in obj:
                if "organism" in doc:
                    otext = doc["organism"]
            else:
                otext = obj["organism"]
            otext = "N/A" if otext in [None, ""] else otext
            otext = otext.lower()
            tax_id_dict = {}
            if otext in otext2taxid_dict:
                for tax_id in otext2taxid_dict[otext]:
                    tax_id_dict[tax_id] =  otext2taxid_dict[otext][tax_id]
            
            if tax_id_dict == {}:
                if doc_id in pubtator_docid2taxid:
                    for tax_id in pubtator_docid2taxid[doc_id]:
                        if tax_id not in glygen_species_dict:
                            continue
                        tax_id_dict[tax_id] = glygen_species_dict[tax_id]
                        tax_name_src = "pubtator"

            if site_text != "" and tax_id_dict != {}:
                seen_site_text[site_text] = True
        
        for site_text in seen_site_text:
            site_list = transform_llm_site(site_text, aa_dict)
            for site in site_list:
                pos = site.split("-")[-1]
                cmb = "%s|%s" % (doc_id, pos)
                tmp_dict[cmb] = True
                #print (cmb, tax_id_dict)

    return tmp_dict





def download_file(url: str, save_path: str, chunk_size: int = 8192) -> bool:
    try:
        with requests.get(url, stream=True, timeout=30) as response:
            response.raise_for_status()
            with open(save_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk: 
                        file.write(chunk)
            return True
    except requests.exceptions.RequestException as e:
        if os.path.exists(save_path):
            os.remove(save_path)
        return False



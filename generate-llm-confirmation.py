import sys,os
import json
import util
import glob
from openai import OpenAI
from optparse import OptionParser
import subprocess



def load_site2genetext(data_dir):

    tmp_dict = {}
    tmp_dict["33608772|P16960-1|268"] = {"CRC":True} 
    file_list = glob.glob(data_dir + "integrated/*.json")
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
                if cmb not in tmp_dict:
                    tmp_dict[cmb] = {}
                tmp_dict[cmb][gene_text] = True


    return tmp_dict




def main():

    home_dir = os.path.expanduser('~')
    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())
    client = OpenAI(api_key=config_obj["llm"]["api_key"])
   
     # create dir if doesn't exist    
    if os.path.isdir(config_obj["data_dir"] + "confirmation/") == False:
        cmd = "mkdir -p " + config_obj["data_dir"] + "confirmation/"
        x = subprocess.getoutput(cmd)


 
    model = config_obj["llm"]["model"]
    role = config_obj["llm"]["role"]
    #res_schema = config_obj["llm"]["schema"]
    res_schema = {
        "type": "object",
        "properties": {
            "answer": {
                "type": "string",
                "enum": ["yes", "no"]
            },
            "evidence":{
                "type": "string"
            }
        }
    }

    in_dict = {}
    species_dict = json.load(open(config_obj["misc_dir"]  + "species_common_name.json"))
    seq_dict, canon2taxid = util.load_seq_dict(config_obj["data_dir"], config_obj["misc_dir"] )
    site2genetext = load_site2genetext(config_obj["data_dir"])
    line_list = open(config_obj["data_dir"]+ "predicted/predicted_tmp.csv", "r").read().split("\n")
    for line in line_list[1:-1]:
        row = line.replace("\"","").strip().split(",")
        doc_id, canon, pos = row[0],row[1],row[2]
        tax_id = canon2taxid[canon]
        tax_name = species_dict[tax_id]
        glygen_status = row[-3]
        llm_flag = row[-1]
        if llm_flag == "no_llm":
            continue
        site = "%s|%s|%s" % (doc_id, canon, pos)
        for gene_text in site2genetext[site]:
            o = {"canon":canon, "tax_name":tax_name, "pos":pos, "gene_text":gene_text, 
                    "glygen_status":glygen_status}
            if doc_id not in in_dict:
                in_dict[doc_id] = []
            in_dict[doc_id].append(o)

    query_tmpl = 'tell me if %s "%s" is glycosylated at position "%s" using only using information provided in the following text'
    #query_tmpl = 'give me answer and evidence if "%s" is glycosylated at position "%s" using only using information provided in the following text'

    for doc_id in in_dict:
        out_file = config_obj["data_dir"]  + "/confirmation/pmid.%s.json" % (doc_id)
        if os.path.isfile(out_file) == True:
            continue
        abstract_file = config_obj["data_dir"] + "medline_abstracts/pmid.%s.txt" % (doc_id)
        abstract_text = open(abstract_file, "r").read()
        abstract_text = abstract_text.strip()
        obj_list = []
        for o in in_dict[doc_id]:
            query = query_tmpl % (o["tax_name"], o["gene_text"], o["pos"])
            #query = query_tmpl % (o["gene_text"], o["pos"])
            prompt = query + ": " + abstract_text
            res = client.chat.completions.create(
                model=model,
                store=True,
                messages=[ {"role": role, "content": prompt} ],
                response_format={
                    "type": "json_schema",
                    "json_schema": {
                        "name": "site_response",
                        "schema": res_schema,
                        "strict": False
                    }
                }
            )
            res_json = json.loads(res.choices[0].message.content)
            o["answer"] = res_json["answer"]
            o["evidence"] = res_json["evidence"]
            obj_list.append(o)

        out_file = config_obj["data_dir"]  + "/confirmation/pmid.%s.json" % (doc_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj_list, indent=4)))


if __name__ == '__main__':
    main()


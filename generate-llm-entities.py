import sys,os
import json
import util
import glob
from openai import OpenAI
from optparse import OptionParser
import subprocess



def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-b","--batch",action="store",dest="batch",help="1/2/3/4 ...")

    (options,args) = parser.parse_args()
    for key in ([options.batch]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    batch = options.batch


    home_dir = os.path.expanduser('~')

    global config_obj
    config_obj = json.loads(open("conf/config.json", "r").read())

    model = config_obj["llm"]["model"]
    role = config_obj["llm"]["role"]
    query = config_obj["llm"]["query"]
    res_schema = config_obj["llm"]["schema"]
   
    doc_id_list = []
    file_list = glob.glob(config_obj["data_dir"] + "/medline_abstracts/pmid.*.txt")
    for in_file in file_list:
        doc_id = in_file.split(".")[-2]
        doc_id_list.append(doc_id)

    out_dir = config_obj["data_dir"] + "llm_entities/%s/" % (batch)
    if os.path.isdir(out_dir) == False:
        cmd = "mkdir -p  " + out_dir
        x = subprocess.getoutput(cmd)
    
    client = OpenAI(api_key=config_obj["llm"]["api_key"])
    
    for doc_id in doc_id_list:
        abstract_file = config_obj["data_dir"] + "medline_abstracts/pmid.%s.txt" % (doc_id)
        if os.path.isfile(abstract_file) == False:
            continue
        out_file = out_dir + "llm.%s.json" % (doc_id)
        if os.path.isfile(out_file) == True:
            continue
        abstract_text = open(abstract_file, "r").read()
        abstract_text = abstract_text.strip()
        prompt = query + ":" + abstract_text

        res = client.chat.completions.create(
            model=model,
            store=True,
            messages=[ {"role": role, "content": prompt} ],
            #response_format= { "type":"json_object" }
            response_format={
                "type": "json_schema",
                "json_schema": {
                    "name": "site_response",
                    "schema": res_schema,
                    "strict": False
                }
            }
        )
        out_doc = json.loads(res.choices[0].message.content)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(out_doc, indent=4)))
    

    



if __name__ == '__main__':
    main()



import os
import sys
import gzip
import json
import glob
import pubmed_parser as pp
import re

import spacy
from optparse import OptionParser
import util
import subprocess


def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-s","--start",action="store",dest="start",help="")
    parser.add_option("-e","--end",action="store",dest="end",help="")

    (options,args) = parser.parse_args()
    for key in ([options.start, options.end]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    global pmid_dict

    config_obj = json.loads(open("conf/config.json", "r").read())

 
    start_idx = int(options.start)
    end_idx = int(options.end)
    xml_dir = config_obj["data_dir"] + "medline_xml/"
    out_dir_one = config_obj["data_dir"] + "medline_extracts/"
    out_dir_two = config_obj["data_dir"] + "medline_abstracts/"
    log_dir = config_obj["data_dir"] + "logs/"
    log_file = log_dir + "medline_abstracts.%s.%s.log" % (start_idx, end_idx)

    DEBUG = False
    #DEBUG = True
    if DEBUG:
        file_list = [xml_dir  + "/pubmed23n0349.xml.gz"]
        debug_doc_id_list = ["10531415"]
        out_dir_one = "tmp/medline_extracts/"
        out_dir_two = "tmp/medline_abstracts/"

    if os.path.isdir(out_dir_one) == False:
        cmd = "mkdir -p " + out_dir_one
        x = subprocess.getoutput(cmd)
    if os.path.isdir(out_dir_two) == False:
        cmd = "mkdir -p " + out_dir_two
        x = subprocess.getoutput(cmd)
    if os.path.isdir(log_dir) == False:
        cmd = "mkdir -p " + log_dir
        x = subprocess.getoutput(cmd)



    # Using EntityRuler
    nlp = spacy.blank("en")
    patterns_json_file = "misc/glyco.json"
    ruler = nlp.add_pipe("entity_ruler")
    pattern_obj_list = json.loads(open(patterns_json_file, "r").read())
    ruler.add_patterns(pattern_obj_list)
    nlp.add_pipe('sentencizer')

    SPACE_PATTERN = re.compile(r'[\s][\s]+')
    with open(log_file, "w") as FL:
        FL.write("Started logging\n")

    file_list = sorted(glob.glob(xml_dir + "/*.xml.gz"))
    aa_dict = util.get_aa_dict("misc/")
    end_idx = len(file_list) if end_idx > len(file_list) else end_idx

    for gz_xml_file in file_list[start_idx-1:end_idx]:
        if not gz_xml_file.endswith('.gz'):
            continue
        file_name = gz_xml_file.split("/")[-1].replace(".xml.gz", "")
        dicts_out = pp.parse_medline_xml(gz_xml_file)
        n_found = 0
        idx = 0
        for obj in dicts_out:
            if "pmid" not in obj:
                continue
            if obj["pmid"] == "":
                continue
            if "abstract" not in obj:
                continue
            if obj["abstract"] == "":
                continue
            doc_id = obj["pmid"]

            if DEBUG:
                if doc_id not in debug_doc_id_list:
                    continue

            if idx%1000 == 0:
                with open(log_file, "a") as FL:
                    FL.write("checked %s PMIDs, %s found from %s\n" % (idx, n_found, gz_xml_file))
            idx += 1
            
            title_text, abstract_text = obj["title"], obj["abstract"]
            mesh = obj["mesh_terms"]
            publication_types = obj["publication_types"]
            text = "{} {}".format(title_text, abstract_text).strip()
            text = re.sub(SPACE_PATTERN, ' ', text)
            doc = {"docId":doc_id, "text":text} 
            if len(title_text) > 0:
                doc["title"] = {"charStart": 0, "charEnd": len(title_text) - 1}
            if len(publication_types) > 0:
                doc["type"] = publication_types
            if len(mesh) > 0:
                doc["mesh"] = mesh
            #sent_list = sent_tokenize(doc["text"])
            tmp_doc = nlp(doc["text"])
            sent_list = [sent.text.strip() for sent in tmp_doc.sents]
            sent_obj_list = []
            doc2lbl = {}
            sent_idx = 0
            for sent in sent_list:
                sent_idx += 1
                sent = util.transform_sent(sent, aa_dict)
                obj = nlp(sent)
                entities = []
                seen_lbl = {}
                for ent in obj.ents:
                    lbl, txt = ent.label_, ent.text
                    entities.append({"start":ent.start_char, "end":ent.end_char, "label":ent.label_, "text":ent.text})
                    doc2lbl[ent.label_] = True
                    seen_lbl[ent.label_] = True
                if "SITE" in seen_lbl or "GLYCOSYLATION" in seen_lbl:
                    sent_obj_list.append({"sentence":sent, "sentence_idx": sent_idx, "entities": entities})
            if len(sent_obj_list) > 0 and "SITE"  in doc2lbl and "GLYCOSYLATION" in doc2lbl:
                n_found += 1
                out_doc = {"doc_id":doc_id, "sent_list":sent_obj_list}
                out_file = out_dir_one + "pmid.%s.json" % (doc_id)
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(out_doc, indent=4)))
                out_file = out_dir_two + "pmid.%s.txt" % (doc_id)
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (doc["text"]))
                if DEBUG:
                    print ("created file:%s\n" % (out_file))
        with open(log_file, "a") as FL:
            FL.write("final %s found from %s\n" % (n_found, gz_xml_file))
       
    
    cmd = "chmod -R 777 reldir/medline_extracts"
    x = subprocess.getoutput(cmd)
    cmd = "chmod -R 777 reldir/medline_abstracts" 
    x = subprocess.getoutput(cmd)
    cmd = "chmod -R 777 logs/*"
    x = subprocess.getoutput(cmd)
 
    return




if __name__ == '__main__':
    main()


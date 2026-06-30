### Setting up venv
After you clone this repo, follow these commands to setup your virutal 
environment.
```
$ python3 -m venv venv
$ source venv/bin/activate
$ pip3 install -r requirements.txt
```


### Step-1: Data download
After cloning this repository, edit the file conf/config.json to 
set your data_dir (for example /data/glycositeminer/). All data used
in this pipeline will be stored in that folder.

The following commands will download data to "$data_dir/medline_xml/",
"$data_dir/pubtator_downloads/", "$data_dir/glygen/", and "$data_dir/gene_info/",
respectively. Run one command at a time and monitor the corresponding download directory
to make sure files are being downloaded.
```
$ nohup python3 download-glygen.py &
$ nohup python3 download-medline.py &
$ nohup python3 download-pubtator.py &
$ nohup python3 download-gene-info.py &
```
The download-medline.py script depends on the file conf/medline.json which needs
to be edited. To adjust the values for the "baseline" section, go to 
https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/ and see the indexes of the 
first and last pubmedxxx.xml.gz files. For example, if the first and last files are
"pubmed25n0001.xml.gz" and "pubmed25n1274.xml.gz", then the base line section
in conf/medline.json should be:
```	
{
	"baseline":{
		"url":"https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/",
		"start":1,
		"end":1274,
     	 "year":26
	}
}
```	
The "updatefiles" section should be updated in the same manner by going to 
https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/

### Step-2: Extract known sites
The following script creates data_dir/glygen/known_sites.csv which
contains known glycosylation sites assuming XX_proteoform_glycosylation_sites?.csv
datasets are created downloaded under data_dir/glygen/.
```
$ nohup python3 extract-known-sites.py &
```


### Step-3: Processing medline downloads
The following script is a wrapper around "extract-medline-data.py" for extracting
data from downloaded medline abstracts in xml format.
```
$ nohup python3 wrap-extract-medline-data.py &

Input:
	data_dir/medline_xml/*
	misc/glyco.json

Output:
	data_dir/medline_extracts/*
	data_dir/medline_abstracts/*
	data_dir/logs/medline_abstracts.%s.%s.log
```


### Step-4: Processing pubtator downloads
```	
$ nohup python3 extract-pubtator-data.py &
  
Input:
	data_dir/pubtator_downloads/bioconcepts2pubtatorcentral.offset.gz

Output:
	data_dir/pubtator_extracts/*
	data_dir/logs/pubtator_extracts.log
```

### Step-5: Making LLM API calls
This step uses openai.com API to extract species, gene, site entities from
abstracts in data_dir/medline_abstracts/ and writes out entities JSON
files under data_dir/llm_entities/$BATCH/ where $BATCH is batch index (1, 2, 3).
You can run this for more than one batches if you think the LLM output will vary. You will 
need to edit the "llm" section of the "conf/config.json" to set your LLM access.
```
$ python3 generate-llm-entities.py  -b $BATCH
Input
	data_dir/medline_abstracts/*
Output
	data_dir/llm_entities/$BATCH/llm.*.json
```

### Step-6 Making entities
```   
$ python3 make-entities.py
Input
	data_dir/medline_extracts/*.json
	data_dir/pubtator_extracts/*.json
	data_dir/glygen/*names*.csv
	data_dir/glygen/*genelocus*.csv
	data_dir/glygen/*geneid*.csv
	data_dir/llm_entities/*/*.json
	data_dir/misc/aa.csv
	data_dir/misc/known_sites.csv
Output
	data_dir/entities/site.*.json
	data_dir/entities/glyco.*.json
	data_dir/entities/species.*.json
	data_dir/entities/gene.*.json
	data_dir/entities/extragene.*.json
```


### Step-7: Integrate entities
```         
$ python3 integrate-entities.py
Input
	data_dir/entities/site.*.json
	data_dir/entities/glyco.*.json
	data_dir/entities/species.*.json
	data_dir/entities/gene.*.json
	data_dir/entities/extragene.*.json
	data_dir/gene_info/All_Data.gene_info
	data_dir/misc/canon_map.csv 
	data_dir/misc/species_map.json
Output
	data_dir/integrated/site.*.json
```

### Step-8: Making match sites
```
$ python3 make-match-sites.py
Input
	data_dir/integrated/site.*.json
	data_dir/misc/known_sites.csv
	data_dir/glygen/*_protein_canonicalsequences.fasta
	data_dir/misc/species_map.json
	data_dir/misc/aa.csv
	data_dir/misc/glyco_aa_list.json
Output
	data_dir/match_sites/match_sites.csv
	data_dir/match_sites/mismatch_sites.csv
```


### Step-9: Making samples/features
```
$ python3 make-samples.py
Input
	data_dir/match_sites/sites.csv
	data_dir/glygen/*protein*canonicalsequences.fasta
	data_dir/misc/species_map.json
	data_dir/pubtator_entities/species.*.json
Output
	data_dir/samples/samples_labeled.csv
	data_dir/samples/samples_all.csv
```

### Step-10: Making predicdtions
```	  		
$ python3 make-predictions.py 
Input
	data_dir/samples/samples_all.csv
	data_dir/glygen/*_protein_canonicalsequences.fasta
	data_dir/misc/known_sites.csv
	data_dir/misc/aa.csv
	data_dir/misc/species_map.json
Output
	data_dir/predicted/predicted_tmp.csv
	data_dir/predicted/stats.txt
```


### Step-11: Generating LLM confirmations
```		
$ python3 generate-llm-confirmation.py
Input
	data_dir/predicted/predicted_tmp.csv
Output
	data_dir/confirmations/*.json
```

### Step-12: Updating predicdtions
```
$ python3 update-predictions.py
Input
	data_dir/predicted/predicted_tmp.csv
	data_dir/confirmations/*.json
Output
	data_dir/predicted/predicted.csv
```

### Step-13: Correcting for species mismatch
```         
$ python3 correct-predictions.py
Input
	data_dir/predicted/predicted.csv
	data_dir/llm_entities/*/*.json
Output
	data_dir/predicted/predicted_corrected.csv
```



	



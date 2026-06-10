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
	



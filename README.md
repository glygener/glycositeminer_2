### Step-1: Data download
After cloning this repository, edit the file conf/config.json to 
set your data_dir (for example /data/glycositeminer/). All data used
in this pipeline will be stored in that folder.

	The following download scripts will download data and place it under
	reldir/glygen/, reldir/medline_xml/, reldir/pubtator_downloads/
	and reldir/gene_info/. This assumes that all protein datasets have been
	created under /data/projects/glygen/generated/datasets/unreviewed/

	The download-medline.py script depends on the file conf/medline.json which needs
	to be edited. To adjust the values for the "baseline" section, go to 
	https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/ and see the indexes of the 
	first and last pubmed*.xml.gz files. For example, if the first and last files are
	"pubmed25n0001.xml.gz" and "pubmed25n1274.xml.gz", then the base line section
	in conf/medline.json should be:
	{
		"baseline":{
			"url":"https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/",
			"start":1,
			"end":1274,
     	 	"year":26
		}
	}
	The "updatefiles" section should be updated in the same manner by going to 
	https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/
	

	$ python3 download-glygen.py
	$ python3 download-medline.py	
	$ python3 download-pubtator.py
	$ python3 download-gene-info.py



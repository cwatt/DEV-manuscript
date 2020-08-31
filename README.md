The scripts in this directory can be used to recreate the results and figures from my paper, **Agricultural Management Affects Root-Associated Microbiome Recruitment over Maize Development (2019)**. Link: https://apsjournals.apsnet.org/doi/10.1094/PBIOMES-03-19-0016-R

Notes:

* The otu count and taxonomies were generated via the hundo pipeline created by Joe Brown (https://github.com/pnnl/hundo) 
* 16S and ITS refer to prokaryotic and fungal data respectively.
* There is an error in USEARCH which caused one of the fungal OTUs in the .txt file to be missing an OTU label. 
See "k__Fungi,p__Basidiomycota,c__Agaricomycetes,o__Agaricales,f__Inocybaceae,g__Crepidotus,s__?" in ITS_OTU.txt.
This needs to be corrected by entering OTU_2583 for that taxa, as it has been in the supplied DEVITS_otu.csv file.
* Soil and Compartment refer to System and Root proximity respectively.

Feel free to email me cjwattenburger@gmail.com if you have any questions.

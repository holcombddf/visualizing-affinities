# visualizing-affinities
This software allows visualization of the potential for the development of antibodies to certain amino acid sequences of a protein in the form of heatmaps and risk populations using your locally installed netMHCIIpan[http://www.cbs.dtu.dk/services/NetMHCIIpan/].  Additionally, by allowing user input, this program also allows for the evaluation of novel protein sequences introduced by the user in the course of protein engineering.  In addition, the software allows the comparison of potential immunogenicity due to the changes made in the protein.  The software is indicative, but not definitive and any developed protein with lower potential for antibody development will need to be tested and reviewed by FDA prior to marketing.  

Please cite :   
Pandey GS, Yanover C, Howard TE, Sauna ZE (2013) Polymorphisms in theF8 Gene and MHC-II Variants as Risk Factors for the Development of Inhibitory Anti-Factor VIII Antibodies during the Treatment of Hemophilia A: A Computational Assessment. PLoS Comput Biol 9(5): e1003066. doi:10.1371/journal.pcbi.1003066

Questions regarding this script to tranform netMHCIIpan results into graphics should be sent to joseph.mcgill@fda.hhs.gov, zuben.sauna@fda.hhs.gov, or holcombddf@gmail.com (for the forked version).

## Affinity_Visualization.py
Contains all methods for the class "Affinity_Visualization". Contains relevant attributes and methods for producing an affinity visualization.

## affinity_viz.py
Using the parameters.py file and importing the class Affinity_Visualization, this script will create an instance of Affinity_Visualization and create the heatmaps and promiscuity plots.

## automate_affinity_viz.py
Separates fasta files containing multiple sequences into component sequences, creating multiple fasta files, then runs affinity_viz.py on all of them (not necessary, just easier to run).

## parameters.py
Contains all the information for this particular run including file locations and titles. This is the only file that should need to be changed for a standard running of these scripts.

## Nworld.csv
This file is not necessary but a file matching the format is and must be referenced in the parameters.py file as allele_csv.

## World.txt
This file is not necessary but a file matching the format is and must be referenced in the parameters.py file as HLA_file_path
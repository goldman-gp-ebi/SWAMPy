Availability of the data used in the publication:

- Nimagen and synthetic ARTIC amplicon sequences are available at the European Nucleotide Archive, BioProject accession PRJEB53222. 
- Real wastewater samples used for parameter estimation, sequenced using the ARTIC v3 protocol are available on request and with permission from the Joint Biosecurity Centre and University of Liverpool.
- Summarised data (e.g. amplicon counts) are available through GitHub at https://github.com/goldman-gp-ebi/SWAMPy/tree/main/supplementary_files.
- The samples used for validation and comparison with ww_simulations are available at ENA BioProject accession PRJNA796340, with the exact samples used listed in the supplementary information file.
- Non-default SWAMPy options used to generate the 73 simulated time point data used in the manuscript:
```  
-ins 0.0002
-del 0.00115
-subs 0.005
-rins 0.0002
-rdel 0.00115
-subs 0.005
--amplicon distribution dirichlet 2
--amplicon pseudocounts 200
```

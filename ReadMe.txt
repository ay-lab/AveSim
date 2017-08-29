AveSim is a pipeline to simulate Hi-C matrices involving CNVs and translocations. 
AveSim can perform the simulation at the restriction site resolution level (RE_level_simulation folder) or at fixed-size bin resolution level (Fixed_size_simulation folder).

RE_level_simulation folder is subdivided into following subfolder
	A. Amplification
		> scripts
		> simulation
	B. Deletion
                > scripts
                > simulation
	C. Translocation
                > scripts
                > simulation
	D. DoubleMinute
                > scripts
                > simulation

Users can generate the simulated paired-end fastq files in the simulation folder. These fastq files can be processed by any standard Hi-C data processing pipeline. 

Fixed_size_simulation folder is subvidived into following subfolder
	A. data
	B. InputFiles
	C. TranslocationScripts
	D. Non-Reciprocal_Balanced_Translocation
	E. Reciprocal_Balanced_Translocation  
	F. Unbalanced_Translocation

Users can perform different translocation simulations (D,E and F) under the respective simulation directories. Fixed-size simulation will create 2D Hi-C matrices.

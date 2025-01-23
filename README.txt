"Systematic identification of allosteric effectors in Escherichia coli metabolism"
This repository holds code for the analysis of progress curves from in vitro enzyme assays on 6520 mass spectrometer (as well as analysis of validation experiments on photometer), and code for the generation of key figures.

I. HOW TO RUN THE MATLAB SCRIPTS
1) Download GitHub repository and move to MATLAB folder.
2) Download preprocessed mass spectrometry data from 10.5281/zenodo.14130718 and copy into folder "=1_input=preprocessed_MSdata_deposited_on_Zenodo"
3) Download MATLAB results files from Zenodo 10.5281/zenodo.14226986 and copy them into folder "=1_output=large_MATfiles_deposited_on_ZenodoO".
4) Inputs and outputs of each script are saved in separate folders, allowing running each one independentely.
4a) For analysis and plotting of MS-based data, open =1=MS_analysis_script.m in MATLAB, change analysis path to match the downloaded repository, enter name of the enzyme of interest and execute the code.
4b) For analysis of photometer-based data, follow same approach using =2=Photometer_analysis_script.m. 
4c) For comparisons to literature and/or between photometer and MS-based measurements, refer to =3=comparative_analysis.
4d) For similarity comparisons and various other subplots, refer to =4=additional_analyses

II. PYTHON SCRIPTS
For the script corresponding to Supplementary Text S1, refer to "=S1=Enzyme_Assay_Simulations.ipynb". Here, we validated our fitting model across different mechanisms by simulating enzymatic progress curves for (i) uni-molecular reversible Haldane kinetics, (ii) irreversible kinetics with “idealized” substrate-binding mechanism, and (iii) “general” irreversible bi-molecular kinetics. Then, we fitted functions based on Michaelis-Menten, the corresponding analytical solutions, and the analytical solution used throughput the manuscript based on "idealized" binding.

III. DATA AVAILABILITY
Note that this script starts from MS data after preprocessing through in-house software, which includes peak picking, annotation and quantification of all measured ions. Raw mass spectrometry files are deposited in the MassIVE depository under MSV000096497.
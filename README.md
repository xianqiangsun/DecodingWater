# DecodingWater

# Information about the entropy calculation for water molecules:

## 1. Center information and coordinates parsing:   

•	a. Read the trajectory file ('.dcd' as file extension), with python script 'readDCDFile.py',you will get the file name 'distance.txt', which includes all the oxygen atom information and coordinates.    
•	b. Cluster all the water oxygen coordinated in the 'distance.txt' file to obtain the centers.   
•	c. Optimizing the centers obtained in the last step to combine centers near each other.   
•	d. Read the center coordinates and save H2O coordinates according to each center. 'getH2OCoordidateAll.py' is used in this step, and the coordinates for H and O atoms were saved into 'WW_allCentre_H2O.txt'.    

## 2. Protein-water correlation entropy calculation:    

•	a. Using 'SW_TranEntr.py' to calculate the protein-water correlation entropy for each water molecule. The output file is saves in 'SW.dat'.   
## 3. Water-Water correlation entropy calculation:    

•	a. 'WW_trans.py' is used to calculate the translational entropy of each water molecules. The output file is saved in 'WW_Trans.dat'.    
•	b.'WW_orientation_R3.py' or 'WW_orientation_R5.py' can be used to save the entropy wo water-water orientational entropy of each water-water correlation. The output file is saved in 'water_Orien.txt'.   

## 4.Water-system interaction energy:   

•	a. NAMDenergy in VMD is used to calculate the interaction energy of each system. 'energy_tcl_generation.py' is used to get the tcl file for namd energy calculation. PLS remind that include the namd directory and add -exe in each calculation.   
•	b. Each interaction energy is extracted with 'extract_energy.py'.   

------------------------------------------------------------------------------
## Ps: ''water_Gww_R1.py' is used to calculate the five angle correlation of bulk water.

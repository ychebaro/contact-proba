# contact-proba
Calculate contact probabilities

Calculate contact probabilities between two chains in a trajectory or a single pdb file. 
Possible selections of atoms for proteins are calpha, heavy atoms of main chain (heavymc) or side chain (heavysc), or all heavy atoms (heavyall). 
Possible selections of atoms for nucleic acids are C5' atoms (C5prime), heavy atoms of sugar and phosphate (heavymc) or nucleotide base (heavybase), or all heavy atoms (heavyall). 
A contact is considered to be formed if at least one heavy atoms between the two chains lies under the specified distance cutoff.
Pymol scripts can be generated where bonds are created between residues (-visu), according to probability intervals (0-0.25,0.25-0.5,0.5-0.75,0.75-1). Use the run option in pymol to execute the scripts. If you have a single pdb
file you would like to use for the visualization, you can specify it in the -pdbvisu option, otherwise the 
script will just create a pdb using the first frame of your trajectory 

Note: residue numbering should be continuous

## Prerequisites
You need :
* Python 2.7 
* NumPy
* MDAnalysis

To install MDAnalysis, you can use this command
pip install mdanalysis

## Description
usage: ./contacts.py [-h] -p PSF_FILE -d DCD_FILE -c1 SEGID1 -c2 SEGID2
                            [-s1 SELECTION1] [-s2 SELECTION2] -co CUTOFF
                            [-o OUTPUT_FILE] [-visu PYMOL]
                            [-pdbvisu PYMOL_PDB]
                            
-h, --help            show this help message and exit

-p PSF_FILE, --psf PSF_FILE
                       topology file used for simulation (pdb, psf)
                        
-d DCD_FILE, --dcd DCD_FILE
                        trajectory (dcd or pdb) or single frame in pdb format
                        
-c1 SEGID1, --chain1 SEGID1
                        segid of first chain for contact calculations
                        
-c2 SEGID2, --chain2 SEGID2
                        segid of second chain for contact calculations (similar as c1 if intrachain contacts needed)
                        
-s1 SELECTION1, --selection1 SELECTION1
                        atom selection for chain 1 (default is all heavy atoms)
                        
-s2 SELECTION2, --selection2 SELECTION2
                        atom selection for chain 2 (default is all heavy atoms)
                        
-co CUTOFF, --cutoff CUTOFF
                        distance cutoff
                        
-o OUTPUT_FILE, --output OUTPUT_FILE
                        name of output file for contact probability
                        
-visu PYMOL, --pymolvisu PYMOL
                        pymol files for visualization [Y/N]
                        
-pdbvisu PYMOL_PDB, --pymol_pdb PYMOL_PDB
                        If you want the pymol scripts to visualize your results on the pdb structure, choose a name for a pdb file otherwise it will just create one :-)
                          

* outputs:
- a contact probability matrix (matrix.csv)
- if specified, four pymol scripts which can be used to visualize the contact probability network (prob-0-0.25.pml, prob-0.25-0.5.pml, prob-0.5-0.75.pml, prob-0.75-1.0.pml) 
and a pdb file if needed named forpymol.pdb

## Example 

* Using a psf file and dcd trajectory, calculate contacts between heavy atoms of chain A and chain B using a 3.5 A cutoff, 
create pymol visualization scripts and use the file forpymol.pdb for the scripts:
./contacts.py -p structurefile.psf -d twoframes.pdb -c1 A -c2 B -co 3.5 -visu Y -s heavy -pdbvisu forpymol.pdb

* If you have no psf file you can still use it with the following options and your pdb file instead of the psf
./contacts.py -p structurefile.pdb -d twoframes.pdb -c1 A -c2 B -co 3.5 -visu Y -s heavy -pdbvisu forpymol.pdb
 
## Plotting with the script plot-contacts.py
If you want to plot the matrix, you can use the script plot-contacts.py (an example of plot is present in the example folder).
For this the arguments needed are of course the matrix and also the pdb file, the script will extract sequence residue name and number for labelling each axis. You can also transpose the matrix if you prefer it this way (the two possibilities are in the example folder).
Example:
Plotting with the transpose option
./plot-contacts.py -m matrix.csv -p forpymol.pdb -c1 A -c2 B -d1 4 -d2 4 -t Y

## Additional tool
If you need to change the name of a chain to another (for instance for making all protein chains into one chain id, if you need to calculate protein dna contacts this could be useful), use the script renumber.psf.


# TODO :
be able to use pdbs with missing residues...

Enjoy! 
Contact: yasmine.chebaro@gmail.com





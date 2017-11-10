load forpymol.pdb, pdbin
hide everything, /pdbin 

color purpleblue, pdbin and ((resi 240 and chain A) or (resi 148 and chain B)) 
bond (pdbin and resi 240 and name ca and chain A), (pdbin and resi 148 and name ca and chain B)
show lines, pdbin and ((resi 240 and chain A) or (resi 148 and chain B)) and name ca

color purpleblue, pdbin and ((resi 243 and chain A) or (resi 141 and chain B)) 
bond (pdbin and resi 243 and name ca and chain A), (pdbin and resi 141 and name ca and chain B)
show lines, pdbin and ((resi 243 and chain A) or (resi 141 and chain B)) and name ca

set_name pdbin, prob_0.25_0.5
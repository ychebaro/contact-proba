load forpymol.pdb, pdbin
hide everything, /pdbin 

color hotpink, pdbin and ((resi 229 and chain A) or (resi 133 and chain B)) 
bond (pdbin and resi 229 and name ca and chain A), (pdbin and resi 133 and name ca and chain B)
show lines, pdbin and ((resi 229 and chain A) or (resi 133 and chain B)) and name ca

color hotpink, pdbin and ((resi 241 and chain A) or (resi 148 and chain B)) 
bond (pdbin and resi 241 and name ca and chain A), (pdbin and resi 148 and name ca and chain B)
show lines, pdbin and ((resi 241 and chain A) or (resi 148 and chain B)) and name ca

color hotpink, pdbin and ((resi 242 and chain A) or (resi 145 and chain B)) 
bond (pdbin and resi 242 and name ca and chain A), (pdbin and resi 145 and name ca and chain B)
show lines, pdbin and ((resi 242 and chain A) or (resi 145 and chain B)) and name ca

set_name pdbin, prob_0.5_0.75
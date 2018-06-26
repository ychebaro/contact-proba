#!/usr/bin/env python
# -*- coding: utf-8 -*-


""" Contact probabilities """

# comment: need to add something for protein-ligand contacts


import subprocess
import string
import sys
import argparse
import numpy as np
from math     import sqrt

from MDAnalysis import Universe
from collections import defaultdict

np.set_printoptions(precision=3, threshold=np.nan)


__author__ = "Yasmine Chebaro"
__email__ = "yasmine.chebaro@gmail.com"




class GetContacts(object):
    """ Contact probabilities
    
    Parameters 
    ----------
    - universe : Universe object from MDAnalysis
    
    """
    
    def __init__(self, universe):
        self.u = universe

    def get_protein_coords(self, chainid, selection, nres, res):
        """ Extract protein coordinates
        
        The options available in this case are  
        Calphas, heavy atoms (side chain, main chain or both)
        
        Parameters
        ----------
        - chain id 
        - selection : calpha, heavyall, heavysc, heavymc
        - number of residues
        - residue in consideration
        
        Returns
        -------
        Coordinates of selected atoms
        """
        # All heavy atoms
        if selection=='heavyall':
                atprot = self.u.select_atoms("segid %s and resid %d and not name H* " %(chainid, nres[res]))
                coordinates = atprot.positions
                
        # Main chain heavy atoms
        if selection=='heavymc':
                atprot = self.u.select_atoms("segid %s and resid %d and backbone" %(chainid, nres[res]))
                coordinates = atprot.positions

        # Side chain heavy atoms
        if selection=='heavysc':
                atprot = self.u.select_atoms("segid %s and resid %d" %(chainid, nres[res]))
                coordinates = atprot.positions
                # Consider only Catoms for GLY residues
                if atprot.residues.resnames == 'GLY':
                    coordinates = atprot.select_atoms("name CA").positions
                else:
                    coordinates = atprot.select_atoms("not name H* and not backbone").positions
        
        # C alpha atoms
        if selection=='calpha':
                atprot = self.u.select_atoms("segid %s and resid %d and name CA" %(chainid, nres[res]))
                coordinates = atprot.positions
                    
        return coordinates

    def get_nucleic_coords(self, chainid, selection, nres, res):
        """ Extract nucleic acid coordinates
        
        The options available in this case are 
        C5' atoms, heavy atoms (nucleic base, sugar+phosphate (i.e. backbone) or both)
        
        Parameters
        ----------
        - chain id
        - selection : C5', heavyall, heavybase, heavymc
        - number of residues
        - residue in consideration
        
        Returns
        -------
        Coordinates of selected atoms
        """
        # All heavy atoms
        if selection=='heavyall':
                atnuc = self.u.select_atoms("segid %s and resid %d and not name H* " %(chainid, nres[res]))
                coordinates = atnuc.positions
                    
        # Main chain (sugar+phosphate) heavy atoms
        if selection=='heavymc':
                atnuc = self.u.select_atoms("segid %s and resid %d and (nucleicsugar or nucleicbackbone or (name O1P or name O2P or name OP1 or name OP2))" %(chainid, nres[res]))
                coordinates = atnuc.positions
            
        # Nucleic base
        if selection=='heavybase':
                atnuc = self.u.select_atoms("segid %s and resid %d and nucleicbase" %(chainid, nres[res]))
                coordinates = atnuc.positions

        # C5' atoms
        if selection=='C5prime':
                atnuc = self.u.select_atoms("segid %s and resid %d and name C5'" %(chainid, nres[res]))
                coordinates = atnuc.positions

        return coordinates

    def assign_biomolecule(self, chainid):
        """ Determine if chain is protein or nucleic
        
        Parameters
        ----------
        - chain id
        
        Returns
        -------
        Protein or Nucleic
        """
        nresids = 0
        chainprot = self.u.select_atoms("segid %s and protein" %(chainid))
        chainnuc = self.u.select_atoms("segid %s and nucleic" %(chainid))
        if chainnuc.n_residues != 0:
            bio = "nucleic"
            nresids = chainnuc.n_residues
        if chainprot.n_residues != 0:
            bio = "protein"
            nresids = chainprot.n_residues
        return bio, nresids
    
    def calc_dist(self, v1, v2):
        """ Calculate distance between two positions defined by 3D coordinates """
        dist = round(sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2)+pow(v1[2]-v2[2],2)),3)
        return dist
    
    def get_list_residues(self, chain, bio):
        """ Get list with number of residues """
        if bio == "protein":
            ca_resids = self.u.select_atoms("segid %s and name CA" %(chain))
            nres = ca_resids.resids
        if bio == "nucleic":
            c5_resids = self.u.select_atoms("segid %s and name C5'" %(chain))
            nres = c5_resids.resids
            
        return nres
    
    def run(self, chain1, chain2, selection1, selection2, cutoff):
        """ Run the contact calculation between chains
        
        Parameters
        ----------
        - chain 1 segid
        - chain 2 segid
        - selection chain 1
        - selection chain 2 
        - distance cutoff
        
        Returns
        -------
        Distance matrix
        """
        bio1, nresids1 = self.assign_biomolecule(chain1)
        bio2, nresids2 = self.assign_biomolecule(chain2)
        nres1 = self.get_list_residues(chain1,bio1)
        nres2 = self.get_list_residues(chain2,bio2)
        COOR1 = dict()
        COOR2 = dict()
        
        contactmat = np.zeros((nresids1, nresids2))
        
        for ts in self.u.trajectory:
            for res in range(len(nres1)):
                # find more elegant way to do this when I have time...
                if bio1 == "protein":
                    COOR1[res+1] = [self.get_protein_coords(chain1, selection1, nres1, res)]
                if bio1 == "nucleic":
                    COOR1[res+1] = [self.get_nucleic_coords(chain1, selection1, nres1, res)]
            for res in range(len(nres2)):                
                if bio2 == "protein":
                    COOR2[res+1] = [self.get_protein_coords(chain2, selection2, nres2, res)]
                if bio2 == "nucleic":
                    COOR2[res+1] = [self.get_nucleic_coords(chain2, selection2, nres2, res)]         
                    
            for key1, value1 in COOR1.iteritems():
                for key2, value2 in COOR2.iteritems():
                    find = False
                    for i in range(len(value1[0])):
                        for j in range(len(value2[0])):
                            dist = self.calc_dist(value1[0][i],value2[0][j])
                            if dist <= float(cutoff):
                                contactmat[key1-1, key2-1] += 1
                                find = True
                                if find:
                                    break
                        if find:
                            break
                        
        return contactmat/len(self.u.trajectory), bio1, bio2
    

def search_for_string(filein,lookup):

    # search for a string in a file and return the line number this string appears for the first time
    list_starts = []
    with open(filein) as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line[20:23] and line[0:4] == 'ATOM' :
                list_starts.append(int(num))

    return list_starts

def pymol_contact_visu(contactprob,pymol_pdb,chain1,chain2,bio1,bio2):

    # first get the first line for each chain 
    lookup1 = search_for_string(pymol_pdb, ' %s ' %chain1)
    lookup2 = search_for_string(pymol_pdb, ' %s ' %chain2)
    if bio1 == 'protein':
        at1 = "ca" 
    else:
        at1 = "c1'"
    if bio2 == 'protein':
        at2 = "ca"
    else:
        at2 = "c1'"

    # next check if the residue number does not start with 1, if not get the shift for each chain, 
    # this is needed for the pymol scripts, so that the row and column elements match properly the residue numbers in the pdb
    shift1 = 0
    shift2 = 0
    with open(pymol_pdb) as myFile:
        lines = myFile.readlines()
        if int(lines[lookup1[0]][23:26]) != 1:
            shift1 = int(lines[lookup1[0]][23:26])-1
        if int(lines[lookup2[0]][23:26]) != 1:
            shift2 = int(lines[lookup2[0]][23:26])-1
    myFile.close()
    
    file_list = ["prob-0-0.25.pml","prob-0.25-0.5.pml","prob-0.5-0.75.pml","prob-0.75-1.0.pml"]
    file_handles = []

    # create pml files for each interval, first thing is header to load the pdb
    try: 
      for name in file_list:
          file_handles.append(open(name,"w")) 
      
      for fh in file_handles:
          fh.write("load %s, pdbin\n" %pymol_pdb)
          fh.write("hide everything, /pdbin \n")
          fh.write("\n")
      
      # append in each file the residue numbers for contacts in each probability interval
      # bond the pairwise residues and color the bonds according to probability (blue to red)
      for row in range(len(contactprob)):
        for col in range(len(contactprob[0])):
          if contactprob[row,col] > 0 and contactprob[row,col] <= 0.25 :
              file_handles[0].write("color blue, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) \n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[0].write("bond (pdbin and resi %s and name %s and chain %s), (pdbin and resi %s and name %s and chain %s)\n" %((shift1+row+1), at1, chain1,(shift2+col+1), at2, chain2))
              file_handles[0].write("show lines, pdbin and ((resi %s and name %s and chain %s) or (resi %s and name %s and chain %s))\n" %((shift1+row+1),at1,chain1,(shift2+col+1),at2,chain2))
              file_handles[0].write( "\n")
          if contactprob[row,col] > 0.25 and contactprob[row,col] <= 0.5 :
              file_handles[1].write("color purpleblue, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) \n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[1].write("bond (pdbin and resi %s and name %s and chain %s), (pdbin and resi %s and name %s and chain %s)\n" %((shift1+row+1),at1,chain1,(shift2+col+1),at2,chain2))
              file_handles[1].write("show lines, pdbin and ((resi %s and name %s and chain %s) or (resi %s and name %s and chain %s))\n" %((shift1+row+1),at1,chain1,(shift2+col+1),at2,chain2))
              file_handles[1].write("\n")
          if contactprob[row,col] > 0.5 and contactprob[row,col] <= 0.75 :
              file_handles[2].write("color hotpink, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) \n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[2].write("bond (pdbin and resi %s and name %s and chain %s), (pdbin and resi %s and name %s and chain %s)\n" %((shift1+row+1),at1,chain1,(shift2+col+1),at2,chain2))
              file_handles[2].write("show lines, pdbin and ((resi %s and name %s and chain %s) or (resi %s and name %s and chain %s))\n" %((shift1+row+1),at1,chain1,(shift2+col+1),at2,chain2))
              file_handles[2].write("\n")
          if contactprob[row,col] > 0.75 and contactprob[row,col] <= 1.0 :
              file_handles[3].write("color red, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) \n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[3].write("bond (pdbin and resi %s and name %s and chain %s), (pdbin and resi %s and name %s and chain %s)\n" %((shift1+row+1),at1,chain1,(shift2+col+1),at2,chain2))
              file_handles[3].write("show lines, pdbin and ((resi %s and name %s and chain %s) or (resi %s and name %s and chain %s))\n" %((shift1+row+1),at1,chain1,(shift2+col+1),at2,chain2))
              file_handles[3].write("\n")
    
      file_handles[0].write("set_name pdbin, prob_0_0.25")
      file_handles[1].write("set_name pdbin, prob_0.25_0.5")
      file_handles[2].write("set_name pdbin, prob_0.5_0.75")
      file_handles[3].write("set_name pdbin, prob_0.75_1.0")


    finally:
      for fh in file_handles:
          fh.close()

    # now if some files are empty (no contacts in this interval) just delete them
    for name in file_list:
      p = subprocess.check_output(['wc %s' %name],shell=True)
      if int(string.split(p)[0]) <= 4:
          subprocess.Popen(['rm',string.split(p)[3]])



def check_pdb(pdbfile):
    # check for weird encoding in pdb file, sometimes happens when trajectory is written with VMD, 
    # MDAnalysis keeps the strange characters in the pdb file of the first frame
    with open(pdbfile) as myFile:
       p = subprocess.check_output(['grep -e REMARK %s | grep -e ^' %pdbfile],shell=True)
       if 'matches' in p:
           # this means that the pdb has something strange in the REMARK part
           # removing the REMARK lines
           subprocess.Popen(['sed -i "/REMARK/d" %s' %pdbfile], shell=True)
       else:
           None
def parse_options():
    parser = argparse.ArgumentParser(description="Calculate contact probabilities between two chains in a trajectory or a single pdb file. Possible selections of atoms for proteins are calpha, heavy atoms of main chain (heavymc) or side chain (heavysc), or all heavy atoms (heavyall). Possible selections of atoms for nucleic acids are C5' atoms (C5prime), heavy atoms of sugar and phosphate (heavymc) or nucleotide base (heavybase), or all heavy atoms (heavyall). Pymol scripts can be generated where bonds are created between residues (-visu), according to probability intervals (0-0.25,0.25-0.5,0.5-0.75,0.75-1).")
    parser.add_argument("-p", "--psf", dest="psf_file", required=True,
                        action="store", type=str,
                        help="topology file used for simulation (pdb, psf)")
    parser.add_argument("-d", "--dcd", dest="dcd_file", required=True,
                        action="store", type=str,
                        help="trajectory (dcd or pdb) or single frame in pdb format")
    parser.add_argument("-c1", "--chain1", dest="segid1", required=True,
                        action="store", type=str,
                        help="segid of first chain for contact calculations")
    parser.add_argument("-c2", "--chain2", dest="segid2", required=True,
                        action="store", type=str, 
                        help="segid of second chain for contact calculations (similar as c1 if intrachain contacts needed)")
    parser.add_argument("-s1", "--selection1", dest="selection1",
                        action="store", type=str,
                        default="heavyall", help="atom selection for chain 1 (default is all heavy atoms)")
    parser.add_argument("-s2", "--selection2", dest="selection2",
                        action="store", type=str,
                        default="heavyall", help="atom selection for chain 2 (default is all heavy atoms)")
    parser.add_argument("-co", "--cutoff", dest="cutoff", required=True,
                        action="store", type=str,
                        help="distance cutoff")
    parser.add_argument("-o", "--output", dest="output_file",
                        action="store", type=str,
                        default="matrix.csv", help="name of output file for contact probability")
    parser.add_argument("-visu", "--pymolvisu", dest="pymol",
                        action="store", type=str,
                        default="N", help="pymol files for visualization [Y/N]")
    parser.add_argument("-pdbvisu", "--pymol_pdb", dest="pymol_pdb",
                        action="store", type=str,
                        help="If you want the pymol scripts to visualize your results on the pdb structure, choose a name for a pdb file otherwise it will just create one :-)")

    return parser.parse_args()
        
       
def main():
    # get options
    options = parse_options()
    psf = options.psf_file
    dcd = options.dcd_file
    chain1 = options.segid1
    chain2 = options.segid2
    selection1 = options.selection1
    selection2 = options.selection2
    co = options.cutoff
    output = options.output_file
    visu = options.pymol
    pdbvisu = options.pymol_pdb
    
    # use MDAnalysis to read trajectory
    u = Universe(psf,dcd)

    # get contact probability 
    cp = GetContacts(u)
    contactprob, bio1, bio2 = cp.run(chain1, chain2, selection1, selection2, co)
    np.savetxt(output, contactprob, fmt='%4.2f', delimiter=" ")

    # generate pymol scripts if needed
    if visu == 'Y':
        # if no pdb file is supplied, write one from trajectory, first frame
        if pdbvisu == None:
            seleforpymol = u.select_atoms("segid %s or segid %s" %(chain1,chain2))
            seleforpymol.write('forpymol.pdb',remarks=None)
            pdbvisu = 'forpymol.pdb'
            # check pdb file format for weird encoding
            check_pdb(pdbvisu)
            
        pymol_contact_visu(contactprob, pdbvisu, chain1, chain2, bio1, bio2)

if __name__ == "__main__":
    main()
         


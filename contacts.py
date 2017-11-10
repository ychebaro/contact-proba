#!/usr/bin/env python
# -*- coding: utf-8 -*-


""" Contact probabilities """


import subprocess
import string
import argparse
import numpy as np
from math     import sqrt

from MDAnalysis import Universe, collection, Timeseries
from collections import defaultdict

np.set_printoptions(precision=3, threshold=np.nan)


__author__ = "Yasmine Chebaro"
__email__ = "yasmine.chebaro@gmail.com"




def getcontacts(u,chain1,chain2,selection,co):
	

    # get number of residues of each chain (intra or inter chain)
    ch1 = u.select_atoms("segid %s and name CA" %(chain1))
    ch2 = u.select_atoms("segid %s and name CA" %(chain2))
    resch1 = ch1.resids
    resch2 = ch2.resids

    COOR1 = dict()
    COOR2 = dict() 
    answer = np.zeros((ch1.n_residues,ch2.n_residues))
    nframes = len(u.trajectory)

    for ts in u.trajectory:

           # extract coordinates of appropriate atoms, for distance calculations (i.e. either Calphas or all heavy atoms)
           for res in range(len(resch1)):
              if selection=='heavy':
                 atprot1 = u.select_atoms("segid %s and not name H* and resid %d" %(chain1, resch1[res]))
                 COOR1[res+1] = [atprot1.positions]
              if selection=='calpha':
                 atprot1 = u.select_atoms("segid %s and name CA and resid %d" %(chain1, resch1[res]))
                 COOR1[res+1] = [atprot1.positions]
           for res in range(len(resch2)):
              if selection=='heavy':
                 atprot2 = u.select_atoms("segid %s and not name H* and resid %d" %(chain2, resch2[res]))
                 COOR2[res+1] = [atprot2.positions]
              if selection=='calpha':
                 atprot2 = u.select_atoms("segid %s and name CA and resid %d" %(chain2, resch2[res]))
                 COOR2[res+1] = [atprot2.positions]

           # now calculate distances, if less than cutoff then fill matrix with value 1
           for key1, value1 in COOR1.iteritems():
              for key2, value2 in COOR2.iteritems():
                  find = False
                  for i in range(len(value1[0])):
                     for j in range(len(value2[0])):
                        dist = round(sqrt(pow(value1[0][i][0]-value2[0][j][0],2)+pow(value1[0][i][1]-value2[0][j][1],2)+pow(value1[0][i][2]-value2[0][j][2],2)),3)
                        if dist <= float(co):
                          # fill matrix, values here are substracted by 1 due to python numbering
                          answer[key1-1,key2-1] += 1
                          find = True
                          if find:
                            break
                     if find:
                       break

    return answer/nframes


def search_for_string(filein,lookup):

    # search for a string in a file and return the line number this string appears for the first time
    list_starts = []
    with open(filein) as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line and line[0:4] == 'ATOM' :
                list_starts.append(int(num))

    return list_starts



def pymol_contact_visu(contactprob,pymol_pdb,chain1,chain2):

    # first get the first line for each chain 
    lookup1 = search_for_string(pymol_pdb, ' %s ' %chain1)
    lookup2 = search_for_string(pymol_pdb, ' %s ' %chain2)

    # next check if the residue number does not start with 1, if not get the shift for each chain, 
    # this is needed for the pymol scripts, so that the row and column elements match properly the residue numbers in the pdb
    with open(pymol_pdb) as myFile:
          lines = myFile.readlines()
          if int(lines[lookup1[0]][23:26]) != 1:
                 shift1 = int(lines[lookup1[0]][23:26])-2
          if int(lines[lookup2[0]][23:26]) != 1:
                 shift2 = int(lines[lookup2[0]][23:26])-2
          else:
                 shift1 = 0    
                 shift2 = 0    
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

      
      # append in each file the residue numbers for contacts in each probability interval
      # bond the pairwise residues and color the bonds according to probability (blue to red)
      for row in range(len(contactprob)):
        for col in range(len(contactprob[0])):
          if contactprob[row,col] > 0 and contactprob[row,col] <= 0.25 :
              file_handles[0].write("color blue, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) \n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[0].write("bond (pdbin and resi %s and name ca and chain %s), (pdbin and resi %s and name ca and chain %s)\n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[0].write("show lines, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) and name ca\n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[0].write( "\n")
          if contactprob[row,col] > 0.25 and contactprob[row,col] <= 0.5 :
              file_handles[1].write("color purpleblue, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) \n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[1].write("bond (pdbin and resi %s and name ca and chain %s), (pdbin and resi %s and name ca and chain %s)\n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[1].write("show lines, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) and name ca\n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[1].write("\n")
          if contactprob[row,col] > 0.5 and contactprob[row,col] <= 0.75 :
              file_handles[2].write("color hotpink, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) \n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[2].write("bond (pdbin and resi %s and name ca and chain %s), (pdbin and resi %s and name ca and chain %s)\n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[2].write("show lines, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) and name ca\n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[2].write("\n")
          if contactprob[row,col] > 0.75 and contactprob[row,col] <= 1.0 :
              file_handles[3].write("color red, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) \n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[3].write("bond (pdbin and resi %s and name ca and chain %s), (pdbin and resi %s and name ca and chain %s)\n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
              file_handles[3].write("show lines, pdbin and ((resi %s and chain %s) or (resi %s and chain %s)) and name ca\n" %((shift1+row+1),chain1,(shift2+col+1),chain2))
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
    parser = argparse.ArgumentParser(description="Calculate contact probabilities between two chains in a trajectory or a single pdb file. Distances are calculated between calpha or heavy atoms, with a user-defined cutoff. Pymol scripts can be generated where bonds are created between residues (-visu), according to probability intervals (0-0.25,0.25-0.5,0.5-0.75,0.75-1).")
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
    parser.add_argument("-s", "--selection", dest="selection",
                        action="store", type=str,
                        default="calpha", help="atom selection calpha or heavy atoms for distance calculation [calpha/heavy] (default is calphas)")
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
    selection = options.selection
    co = options.cutoff
    output = options.output_file
    visu = options.pymol
    pdbvisu = options.pymol_pdb

    # use MDAnalysis to read trajectory
    u = Universe(psf,dcd)

    # get contact probability 
    contactprob=getcontacts(u, chain1, chain2, selection, co)
    np.savetxt(output, contactprob,fmt='%4.2f',delimiter=" ")

    # generate pymol scripts if needed
    if visu == 'Y':
      # if no pdb file is supplied, write one from trajectory, first frame
      if pdbvisu == None:
        seleforpymol = u.select_atoms("segid %s or segid %s" %(chain1,chain2))
        seleforpymol.write('forpymol.pdb',remarks=None)
        pdbvisu = 'forpymol.pdb'     
        # check pdb file format for weird encoding
        check_pdb(pdbvisu)

      pymol_contact_visu(contactprob,pdbvisu,chain1,chain2)



if __name__ == "__main__":
    main()
         


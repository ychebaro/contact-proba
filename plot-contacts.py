#!/usr/bin/env python
import numpy as np
import argparse
import sys
from string import split
from numpy import genfromtxt
from scipy import stats
from matplotlib import rc
from matplotlib import colors
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MultipleLocator 
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib as mpl
import matplotlib.gridspec as gridspec



def getonelettercode(resname):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
     'CYR':'C','HSE':'H','HSD':'H','HSP':'H', 'ADE': 'A', 'CYT': 'C', 'GUA': 'G', 'THY': 'T'}
    if resname in d.keys():
        return d[resname]

def getseq(filein,chainid,freq):
    seqall = []
    with open(filein) as myFile:
        lines = myFile.readlines()
        for line in lines:
            if line[0:4]=='ATOM' and line[21]==str(chainid):
                if line[13:15] == "CA":
                    seqall.append(getonelettercode(line[17:20])+split(line[22:26])[0])
                if line[13:16] == "C5'":      
                    seqall.append(getonelettercode(line[17:20])+split(line[22:26])[0])
    seqfinale = []
    for i in range(0,len(seqall),freq):
        seqfinale.append(seqall[i]) 
    
    return seqfinale


def plotmatrix(filematrix,pdbfile,chain1,chain2,del1,del2,output,transpose):
   # read list of contact matrices to plot and pdb files for each (for axis legends)
   matrix = genfromtxt(filematrix,delimiter=' ')
   
   cmap = colors.ListedColormap(['white','#99FFFF','#66CCFF','#3399FF','#3366FF','#3333FF','#3300FF','#000099'])
   bounds= [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1]
   norm= colors.BoundaryNorm(bounds,cmap.N)
   
   # define figure size and x,y tick label size USER: adjust as you see fit
   fig = plt.figure(1, figsize=(10,10))
   gs = gridspec.GridSpec(1, 1)
   ax1 = fig.add_subplot(111)
  
   majorLocatorc1 = MultipleLocator(del1)
   majorLocatorc2 = MultipleLocator(del2)
   minorLocatorc1 = MultipleLocator(1)
   minorLocatorc2 = MultipleLocator(1)

   # plot matrice
   # if user wants to transpose matrix for clarity
   if transpose == "Y":
       im = ax1.matshow(matrix.T,cmap=cmap,norm=norm)

       # now get sequence info from pdb files, with the third argument of the function being the frequency of numbering in the plot for each chain
       chain1lab = getseq(pdbfile,chain1,del1)
       chain2lab = getseq(pdbfile,chain2,del2)

       # and now label axes accordingly and dont forget to adjust the frequency accordingly 
       ax1.xaxis.set_major_locator(majorLocatorc1)
       ax1.xaxis.set_minor_locator(minorLocatorc1)
       ax1.yaxis.set_major_locator(majorLocatorc2)
       ax1.yaxis.set_minor_locator(minorLocatorc2)
       ax1.set_xticklabels(['']+chain1lab, rotation='vertical',size=10)
       ax1.set_yticklabels(['']+chain2lab, rotation='horizontal',size=10)

   else: 
       im = ax1.matshow(matrix,cmap=cmap,norm=norm)

       # now get sequence info from pdb files, with the third argument of the function being the frequency of numbering in the plot for each chain
       chain1lab = getseq(pdbfile,chain1,del1)
       chain2lab = getseq(pdbfile,chain2,del2)

       # and now label axes accordingly and dont forget to adjust the frequency accordingly 
       ax1.xaxis.set_major_locator(majorLocatorc2)
       ax1.xaxis.set_minor_locator(minorLocatorc2)
       ax1.yaxis.set_major_locator(majorLocatorc1)
       ax1.yaxis.set_minor_locator(minorLocatorc1)
       ax1.set_xticklabels(['']+chain2lab, rotation='vertical',size=10)
       ax1.set_yticklabels(['']+chain1lab, rotation='horizontal',size=10)



   # add more global labels
   ax1.set_ylabel('Residue number',fontsize=12)
   ax1.set_xlabel('Residue number',fontsize=12)

   # and the colorbar
   fig.colorbar(im,fraction=0.03)

   plt.savefig(output)
#   plt.show()


def parse_options():
    parser = argparse.ArgumentParser(description="Plot contact probabilities")
    parser.add_argument("-m", "--matrix", dest="matrix_file", required=True,
                        action="store", type=str,
                        help="Matrix of contact probabilities")
    parser.add_argument("-p", "--pdb", dest="pdb_file", required=True,
                        action="store", type=str,
                        help="Pdb file to extract residues names and numbers")
    parser.add_argument("-c1", "--chain1", dest="segid1", required=True,
                        action="store", type=str,
                        help="segid of first chain for contact calculations")
    parser.add_argument("-c2", "--chain2", dest="segid2", required=True,
                        action="store", type=str,
                        help="segid of second chain for contact calculations (similar as c1 if intrachain contacts needed)")
    parser.add_argument("-d1", "--del1", dest="del1",
                        action="store", type=int, 
                        default="2", help="parameter for frequency of x axis labelling, corresponding to chain1")
    parser.add_argument("-d2", "--del2", dest="del2",
                        action="store", type=int,
                        default="2", help="parameter for frequency of y axis labelling, corresponding to chain2")
    parser.add_argument("-o", "--output", dest="output_file",
                        action="store", type=str,
                        default="matrix.png", help="name of output plot file")
    parser.add_argument("-t","--transpose", dest="transpose",
                        action="store", type=str,
                        default="N", help="[Y/N] if you want to transpose matrix for plotting, if chain2 if much smaller than chain1 this could be usefull")

    return parser.parse_args()



def main():
    # get options
    options = parse_options()

    matrix = options.matrix_file
    pdb    = options.pdb_file
    chain1 = options.segid1
    chain2 = options.segid2
    del1   = options.del1
    del2   = options.del2
    output = options.output_file
    transp = options.transpose

    # Plot matrix
    plotmatrix(matrix,pdb,chain1,chain2,del1,del2,output,transp)

if __name__ == "__main__":
    main()



#!/usr/bin/python2.7
# -*- coding: utf-8 -*-


import sys
import os
import string
import numpy as np
import argparse


def getlastresidnumber(filein, chainid, ext):
    with open(filein) as myFile:
        # first line is skipped to avoid problems extract residue number 
        lines = myFile.readlines()[1:]
        for line in lines:
            if ext == "ext" and line[11:12] == chainid:
                last_resid = int(line[20:25])
            if ext == "psf" and line[9:10] == chainid:
                last_resid = int(line[14:19])
            if ext == "pdb" and line[21:22] == chainid and line[0:4] == 'ATOM':
                last_resid = int(line[22:26])
        return last_resid
                    
def writenewfile(filein, chain_k, chain_r, lastresid, ext, fileoutname):
    fileout = open(fileoutname,'w')
    with open(filein) as myFile:
        lines = myFile.readlines()
        for line in lines:
            if ext == "ext":
                if line[11:12] == chain_r:
                    space = len(line[20:25]) - len(str(int(line[20:25])+lastresid))
                    fileout.write('%s %s%s%s%s%s' %(str(line[0:10]), chain_k, line[12:20], str(int(line[20:25])+lastresid), space*str(' '), line[25:]))
                else:
                    fileout.write('%s' %line)
            if ext == "psf":
                if line[9:10] == chain_r:
                    space = len(line[14:19]) - len(str(int(line[14:19])+lastresid))
                    fileout.write('%s %s%s%s%s%s' %(str(line[0:8]), chain_k, line[10:14], str(int(line[14:19])+lastresid), space*str(' '), line[19:]))
                else:
                    fileout.write('%s' %line)
            if ext == "pdb":
                if line[21:22] == chain_r:
                    space = len(line[22:26]) - len(str(int(line[22:27])+lastresid))
                    fileout.write('%s %s%s%s%s%s%s' %(str(line[0:20]), chain_k, space*str(' '), str(int(line[22:27])+lastresid), line[26:72], chain_k, line[73:]))
                else:
                    fileout.write('%s' %line)
    return

def parse_options():
    parser = argparse.ArgumentParser(description="Modify chain id and residue numbering, useful to merge two chains")
    parser.add_argument("-f", "--file", dest="inputfile", required=True,
                        action="store", type=str,
                        help="Input file in psf or pdb format")
    parser.add_argument("-c1", "--chain1", dest="chain1", required=True,
                        action="store", type=str,
                        help="Segid of the chain needing to be replaced")
    parser.add_argument("-c2", "--chain2", dest="chain2", required=True,
                        action="store", type=str,
                        help="Segid of the chain replacing chain1")
    parser.add_argument("-fo", "--formatfile", dest="formatfile", required=True,
                        action="store", type=str,
                        help="Extension of input file: psf, ext (extended psf) or pdb")
    parser.add_argument("-o", "--output", dest="output_file", required=True,
                        action="store", type=str,
                        help="Output psf or pdb file")
    
    return parser.parse_args()
    
       
def main():
    # get options
    options = parse_options()
    inputfile = options.inputfile
    chain1 = options.chain1 #replaced
    chain2 = options.chain2 #kept
    ext = options.formatfile
    output = options.output_file
    
    # first get last residue of the chain to keep
    lastresid = getlastresidnumber(inputfile, chain2, ext)
    writenewfile(inputfile, chain2, chain1, lastresid, ext, output)

if __name__ == "__main__":
    main()
         



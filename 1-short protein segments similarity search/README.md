

DESCRIPTION
=========================
This is a short protein segments similarity search code Written in Perl. 
This code takes PDB file [3D structure], run dssp file and convert it to angles encoded to character [1D Structural Alphabet(SA) sequence].
The encodes 1D Structural Alphabet used as query input and compared to the 1D Structural Alphabet in th database of C-alpha structural datadabse. 


Requirements
=========================
In order to use this code the former installation of the DSSP package is required. 
the current version of DSSP available here: https://swift.cmbi.umcn.nl/gv/dssp/index.html

This version is suitable for running on a computer with a Linux operating system.


 Contents of the package
=========================

  1.  code.pl					    -  The source code of short protein segments similarity search
  2.  database/database.txt			-  Phi and Psi databse
  3.  examples/1l4w.pdb				-  PDB file used to demonstrate the running of this code
  4.  output						-  A folder where the software output files are stored
  5.  README						-  This document
 

  NOTE:
  -----
  While the file 'database.txt' contain all the PDBid and their 1D Structural Alphabet (updated to Apr 2018), 
  the DSSP files themselves, which are necessary to more calculation, appear in the folder ../databases/dssp_files partially,
  so that only the DSSP files needed to run the below example are in this folder.
  In order to update the databse, read the file README in the '../databases' folder.
 

	SYNOPSIS
=========================

	perl code.pl -i queryfile [options]

	options:
	queryfile 				-i	<PDB file> [String] 
	chain					-c	<Chain> [String] (default = A)
	firstres				-f	<First amino acid position> [Real] (default = '')
	lastres					-l	<Last amino acid position> [Real] (default = '')
	seqmotif (amino acid)	-a	<Amino acid sequence motif> [String] (default = '')
	threshold				-t	<Minimum score to add a word to the BLAST lookup table> [Real] (default = 11)
	
	While the queryfile option is necessary, the other options are not.

         
	Example
=========================
	
	Example to demostrate how to use pdb file as query to search against Phi and Psi angles database.
	
	Run code.pl:
	perl code.pl -i examples/1l4w.pdb -c A -t 8 -f 7 -l 10

    The results of similiar strucures are saved in the output directory.

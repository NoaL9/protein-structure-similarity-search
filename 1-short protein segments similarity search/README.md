

DESCRIPTION
=========================
This is a short protein segments similarity search code Written in Perl. </br>
This code takes PDB file [3D structure], run dssp file and convert it to angles encoded to character [1D Structural Alphabet(SA) sequence].</br>
The encodes 1D Structural Alphabet used as query input and compared to the 1D Structural Alphabet in th database of C-alpha structural datadabse. </br>


Requirements
=========================
In order to use this code the former installation of the DSSP package is required. </br>
the current version of DSSP available here: https://swift.cmbi.umcn.nl/gv/dssp/index.html</br>

This version is suitable for running on a computer with a Linux operating system.</br>


 Contents of the package
=========================

  1.  code.pl					    	-  The source code of short protein segments similarity search.
  2.  database/database.txt				-  Phi and Psi databse. Create directory with the name 'database'. Download to this directory the file from this link: https://drive.google.com/file/d/12CWurBNHcxKpswsBWLgHR1Q9ML4FvBTj/view?usp=sharing and unzip the file.
  3.  examples/1l4w.pdb					-  PDB file used to demonstrate the running of this code.
  4.  output						-  Create folder with the name 'output'. This folder is where the software output files are stored.
  5.  README						-  This document.
 

  NOTE:
  -----
  While the file 'database.txt' contain all the PDBid and their 1D Structural Alphabet (updated to Apr 2018), </br>
  the DSSP files themselves, which are necessary to more calculation, appear in the folder ../databases/dssp_files partially,</br>
  so that only the DSSP files needed to run the below example are in this folder.</br>
  In order to update the databse, read the file README in the '../databases' folder.
 </br>

 SYNOPSIS
=========================

	perl code.pl -i queryfile [options]

	options:
	queryfile 				-i	<PDB file> [String] 
	chain					-c	<Chain> [String] (default = A)
	firstres				-f	<First amino acid position> [Real] (default = '')
	lastres					-l	<Last amino acid position> [Real] (default = '')
	seqmotif (amino acid)			-a	<Amino acid sequence motif> [String] (default = '')
	threshold				-t	<Minimum score to add a word to the BLAST lookup table> [Real] (default = 11)
	
	While the queryfile option is necessary, the other options are not.

         
 Example
=========================
	
	Example to demostrate how to use pdb file as query to search against Phi and Psi angles database.
	
	Run code.pl:
	perl code.pl -i examples/1l4w.pdb -c A -t 8 -f 7 -l 10

    The results of similiar strucures are saved in the output directory.

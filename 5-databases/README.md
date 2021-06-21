  NOTE:
  -----
  AngLAST databases in FASTA format and in NCBI binary format are allready exist in this local directory (Updated to Apr 2018).
  The bellow description and instructions are suitble for establishing new AngLAST database.

  Requirements
=========================

  In order to establishing angblast database the former installation of the DSSP package and ncbi-blast package is required. 
  the current version of DSSP is available here: https://swift.cmbi.umcn.nl/gv/dssp/index.html
  the ncbi package version 2.2.7 is available here: https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.7/



 Contents of the package
=========================

  1.  database_mfa.pl				-  Script to use in order to turn PDB files to C-alpha structural datadabse 
  2.  pdb_rsync.pl					-  Update C-alpha structural datadabse (optionally) 
  3.  pdb_rsync						-  Directory for decompresses PDB files for rsync (optionally)     
  4.  pdb_files						-  Directory of all PDB files (The files are necessary for the creation of the AngLAST database and for additional calculations during the running of AngLAST)
  5.  dssp_files					-  Directory of all DSSP files for generating database (The files are necessary for the creation of the AngLAST database and for additional calculations during the running of AngLAST)
  6.  anglastDB.mfa					-  C-alpha structural datadabse  - FASTA format
  7.  ncbi_anglast_format			-  Directory of C-alpha structural datadabse - NCBI binary format
  8.  README						-  This document
 

  Instructions
=========================
   
	Before using Anglast, user needs to generate the structural C-alpha database. There are two ways:
 
	Generating Angblast Database using the script database_mfa.pl.
	
	database_mfa.pl:
	DESCRIPTION
	Script to use in order to turn PDB files to C-alpha structural datadabse. 
	At the first stage, this script reads a list of pdb-style and dssp-style files, which inside the pdb_files and dssp_files directories, respectively 
	and translates all the PDB structures into Ca angles encoded to character [1D Structural Alphabet(SA) sequence]. 
	The 1D Structural Alphabet(SA) printed to file in FASTA format (anglastDB.mfa). 
	Then, database_mfa.pl runs 'makeblastdb' on the anglastDB.mfa file (AngLAST database in FASTA format).
    (The makeblastdb command: 'makeblastdb -in angblastDB.mfa -out ncbi_anglast_format/anglastDB') 	
	The output files are AngLAST database in NCBI binary format with three different suffix:phr, pin and psq.
	These files are generated in the 'ncbi_anglast_format' directory.
	
	SYNOPSIS
	--------
	perl database_mfa.pl 
	
		
	pdb_rsync.pl:
	This script allow to update the anglastDB by update the PDB files and calculate DSSP for only the new PDB structures.
	
	
	SYNOPSIS
	--------	
    perl database_mfa.pl 
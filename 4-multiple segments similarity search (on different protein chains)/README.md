DESCRIPTION
=========================
This code is a perl script similiar to the code of protein domain similarity search with a little modification.
While protein domain similarity search can get one fragmet/chains as a query, this code can get multiple fragments/chains.

The instructions and requirements are similiar to those of protein domain similarity search and the difference is reflected only in the command line.
Unlike protein domain similarity search command line, in this code command line can conculde multiple chains, multiple first amino acid position and multiple last amino acid position.

For more details, read the README of protein domain similarity search.


 Contents of the package
=========================

  1.  code.pl			        	-  The source code 
  2.  database/anglastDB.mfa		-  The databse
  3.  lib/aa_matrix.pm				-  BLOSUM62 amino acid substitution matrix
  4.  lib/BLOSUM80					-  C-alpha angles ansubstitution matrix (for long matches)
  5.  lib/BLOSUM45					-  C-alpha angles ansubstitution matrix (for short matches)
  6.  examples/1avo.pdb			    -  PDB file used to demonstrate the running of this code
  7.  output						-  A folder where the software output files are stored
  8.  README						-  This document  
  
        
	Example
=========================
	
	Example to demostrate how to use pdb file as query to search against C-alpha angles database.
	
	Run code.pl:
	perl code.pl -i example/1avo.pdb -c A -c B -w 3 -g 11 -e 1 -t 11 -p 1 -d 20 -h 10000
	  

	In the example command line below,  both chain A and chain B of 1avo.pdb structure are searched against C-alpha structural database. 
	The results are saved in the output directory.
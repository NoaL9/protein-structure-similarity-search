DESCRIPTION</br>
=========================
</br>
This script is a perl wrapper for BlastP which allow to </br>
Implement new strategy for searching protein structure similarity using the </br>
structural information embedded in C-alpha angles.</br>
</br>
This code takes PDB file [3D structure], create its dssp file and convert it to C-alpha angles encoded to character [1D Structural Alphabet(SA) sequence] </br>
and than printed the 1D Structural Alphabet to FASTA file which used as query input. </br>
This allows to compare the C-alpha angles query sequence with a database of C-alpha 1D Structural Alphabet datadabse</br>
by runing blastp using the C-alpha angles ansubstitution matrix (CASM).</br>
 </br>
This package aim to be similar to NCBI-Blast in terms of usage. </br>
This means that you can pass to the script the majority of the options </br>
that you normally would use for NCBI-Blast. By default the script call Blast </br>
with the same gap opening, gap extension and filtering settings used </br>
in our benchmark.
</br>
</br>
Requirements</br>
=========================
</br>
In order to use this package the former installation of the DSSP package and ncbi-blast package is required. </br>
the current version of DSSP is available here: https://swift.cmbi.umcn.nl/gv/dssp/index.html</br>
the ncbi package version 2.2.7 is available here: https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.7/</br>
The path to BLAST matrices should be changed in the .ncbirc file.</br>
You should specified the location of CASM matrices. i.e:</br>
BLASTMAT=../lib</br>
lib is a local directory in this package that contain two matrices files:</br>
BLOSUM80- contain CASM1 matrix.</br>
BLOSUM45- contain CASM1 matrix.</br>
In addition,the following perl moduls should be installed: Parallel::Loops and Bio::PDB::Structure </br>
this package itself doesn't require any installation process. Just run the perl scripts which in turn will call NCBI-Blast on your computer.</br>
</br>
This version is suitable for running on a computer with a Linux operating system.</br>
</br>
</br>
 Contents of the package</br>
=========================
</br>
  1.  code.pl				 	    	-  The source code </br>
  2.  database/anglastDB.mfa		                -  The C alpha angles databse (unzip the database.zip file)</br>
  3.  lib/aa_matrix.pm					-  BLOSUM62 amino acid substitution matrix</br>
  4.  lib/BLOSUM80					-  C-alpha angles ansubstitution matrix (for long matches)</br>
  5.  lib/BLOSUM45					-  C-alpha angles ansubstitution matrix (for short matches)</br>
  6.  examples/2i0l.pdb			    		-  PDB file used to demonstrate the running of this code</br>
  7.  output						-  You shuld create directory with the name 'output'. This is A folder where the software output files are stored.</br>
  8.  README						-  This document</br>  
  </br>
</br>
  NOTE:</br>
  -----
 </br>
  While the file 'anglastDB.mfa' contain all the PDBid and their 1D Structural Alphabet (updated to Apr 2018), 
  the DSSP files and PDB files themselves, which are necessary to more calculation, appear in the folder ../databases/dssp_files and ../databases/pdb_files partially,
  so that only the DSSP files and PDB files needed to run the below example are in this folder.
  In order to update the anglastDB.mfa databse, read the file README in the '../databases' folder. 
</br>
</br>
  SYNOPSIS</br>
=========================
</br>
	perl code.pl -i queryfile [options]

	options:
	filename 				-i	<PDB file> [String] 
	chain_q					-c	<Chain> [String] (default = A)
	dssp					-b	<DSSP files path> [String] (default = "../databases/dssp/")
	pdb					-j	<PDB files path> [String] (default = "../databases/pdb/")
	firstres				-f	<First amino acid position> [Real] (default = '')
	lastres					-l	<Last amino acid position> [Real] (default = '')
	conserved_aa 				-s	<Conserved Amino acid / sequence motif> [String] (default = '')
	word_size				-w	<Length of initial exact match> [Real] (default = 3)
	window_size				-d	<Multiple hits window size> [Real] (default = 40)
	gap_score				-g	<Cost to open a gap> [Real] (default = 11)
	ext_score				-e	<Cost to extend a gap> [Real] (default = 1)
	threshold				-t	<Minimum score to add a word to the BLAST lookup table> [Real] (default = 11)
	evalue					-v	<E-value threshold> [Real] (default = 10.0)
	number of hits 				-h	<Number of sequences to show alignments> [Integer] (default = 1000)
	maxprocs 				-p	<Number of threads (CPUs) to use in blast search> [Integer] (default = 4)
	matrix					-m	<C-alpha angles ansubstitution matrix> [String] (default = CASM1 (for longmatch)) (alternative = CASM2 (for short_match))

	

	The default values : -c A -w 3 -d 40 -g 11 -e 0.01 -t 11 -v 1 -h 1000 -p 1.
	However, for short fragments (less than 30 aa) the reccomented values are : c A -w 2 -d 15 -g 9 -e 1 -t 6 -v 10000 -h 1000 -p 1.
</br>
 </br>       
   Example</br>
=========================
</br>	
	Example to demostrate how to use pdb file as query to search against C-alpha angles database.
	
	Run code.pl:
	perl code.pl -i example/2i0l.pdb -c A -s GLGF -w 3 -g 11 -e 1 -t 20 -p 1

	The results are saved in the output directory.
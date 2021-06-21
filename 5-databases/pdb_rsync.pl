#!/usr/bin/perl

# Update_PDB_Files_and_AngLAST_Database
#
# Script to download PDB Files from the RCSB ftp and create formatted BLAST databases.
#
# This script will download multiple tar files of PDB structures. 
# The command rsync allow to compare the local copy of the database tar file(s) 
# to the newer version of the database and only download new tar files. 
#
# This script can be run without any arguments.


#######################
##PDB- download files##
#######################

# download PDB files
system('rsync -rlpt -v -z --delete --port=33444 \rsync.rcsb.org::ftp_data/structures/divided/pdb/pdb_rsync');

# Decompresses the files
system ('mkdir pdb_gzfiles');

system('find pdb_rsync/ -type f -print0 | xargs -0 cp -t pdb_gzfiles');

system ('ls pdb_gzfiles > pdb_gzfiles.txt'); 

open IN, "pdb_gzfiles.txt";
while (<IN>){
	chomp ($_);
	system("gunzip pdb_gzfiles/$_");
}
close IN;
system ('rm  pdb_gzfiles.txt');

# Extracts the pdbid for new file name- pdbid.pdb

system ('ls pdb_gzfiles > pdb_gzfiles.txt');
open IN, "pdb_gzfiles.txt";
while (<IN>){
	chomp ($_);
	my $pdbid = substr($_,3,4);
	my $new_name = "pdb_gzfiles/".$pdbid.".pdb";
	my $old_name = "pdb_gzfiles/$_";
	rename $old_name, $new_name;
}
close IN;

system ('rm  pdb_gzfiles.txt');


####################
##DSSP_Calculation##
####################

# Obtains the list of new pdb files

system ('ls pdb_files > pdb_files_list1.txt'); 
system ('ls pdb_gzfiles > pdb_files_list2.txt');

my $first_file  = shift || 'pdb_files_list1.txt';
my $second_file = shift || 'pdb_files_list2.txt';

open my $a_fh, '<', $first_file  or die "$first_file: $!";
open my $b_fh, '<', $second_file or die "$second_file: $!";

my %second_file;
@second_file{map { unpack 'A*', $_ } <$a_fh>} = ();


open OUT, ">compare.txt";
while (<$b_fh>) {
    print OUT unless exists $second_file{unpack 'A*', $_};
}
close OUT;


# DSSP calculation for new PDB files 

open IN, "compare.txt";

while (my $line=<IN>){
	chomp ($line);
	my $pdbid=substr($line,0,5);
	my $dsspfile = $pdbid."dssp";
	system ( "dssp pdb_gzfiles/$line dssp_files/$dsspfile");
}
close IN;


system ('ls pdb_gzfiles > pdb.txt'); 
system ('ls dssp_files > /home/samsona/public_html/dssp.txt');


my $pdb_list_file = "pdb.txt";
open (PDB_LIST,$pdb_list_file);

my %pdb;
while (<PDB_LIST>){
	chomp($_);
	my $pdbid = substr($_,0,4);
	$pdb{$pdbid}++;
}

my $dssp_list_file = "dssp.txt";
open (DSSP_LIST,$dssp_list_file);

while (<DSSP_LIST>){
	chomp($_);
	my $pdbid = substr($_,0,4);
	if (!$pdb{$pdbid}){
		my $dssp_file = "dssp/".$pdbid.".dssp";
		system ("rm $dssp_file");
	}
}

system ('ls dssp_files > dssp.txt');


system ('rm -r pdb_files');

system ('mv pdb_gzfiles pdb_files');



#################################
##Create_TXT_Database_of_Angles##
#################################

my @code=("Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","R","R","R","R","R","R","R","R","R","R","R","R","R","R","R","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","I","I","I","I","I","I","I","I","I","I","I","I","I","I","I","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","K","K","K","K","K","K","K","K","K","K","K","K","K","K","K","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","X","X","X","X","X","X","X","X","X","X","X","X","X","X","X","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","V","V","V","V","V","V","V","V","V","V","V","V","V","V","V","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M");


open (IN, "dssp.txt");
open (OUT, ">anglastDB.mfa");


while ($file=<IN>){

	chomp($file);

	my $pdbid = substr($file,0,4);
	my $pdbfile = $pdbid.".pdb";
	my $dsspfile = $pdbid.".dssp";
	open (PDB, "pdb_files/$pdbfile") or die "cannot open $pdbfile:$!\n";
	my $title;
	while (<PDB>){
		chomp($_);
		if ($_=~/TITLE/){
			$title=substr($_,10);
			last;
		}
	}
	$title =~ s/^\s+|\s+$//g;
	$title =~ s/ /_/g;
	close PDB;
		
	open (DSSP, "dssp_files/$dsspfile")or die "cannot open $dsspfile:$!\n";;

	my @chains;
	my $ref = {}; my $countq= 1; 

	while (<DSSP>){
		next if (($_=~/\.$/) or ($_=~/\#/));
		my $chain_here = substr ($_, 11, 1);
		my $chain_near = substr ($_, 13, 1);		
		if (($chain_here eq ' ') and ($chain_near ne '!')){$chain_here='A';}
		next if $_=~/!/;
									
		push @chains, $chain_here;
		my $new_ca = int(substr($_, 97, 6));			
		if ($new_ca==360) {
			$new_ca=0; 
		}	
		next if !$new_ca;
		
		my $res=substr ($_,13,1);
        if($res=~/[a-z]/){
              $res="C";
        }		
		next if !$an_cnvd{$res};
		
		$ref -> {$chain_here} -> {ca} .= $code[$new_ca];
	}
	close DSSP;
	
	my @chains_filtered = uniq(@chains);
	
	foreach (@chains_filtered){
		my $chain=$_;
		my $str_ca=$ref -> {$chain} -> {ca};
		my $resnum=$ref -> {$chain} -> {resnum};
		my @str=split(//,$str_ca);
		my $c=1;
		my $ca_str;
		foreach(@str){
			if ($c==60) {
				$ca_str.= $_."\n"; 
				$c = 0;
			}
			else{
				$ca_str.= $_;
			}
			$c++;
		}
		$ref -> {$chain} -> {count}=@str;
		next if !$ca_str;
		chomp($ca_str);
		my $output = $pdbid.$chain."/".$ref -> {$chain} -> {count}."/".$title."\n".$ca_str."\n";
		print OUT ">$output";
		
	}
}
close OUT;
close IN;



sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
}


########################################
##Convert_TXT_Database_to_BLAST_Format##
########################################

system ("makeblastdb -in anglastDB.mfa -out ncbi_anglast_format/anglastDB");

#!/usr/bin/perl

#######################
##Encoding_Structures##
#######################

my @code=("Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","R","R","R","R","R","R","R","R","R","R","R","R","R","R","R","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","I","I","I","I","I","I","I","I","I","I","I","I","I","I","I","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","K","K","K","K","K","K","K","K","K","K","K","K","K","K","K","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","X","X","X","X","X","X","X","X","X","X","X","X","X","X","X","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","V","V","V","V","V","V","V","V","V","V","V","V","V","V","V","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M");

#C-alpha structural datadabse - FASTA format
open (OUT, ">anglastDB.mfa");
system ("ls pdb_files/  > files.txt");

open (FILES, "files.txt");

my $numf=0;
while ($file=<FILES>){
	
	$numf++;
	chomp($file);
	next if ($file !~ /dssp$/);

	$pdbid=substr($file,0,4);
	$pdbfile=$pdbid.".pdb";
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
		
	open (DSSP,"dssp_files/$dsspfile")or die "cannot open $dsspfile:$!\n";;

	my @chains;
	my $ref = {}; my $countq= 1; 

	while (<DSSP>){
		next if (($_=~/\.$/) or ($_=~/\#/));
		my $chain_here = substr ($_, 11, 1);
		my $chain_near = substr ($_, 13, 1);
		if ($chain_here eq ' '){$chain_here='A';}	
		next if ($chain_near eq "!");
							
		my $num_residue = substr ($_,6,4);	
		$ref -> {$chain_here} -> {resnum} .= $num_residue.",";
		
		push @chains, $chain_here;
		my $new_ca = int(substr($_, 97, 6));			
		if ($new_ca==360) {
			$new_ca="0"; 
		}	

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
		my @resnum=split(/,/,$resnum); 
		my $last_one = pop @resnum;
		$last_one =~ s/^\s+|\s+$//g;
		my $first = shift @resnum;
		$first =~ s/^\s+|\s+$//g;
		my $output = $pdbid.$chain."/".$first."-".$last_one."/".$title."\n".$ca_str."\n";
		print OUT ">$output";
		
	}
}
close OUT;
close FILES;

sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
}

########################################
##Convert_TXT_Database_to_BLAST_Format##
########################################

#C-alpha structural datadabse - NCBI binary format
system ("makeblastdb -in anglastDB.mfa -out ncbi_anglast_format/anglastDB");
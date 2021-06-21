#!/usr/bin/perl -w


# ----------------------------------- #
#       Open source libraries         #
# ----------------------------------- #

$| = 1;

use File::Basename;
use List::Util qw[min max];
use POSIX;
use Time::HiRes qw/gettimeofday tv_interval/;
use Parallel::Loops;
use Bio::PDB::Structure;

# ----------------------- #
#    Custom libraries     #
# ----------------------- #

use lib './lib';
use aa_matrix;

# ----------------------- #
#	input_parameters	  #
# ----------------------- #

my $start_time = [gettimeofday];


#Input_variable_initialization

#<PDB file> [String] 
my $filename = "";
#<Chain> [String] (default = A)
my $chain_q = "A";
#<First amino acid position> [Real] (default = '')
my $firstres = '';
#<Last amino acid position> [Real] (default = '')
my $lastres = '';
#<Amino acid sequence motif> [String] (default = '')
my $conserved_aa = '';
#<Length of initial exact match> [Real] (default = 3)
my $word_size = 3;
#<Multiple hits window size> [Real] (default = 40)
my $window_size = 40;
#<Cost to open a gap> [Real] (default = 11)
my $gap_score = 11;
#<Cost to open a gap> [Real] (default = 11)
my $ext_score = 0.01;
#<Minimum score to add a word to the BLAST lookup table> [Real] (default = 11)
my $threshold = 11;
#<E-value threshold> [Real] (default = 10.0)
my $evalue = 10;
#<Number of sequences to show alignments> [Integer] (default = 1000)
my $number_of_hits = 1000;
#<Ca angles ansubstitution matrix>
my $matrix = 'BLOSUM80';
#<Number of threads (CPUs) to use in blast search> [Integer] (default = 1)
my $maxProcs = 1;
my $pl = Parallel::Loops->new($maxProcs);

my $job = $$;

#Get input from user
my @arrIn;
while(@ARGV){
	my $i = shift;
	push @arrIn, $i;
}
$l = @arrIn;

my %hash;
for ($j=0;$j<$l-1;$j=$j+2){
	$hash{$arrIn[$j]}=$arrIn[$j+1];
}

for my $outer_key ( keys %hash ){
	if ($outer_key eq "-i" or $outer_key eq "-c" or $outer_key eq "-f" or $outer_key eq "-l" or $outer_key eq "-s" or $outer_key eq "-w" or $outer_key eq "-d" or $outer_key eq "-g" or $outer_key eq "-e" or $outer_key eq "-t" or $outer_key eq "-v" or $outer_key eq "-h" or $outer_key eq "-p" or $outer_key eq "-m" or $outer_key eq "-o"){
		if ($outer_key eq "-i"){
			$filename=$hash{-i};
		}
		elsif ($outer_key eq "-c"){
			$chain_q=$hash{-c};
		}
		elsif ($outer_key eq "-f"){
			$firstres=$hash{-f};
		}
		elsif ($outer_key eq "-l"){
			$lastres=$hash{-l};
		}
		elsif ($outer_key eq "-s"){
			$conserved_aa=$hash{-s};
			$conserved_aa = uc $conserved_aa;
		}
		elsif ($outer_key eq "-w"){
			$word_size=$hash{-w};
		}
		elsif ($outer_key eq "-d"){
			$window_size=$hash{-d};
		}
		elsif ($outer_key eq "-g"){
			$gap_score=$hash{-g};
		}
		elsif ($outer_key eq "-e"){
			$ext_score=$hash{-e};
		}
		elsif ($outer_key eq "-t"){
			$threshold=$hash{-t};
		}
		elsif ($outer_key eq "-v"){
			$evalue=$hash{-v};
		}
		elsif ($outer_key eq "-h"){
			$number_of_hits=$hash{-h};
		}
		elsif ($outer_key eq "-m"){
			$matrix=$hash{-m};
			if ($matrix eq "CASM1"){
				$matrix = "BLOSUM80";
			}
			elsif ($matrix eq "CASM2"){
				$matrix = "BLOSUM45";
			}
			else{
				print "Error! Please enter a valid matrix name (CASM1 or CASM2)";
				exit;
			}
		}
		elsif ($outer_key eq "-p"){
			$maxprocs=$hash{-p};
		}		
	}
	else{
		print "Error! Please check your command.\n"; 
		exit;
	}
}

	
#--------------------------------------------
#Validating variables


if($filename){	
	unless (-e $filename){
		print "The requested file $filename was not found\n";
		exit;
	}
}		
else{
	print "The query file doesn't exist.";
	exit;
}

$chain_q =~ tr/[a-z]/[A-Z]/;
if ($chain_q eq ''){
		$chain_q = "A";
}
elsif ($chain_q !~ /\w/){
		print "Error! Please enter an alphabetical value for chain or leave blank\n";
		exit;
}

if (($firstres !~ /\d+/) and ($firstres ne '')){
		print "Please enter a numerical value for 'First Residue or leave blank'\n";
		exit;
}

if (($lastres !~ /\d+/) and ($lastres ne '')){
		print "Please enter a numerical value for 'Last Residue' or leave blank\n";
		exit;
}

if (($conserved_aa!~/[QWERTYIOPASDFGHKLCVNM]/)and($conserved_aa ne '')){
  print "Error! Please enter a valid amino acid (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, and Y) or leave blank\n"; 
  exit;
}

if ($filename=~/\./){
	@path=split(/\//,$filename);
	$pdb=$path[0];
}
else{
	$pdb=$filename;
}

my $pdbidq=substr($pdb,0,4).$chain_q;

# DSSP calculation - QUERY FILE
#--------------------

my $dsspfile = $file.".dssp";
system ("dssp $input $dsspfile");



#ENCODING QUERY STRUCTURE	
my @code=("Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","Q","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","R","R","R","R","R","R","R","R","R","R","R","R","R","R","R","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","I","I","I","I","I","I","I","I","I","I","I","I","I","I","I","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","J","J","J","J","J","J","J","J","J","J","J","J","J","J","J","K","K","K","K","K","K","K","K","K","K","K","K","K","K","K","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","X","X","X","X","X","X","X","X","X","X","X","X","X","X","X","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","V","V","V","V","V","V","V","V","V","V","V","V","V","V","V","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M");


# DSSPFileExtract - Encoding C-alpha angles to characters(query)
my ($ca_str_query,$aa_query_all,$refcaq,$refaaq,$refaccq) = DSSPFileExtract(1,\@code,$dsspfile, $chain_q, $firstres, $lastres,1);

my %refcaq = %{$refcaq};
my %refaaq = %{$refaaq};
my %refaccq = %{$refaccq};

# Sequence motif  
if (($conserved_aa) and ($aa_query_all !~ /$conserved_aa/i)){
	print "Please enter a right regular expression value for 'Motif Residues' or leave blank\n";
	exit;
}

# Print C-alpha characters to new file -> this file is the input file of BLAST program
my $TranslatedQuery="output/query_".$job;
open (QUERY, ">$TranslatedQuery") or die "cannot open $TranslatedQuery:$!\n";
print QUERY ">$filename/query\n$ca_str_query\n";
close QUERY;
	
	

# ----------------------- #
#       Main Program      #
# ----------------------- #

#--------------------------------------------
# BLAST 

my $blastout = "output/blastout_".$job;
	
system("blastp -i $filename -M $matrix -d database/anglastDB.mfa -o $blastout -W $word_size -G $gap_score -E $ext_score -f $threshold -e $evalue -b $number_of_hits -A $window_size -p $maxprocs -V T -v 0 -C 0 -F \"\" -a 1");


#--------------------------------------------
# EDITING BLAST OUTPUT FILE

my ($pdbid_list,$expv,$BLASTOUT) = BLASTOUT_EXTRACT($blastout);  # hash of hashes

my %pdbid = %{$pdbid_list};
my @pdbid = keys %pdbid;
my $hits = @pdbid;
%pdbid ={};
my %expv = %{$expv};
my %BLASTOUT = %{$BLASTOUT};
	
			
my $nun_of_word=@pdbid;
my $part= int($nun_of_word/4);
my $m = $nun_of_word%4;

my @one = @pdbid[0..$part-1];
my @two = @pdbid[$part..$part+$part-1];
my @three = @pdbid[2*$part..3*$part-1];
my @four = @pdbid[3*$part..4*$part+$m-1];

push @one, "one";
push @two, "two";
push @three, "three";
push @four, "four";

my @i=( "@one", "@two", "@three", "@four");


my %returnValues;
$pl->share( \%returnValues);
		
$pl->foreach( \@i, sub {
		my @part = split(" ",$_);
		my $last_one = pop @part;
		while(@part) {			
			my $pdbid = shift @part;
			my $pdbidt = substr($pdbid,0,4);	
			my $chaint = substr($pdbid,4,1);
			my $dsspfile_target ="../databases/dssp_files/$pdbidt.dssp";	
			
			# Query	
									
			my $lenq = $BLASTOUT{$pdbid}{lposq} - $BLASTOUT{$pdbid}{fposq} + 1;
			my $fpos = $BLASTOUT{$pdbid}{fposq}-1;
			my $lpos = $BLASTOUT{$pdbid}{lposq}-1;	
			my $aa_query = substr($aa_query_all,$fpos,$lenq);
					
			next if ($conserved_aa and $aa_query !~ /$conserved_aa/i);
			
			my $exposedq;			
			for (my $i=$BLASTOUT{$pdbid}{fposq}; $i<=$BLASTOUT{$pdbid}{lposq};$i++){
				if($refaccq{1}{$i}{acc}>=0.1){
					$exposedq++;
				}			
			}
			my $percent_exposed_q=($exposedq/$lenq)*100;
			
			
			# Target
					
			my ($ca_str_target,$aa_target, $refcat,$refaat,$refacct) = DSSPFileExtract($pdbid,\@code,$dsspfile_target, $chaint, $BLASTOUT{$pdbid}{fpost}, $BLASTOUT{$pdbid}{lpost},0);
			
			next if ($conserved_aa and $aa_target !~ /$conserved_aa/i);	
			
			my %refcat = %{$refcat};
			my %refaat = %{$refaat};
			my %refacct = %{$refacct};

			my $lent = $BLASTOUT{$pdbid}{lpost} - $BLASTOUT{$pdbid}{fpost} + 1;
					
			my $exposedt;
			for (my $i=$BLASTOUT{$pdbid}{fpost}; $i<=$BLASTOUT{$pdbid}{lpost};$i++){
				if($refacct{$pdbid}{$i}{acc}>=0.1){
					$exposedt++;
				}			
			}
			my $percent_exposed_t=($exposedt/$lent)*100;

			
			# Convert angles char to amino acid

			my @aa_query = ANGCHAR_TO_AA(1,$BLASTOUT{$pdbid}{query}, $BLASTOUT{$pdbid}{fposq}, \%refaaq);
						
			my @aa_target = ANGCHAR_TO_AA($pdbid,$BLASTOUT{$pdbid}{target}, $BLASTOUT{$pdbid}{fpost}, \%refaat);
			

			# Calculate amino acid similarity
			
			$aa_query[0] =~ s/^\s+|\s+$//g;
			$aa_target[0] =~ s/^\s+|\s+$//g;
			$aa_query[0]=~s/\[a\-z\]/C/g;
			$aa_target[0]=~s/\[a\-z\]/C/g;
			
			
			my @aa_similarity = AA_SIMILARITY($aa_query[0],$aa_target[0]);

		
			my $start=-100;
			my $end=-100;	
			if ($conserved_aa){
				my $startq;
				my $endq;
				my $startt;
				my $endt;
				while ($aa_similarity[0] =~ /$conserved_aa/pgi) {
						$startq = $-[0];
						$endq = $+[0];
					while ($aa_similarity[1] =~ /$conserved_aa/pgi){
						$startt = $-[0];
						$endt = $+[0];
						if ($startt==$startq){
							$start=$startt;
							$end=$endt;
							last;
						}
					}
				}	
			}
			next if ($conserved_aa and $start==-100);
			
			# RMSD and SAS calculation -angles-Ca
			my $rmsd = "rmsd".$pdbid;
			my $sas = "sas".$pdbid;		
			my $rmsdca = "rmsdca".$pdbid;
			my $sasca = "sasca".$pdbid;
		
			my $pdbq_path="output/".$q_id."_".$job;
			my $rmsd_qfile = "output/query_sa".$pdbid."_".$job;			
			my $pdbq = PDBFileExtract($pdbq_path,$chain_q,$BLASTOUT{$pdbid}{fposq}, $BLASTOUT{$pdbid}{lposq});
			
			my $pdbt_path="../databases/pdb_files/$pdbidt.pdb";	
			my $rmsd_tfile = "output/target_sa".$pdbid."_".$job;
			my $pdbt = PDBFileExtract($pdbt_path,$chaint,$BLASTOUT{$pdbid}{fpost}, $BLASTOUT{$pdbid}{lpost});		

				
			my $align_len;
			($returnValues{$rmsd}, $returnValues{$sas},$align_len) = RMSD_SAS_CALC($pdbid,$BLASTOUT{$pdbid}{query},$BLASTOUT{$pdbid}{target}, $BLASTOUT{$pdbid}{fposq},$BLASTOUT{$pdbid}{fpost},\%refcaq, \%refcat,$rmsd_qfile,$rmsd_tfile,$pdbq,$pdbt); 

			
			my $mol1= Bio::PDB::Structure::Molecule -> new;
			my $mol2= Bio::PDB::Structure::Molecule -> new;
			
			$mol1 -> read($rmsd_qfile);         #read the first model
			$mol2 -> read($rmsd_tfile);         #read the second model
			my $mol1b = $mol1 -> backbone;              #create a list with the backbone of mol1
			my $mol2b = $mol2 -> backbone;              #create a list with the backbone of mol2
			my @transform = $mol2b ->superpose($mol1b); #compute alignment of mol2 to mol1
			$mol2 ->rotate_translate(@transform);    #rotate and translate mol2
			my $rmsdcav = $mol2 -> rmsd($mol1);            #compute the rmsd between mol2 and mol1
			
			$returnValues{$rmsdca} = sprintf("%.2f",$rmsdcav);
			$returnValues{$sasca} = sprintf("%.2f",(100*$rmsdcav)/$align_len);
			
			#solvent-accessible surface area
			my $dsas = "dsas".$pdbid;
			my $percent_exposed = abs($percent_exposed_q-$percent_exposed_t);
			$returnValues{$dsas} = sprintf("%.2f",$percent_exposed);
			
			
			# Re-edit - structural alignment
			my $sa = "sa".$pdbid;
			$returnValues{$sa} = SA_REDIT($aa_similarity[0],$aa_similarity[1],$aa_similarity[2], $BLASTOUT{$pdbid}{fposq},$BLASTOUT{$pdbid}{fpost},$BLASTOUT{$pdbid}{lposq},$BLASTOUT{$pdbid}{lpost},$start,$end);	
		}# while
});

	my $end_time = [gettimeofday];
	my $time_interval = tv_interval($start_time, $end_time);
	my ($time_interval2) = $time_interval =~ /(\d+\.\d{1})\d+/;
	
	Print_Results($q_id, $pdbidq, $chain_q ,$blastout, $time_interval2, \%BLASTOUT, \%returnValues,\%expv,$job,$hits);
	system ("rm $dsspfile");


# ----------------------- #
#        FUNCTIONS        #
# ----------------------- #

sub DSSPFileExtract {
	my ($i,$code,$dsspfile, $chain, $firstres, $lastres,$ca) = @_;
	my $c =0; my %refca; my %refaa; my %refacc; my $aa = ""; my $ca_str = "";
	my @code = @$code; my $exposed=0; my $count_res=0; my $count=0;
	my $all_aa_str = "CSTPAGNDEQHRKMILVFYW";
	my @all_aas = split(//,$all_aa_str);
	my %an_cnvd;
	foreach (@all_aas){
		$an_cnvd{$_}++;
	}
	open (DSSP, "$dsspfile")or die "cannot open $dsspfile :$!\n";	
	while (<DSSP>){	
		next if (($_=~/\.$/) or ($_=~/\#/));
		my $chain_here = substr ($_, 11, 1);
		my $chain_near = substr ($_, 13, 1);
		next if $_=~/!/;
		if (($chain_here eq ' ') and ($chain_near ne '!')){$chain_here='A';}
		next if ($chain_here!~/$chain/);		
		
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
		
		$count_res++;
		
		my $num_residue = substr ($_,5,5);
		$num_residue =~ s/^\s+|\s+$//g;

		if (!$ca){
			$num_residue = $count_res;
		}

		if ($firstres){
			next if (($num_residue < $firstres)and($firstres =~/\w/));
		}
		if ($lastres){
			last if (($num_residue > $lastres)and($lastres=~/\w/));
		}
		$count++;		
						
        my $acc=substr ($_,35,3);
        $acc=~s/\s//g;
                
        my $relative_acc= get_relative_acc($acc,$res);
		
		$aa.=$res;
		$refca{$i}{$count_res}{ca} = $new_ca;											
		$refaa{$i}{$count_res}{aa} = $res;
		$refacc{$i}{$count_res}{acc} = $relative_acc;
                					
		next unless $ca;

		$c++;
		$ca_str .= $code[$new_ca]; 
		if ($c==60) {
		  $ca_str = $ca_str."\n"; 
		  $c = 0;
		}
	}
	close DSSP;
	return ($ca_str, $aa,\%refca,\%refaa,\%refacc);
}

sub get_relative_acc{
    
    my ($acc,$res)=@_;
    
    my %amino_acid = qw(
		A  	115  
		R 	225  
		D 	150  
		N 	160  
		C 	135
		E 	190  
		Q 	180  
		G 	75 	 
		H 	195  
		I 	175  
		L 	170  
		K 	200  
		M 	185  
		F 	210  
		P 	145  
		S 	115  
		T 	140  
		W 	255  
		Y 	230  
		V 	155
        X   145
	);
    
    my $relative_acc=$acc/$amino_acid{$res};
    return $relative_acc;
    
}


sub BLASTOUT_EXTRACT {
	my @blastoutfile =@_;
	my %hash;
	my %pdbid;
	my %expv;	
	open (BLASTOUT, "$blastoutfile[0]")or die "cannot open $blastoutfile[0] :$!\n";
	while (<BLASTOUT>){
		chomp $_;
		next unless $_=~/^>/;
		my @title = split (/\//,$_);
		my $pdbid = substr($title[0],1,6);
		$pdbid =~ s/^\s+|\s+$//g;
		$pdbid{$pdbid}++;
		$hash{$pdbid}{pdbid}=$pdbid;
		my $last_pos = $title[1];
		$last_pos=~ s/^\s+|\s+$//g;
		$hash{$pdbid}{title} = $title[2];
		my $line=<BLASTOUT>;
		chomp $line;
		$line =~ s/^\s+|\s+$//g;
		if ($line !~/Length/){
			$hash{$pdbid}{title}.= $line;
			$line=<BLASTOUT>;
			chomp $line;
		}

		my @len=split(/=/,$line);
		$hash{$pdbid}{length}=$len[1];
		$line=<BLASTOUT>;
		$line=<BLASTOUT>;
		chomp $line;
		my @paraf=split(/,/,$line);
		my @score =split(/=/,$paraf[0]);
		my @score2 =split(/\s+/,$score[1]);
		$hash{$pdbid}{score} = $score2[1];
		my @exp =split(/=/,$paraf[1]);
		$hash{$pdbid}{expv} =sprintf("%.1000f", $exp[1]);
		if ($hash{$pdbid}{expv}>1){
			$hash{$pdbid}{exp}= ">1";
		}
		else{
			$hash{$pdbid}{exp}= $exp[1];
		}
		$expv{$pdbid}=$hash{$pdbid}{exp};
		#$expv{$pdbid}='e-1000000' if $exp[1] eq '0.0';
		$line=<BLASTOUT>;
		chomp $line;
		my @paras=split(/,/,$line);
		my @positive =split(/=/,$paras[1]);
		my @positive2 =split(/\s+/,$positive[1]);
		my @positive3= split(/%/,$positive2[2]); 
		$hash{$pdbid}{positive} = substr($positive3[0],1);
		$line=<BLASTOUT>;

		my $sa;
		my $a=0;
		my $query;
		my $target;	
		while (<BLASTOUT>){
			chomp($_);
			my @line;
			$a++;
			if (($a % 4) and ($_ eq "")){
				last;
			}
			else{
				@line=split(/\s+/,$_);
				if ($a==1){
					$hash{$pdbid}{fposq} = $line[1];
				}
				elsif ($a==3){
					$hash{$pdbid}{fpost} = $line[1];
				}
				if ((defined($line[0])) && ($line[0] =~/^Query/)){
					$query.=$line[2];				
				}
				if ((defined($line[0])) && ($line[0] =~/^Sbjct/)){
					$target.=$line[2];				
				}
			}
			my $b=4*int($a/4)+1;
			if ($b==$a){
				$hash{$pdbid}{lposq} = $line[-1];
			}
			elsif ($_){
				$hash{$pdbid}{lpost} = $line[-1];
			}
		}

	$hash{$pdbid}{query}=$query;
	$hash{$pdbid}{target}=$target;
	$hash{$pdbid}{query}=~ s/^\s+|\s+$//g;
	$hash{$pdbid}{target}=~ s/^\s+|\s+$//g;
	}
	close BLASTOUT;
	return (\%pdbid,\%expv,\%hash);
}



sub PDBFileExtract {
	my ($pdbfile, $chain, $firstres, $lastres) = @_;
	my $count_res=0;
    my %pdb;
	open (PDB, "$pdbfile")or die "cannot open $pdbfile :$!\n";	

	while (<PDB>){
	
		next if !/^ATOM/;

		my $atom = substr ($_,13,4);
				
		next unless $atom =~ /CA/;
		
		my $chain_here = substr ($_,21,1);

		if ($chain_here eq ' '){$chain_here='A';}
		
		next if ($chain_here!~/$chain/);
		
		$count_res++;
		next if (($count_res < $firstres)and($firstres =~/\w/));
		last if (($count_res > $lastres)and($lastres=~/\w/));
		
		$pdb{$count_res}=$_;		
	}
	close PDB;
	return (\%pdb);
}


sub ANGCHAR_TO_AA{
	my ($i,$angchar, $fpos, $ref) = @_;
	my @alignseq=split(//,$angchar);
	my %ref = %{$ref};
	
	my $aa_align;
	foreach(@alignseq){
		if ($_ eq "-"){
			$aa_align.="-";
		}
		else{
			$aa_align.=$ref{$i}{$fpos}{aa};
			$fpos++;
		}
		
	}
	return ($aa_align);
}


sub AA_SIMILARITY{
	my ($aa_query,$aa_target) = @_;
	
	# IMPORT AMINO-ACID SCORE MATRIX 
	my $ref_score_matrix = &return_blosum62;
	my @score_matrix = @{$ref_score_matrix};
	
	my @alignq=split(//,$aa_query);
	my @alignt=split(//,$aa_target);
	
	# Hash for conversion of amino acid to number
	my $all_aa_str = "CSTPAGNDEQHRKMILVFYW-X";
	my @all_aas = split(//,$all_aa_str);
	my %an_cnvd;
	my $num = 0;
	foreach (@all_aas){
		$an_cnvd{$_} = $num;
		$num++;
	}
	
	my $len = @alignq;
	for (my $i; $i<$len;$i++){
		if (!$an_cnvd{$alignq[$i]}){
			$alignq[$i]="X";
		}
	}

	my $lent = @alignt;
	for (my $i; $i<$lent;$i++){
		if (!$an_cnvd{$alignt[$i]}){
			$alignt[$i]="X";
		}
	}
	
	# Convert amino acid alphapet to number
	$_ = $an_cnvd{$_} foreach (@alignq);
	$_ = $an_cnvd{$_} foreach (@alignt);

	my @homol;
	my $match_num = 0;
	my $gap_num = 0;
	my $sim_num = 0;
	foreach (0..$#alignq){
		if ($alignq[$_] eq $alignt[$_]){
			push @homol,'|';
			$match_num++;
		}elsif($alignq[$_] eq '20' || $alignt[$_] eq '20'){
			push @homol,' ';
			$gap_num++;
		}elsif($alignq[$_] ne $alignt[$_]){
			my $scorem = $score_matrix[$alignq[$_]][$alignt[$_]];
			if ($scorem == 2 || $scorem == 3){
				push @homol,':';
				$sim_num++;
			}elsif($scorem == 0 || $scorem == 1){
				push @homol,'.';
				$sim_num++;
			}else{
				push @homol,' ';
			}
		}
	}
	
	# Convert number back to amino acid or "-"	
	foreach (@alignq){
		$_ = $_ eq '-'?'-':$all_aas[$_];
	}
	my $alignq = join('',@alignq);
	
	foreach (@alignt){
		$_ = $_ eq '-'?'-':$all_aas[$_];
	}
	my $alignt = join('',@alignt);
	
	my $homol = join('',@homol);
	
	return ($alignq, $alignt,$homol);
}

sub RMSD_SAS_CALC{
	
	my ($i,$query, $target, $fposq, $fpost, $refq, $reft,$rmsd_qfile,$rmsd_tfile,$pdbq,$pdbt) = @_;
		
	my @alignq=split(//,$query);	
	my @alignt=split(//,$target);
	my %refq = %{$refq};
	my %reft = %{$reft};
	my %pdbq = %{$pdbq};
	my %pdbt = %{$pdbt};

	
	open (RMSDQ, ">$rmsd_qfile")or die "cannot open $rmsd_qfile :$!\n";
	open (RMSDT, ">$rmsd_tfile")or die "cannot open $rmsd_tfile :$!\n";	
	
	my $align_len = @alignq;
	my $num = 0;
	my $diff = 0;
	for (my $j=0; $j<$align_len; $j++){
		if (($alignt[$j] eq "-") && !($alignq[$j] eq "-")){
			$fposq++;
			next;
		}
		elsif (($alignq[$j] eq "-") && !($alignt[$j] eq "-")){
			$fpost++;
			next;
		}
		elsif (($alignq[$j] eq "-") && ($alignt[$j] eq "-")){
			next;
		}
		
		$diff += (($refq{1}{$fposq}{ca} - $reft{$i}{$fpost}{ca})**2);
		print RMSDQ "$pdbq{$fposq}";
		print RMSDT "$pdbt{$fpost}";
		$fposq++;
		$fpost++;
		$num++;

	}
	
	my $rmsd = sprintf("%.2f",sqrt($diff/$num));
	my $sas = sprintf("%.2f",(100*$rmsd)/$align_len);
	
	close RMSDQ;
	close RMSDT;
	return ($rmsd, $sas,$align_len);

}	

sub SA_REDIT{
	my ($alignq,$alignt,$homol,$fposq,$fpost,$lposq,$lpost,$start,$end)=@_;
	my @alignq=split(//,$alignq);
	my @alignt=split(//,$alignt);	
	my @homol=split(//,$homol);
	
	my $pos = 0;		
	my $k=1;
	my $q_str;
	my $o_str;
	my $d_str;
	my $out_str;
	my $num=0;
	
	my $a= $fposq;
	$a =~ s/^\s+|\s+$//g;
	my $aq = 5 - length ($a);
	$q_str.= "Query: $a";
	$q_str .= " " for (1..$aq);
	$o_str .= "            "; 
	my $e= $fpost;
	$e =~ s/^\s+|\s+$//g;
	my $bs = 5 - length ($e);
	$d_str.= "Sbjct: $e";
	$d_str .= " " for (1..$bs);	
	foreach(@alignq){
			if ($pos==$start) {
				$q_str.= "<b>";
				$d_str.= "<b>";
			}
			if ($pos==$end) {
				$q_str.= "</b>";
				$d_str.= "</b>";
			}
			if ($k==60) {				
				$q_str.= $_;
				my $b=$a+59;
				$q_str.= "   $b\n";
				$o_str .= shift @homol;
				$o_str.= "\n";
				$d_str.= shift @alignt;
				my $f=$e+59;
				$d_str.= "   $f\n\n";
				$out_str.=$q_str.$o_str.$d_str;
				$q_str=();
				$d_str=();
				$o_str=();			
				$k = 0;
				$num=$num+60;
				$a= $a+60;
				$aq = 5 - length ($a);
				$q_str.= "Query: $a";
				$q_str.= " " for (1..$aq);			 						
				$o_str .= "            "; 						 
				$e=$e+60;
				$bs = 5 - length ($e);
				$d_str.= "Sbjct: $e";
				$d_str.= " " for (1..$bs);				
			}
			else{
				$q_str.= $_;
				$o_str .= shift @homol;
				$d_str .= shift @alignt;
			}
			$k++;
			$pos++;
	}
	$q_str.= "   $lposq\n";
	$o_str.= "\n";
	$d_str.= "   $lpost\n\n";
	$out_str.=$q_str.$o_str.$d_str;	

	return ($out_str);
}


sub Print_Results{

	my ($q_id, $pdbidq, $chain_q ,$blastoutf, $time_interval2,$BLASTOUT, $returnValues, $expv,$job,$hits) = @_;
	
	my %BLASTOUT = %{$BLASTOUT};
	my %returnValues = %{$returnValues};
	my %expv = %{$expv};

	
	open (OUT, ">output/results_$job");

	my $numseqdb;
	open (BLASTOUT, "$blastoutf") or die "error\n";
	while (<BLASTOUT>){
		chomp $_;
		if ($_=~/sequences;/){
			my @s= split(/\s+/,$_);
			$numseqdb = $s[1];
			last;
		}
	}
	close BLASTOUT;

	$pdbidq=substr($q_id,0,4).$chain_q;
	
	print OUT "Results-$pdbidq\n\nNumber of protein structures in database: $numseqdb.	Number of hits: $hits.	Process time: $time_interval2 sec.\n\nSummery\n\n";
	print OUT "No\tPDBID\tLength\tScore\tE-value\tPositive\tcangRMSD\tcangSAS\tRMSD\tSAS\tDSASA\tTitle\n";
	
	my $numc = 0;
	foreach my $pdbid (sort { $expv{$a} <=> $expv{$b} } keys %expv) {
		my $rmsd = "rmsd".$pdbid;
		my $sas = "sas".$pdbid;
		my $rmsdca = "rmsdca".$pdbid;
		my $sasca = "sasca".$pdbid;	
		my $dsas = "dsas".$pdbid;	
		next if !$returnValues{$sas};
		$numc++;
		
		print OUT "$numc\t$BLASTOUT{$pdbid}{pdbid}\t$BLASTOUT{$pdbid}{length}\t$BLASTOUT{$pdbid}{score}\t$BLASTOUT{$pdbid}{exp}\t$BLASTOUT{$pdbid}{positive}\t$returnValues{$rmsdca}\t$returnValues{$sasca}\t$returnValues{$rmsd}\t$returnValues{$sas}\t$returnValues{$dsas}\t$BLASTOUT{$pdbid}{title}\n";
	}


	print OUT "\nPairwise Structural Alignment\n\n";
	
	my $numc2=0;
	foreach my $pdbid (sort { $expv{$a} <=> $expv{$b} } keys %expv) {
		my $rmsd = "rmsd".$pdbid;
		my $sas = "sas".$pdbid;		
		my $sa = "sa".$pdbid;
		my $rmsdca = "rmsdca".$pdbid;
		my $sasca = "sasca".$pdbid;
		my $dsas = "dsas".$pdbid;			
		next if !$returnValues{$sas};
		$numc2++;
	
		print OUT "$numc2. Query: $pdbidq Sbjct: $BLASTOUT{$pdbid}{pdbid} Score = $BLASTOUT{$pdbid}{score} E-value = $BLASTOUT{$pdbid}{exp} cangRMSD = $returnValues{$rmsd} cangSAS = $returnValues{$sas} RMSD = $returnValues{$rmsdca} SAS = $returnValues{$sasca} DSASA = $returnValues{$dsas}\n\n$returnValues{$sa}\n";
	}

	close OUT;
	
}
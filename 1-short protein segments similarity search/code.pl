#!/usr/bin/perl


$sttime = time;


###input parameters

my $input="";
my $chain="";
my $tolerance; 
my $aa_seq='';
my $save_aa= 0;
my $firstres;
my $lastres;

my $job = $$;

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
	if ($outer_key eq "-i" or $outer_key eq "-c" or $outer_key eq "-t" or $outer_key eq "-a" or $outer_key eq "-f" or $outer_key eq "-l"){
		if ($outer_key eq "-i"){
			$input=$hash{-i};
		}
		elsif ($outer_key eq "-c"){
			$chain=$hash{-c};
		}
		elsif ($outer_key eq "-t"){
			$tolerance=$hash{-t};
		}
		elsif ($outer_key eq "-a"){
			$aa_seq=$hash{-a};
			$aa_seq = uc $aa_seq;
		}
		elsif ($outer_key eq "-f"){
			$firstres=$hash{-f};
		}
		elsif ($outer_key eq "-l"){
			$lastres=$hash{-l};
		}
	}
	else{
		print "Error! Please check your command.\n"; 
		exit;
	}
}


$h = '<html><p>';


if (($aa_seq!~/[QWERTYIOPASDFGHKLCVNM]/)and($aa_seq ne '')){
  print "Error! Please enter a valid amino acid (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, and Y) or leave blank\n"; 
  exit;
}

if (($chain!~/\w/)and ($chain ne '')){
  print "Error! Please enter an alphabetical value for 'Chain' or leave blank\n"; 
  exit;
}

if (($firstres!~/\d+/)and ($firstres ne '')){
  print "Please enter a numerical value for 'First Residue or leave blank'\n"; 
  exit;
}


if (($lastres!~/\d+/)and ($lastres ne '')){
  print "Please enter a numerical value for 'Last Residue' or leave blank\n";  
  exit;
}


my @results = ();

if ($input=~/\//){
	@path=split(/\//,$input);
	$len=@path;
	$file=$path[$len-1];
}
else{
	$file=$input;
}

$input_dssp = $file.".dssp";
system ("dssp $input $input_dssp");


####Encoding input
@code=("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","B0","B1","B2","B3","B4","B5","B6","B7","B8","B9","C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","D0","D1","D2","D3","D4","D5","D6","D7","D8","D9","E0","E1","E2","E3","E4","E5","E6","E7","E8","E9","F0","F1","F2","F3","F4","F5","F6","F7","F8","F9","G0","G1","G2","G3","G4","G5","G6","G7","G8","G9","H0","H1","H2","H3","H4","H5","H6","H7","H8","H9","I0","I1","I2","I3","I4","I5","I6","I7","I8","I9","J0","J1","J2","J3","J4","J5","J6","J7","J8","J9","K0","K1","K2","K3","K4","K5","K6","K7","K8","K9","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9","M0","M1","M2","M3","M4","M5","M6","M7","M8","M9","N0","N1","N2","N3","N4","N5","N6","N7","N8","N9","O1","O2","O3","O4","O5","O6","O7","O8","O9","P0","P1","P2","P3","P4","P5","P6","P7","P8","P9","Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","R0","R1","R2","R3","R4","R5","R6","R7","R8","R9","a0","a1","a2","a3","a4","a5","a6","a7","a8","a9","b0","b1","b2","b3","b4","b5","b6","b7","b8","b9","c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","d0","d1","d2","d3","d4","d5","d6","d7","d8","d9","e0","e1","e2","e3","e4","e5","e6","e7","e8","e9","f0","f1","f2","f3","f4","f5","f6","f7","f8","f9","g0","g1","g2","g3","g4","g5","g6","g7","g8","g9","h0","h1","h2","h3","h4","h5","h6","h7","h8","h9","i0","i1","i2","i3","i4","i5","i6","i7","i8","i9","j0","j1","j2","j3","j4","j5","j6","j7","j8","j9","k0","k1","k2","k3","k4","k5","k6","k7","k8","k9","l0","l1","l2","l3","l4","l5","l6","l7","l8","l9","m0","m1","m2","m3","m4","m5","m6","m7","m8","m9","n0","n1","n2","n3","n4","n5","n6","n7","n8","n9","o1","o2","o3","o4","o5","o6","o7","o8","o9","p0","p1","p2","p3","p4","p5","p6","p7","p8","p9","q0","q1","q2","q3","q4","q5","q6","q7","q8","q9","r0","r1","r2","r3","r4","r5","r6","r7","r8","r9","s0");

open (DSSP, "$input_dssp")||die "cannot open input dssp:$!\n";


$aa = "";  $phi = "";  $psi = "";


my $refi = {};
my $counti= 1;


while (<DSSP>) { 

  next if (/\.$/||/\#/);

  $residue = substr ($_,6,4); $residue+=0; 
  next if (($residue < $firstres)and($firstres=~/\w/));
  next if (($residue > $lastres)and($lastres=~/\w/));
  $chain_here = substr ($_, 11, 1);
  if ($chain_here eq ' '){$remember=1}
  next if ($chain_here!~/$chain/);
  if ($remember == 1){$aa.=".";$phi.="..";$psi.="..";$remember=0;}
  $aa= $aa. substr($_,13,1);

  $new_phi =  substr ($_, 103, 6);
  if ($new_phi==360) {
     $new_phi="0"; 
  }
  $new_phi =int($new_phi);
  
  $refi -> {$counti} -> {phi} = $new_phi;
  
  $phi = $phi . $code[$new_phi] if ($tolerance eq 0);
  if ($tolerance  ne 0){
    $phi_tol= "(";   
    for ($i=-$tolerance;$i<=$tolerance;$i++){
       $phi_tol.=$code[$new_phi+$i]."|";
    }
    chop ($phi_tol);      $phi_tol.=")";
    $phi = $phi . $phi_tol;
  }

  $new_psi =  substr ($_, 109, 6);
  if ($new_psi==360) {
     $new_psi="0";
  }
  $new_psi=int($new_psi);
  
  $refi -> {$counti} -> {psi} = $new_psi;

  $psi = $psi . $code[$new_psi] if ($tolerance eq 0);
  if ($tolerance  ne 0){
    $psi_tol= "(";   
    for ($i=-$tolerance;$i<=$tolerance;$i++){
       $psi_tol.=$code[$new_psi+$i]."|";
    }
    chop ($psi_tol);      $psi_tol.=")";
    $psi = $psi . $psi_tol;
  }
  
  $counti ++;
}

if ($aa_seq){
  $save_aa=1;
}

$aa=~s/[a-z]/\[a-z\]/g;  #cysteines are in lower case.
$phi=~s/A0/\.\./g; #0 angles are disregarded
$psi=~s/A0/\.\./g; #0 angles are disregarded
if ($aa=~/^\./){$aa=substr($aa,1);$phi=substr($phi,2);$psi=substr($psi,2);}



unless ($phi){
  print "Chain $chain could not be found for PDB $pdbid\n"; 
  exit
}

#open output file

open (OUT, ">output/output_".$job.".txt");	



########################################
########################################
$counter=0;
open (DATA,"database/database.txt")||die "cannot open database:$!\n";
########if conserved sequence###########
if ($save_aa == 1){

while (<DATA>){

  #next unless  /$aa.*$phi.*$psi/;
  next unless  /$phi.*$psi/;
  
  next unless $aa=~/$aa_seq/;
  
  my $refj = {};
  
  my $hM = "";
  @line =split ('\|',$_);
  $line[2]=~/$phi/; $posphi=$-[0]/2;  
  $line[3]=~/$psi/; $pospsi=$-[0]/2;  
  $line[1]=~/^.{$posphi}($aa)/; $seq_match=$1; $posseq=$-[1];

  $as = substr($line[1],$posseq,length($aa));  
  $aa=~s/\[a\-z\]/C/g; 
  $as=~s/\[a\-z\]/C/g; 
  $pipes="";

  if ($pospsi eq $posphi){
  
	next unless $as=~/$aa_seq/;
	
	$as=~/$aa_seq/;$posas=$-[0];  
	$aa=~/$aa_seq/;$posaa=$-[0];

	next unless	$posas==$posaa;
	
	$seq_match=$as;
	
    $counter++;
    
	$f ="../databases/dssp_files/".$line[0].".dssp";
    #open the hit dssp file.
    open (DSS, $f) || print "cannot open DSSP reference file:$!\n";
		
    my $pos_end = $posphi+(length($phi)/2);
    my $pos_init = $posphi+1;
        
    while (<DSS>){
      next if (/\.$/||/\#/);
      $id2 = substr ($_, 0, 5); 
      $id2-=0;

      if ($posphi+1 == $id2 && $pos_init <= $pos_end){
      	setAnglePos($refj,$id2,'phi',int(substr ($_, 103, 6)));
      	setAnglePos($refj,$id2,'psi',int(substr ($_, 109, 6)));
      	
        $chain_type= substr ($_, 11, 1);
        $pos_init ++;
      }
    }
	close (DSS);
   
    my $rmsd = calcRMSD($refi,$refj);
    
    $pinc = $posseq+1;
    $plength = $posseq+length($aa);
		
	
    $hM .= "\nPDB: $line[0]";
    $hM .= " (Chain $chain_type) " if $chain_type;
	$a=$posseq+1;
	$b=$posseq+length($aa);
    $hM .= "position: ". $a .":". $b ."\tRMSD:  $rmsd \n";
    $aa=~s/\[a\-z\]/C/g; $seq_match=~s/[a-z]/C/g;
    @aas=split(//,$aa); @seq_matches=split(//,$seq_match);
    foreach ($i=0;$i<=$#aas;$i++){
      $pipes.=" " if ($aas[$i] ne $seq_matches[$i]); 
      $pipes.="|" if ($aas[$i] eq $seq_matches[$i]);
    } 
    $hM .=  "Query:	$aa\n		$pipes\nResult:	$seq_match";

    $h .= $hM;

    push (@results, $count_aa."+".$hM);		

  }
}
}


#########if not conserved sequence##########
if ($save_aa == 0){

while (<DATA>){

  next unless  /$phi.*$psi/;
  my $refj = {};
  
   my $hM = "";

  @line =split ('\|',$_); 
  
  $line[2]=~/$phi/; $posphi=$-[0]/2;
  $line[3]=~/$psi/; $pospsi=$-[0]/2; 
  $aa=~s/\[a\-z\]/C/g; 

  $line[1]=substr($line[1],$posphi,length($aa));
  $pipes="";
  
  if ($pospsi eq $posphi){
    $counter++;
	
	$f ="../databases/dssp_files/".$line[0].".dssp";
    #open the hit dssp file.
    open (DSS, $f) || print "cannot open DSSP reference file:$!\n";
       
    my $pos_end = $posphi+(length($phi)/2);
    my $pos_init = $posphi+1;
    
    while (<DSS>){
      next if (/\.$/||/\#/);
      
      $id2 = substr ($_, 0, 5); $id2-=0;
      if ($posphi+1 == $id2) { $pos_x = substr ($_, 5, 5); $pos_x-=0;}           
      if ($posphi+1 <= $id2 && $pos_init <= $pos_end){
      	
      	setAnglePos($refj,$id2,'phi',int(substr ($_, 103, 6)));
      	setAnglePos($refj,$id2,'psi',int(substr ($_, 109, 6)));
      	
        $chain_type= substr ($_, 11, 1);
        $pos_init ++;
        
      }
    }
	close (DSS);
    
    my $rmsd = calcRMSD($refi,$refj);
        
    $hM .= "\nPDB: $line[0]";
    $hM .= " (Chain $chain_type)" if $chain_type;
 
    $pinc = $pos_x;
    $plength = $pos_x-1+length($aa);

    $hM .= " position: $pinc:$plength	RMSD:  $rmsd \n";

    $aa=~s/\[a\-z\]/C/g; $line[1]=~s/[a-z]/C/g;
    @aas=split(//,$aa); @seq_matches=split(//,$line[1]);

    my $count_aa = 0;
	
    foreach ($i=0;$i<=$#aas;$i++){
      $pipes.=" " if ($aas[$i] ne $seq_matches[$i]);

      if ($aas[$i] eq $seq_matches[$i]){
		 $pipes.="|";
	 	$count_aa ++; 	
      }
    }

    $hM .= "Query:	$aa\n		$pipes\nResult:	$line[1]";
    $h .= $hM;
	
    push (@results, $rmsd."+".$hM);

  }
}
}

system ("rm $input $input_dssp");


print "There were no matches for your query" if ($h eq '<html><p>');
$h .= "There were no matches for your query" if ($h eq '<html><p>');

$sttime2 = time;
$time= $sttime2-$sttime;
$t="Time:". $time." Sec\n";
print OUT "$t";
printResultSortAA(sort { $a <=> $b } @results);


#--

sub calcRMSD{
	my ($refi,$refj) = @_;

   my @ri=  keys (%{$refi});
   my @rj=  keys (%{$refj});
   
   @rj = sort { $a <=> $b } @rj;
   @ri = sort { $a <=> $b } @ri;
    
   
   my $phi_diff_sqr = 0;
   my $psi_diff_sqr = 0;
      
	for (my $i=0; $i < scalar(@ri); $i++){
		my $phi_i = abs($refi -> {$ri[$i]} -> {phi});
		my $phi_j = abs($refj -> {$rj[$i]} -> {phi});
		my $psi_i = abs($refi -> {$ri[$i]} -> {psi});
		my $psi_j = abs($refj -> {$rj[$i]} -> {psi});
		
		if ($phi_i != 0 && $psi_i != 0 && $phi_j != 0 && $psi_j != 0){
			$phi_diff_sqr += ($phi_i - $phi_j)**2;
			$psi_diff_sqr += ($psi_i - $psi_j)**2;
		}
	}
	
	my $rmsd = sqrt($phi_diff_sqr + $psi_diff_sqr);
	$rmsd = sprintf("%.3f", $rmsd);
	
	return $rmsd;
}

sub setAnglePos{
	my ($ref,$id,$type,$val) = @_;
	
	 if ($val==360) {
     	    $val="0"; 
  	 }
  
        $ref -> {$id} -> {$type} = $val;
}

sub printResultSortAA{
	my (@r) = @_;
	my $cc = 1;
	
	print OUT "Results";
	
	foreach my $r(@r){
		my @a = split(/\+/,$r); 
		print OUT "$cc. $a[1]\n";
		
		$cc++;
	}
}
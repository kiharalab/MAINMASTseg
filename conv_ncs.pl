#!/usr/bin/perl
if(@ARGV!=1){
 print "$0 [file]\n";
 exit;
}

($file,$path)=@ARGV;


@B=firstfile("$file");

$cnt=0;
foreach my $l (@B){
 @MTX=() if($$l[0] eq "new_operator");

 if($$l[0] eq "rota_matrix"){
  push @MTX,($$l[1],$$l[2],$$l[3]);
  #print "@MTX\n";
 }
 if($$l[0] eq "tran_orth"){
	 $cnt++;
  printf("REMARK 350   BIOMT1%4d%10.6f%10.6f%10.6f%15.5f\n",$cnt,$MTX[0],$MTX[1],$MTX[2],$$l[1]);
  printf("REMARK 350   BIOMT2%4d%10.6f%10.6f%10.6f%15.5f\n",$cnt,$MTX[3],$MTX[4],$MTX[5],$$l[2]);
  printf("REMARK 350   BIOMT3%4d%10.6f%10.6f%10.6f%15.5f\n",$cnt,$MTX[6],$MTX[7],$MTX[8],$$l[3]);
 }
}


#@key = sort { $hash{$a} <=> $hash{$b} || $a <=> $b} keys %hash;


sub firstfile{
my $cnt=0;
my @A;
open(IN,$_[0]) or die;
while(<IN>){
  next if(/^#/);
 chomp;
 my $item;
 @{$item}=split(/[\s\t]+/,$_);
 push @A, $item
}
close(IN);
return @A;
}
sub onetothree{
 %amin123 = ("W"=>"TRP","F"=>"PHE","Y"=>"TYR","L"=>"LEU","I"=>"ILE","V"=>"VAL","M"=>"MET","A"=>"ALA","G"=>"GLY","P"=>"PRO","C"=>"CYS","T"=>"THR","S"=>"SER","Q"=>"GLN","N"=>"ASN","E"=>"GLU","D"=>"ASP","H"=>"HIS","K"=>"LYS","R"=>"ARG");
 %amin321 = ("TRP"=>"W","PHE"=>"F","TYR"=>"Y","LEU"=>"L","ILE"=>"I","VAL"=>"V","MET"=>"M","ALA"=>"A","GLY"=>"G","PRO"=>"P","CYS"=>"C","THR"=>"T","SER"=>"S","GLN"=>"Q","ASN"=>"N","GLU"=>"E","ASP"=>"D","HIS"=>"H","LYS"=>"K","ARG"=>"R");
}
sub firstfile_line{
my $cnt=0;
my @A;
open(IN,$_[0]) or die;
while(<IN>){
  next if(/^#/);
 chomp;
 push @A, $_;
}
close(IN);
return @A;
}

sub readpdb{
 my $cnt=0;
 my @A;
 my ($x,$y,$z);
 my ($file)=@_;

 open(IN,$file) or die;
 while(<IN>){
  next unless(/^ATOM/);
  chomp;

  $x=substr($_,30,8);
  $y=substr($_,38,8);
  $z=substr($_,46,8);

  my $atm=substr($_,13,3);
  my $res=substr($_,17,3);
  my $rnum=substr($_,22,4);
  #my $m_tag=substr($_,17,9);
  my $m_tag=substr($_,13,13);

  my $item;
  @{$item}=($res,$atm,$x,$y,$z,$rnum,$m_tag);
  push @A, $item;
 }
 close(IN);
 return @A;
}


sub Ave{
 my $inp=$_[0];
 my $n=0;
 my $sum=0;
 foreach $l (@{$inp}){

  $sum+=$l;
  $n++;
 }
 return ($sum/$n)
}

sub Std{
 my $inp=$_[0];
 my $n=0;
 my $sum=0;
 my $ave=Ave($inp);
 foreach $l (@{$inp}){
  
  $sum+=($l-$ave)*($l-$ave);
  $n++;
 }
 return sqrt($sum/$n)
}


#!/usr/bin/perl -w
use strict;
use threads; 
use Config;
my $Avr_read_lengh=$ARGV[1];#read_length
my $window=2*$Avr_read_lengh-16;
my %chongfu;
my %chongfu2;
my $read_count_outfile=$ARGV[3];
$read_count_outfile=~s/\/$//gi;
unless (-e $read_count_outfile) {
    `mkdir $read_count_outfile`;
}
my $SRR_acc=$ARGV[3];
open(IN,"$ARGV[0]") or die "$ARGV[0]";
my $junction_file1=$read_count_outfile."/junction_sequence_length_".$Avr_read_lengh.".txt";
my $junction_file2=$read_count_outfile."/junction_sequence_uniq_length_".$Avr_read_lengh.".txt";
open OUT,">$junction_file1";
open OUT2,">$junction_file2";
print OUT "#Transcript\tExon1Start_exon2End_ch_intronStart_intronEnd\tE1I_junction\tIntron\tIE2_junction\tE1E2_junction\n";
print OUT2 "#Junction\tLength\tExon1Start_IntronStart_IntronEnd_Exon2_end\n";
my %intron_position;
while (<IN>) {
    chomp(my $line=$_);
    $line=~s/B73V4_ctg/B73V4ctg/g;
    my ($tr,$ch,$strand,$exon1_start,$exon2_end,$tr2,$intron_start,$intron_end,$tag,$position)=split /\t/, $line;
    if ($tag eq "-") {
        my $exon1_end=$intron_start-1;
        my $exon2_start=$intron_end+1;
        my $length_e1=$intron_start-$exon1_start;
        my $length_e2=$exon2_end-$intron_end;
        my $length_i=$intron_end-$intron_start+1;
        if (($length_e1>=8)&&($length_e2>=8)&&(($length_i>=8))) {
            my $intron=$ch."_".$intron_start."_".$intron_end;
            $intron_position{$intron}=$position;
            my $e1i_junction;
            my $e1i_length;
            my $window_half=$Avr_read_lengh-8;
            if ((($length_e1-$window_half)>=0)&&(($length_i-$window_half)>=0)) {
                my $e1i_start=$intron_start-($Avr_read_lengh-8);
                my $e1i_end=$intron_start+($Avr_read_lengh-8)-1;
                $e1i_junction=$ch."_".$e1i_start."_".$e1i_end;
                $e1i_length=$window
            }elsif((($length_e1-$window_half)>=0)&&(($length_i-$window_half)<0)){
                my $e1i_start=$intron_start-($Avr_read_lengh-8);
                my $e1i_end=$intron_end;
                $e1i_junction=$ch."_".$e1i_start."_".$e1i_end;
                $e1i_length=$Avr_read_lengh-8+$length_i
            }elsif((($length_e1-$window_half)<0)&&(($length_i-$window_half)>=0)){
                my $e1i_start=$exon1_start;
                my $e1i_end=$intron_start+($Avr_read_lengh-8)-1;
                $e1i_junction=$ch."_".$e1i_start."_".$e1i_end;
                $e1i_length=$Avr_read_lengh-8+$length_e1
            }elsif((($length_e1-$window_half)<0)&&(($length_i-$window_half)<0)){
                  $e1i_junction=$ch."_".$exon1_start."_".$intron_end;
                  $e1i_length=$intron_end-$exon1_start+1;
            }
            my $e1e2_junction;
            my $e1e2_length;
            if ((($length_e1-$window_half)>=0)&&(($length_e2-$window_half)>=0)) {
                my $e1e2_start=$intron_start-($Avr_read_lengh-8);
                my $e1e2_end=$intron_end+($Avr_read_lengh-8);
                $e1e2_junction=$ch."_".$e1e2_start."_".$e1e2_end;
                $e1e2_length=$window
            }elsif((($length_e1-$window_half)>=0)&&(($length_e2-$window_half)<0)){
                my $e1e2_start=$intron_start-($Avr_read_lengh-8);
                my $e1e2_end=$exon2_end;
                $e1e2_junction=$ch."_".$e1e2_start."_".$e1e2_end;
                $e1e2_length=$Avr_read_lengh-8+$length_e2
            }elsif((($length_e1-$window_half)<0)&&(($length_e2-$window_half)>=0)){
                my $e1e2_start=$exon1_start;
                my $e1e2_end=$intron_end+($Avr_read_lengh-8);
                $e1e2_junction=$ch."_".$e1e2_start."_".$e1e2_end;
                $e1e2_length=$Avr_read_lengh-8+$length_e1
            }elsif((($length_e1-$window_half)<0)&&(($length_e2-$window_half)<0)){
                  $e1e2_junction=$ch."_".$exon1_start."_".$exon2_end;
                  $e1e2_length=$length_e1+$length_e2;
            }
            my $ie2_junction;
            my $ie2_length;
            if ((($length_i-$window_half)>=0)&&(($length_e2-$window_half)>=0)) {
                my $ie2_start=$intron_end-($Avr_read_lengh-8)+1;
                my $ie2_end=$intron_end+$Avr_read_lengh-8;
                $ie2_junction=$ch."_".$ie2_start."_".$ie2_end;
                $ie2_length=$window;
            }elsif((($length_i-$window_half)>=0)&&(($length_e2-$window_half)<0)){
                my $ie2_start=$intron_end-($Avr_read_lengh-8)+1;
                my $ie2_end=$exon2_end;
                $ie2_junction=$ch."_".$ie2_start."_".$ie2_end;
                $ie2_length=$Avr_read_lengh-8+$length_e2;
            }elsif((($length_i-$window_half)<0)&&(($length_e2-$window_half)>=0)){
                my $ie2_start=$intron_start;
                my $ie2_end=$intron_end+$Avr_read_lengh-8;
                $ie2_junction=$ch."_".$ie2_start."_".$ie2_end;
                $ie2_length=$Avr_read_lengh-8+$length_i;
            }elsif((($length_i-$window_half)<0)&&(($length_e2-$window_half)<0)){
                $ie2_junction=$ch."_".$intron_start."_".$exon2_end;
                $ie2_length=$length_i+$length_e2;
            }
            my $i_junction=$ch."_".$intron_start."_".$intron_end;
            my $i_length=$intron_end-$intron_start+1;
            my $exon_intron=$exon1_start."_".$exon2_end."_".$intron;
            unless (exists $chongfu{$exon_intron}){
                    print OUT "$tr\t$exon_intron\t$e1i_junction\t$i_junction\t$ie2_junction\t$e1e2_junction\n";
                    $chongfu{$exon_intron}++;
            }
            my $string1=$exon1_start."_".$intron_start."_".$intron_end."_".$exon2_end;
            my $intron_temp="i-".$i_junction;
            my $long_line="$intron_temp\t$i_length\t$string1";
            unless (exists $chongfu2{$long_line}){
                print OUT2 "$long_line\n";
                $chongfu2{$long_line}++
            }
            my $e1i_junction_temp="e1i-".$e1i_junction;
            $long_line="$e1i_junction_temp\t$e1i_length\t$string1";
            unless (exists $chongfu2{$long_line}){
                print OUT2 "$long_line\n";
                $chongfu2{$long_line}++
            }
            my $e1e2_junction_temp="e1e2-".$e1e2_junction;
            $long_line="$e1e2_junction_temp\t$e1e2_length\t$string1";
            unless (exists $chongfu2{$long_line}){
                print OUT2 "$long_line\n";
                $chongfu2{$long_line}++
            }   
            my $ie2_junction_temp="ie2-".$ie2_junction;
            $long_line="$ie2_junction_temp\t$ie2_length\t$string1";
            unless (exists $chongfu2{$long_line}){
                print OUT2 "$long_line\n";
                $chongfu2{$long_line}++
            }
        }    
    }
}
close IN;
close OUT;
close OUT2;
my $limit_length=$Avr_read_lengh-8;
my $sam_file=$ARGV[2];
my %IR_postion;
my %long_junction;
my %tag_position;
open(IN,"$junction_file2") or die "$junction_file2";
while (<IN>) {
    chomp(my $line=$_);
    unless ($line=~m/^#/){
        my ($junction,$length,$tag_pos)=(split /\t/, $line)[0,1,2,3];
        $junction=$junction."_".$tag_pos;
        $long_junction{$junction}=$length;
        my $important=int((split /_/,$junction)[2]/1000);
        my $ch=(split /_/,$junction)[0];
        $ch=(split /-/,$ch)[1];
        my $tag=$ch."_".$important;
        my $tag1="long";
        if ($length<$window) {
            $tag1="short";
        }
        if (exists  $IR_postion{$tag}{$tag1}) {
            my $temp=$IR_postion{$tag}{$tag1};
            $IR_postion{$tag}{$tag1}="$temp,$junction"
        }else{
            $IR_postion{$tag}{$tag1}=$junction;
        }
    }
}
close IN;
my %type_value1;
my %type_value2;
my %read_pos;
# read_count
my %types_count;
my %junction_count;
my %windows;
my @canshu=@_;
#read count
    open(IN,"$sam_file") or die "can't open $sam_file";
    while (<IN>){
        chomp(my $line=$_);
        unless($line=~m/^@/){
            if ($line=~m/NH:i:1/){
                my $mismatch=0;
                if ($line=~m/XM:i:[0-9]+/){
                    $line=~m/(XM:i:[0-9]+)/;
                    $mismatch=$1;
                    $mismatch=(split /:/,$mismatch)[2];
                }
                if ($mismatch<=2) {
                    my ($read_name,$ch,$read_start,$code)=(split /\t/, $line)[0,2,3,5];
                    if ($code=~m/[0-9]+M/gi) {
                        $code=~s/D/M/gi;
                        $code=~s/[0-9]+I//gi;
                        for my $i(1...5){
                            if ($code=~m/([0-9]+M[0-9]+M)/) {
                                my $temp=$1;
                                my @number_temp=split /M/,$temp;
                                my $number_temp=$number_temp[0]+$number_temp[1];
                                $number_temp=$number_temp."M";
                                $code=~s/$temp/$number_temp/;
                            }else{
                                last
                            }
                        }
                        for my $i(1...5){
                            if ($code=~m/([0-9]+M[0-9]+M)/) {
                                my $temp=$1;
                                my @number_temp=split /M/,$temp;
                                my $number_temp=$number_temp[0]+$number_temp[1];
                                $number_temp=$number_temp."M";
                                $code=~s/$temp/$number_temp/;
                            }else{
                                last
                            }
                        }
                        $code=~s/[0-9]+[S|H]//gi;
                        my $code_temp=$code;
                        $code_temp=~s/M//gi;
                        $code_temp=~s/N//gi;
                        if ($code_temp=~m/^\d+$/g){
                            if ($code=~m/N/) {
                                my $min_value=int($read_start/1000)-5;
                                my $max_value=int($read_start/1000)+5;
                                for (my $important=$min_value;$important<=$max_value;$important++){
                                    my $tag=$ch."_".$important;
                                    my @n_cishu=($code=~m/([0-9]+N[0-9]+)/g);
                                    my $n_cishu=$#n_cishu+1;
                                    #只有一个N
                                    if ($n_cishu==1) {
                                        my @n_cishu=($code=~m/([0-9]+[A-Z])/gi);
                                        my $read_e1_length=$n_cishu[0];
                                        my $read_intron_length=$n_cishu[1];
                                        my $read_e2_length=$n_cishu[2];
                                        $read_e1_length=~s/[A-Z]//gi;
                                        $read_intron_length=~s/[A-Z]//gi;
                                        $read_e2_length=~s/[A-Z]//gi;
                                        my $read_intron_start=$read_start+$read_e1_length;
                                        my $read_exon1_end=$read_intron_start-1;
                                        my $read_intron_end=$read_start+$read_e1_length+$read_intron_length-1;
                                        my $read_exon2_end=$read_start+$read_e1_length+$read_intron_length-1+$read_e2_length;
                                        my $read_length=$read_e1_length+$read_e2_length;
                                        if (exists $IR_postion{$tag}{long}){
                                            my $junctions=$IR_postion{$tag}{long};
                                            my @junctions=split /,/,$junctions;
                                            foreach my $junction(@junctions){ 
                                                my ($junction_type,$junction_start,$junction_end)=(split /_/,$junction)[0,1,2];
                                                $junction_type=(split /-/,$junction_type)[0];
                                                my ($exon_start,$intron_start,$intron_end,$exon_end)=(split /_/,$junction)[3,4,5,6];
                                                if ($junction_type eq "e1e2"){
                                                    if (($read_start<=($intron_start-8))&&($read_exon2_end>=($intron_end+8))&&($intron_start==$read_intron_start)&&($intron_end==$read_intron_end)) {
                                                        if ($read_length == $Avr_read_lengh) {
                                                            $types_count{$junction}{$read_start}++;
                                                        }else{
                                                            $junction_count{$junction}++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if (exists $IR_postion{$tag}{short}){
                                            my $junctions=$IR_postion{$tag}{short};
                                            my @junctions=split /,/,$junctions;
                                            foreach my $junction(@junctions){
                                                my ($junction_type,$junction_start,$junction_end)=(split /_/,$junction)[0,1,2];
                                                $junction_type=(split /-/,$junction_type)[0];
                                                my ($exon_start,$intron_start,$intron_end,$exon_end)=(split /_/,$junction)[3,4,5,6];
                                                my $read_exon2_start=$read_intron_end+1;
                                                my $read_exon1_start=$read_start;
                                                if ($junction_type eq "i"){
                                                    if (($intron_end-$intron_start+1)<$Avr_read_lengh) {
                                                        if ((($read_start<=$intron_start)&&($read_exon1_end==$exon_end))||(($read_exon2_start==$exon_start)&&($read_exon2_end>=$intron_end))) {
                                                        $junction_count{$junction}++;
                                                        }elsif((($read_exon1_start==$exon_start)&&($read_exon1_end==$exon_end))||(($read_exon2_start==$exon_start)&&($read_exon2_end==$exon_end))){
                                                            $junction_count{$junction}++;
                                                        }
                                                    }
                                                }
                                                if ($junction_type eq "e1e2") {
                                                    $read_exon2_start=$read_intron_end+1;
                                                    if (($read_intron_start==$intron_start)&&($read_intron_end==$intron_end)&&($read_start<=($intron_start-8))&&($read_exon2_end>=($intron_end+8))) {
                                                        $junction_count{$junction}++;
                                                    }
                                                }
                                                if ($junction_type eq "e1i"){
                                                    my $intron_e1i=$ch."_".$intron_start."_".$intron_end;
                                                    my $intron_e1i_position=$intron_position{$intron_e1i};
                                                    if ($intron_e1i_position eq "NotFirst") {
                                                        if (($read_exon2_end>=($intron_start+7))&&($read_intron_end==($exon_start-1))&&($read_exon2_end<=$intron_end)) {
                                                            $junction_count{$junction}++;
                                                        }
                                                    }
                                                }
                                                if ($junction_type eq "ie2"){
                                                    my $intron_ie2=$ch."_".$intron_start."_".$intron_end;
                                                    my $intron_ie2_position=$intron_position{$intron_ie2};
                                                    if ($intron_ie2_position eq "NotLast") {
                                                        if (($read_start<=($intron_end-7))&&($read_exon1_end==$exon_end)&&($read_start>=$intron_start)) {
                                                            $junction_count{$junction}++;
                                                        }  
                                                    }
                                                }    
                                            }
                                        }
                                    }else{
                                        my @zuhe;
                                        my $zuhe_type;
                                        my $zuhe_index=-1;
                                        my $code_temp=$code;
                                        for(my $i=0;$i<=10;$i++){
                                            if ($code_temp=~m/N/) {
                                                $code_temp=~m/([0-9]+M[0-9]+N[0-9]+M)/;
                                                my $pattern=$1;
                                                $zuhe_index++;
                                                $zuhe[$zuhe_index]=$pattern;
                                                $pattern=~m/([0-9]+M[0-9]+N)/;
                                                $code_temp=~s/^$1//;
                                            }else{
                                                last
                                            }
                                        }
                                        my $read_start1=$read_start;
                                        my $zuhe_num=$#zuhe+1;
                                        my $zuhe_order=0;
                                        foreach my $code1(@zuhe){
                                            $zuhe_order+=1;
                                            if ($zuhe_order==1) {
                                                $zuhe_type="First"
                                            }elsif($zuhe_order==$zuhe_num){
                                                $zuhe_type="Last"
                                            }else{
                                                $zuhe_type="Other"
                                            }
                                            my @n_cishu=($code1=~m/([0-9]+[A-Z])/gi);
                                            my $read_e1_length=$n_cishu[0];
                                            my $read_intron_length=$n_cishu[1];
                                            my $read_e2_length=$n_cishu[2];
                                            $read_e1_length=~s/[A-Z]//gi;
                                            $read_intron_length=~s/[A-Z]//gi;
                                            $read_e2_length=~s/[A-Z]//gi;
                                            my $read_intron_start=$read_start1+$read_e1_length;
                                            my $read_exon1_start=$read_start1;
                                            my $read_exon1_end=$read_intron_start-1;
                                            my $read_intron_end=$read_start1+$read_e1_length+$read_intron_length-1;
                                            my $read_exon2_end=$read_start1+$read_e1_length+$read_intron_length-1+$read_e2_length;
                                            my $read_length=$read_e1_length+$read_e2_length;
                                            $read_start1=$read_intron_end+1;
                                            if (exists $IR_postion{$tag}{short}){
                                                my $junctions=$IR_postion{$tag}{short};
                                                my @junctions=split /,/,$junctions;
                                                foreach my $junction(@junctions){
                                                    my ($junction_type,$junction_start,$junction_end)=(split /_/,$junction)[0,1,2];
                                                    $junction_type=(split /-/,$junction_type)[0];
                                                    my ($exon_start,$intron_start,$intron_end,$exon_end)=(split /_/,$junction)[3,4,5,6];
                                                    my $read_exon2_start=$read_intron_end+1;
                                                    if ($junction_type eq "i"){
                                                        if (($intron_end-$intron_start+1)<$Avr_read_lengh){
                                                            if ($zuhe_type eq "First"){
                                                                if (($read_start<=$intron_start)&&($read_exon1_end==$exon_end)){
                                                                    $junction_count{$junction}++;
                                                                }
                                                            }
                                                            if ($zuhe_type eq "Last"){
                                                                if (($read_exon2_start==$exon_start)&&($read_exon2_end>=$intron_end)) {
                                                                    $junction_count{$junction}++;
                                                                }
                                                            }
                                                            if ($zuhe_type eq "Other"){
                                                                if ((($read_exon1_start==$exon_start)&&($read_exon1_end==$exon_end))||(($read_exon2_start==$exon_start)&&($read_exon2_end==$exon_end))) {
                                                                    $junction_count{$junction}++;
                                                                }
                                                            }
                                                        }
                                                    }
                                                    if ($junction_type eq "e1e2") {
                                                        $read_exon2_start=$read_intron_end+1;
                                                        if ($zuhe_type eq "First"){
                                                            if (($read_intron_start==$intron_start)&&($read_intron_end==$intron_end)&&($read_exon1_start<=($intron_start-8))&&($read_exon2_end==$exon_end)) {
                                                                $junction_count{$junction}++;
                                                            }
                                                        }
                                                        if ($zuhe_type eq "Last"){
                                                            if (($read_intron_start==$intron_start)&&($read_intron_end==$intron_end)&&($read_exon1_start==$exon_start)&&($read_exon2_end>=($intron_end+8))) {
                                                                $junction_count{$junction}++;
                                                            }
                                                        }
                                                        if ($zuhe_type eq "Other"){
                                                            if (($read_intron_start==$intron_start)&&($read_intron_end==$intron_end)&&($read_exon1_start==$exon_start)&&($read_exon2_end==$exon_end)) {
                                                                $junction_count{$junction}++;
                                                            }
                                                        }
                                                    }
                                                    if ($junction_type eq "e1i") {
                                                        if ($zuhe_type eq "Last") {
                                                            my $intron_e1i=$ch."_".$intron_start."_".$intron_end;
                                                            my $intron_e1i_position=$intron_position{$intron_e1i};
                                                            if ($intron_e1i_position eq "NotFirst") {
                                                                if (($read_exon2_end>=($intron_start+7))&&($read_intron_end==($exon_start-1))&&($read_exon2_end<=$intron_end)) {
                                                                    $junction_count{$junction}++;
                                                                }
                                                            }
                                                        }  
                                                    }
                                                    if ($junction_type eq "ie2"){
                                                        if ($zuhe_type eq "First") {
                                                            my $intron_e1i=$ch."_".$intron_start."_".$intron_end;
                                                            my $intron_e1i_position=$intron_position{$intron_e1i};
                                                            if ($intron_e1i_position eq "NotLast") {
                                                                if (($read_start<=($intron_end-7))&&($read_exon1_end==$exon_end)&&($read_start>=$intron_start)) {
                                                                    $junction_count{$junction}++;
                                                                }
                                                            }
                                                        }
                                                    }  
                                                }
                                            }
                                        }
                                    }
                                }
                            }else{
                                $code=~m/([0-9]+)M/;
                                my $read_length=$1;
                                my $read_end=$read_start+$read_length-1;
                                my $min_value=int($read_start/1000)-5;
                                my $max_value=int($read_start/1000)+5;
                                for (my $important=$min_value;$important<=$max_value;$important++){
                                    my $tag=$ch."_".$important;
                                    if (exists $IR_postion{$tag}{long}) {
                                        my $junctions=$IR_postion{$tag}{long};
                                        my @junctions=split /,/,$junctions;
                                        foreach my $junction(@junctions){
                                            my ($junction_type,$junction_start,$junction_end)=(split /_/,$junction)[0,1,2];
                                            $junction_type=(split /-/,$junction_type)[0];
                                            unless ($junction_type eq "e1e2"){
                                                if (($read_start>=$junction_start)&&($read_end<=$junction_end)) {
                                                    if ($read_length == $Avr_read_lengh) {
                                                        $types_count{$junction}{$read_start}++;
                                                    }else{
                                                        my ($exon_start,$intron_start,$intron_end,$exon_end)=(split /_/,$junction)[3,4,5,6];
                                                        if ($junction_type eq "e1i") {
                                                            if (($read_end<=$intron_end)&&($read_end>=($intron_start+7))&&($read_start<=($intron_start-8))) {
                                                                $junction_count{$junction}++;
                                                            }
                                                        }elsif($junction_type eq "ie2"){
                                                            if (($read_start>=$intron_start)&&($read_start<=($intron_end-7))&&($read_end>=($intron_end+8))) {
                                                                $junction_count{$junction}++;
                                                            }
                                                        }elsif($junction_type eq "i"){
                                                            if (($read_start<=$intron_start)&&($read_end>=$intron_end)) {
                                                                $junction_count{$junction}++;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if (exists $IR_postion{$tag}{short}) {
                                        my $junctions=$IR_postion{$tag}{short};
                                        my @junctions=split /,/,$junctions;
                                        foreach my $junction(@junctions){
                                            my ($junction_type,$junction_start,$junction_end)=(split /_/,$junction)[0,1,2];
                                            $junction_type=(split /-/,$junction_type)[0];
                                            my ($intron_start,$intron_end)=(split /_/,$junction)[4,5];
                                            if ($junction_type eq "e1i") {
                                                if (($read_end<=$intron_end)&&($read_end>=($intron_start+7))&&($read_start<=($intron_start-8))) {
                                                    $junction_count{$junction}++;
                                                }
                                            }elsif($junction_type eq "ie2"){
                                                if (($read_start>=$intron_start)&&($read_start<=($intron_end-7))&&($read_end>=($intron_end+8))) {
                                                    $junction_count{$junction}++;
                                                }
                                            }elsif($junction_type eq "i"){
                                                if (($intron_end-$intron_start+1)>=$Avr_read_lengh) {
                                                    if (($read_start>=$intron_start)&&($read_end<=$intron_end)) {
                                                        $junction_count{$junction}++;
                                                    }
                                                }else{
                                                    if (($read_start<=$intron_start)&&($read_end>=$intron_end)) {
                                                        $junction_count{$junction}++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    foreach my $junction(keys %long_junction){
        my ($type1,$type2)=(0,0);
        if (exists $junction_count{$junction}) {
            $type1=$junction_count{$junction};
        }
        my $read_positions="";
        if (exists $types_count{$junction}) {
            my @keys =keys %{$types_count{$junction}};
            $read_positions= join ",",@keys;
            my $key_num=$#keys+1;
            foreach my $read_start(@keys){
                $type2+=$types_count{$junction}{$read_start};#type2 raw count
            }
        }
        $type_value1{$junction}+=$type1;
        $type_value2{$junction}+=$type2;
        my @read_positions=split /,/,$read_positions;
        foreach my $read_p(@read_positions){
            $read_pos{$junction}{$read_p}++;
        }
    }
    
my $out_result_count=$SRR_acc."_read_count\.txt";
open OUT,">$read_count_outfile/$out_result_count";
print OUT "#Info_information\tNo. of Incomplete reads\tNo. of complete reads\tNo. of normalized complete reads\tNo. of all reads\tNo. of junction-length-normalized reads\n";
my %junctions;
my %junctions_normaliztion;
my %junctions2;
my %junctions2_normaliztion;
foreach my $info(keys %long_junction){
    my $junc_length=$long_junction{$info}; 
    if (exists $type_value1{$info}){
        my $junction_temp=(split /-/,$info)[1];
        my $value1=$type_value1{$info};
        my $value2=$type_value2{$info};
        my $effect=0;
        my $max=$long_junction{$info}-$Avr_read_lengh+1;
        my $value2_normal=0;
        if (exists $read_pos{$info}) {
            my @keys= keys %{$read_pos{$info}};
            $effect=$#keys+1;
            $value2_normal=$value2*$max/$effect;
        }
        my $count_nor=$value1+$value2_normal;
        my $count=$value1+$value2;
        my $count2=$count/$junc_length;
        my $count2_nor=$count_nor/$junc_length;
        print OUT "$info\t$value1\t$value2\t$value2_normal\t$count\t$count2\n";
        $junctions{$junction_temp}=$count;
        $junctions_normaliztion{$junction_temp}=$count_nor;
        $junctions2{$junction_temp}=$count2;
        $junctions2_normaliztion{$junction_temp}=$count2_nor;
    }else{
        print OUT "$info\t0\t0\t0\t0\t0\n";
    }
}
close OUT;
my %read_nums;
my %pirs;
my %warnings;
open(IN,"$junction_file1") or die "$junction_file1";
my $out_result_pir_temp=$SRR_acc."_psi_temp\.txt";
my $out_result_pir_flter=$SRR_acc."_psi_filter\.txt";
while (<IN>) {
    chomp(my $line=$_);
    unless ($line=~m/^#/){
        my ($tr,$intron_position,$junction1,$junction2,$junction3,$junction4)=(split /\t/,$line)[0,1,2,3,4,5];
        #maize
        $line=~s/B73V4_ctg/B73V4ctg/g;
        my ($exon_start,$exon_end,$ch,$intron_start,$intron_end)=split /_/,$intron_position;
        my $string1=$exon_start."_".$intron_start."_".$intron_end."_".$exon_end;
        $junction1=$junction1."_".$string1;
        $junction2=$junction2."_".$string1;
        $junction3=$junction3."_".$string1;
        $junction4=$junction4."_".$string1;
        my $junction4_temp="e1e2-".$junction4;
        my $length_e1e2=$long_junction{$junction4_temp};
        my $e1i_junction=$junctions2{$junction1};
        my $intron=$junctions2{$junction2};
        my $ie2_junction=$junctions2{$junction3};
        my $e1e2_junction=$junctions2{$junction4};  
        my $e1i_junction2=$junctions{$junction1};
        my $intron2=$junctions{$junction2};
        my $ie2_junction2=$junctions{$junction3};
        my $e1e2_junction2=$junctions{$junction4};
        my $intron_length=$intron_end-$intron_start+1;
        my ($e1i_s,$e1i_e)=(split /_/,$junction1)[1,2];
        my $length_e1i=$e1i_e-$e1i_s+1;
        my ($ie2_s,$ie2_e)=(split /_/,$junction3)[1,2];
        my $length_ie2=$ie2_e-$ie2_s+1;
        my $sum=median($e1i_junction2,$intron2,$ie2_junction2)+$e1e2_junction2;
        my $m=median($e1i_junction2,$ie2_junction2,$intron2); 
        my $test=($e1i_junction+$ie2_junction)/2+$e1e2_junction;
        my $warning="-";
        my $gene=(split /\./,$tr)[0];
        my $gene_intron=$gene."_".$intron_start."_".$intron_end;
        $read_nums{$gene_intron}=$sum;
        if ($test>0){
            my $pir=($e1i_junction+$ie2_junction)/2/(($e1i_junction+$ie2_junction)/2+$e1e2_junction);
            if (exists $pirs{$gene_intron}) {
                my $temp=$pirs{$gene_intron};
                my $temp2=($temp+$pir)/2;
                $pirs{$gene_intron}=$temp2
            }else{
                $pirs{$gene_intron}=$pir  
            }  
        }
        if ($sum>=$ARGV[4]) {
            unless ($m>=1){
                $warning="Imbalance";
            }
        }else{
            if ($m>=1) {
                $warning="LowReads";
            }else{
                $warning="LowReads;Imbalance";
            }
        }
        $warnings{$gene_intron}=$warning;
    }
}
close IN;
open OUT,">$read_count_outfile/$out_result_pir_temp";
open OUT2,">$read_count_outfile/$out_result_pir_flter";
print OUT "Gene\tIntron\tPSI\tWarning\tReadNum\n";
print OUT2 "Gene\tIntron\tPSI\n";
foreach my $gene_intron(sort keys %pirs){
    my $pir=$pirs{$gene_intron};
    my ($gene,$start,$end)=split /_/,$gene_intron;
    my $intron=$start."_".$end;
    my $warning=$warnings{$gene_intron};
    my $read_num=$read_nums{$gene_intron};
    if (($warning eq "-")&&($pir>=$ARGV[5])) {
        print OUT2 "$gene\t$intron\t$pir\n"; 
    }
    print OUT "$gene\t$intron\t$pir\t$warning\t$read_num\n";
}
close OUT;
close OUT2;
sub min{
    my @numbers=@_;
    my $min=$numbers[0];
    foreach my $number (@numbers){
        if ($number < $min) {
            $min=$number
        }
    }
    return $min;
}
sub max{
    my @numbers=@_;
    my $max=$numbers[0];
    foreach my $number (@numbers){
        if ($number > $max) {
            $max=$number
        }
    }
    return $max
}
sub median{
    my @numbers=@_;
    my @number_new=sort {$a<=>$b} @numbers;
    my $median=$number_new[1];
    return $median;
}


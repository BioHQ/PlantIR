#!/usr/bin/perl -w
use strict;
my %number;
my %exon_site;
my %exon_strand;
my %exon_ch;
my $transcript;
my %transcript_gene;
my %known_exons;
my %known_coding_gene;
my %tr_exon;
my $file_type= (split /\./,$ARGV[0])[-1];
if ($file_type eq "gff3"){
    open(IN,"$ARGV[0]") or die "can't open $ARGV[0]";
    my $number=0;
    while (<IN>) {
        chomp(my $line=$_);
        unless ($line=~m/^#/){
            my ($ch,$tag,$start,$end,$strand,$descript)=(split /\t/,$line)[0,2,3,4,6,8];
            $ch=~s/chr|ch//gi;
            if($tag=~m/exon/gi){
                $number++;
                my $tr=(split /\t/,$line)[8];
                $tr=(split /;/,$tr)[0];
                $tr=(split /:/,$tr)[1];
                my @tr=split /,/, $tr;
                my $exon=$tr.",,".$ch.",,".$start.",,".$end;
                my $tag1 =int($start/1000);
                my $tag2 =int($end/1000);
                foreach (my $tag=$tag1;$tag<=$tag2;$tag++){
                    $known_exons{$tag}{$exon}++;
                }
                $tr_exon{$tr}{$exon}++;
                $exon_strand{$tr}=$strand;
                $exon_ch{$tr}=$ch;
            }elsif($tag=~m/gene/g){
                my $gene=(split /\t/,$line)[8];
                if ($gene=~m/nontranslating_CDS|protein_coding/gi) {
                    $gene=(split /;/,$gene)[0];
                    $gene=(split /:/,$gene)[1];
                    $known_coding_gene{$gene}++;
                }
            }elsif(($descript=~m/^ID/)&&($descript=~m/Parent/)) {
                my ($transcript,$gene);
                my @descript=split /;/,$descript;
                foreach my $des(@descript){
                    if ($des=~m/^ID/gi) {
                        $transcript=(split /=|:/,$des)[-1];
                    }elsif($des=~m/Parent/gi){
                        $gene=(split /=|:/,$des)[-1];
                    }
                }
                if ($tag=~m/gene|RNA|transcript/gi) {
                    $transcript_gene{$transcript}=$gene;
                }
            }
        }
    }
    close IN;
}elsif ($file_type eq "gtf"){
    open(IN,"$ARGV[0]") or die "can't open $ARGV[0]";
    while (<IN>) {
        chomp(my $line=$_);
        unless ($line=~m/^#/){
            my ($ch,$tag,$start,$end,$strand)=(split /\t/,$line)[0,2,3,4,6];
            $ch=~s/chr|ch//gi;
            if ($tag=~m/exon/gi) {
                my $descrip=(split /\t/,$line)[8];
                $descrip=~s/"//gi;
                my @descript=split /; /,$descrip;
                my ($transcript,$gene);
                foreach my $des(@descript){
                    if ($des=~m/^transcript_id/gi) {
                        $transcript=(split / /,$des)[1];
                    }elsif($des=~m/gene_id/gi){
                        $gene=(split / /,$des)[1];
                    }
                }
                $transcript_gene{$transcript}=$gene;
                my $exon=$transcript.",,".$ch.",,".$start.",,".$end;
                my $tag1 =int($start/1000);
                my $tag2 =int($end/1000);
                foreach (my $tag=$tag1;$tag<=$tag2;$tag++){
                    $known_exons{$tag}{$exon}++;
                }
                $tr_exon{$transcript}{$exon}++;
                $exon_strand{$transcript}=$strand;
                $exon_ch{$transcript}=$ch;
            }elsif($tag=~m/gene/g){
                my $gene=(split /\t/,$line)[8];
                $gene=~s/"//gi;
                if ($gene=~m/nontranslating_CDS|protein_coding/gi) {
                    $gene=(split /; /,$gene)[0];
                    $gene=(split / /,$gene)[1];
                    $known_coding_gene{$gene}++;
                }
            }
        }
    }
    close IN;
}else{
    
    print "The format of gene annotation file should be \*\.gff3 of \*\.gtf!\n"
}

if (($file_type eq "gtf")||($file_type eq "gff3") ) {
    my %introns;
    my %introns2;
    foreach my $transcript(keys %tr_exon){
        my %qh_sites;
        for my $qh_exon(keys %{$tr_exon{$transcript}}){
            my ($start,$end)=(split /,,/,$qh_exon)[2,3];
            my $new_exon=$start."_".$end;
            $qh_sites{$new_exon}++;
        }
        my $qh_index=-1;
        my @number;
        foreach my $qh_tmp( keys %qh_sites){
            my ($q_start,$q_end)=split /_/,$qh_tmp;
            $qh_index++;
            $number[$qh_index]=$q_start;
            $qh_index++;
            $number[$qh_index]=$q_end;
        }
        @number= sort { $a <=> $b } @number;
        my $exon_site=$#number+1;
        my $strand=$exon_strand{$transcript};
        if ($exon_site>=4) {
            for(my $i=1;$i<$exon_site-2;$i+=2){
                my $start=$number[$i]+1;
                my $end=$number[$i+1]-1;
                my $ch=$exon_ch{$transcript};
                $ch=~s/chr|ch//gi;
                my $gene=$transcript_gene{$transcript};          
                if (exists $introns2{$gene}{$start}) {
                    my $temp=$introns2{$gene}{$start};
                    $introns2{$gene}{$start}=$temp.",".$end;
                }else{
                    $introns2{$gene}{$start}=$end;
                }
                if (exists $introns{$gene}{$strand}) {
                    my $temp=$introns{$gene}{$strand};
                    $introns{$gene}{$strand}="$temp\t".$start.",,".$end.",,".$transcript;
                }else{
                    $introns{$gene}{$strand}=$start.",,".$end.",,".$transcript;
                }   
            }
        }
    }
    my %intron_position;
    for my $gene(keys %introns2){
        my $intron_sort=0;
        my @strand= keys %{$introns{$gene}};
        my $strand=$strand[0];
        if ($strand eq "+") {
            for my $start (sort {$a <=>$b} keys %{$introns2{$gene}}){
                my $ends=$introns2{$gene}{$start};
                my @end=split /,/,$ends;
                my %end_temps;
                foreach my $end_tmp(@end){
                    $end_temps{$end_tmp}++  
                }
                my @end2=keys %end_temps;
                foreach my $end(sort {$a <=>$b} @end2){
                    my $intron=$start."_".$end;
                    $intron_position{$gene}{$intron}=$intron_sort;
                    $intron_sort++;
                }
            }  
        }else{
            for my $start (sort {$b <=>$a} keys %{$introns2{$gene}}){
                my $ends=$introns2{$gene}{$start};
                my @end=split /,/,$ends;
                my %end_temps;
                foreach my $end_tmp(@end){
                    $end_temps{$end_tmp}++  
                }
                my @end2=keys %end_temps;
                foreach my $end(sort {$b <=>$a} @end2){
                    my $intron=$start."_".$end;
                    $intron_position{$gene}{$intron}=$intron_sort;
                    $intron_sort++;
                }
            } 
        }
    }
    open OUT,">$ARGV[1]";
    my %intron_pairs;
    foreach my $tr(keys %tr_exon){
        $transcript=$tr;
        my $gene=$transcript_gene{$transcript};
        my $strand=$exon_strand{$tr};
        foreach my $exon(keys %{$tr_exon{$tr}}){
            my ($start,$end)=(split /,,/,$exon)[2,3];
            if (exists $introns{$gene}{$strand}) {
                my $intron=$introns{$gene}{$strand};
                my @intron=split /\t/,$intron;
                my $code;
                foreach $intron (@intron){
                    my ($db_start,$db_end,$transcript2)=split /,,/,$intron;
                    $code=&IR($start,$end,$db_start,$db_end);
                    if ($code == 2) {
                        my $intron=$transcript2.",,".$db_start.",,".$db_end;
                        if (exists $intron_pairs{$intron}) {
                            my $temp=$intron_pairs{$intron};
                            $intron_pairs{$intron}="$temp\t$transcript";
                        }else{
                            $intron_pairs{$intron}=$transcript;
                        }
                    }
                }
            }
        }
    }
    foreach my $tr(keys %tr_exon){
        $transcript=$tr;
        my $gene=$transcript_gene{$transcript};
        my $strand=$exon_strand{$tr};
        my $ch=$exon_ch{$tr};
        foreach my $exon(keys %{$tr_exon{$tr}}){
            my ($start,$end)=(split /,,/,$exon)[2,3];
            if (exists $introns{$gene}{$strand}) {
                my $intron=$introns{$gene}{$strand};
                my @intron=split /\t/,$intron;
                my $code;
                foreach $intron (@intron){
                    my ($db_start,$db_end,$transcript2)=split /,,/,$intron;
                    $code=&IR($start,$end,$db_start,$db_end);
                    if ($code == 2) {
                        my $intron=$transcript2.",,".$db_start.",,".$db_end;
                        my $tr=$intron_pairs{$intron};
                        my $tag1=int($db_start/1000);
                        my $tag2=int($db_end/1000);
                        my $overlap="-";
                        my %exons_temp;
                        for (my $i=$tag1;$i<=$tag2;$i++){
                            if (exists $known_exons{$i}) {
                                foreach my $exon(keys %{$known_exons{$i}} ){
                                    if (exists $exons_temp{$exon}) {
                                        next
                                    }
                                    $exons_temp{$exon}++;
                                    my ($tr_exon,$ch_exon,$start,$end);
                                    if ($exon=~m/Zm/gi) {
                                        ($tr_exon,$ch_exon,$start,$end)=split /,,/,$exon;
                                    }else{
                                        ($tr_exon,$ch_exon,$start,$end)=split /,,/,$exon;
                                    }
                                    if ($ch_exon eq $ch) {
                                        unless ($tr=~m/$tr_exon/){
                                            if ((($start<$db_start)&&($end>$db_start))||(($start>$db_start)&&($start<$db_end))) {
                                                $overlap="Overlap with known exon";
                                                last;
                                            }
                                        }
                                    }
                                }
                            }
                            if ($overlap=~m/Overlap with known exon/gi) {
                                last;
                            }
                        }
                        my $gene=$transcript_gene{$transcript};
                        $intron=$db_start."_".$db_end;
                        my @introns=keys %{$intron_position{$gene}};
                        my $intron_num=$#introns+1;
                        my $type="NA";
                        if ($intron_num>2) {
                            my $position= $intron_position{$gene}{$intron};
                            if ($position>0) {
                                $type="NotFirst";
                            }
                            if ($position<$#introns){
                                $type="NotLast";
                            }
                        }
                        if (exists $known_coding_gene{$gene}) {
                            print OUT "$transcript\t$ch\t$strand\t$start\t$end\t$transcript2\t$db_start\t$db_end\t$overlap\t$type\n";
                        }
                    }
                }
            }
        }
    }
    close OUT;
}
sub IR{
    my ($query_start,$query_end,$db_start,$db_end)=@_;
    my $code=0;
        if (($query_start < $db_start) && ($query_end > $db_end)) {
            $code=2
        } else{
            $code=0
        }
    return $code;
}

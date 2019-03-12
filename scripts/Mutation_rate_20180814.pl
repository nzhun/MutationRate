###
 ## Author:    Na Zhu
 ## Created:   03.08.2019
 ## 
 ## (c) Copyright by ShenLab.
###
 
 
#!/user/bin/perl

use strict;
use warnings;


print "Input perl Mutation_rate_20180814.pl transcript_file outname flag_filter_with_mappability\n";
print "#filter_with_mappability: 0 or 1, 0 means no mappablity filter, otherwise filter the position with mappablity\n";

our %rate=();
our %gene_rate=();
our $cutoff=0.025;
our @nucleotides=("A","T","C","G");

our %trichanges=();
our %amino_tables=();

my $path="~/MutationRate";
my $file_rate="$path/data/3mer_table.txt";
## read amino acid change table
my $file_aminoacid="$path/data/amino_acid_codon.txt";
my $file_map="$path/test/test.hg19.152bp.mappablity.bed.gz";
my $fmeta="$path/test/test.hg19_dbnsfp35a.bed.gz";  ##exon_transcript.dbnsfp30a.gene.bed.gz";


my %complement=(
"A" =>"T",
"T" =>"A",
"C" => "G",
"G" => "C"	
);






sub main {
		#my $dbnsfp="/home/local/ARCS/nz2274/PSAP/Data/exon_transcript.dbnsfp30a.bed.gz";
		## read 3mer
		#my $path="/Users/nazhu/server";
		my $gene="";
		my $range="";
		my $ftranscript=$_[0];
		my $fout=$_[1];
		my $QMAP=$_[2];
		print "input: $ftranscript\t  output to $fout\n";
		my %gene_rate=();
		open OUT, ">$fout";
		print OUT "##MAP_filter=".$QMAP."\n";	
			print OUT "#Reference\tGeneVersion\tTranscript\tTranscriptVersion\tExonCount\tChr\tStartPos\tEndPos\t";		
		    print OUT  "p_synonymous\t".
			"p_misense\t".
			"p_nonsense\t".
			"p_splicing\t".
			"p_frameshift\t".
			"p_LGD\t".
			"p_revel\t".
			"p_revel5_or_cadd25\t".
			"p_mcap\t".
			"p_metaSVM\t".
			"p_metalr\t".
			"p_pdmis\t".
			"p_stoploss\t".
			"p_phvar_dmis\t".
			"p_cadd10\t".
			"p_cadd15\t".
			"p_cadd20\t".
			"p_cadd25\t".
			"p_cadd30\t".
			"p_cadd35\t".
			"p_pp2Hvar\t".
			"p_fahmm\t".
			"peigen_pred10\t".
			"peigen_pred15\t".
			"peigen_pc10\t".
			"peigen_pc15".
			#join("\t", @p_cnn_nms).
			"\n";
		
		
	
	open IN, "zless $ftranscript|" or die $ftranscript." cannot find!\n";
	my @lastgenes=();
		my @lastpos=();
		#my @lastintervals=();
		my $line=<IN>;
	    my $strseq="";
		my @exon_starts=();
        my @exon_ends=();
		my $strand="+";
		do{
			if(!$line){last;}
			#print $line."\t";
		#	exit;
			$line=~ s/range=chr//g;
			$line=~ s/strand=//g;
			my @sets=split(/\s+/,$line);
			if(@sets<6){die "uncomplete info in $line ";}
			my @pos_info=split(/:|-/,$sets[1]);
			my @gene_info=split(/_|\./,$sets[0]);
			
			if(@gene_info==4){$gene_info[3]=$gene_info[3]."\t1";}
	
			if($strseq ne "" && $gene_info[2] ne $lastgenes[2]){
				my $rs=callinfo($strseq,$lastgenes[2],\@lastpos, \@exon_starts,\@exon_ends,$strand,$QMAP);
				if($rs){
					print OUT join("\t",@lastgenes)."\t".join("\t",@lastpos)."\t".$rs."\n"; #.join("\t",@rdmis)."\n";
				}
				
				$strseq="";
				@lastgenes=();
				@lastpos=();
			    @exon_starts=();
				@exon_ends=();
				push(@exon_starts,$pos_info[1]);
				push(@exon_ends,$pos_info[2]);
			} else{
				push(@exon_starts,$pos_info[1]);
				push(@exon_ends,$pos_info[2]);
				if($strseq ne ""){
					if($gene_info[1] =~ /\d+/ && $gene_info[2] =~ /\d+/ && $gene_info[4] =~ /\d/){
						$pos_info[1]=$pos_info[1] > $lastpos[1] ? $lastpos[1]:$pos_info[1];
						$pos_info[2]=$pos_info[2] < $lastpos[2] ? $lastpos[2]:$pos_info[2];
						$gene_info[4]=$gene_info[4] < $lastgenes[4] ? $lastgenes[4]:$gene_info[4];
			
				    }
		    }
		 }
		
			$strand=$sets[4];   	
			## read next line until start with :\>:
			$line=<IN>;
			my $exonseq="";
			while($line&& !($line =~ /^>/)){
				chomp($line);
				$exonseq=$exonseq."".$line;
				$line=<IN>;
			}
			$strseq.=$exonseq;
	
			@lastgenes=@gene_info;
			@lastpos=@pos_info;
			if($lastpos[1]>$lastpos[2]){ 
				my $temp=shift(@lastpos); 
				push(@lastpos,$temp); 
			}
			#push(@exon_starts,$pos_info[1]);
			#push(@exon_ends,$pos_info[2]);
		#	exit;
		} while($line);
		my $rs=callinfo($strseq,$lastgenes[2],\@lastpos, \@exon_starts,\@exon_ends,$strand,$QMAP);
		if($rs){
			print OUT join("\t",@lastgenes)."\t".join("\t",@lastpos)."\t".$rs."\n";
		}else{
			print join("\t",@lastgenes)."\t".join("\t",@lastpos)."\n";
		} 
		close IN;
	    close OUT;
		
	}
	
	
	
	sub load_dmis {
		
		my ($region)=@_;
		my %mcap=();
		my %mrevel=();
		my %meta=();
		my %pdmis=();
		my %phdmis=();
		my %pcadd=();
		my %pp2hvar=();
		my %lr=();
		my %fathmm=();
		my %eigen_phred=();
		my %eigen_pc_phred=();
		my %cnn_map=();
		
		### change to position=> mcap
		if (! -e $fmeta) {print "$fmeta cannot find\t";exit;}		
		open INM, "tabix $fmeta $region|sort -k1,1n -k2,2n -k3,3n -k4,4d -k5,5d -u|";
		
	#	print "CADD 47 |  METASVM 34 | METALR 37  | PP2HVAR 14 | PP2HDIV 11 | MCAP 38 | REVEL 41| FATHMM 51| EIGN 54, 55\n";
		while(my $line=<INM>){
		    my @sets=split(/\t/,$line);
		    my $key=join(":",@sets[(1,3,4)]);
			$pcadd{$key}=$sets[46];
		#	print "CADD ".$sets[46]."\n";
			#print $key."\t".$sets[23]."\n";
			if($sets[33] ne "T"){
		#		print "Meta ".$sets[33]."\n";
				$meta{$key}=$sets[33];  ## >0
			}
			if($sets[36] ne "T"){

		#		print "MetaLR ".$sets[36]."\n";
			    $lr{$key}=$sets[36];  ## >0
			}
			if($sets[13] eq "D" && $sets[46]>=15 && $sets[33] eq "T" ){
			   $pdmis{$key}=1;
			}
			if($sets[10] eq "D" && $sets[46]>=15 && $sets[33] eq "T" ){
			   $phdmis{$key}=1;
			}
			if($sets[13] eq "D" ){
			   $pp2hvar{$key}=1;
			}
			if($sets[51] eq "D"){
				$fathmm{$key}=1;
			}
			$mcap{$key}=$sets[37];
			$mrevel{$key}=$sets[40];
			 # <-1.61
			#$eigen_phred{$key}=$sets[53];  ## >10 >15
			#$eigen_pc_phred{$key}=$sets[54]; ## >10
			
			
		}
		close INM;
		return (\%mrevel,\%mcap,\%meta,\%lr,\%pdmis,\%phdmis,\%pcadd,\%pp2hvar,\%fathmm,\%eigen_phred,\%eigen_pc_phred,\%cnn_map);
	}
#		if (! -e $feigen) {print "$feigen cannot find\t";exit;}	
		### start eigen ###
#		open INE, "tabix $feigen $region|sort -k1,1n -k2,2n -k3,3n -k4,4d -k5,5d -u|";
#		while(my $line=<INE>){
#			my @sets=split(/\s+/,$line);
#			if($sets[3] eq ""){print $line."\n"};
		
		    
#		}
#		close INE;
		#open INM, "tabix $fmeta $region|awk '{OFS=\"\t\";}{if(\$29 !=\"T\"){print}}'|sort -k1,1n -k2,2n -k3,3n -k4,4d -k5,5d -u|";
		
		
		#  if (! -e $fmcap) {print "$fmcap cannot find\t";exit;}		
	#	open INC, "tabix $fmcap $region|sort -k1,1n -k2,2n -k3,3n -k4,4d -k5,5d -u|";
	#	while(my $line=<INM>){
#			my @sets=split(/\s+/,$line);
#			if($sets[3] eq ""){print $line."\n"};
#			$mcap{join(":",@sets[(1,3,4)])}=$sets[5]; ### >0.025
		    
			#if($sets[8] eq "1"){
		#		$mcap_correct{join(":",@sets[(1,3,4)])}=$6;
		#	}
#		}
#		close INC;
### start cnn ###		
	#	if (! -e $fcnn) {print "$fcnn cannot find\t";exit;}	
#		open INE, "tabix $fcnn $region|sort -k1,1n -k2,2n -k3,3n -k4,4d -k5,5d -u|";
#		while(my $line=<INE>){
#			my @sets=split(/\s+/,$line);
#			if($sets[3] eq ""){print $line."\n"};
#			$cnn_map{join(":",@sets[(1,2,3)])}=$sets[4];  ## >0.05..
		    
#		}
#		close INE;
			
		#### ##############

	
	
	
	sub load_info {
				
		### transcript 
		
		## load amino acide change
		if (! -e $file_rate) {print "$file_rate cannot find\t";exit;}	
		open IN, $file_rate or die $file_rate." cannot find!\n";
		while(my $line=<IN>){
			if(!($line =~ /^[ACGT]/)){next;}
			chomp($line);
			my @sets=split(/\s+/,$line);
			if(@sets<3){die "incomplete information in $line ";}
			$trichanges{$sets[0].">".$sets[1]}=$sets[2];
		}
		close IN;
		
		open IN, $file_aminoacid or die $file_rate." cannot find!\n";
		while(my $line=<IN>){
			if(!($line =~ /^[ACGT]{3}/)) {next;}
			chomp($line);
			my @sets=split(/\s+/,$line);
			if(@sets<4){die "uncomplete line : $line "}
			$amino_tables{$sets[0]}=$sets[3];
		}
		close IN;
	}	
		
	sub splicing_core {
		my @bases=@{$_[0]};
		my $flag=$_[1];
		my $rate=0;
		if((uc join("",@bases[(1,2)])) eq $flag ){
			for(my $i=1;$i<3;$i++){
				foreach my $nc(@nucleotides){
						$bases[$i]=uc $bases[$i];
		 		       if($nc eq $bases[$i]){
		 		       	   next;
		 		       }else{
		 				   my $crate=0;
						   my $key=(uc join("",@bases[(($i-1)..($i+1))])).">".(uc $bases[$i-1]).$nc.(uc $bases[$i+1]);
						#   print "try ".$key."\n";
		 				   if(exists($trichanges{$key})){
		 					    $rate+=$trichanges{$key};
							#	print $key."\t".$rate."\n";
		 				   }else{
		 					   print "Error:". $key." does not exist!\n";
		 				       return $rate;
		 				   }
					   }
				}
			}
		}
		return $rate;
	}	
	
	
	sub splicing_check {
		my @intron=@{$_[0]};
		my $rate=0;
	#	print join(" ", @intron)."\n";
		my @query=@intron[(0..3)];
		$rate+=splicing_core(\@query,"GT");
		@query=@intron[(4..7)];
		$rate+=splicing_core(\@query,"AG");
		
		return $rate;	
	}
		
	sub callinfo {
		### $strseq,$lastgenes[2],\@lastpos, \@exon_starts,\@exon_ends
		##############################
		my @sequence=split("",$_[0]);
		my $gene=$_[1];
		my @posinfo=@{$_[2]};
		my @exon_starts=@{$_[3]};
		my @exon_ends=@{$_[4]};
		my $strand=$_[5];
		my $QMAP=$_[6];
		my $region=$posinfo[0].":".$posinfo[1]."-".$posinfo[2];
	    if($sequence[0] eq "n"){return "";}
		my $pre="";
	    my $cur="";
		my $post="";
		my $p_non=0;
		my $p_mis=0;
		my $p_syn=0;
		my $p_revel=0;
		my $p_revel5_or_cadd25=0;
		my $p_mcap=0;
		my $p_mcap_005=0;
		my $p_meta=0;
		my $p_lr=0;
		my $p_pdmis=0;
		my $p_phdmis=0;
		my $p_pcadd10=0;
		my $p_pcadd15=0;
		my $p_pcadd20=0;
		my $p_pcadd25=0;
		my $p_pcadd30=0;
		my $p_pcadd35=0;
		my $p_phvar=0;
		my $p_fahmm=0;
		my $peigen_pred10=0;
		my $peigen_pred15=0;
		my $peigen_pc10=0;
		my $peigen_pc15=0;
		my @p_cnn=();
		for(my $j=0;$j<20;$j++){
		   push(@p_cnn,0);
		}
		my $p_stoploss=0;
		my $exon_cnt=0;
		my $p_splice=0;
		my $protein="";
		my $p_frame=0;
		my $coordinate=$posinfo[1];
		if($strand eq "-"){
			$coordinate=$posinfo[2];
		}
		my @intron=();
		
		my($addr_revel,$addr_mcap,$addr_meta,$addr_lr,$addr_pdmis,$addr_phdmis,$addr_cadd,$addr_pp2hvar,
		$addr_fathmm,$addr_eigen_phred,$addr_eigen_pc_phred,$addr_cnn)=load_dmis($region);
		my %revel=%{$addr_revel};
		my %mcap=%{$addr_mcap};
		my %meta=%{$addr_meta};
		my %lr=%{$addr_lr};
		my %pdmis=%{$addr_pdmis};
		my %phdmis=%{$addr_phdmis};
		my %pcadd=%{$addr_cadd};
		my %phvar=%{$addr_pp2hvar};
		my %pfathmm=%{$addr_fathmm};
		my %peigen=%{$addr_eigen_phred};
		my %peigen_pc=%{$addr_eigen_pc_phred};
		my %pcnn=%{$addr_cnn};
		
		my $pointer_exon=0;
		
		### check whether the mappability is 1 in this region
		## if return several lines, check whether all of them are 1,
		## if yes.....
		## otherwise: check position wheher in the low mappability region ....
		my @starts=();
		my @ends=();
		if(-e $file_map){
			if($QMAP!=0){
		
		    	open MAP, "tabix $file_map $region|awk '{if(\$4!=1){print}}'| ";
		   
				while(my $map_line=<MAP>){
					my @sets=split(/\s+/,$map_line);
			
					push(@starts,$sets[1]);
					push(@ends,$sets[2]);
				#	print $map_line."\n";
					#print 
					
			    }
			    close MAP;
		    }
		}
		##
		my $BIM=0;
		for(my $i=1;$i<@sequence-3;$i++){
			if($strand eq "-"){
				$coordinate-=1;
				if($coordinate<$exon_starts[$pointer_exon]){
					if($pointer_exon < @exon_starts){
						$pointer_exon+=1;
						$coordinate=$exon_ends[$pointer_exon];
					}
				}
			}else{
			   $coordinate+=1;
			   if($coordinate>$exon_ends[$pointer_exon]){
					if($pointer_exon < @exon_starts){
				   	$pointer_exon+=1;
				   	$coordinate=$exon_starts[$pointer_exon];
					}
				#   print "Skip  ".$gene."\t".join("\t",@posinfo)."\t".$pointer_exon."\t".join(":",@exon_starts)."\t".$exon_starts[$pointer_exon]."\t".$coordinate."\t"."\n";
			   }	
			}
			my $flag_map=0;
		    for(my $im=$BIM;$im<@starts;$im++){
				#print "lost: ". @starts."_".$im." BIM ".($BIM)."\t".$coordinate."\t".$starts[$im]."\t".$ends[$im]."\t".qx(tabix $file_map $posinfo[0]:$coordinate-$coordinate)."\t\n";
			
				if($coordinate < $starts[$im] && $strand eq "+"){$flag_map=0;	$BIM=$im;last;}
		    	if($coordinate>$starts[$im]-1 && $coordinate<$ends[$im]){$flag_map=1;last;}
		    }
			if($flag_map){ next;} #print "missed:".$coordinate."\t".qx(tabix $file_map $posinfo[0]:$coordinate-$coordinate)."\t\n";
	#		print "take: $region\t".$coordinate."\t".qx(tabix $file_map $posinfo[0]:$coordinate-$coordinate)."\t\n";
			
			### only check the uppercase, as those bases are in exon.
			if($sequence[$i] eq (uc $sequence[$i])  ){
				#print "$coordinate\t".$sequence[$i]."\textron\n";
				$exon_cnt=$exon_cnt+1;}
			else{
	#			print "$coordinate\t".$sequence[$i]."\tintron\n";
			### get the bases in the introns, it consists of donor and acceptor of splicing sites.
				if($i>3){
					if(@intron ==0){
						push(@intron,$sequence[$i-1]);
				    }
					push(@intron,$sequence[$i]);
					if(@intron==7){
						push(@intron,$sequence[$i+1]);
						$p_splice+=splicing_check(\@intron);
		#				print "$gene\t $coordinate\t".join(" ",@intron)."\t".$p_splice."\n";
						@intron=();
					#	exit;
					}
				}
				next;
			}
			
			if($exon_cnt%3==1){
				$protein="";
				my $k=$i;
				while(length($protein)<3 && $k<@sequence){
					if($sequence[$k] eq (uc $sequence[$k])){
				    	$protein.=$sequence[$k];
					}
					$k=$k+1;
			    }
			}
			$pre=uc $sequence[$i-1];
			$cur=$sequence[$i];
			$post=uc $sequence[$i+1];
			
			
			foreach my $nc(@nucleotides){
		       if($nc eq $sequence[$i]){
		       	   next;
		       }else{
				   my $crate=0;
				   if(exists($trichanges{$pre.$cur.$post.">".$pre.$nc.$post})){
					    $crate=$trichanges{$pre.$cur.$post.">".$pre.$nc.$post};
				   }else{
					   print "Error: $gene\t".$pre.$cur.$post ." > ". $pre.$nc.$post." does not exist!\n";
				       return ;
				   }
			#	   print $coordinate."\t".$nc." < ".$cur.":\t".$pre.$cur.$post.">".$pre.$nc.$post."\t".$crate."\n";
				   my $cprotein=uc $nc."".(substr $protein,1);
				  
				   if($exon_cnt%3==0){
				   		$cprotein=(substr $protein,0,2).$nc;
				   }elsif($exon_cnt%3==2){
				   		$cprotein=(substr $protein,0,1).$nc.(substr $protein,2);
				   }
				   if(length($protein)<3){next;}
				   if(!exists($amino_tables{$protein})||!exists($amino_tables{$cprotein})){
			#		   print  "$gene\t$region\t$i\t$protein\t $cprotein\n";
					   return ;
				   }
				  # print  "Amino Acid change \t$protein\t $cprotein\n";
				 #  print $amino_tables{$protein}." > > > ".$amino_tables{$cprotein}."\n";
				   if($amino_tables{$protein} ne $amino_tables{$cprotein}){
				   	   if($amino_tables{$protein} eq "X"){
						   $p_stoploss=$p_stoploss+$crate;
			#			   print "stoploss\n";
					   }
				   	   if($amino_tables{$cprotein} eq "X"){
						   ## nonsense
						   $p_non=$p_non+$crate;
						   $p_frame=$p_frame+$crate*1.25;
		#				   print "Nonsense\n";
				   	   }else{
						   if($amino_tables{$protein} ne "X"){
						   ## missense 
						  	 $p_mis=$p_mis+$crate;
							 ### check whether the nucletides change in mcap or meta list
						     my $temp_ref=$cur;
							 my $temp_alt=$nc;
							 if($strand eq "-"){
							  	$temp_ref=$complement{$cur};
							  	$temp_alt=$complement{$nc};
						  	 }	
							 my $key=$coordinate.":".$temp_ref.":".$temp_alt;
							 my $union_flag=0;
							 if(exists($revel{$key})){
							 	 my $vrevel=$revel{$key};
								 if($vrevel ne "." && $vrevel >=0.5){
						 			 $p_revel+=$crate;
									 $p_revel5_or_cadd25+=$crate;
									 $union_flag=1;
								# print "$key uninon revel $vrevel adding\n";
							 	}
									 
							 }
							 ### start query
						  #  my $flag_revel_or_cadd=0;
							 if(exists($mcap{$key})){
								#print $key."\t".scalar(keys %mcap)."\n";     
							     my $value=$mcap{$key};
#	print $key."\t".$value."\n";
								  if($value ne ".") {
							     	if($value>0){$p_mcap+=$crate;}
							     	if($value>0.05){$p_mcap_005+=$crate;}
								  }							 
							 }
							 
						### mcap
							 if(exists($lr{$key})){
								 $p_mcap+=$crate;		
								# print "Yes, $region\t $key Mcap\n";							 
							 }
						 
							 if(exists($meta{$key})){
							 	 $p_meta+=$crate;
								#  print "Yes,$region $key Meta\n";	 	 
							 }  
							 if(exists($pdmis{$key})){
							 	 $p_pdmis+=$crate;
								#  print "Yes,$region $key pdmis\n";	 	 
							 }
							 if(exists($phdmis{$key})){
							 	 $p_phdmis+=$crate;
								#  print "Yes,$region $key pdmis\n";	 	 
							 }
							 
							 if(exists($phvar{$key})){
							 	 $p_phvar+=$crate;
								#  print "Yes,$region $key pdmis\n";	 	 
							 }
							 
							 if(exists($pfathmm{$key})){
							 #	 my $value=$pfathmm{$key};		 
							 	# if($value ne "." && $value< -1.61){
							 	 	 $p_fahmm+=$crate;
							#	 }
							 }
							
							 if(exists($peigen{$key})){
							 	 my $value=$peigen{$key};
								 if($value ne "."){
								 	 if($value>10){
								 	 	 $peigen_pred10+=$crate;
									 }
									 if($value>15){
								 	 	 $peigen_pred15+=$crate;
									 }
							 	}
							 }
							
							
							 if(exists($peigen_pc{$key})){
							 	 my $value=$peigen_pc{$key};
								 if($value ne "."){
								 	 if($value>10){
								 	 	 $peigen_pc10+=$crate;
									 }
									 if($value>15){
								 	 	 $peigen_pc15+=$crate;
									 }
								 }
								 
							 }
							 
							 if(exists($pcnn{$key})){
							     my $value=$pcnn{$key};
								  if($value ne "."){
							     my $best=0.56;
							     for(my $t=0;$t<19;$t++){
							        if($value>0.05*($t+1)){
							          $p_cnn[$t]+=$crate;
							        }
							     }
							     
							     if($value>$best){
							        $p_cnn[19]+=$crate;
							     }
							  }
							    
							 }
							  
							 
							 
							 if(exists($pcadd{$key})){
							 	 my $value=$pcadd{$key};
					#		 	 print $key."\t".$value."\t".$crate."\n";
								if($value ne "."){
									 if( $value>=10){
								 	 		$p_pcadd10 +=$crate;
								 	 }
							 	 
								 	 if($value>=15){
								 	 		$p_pcadd15 +=$crate;
								 	 }
								 	 if($value>=20){
								 	 		$p_pcadd20 +=$crate;
								 	 }
							 	 							 	 
								 	 if($value>=25){
								 	 		$p_pcadd25 +=$crate;
											if($union_flag==0){
												 $p_revel5_or_cadd25=$p_revel5_or_cadd25+$crate;
									#			 print "$key union $value adding\n";
											}
										
								 	 }
							 	 
								 	 if($value>=30){
								 	 		$p_pcadd30 +=$crate;
								 	 }
							 	 
								 	 if($value>=35){
								 	 		$p_pcadd35 +=$crate;
								 	 }
							 	 
							 	 
							 	 
								#  print "Yes,$region $key pdmis\n";	 	 
							 }
						 }
							 
							   
					#	   	print "Missense\n";
					      }
				   	   }
				   }else{
					   ## syn
					   $p_syn=$p_syn+$crate;
					#   print "Synonymous: ";
					#   print $coordinate."\t".$nc." < ".$cur.":\t".$pre.$cur.$post.">".$pre.$nc.$post."\t".$crate."\n";
				   }
		     #  print $crate."\t".$cprotein."\t".$protein."\t".$pre.$cur.$post.">".$pre.$nc.$post."\n";
			   }
			  
			}
		
		
		}
	#	print "final  $p_revel	Union  $p_pcadd25 	= $p_revel5_or_cadd25\n";
		%mcap=();
		%meta=();
		my $rs=join("\t",(
			$p_syn,
                        $p_mis,
                        $p_non,
                        $p_splice,
                        $p_frame,
                        $p_frame+$p_splice+$p_non,
								$p_revel,
								$p_revel5_or_cadd25,
                        $p_mcap,
                        $p_meta,
                        $p_lr,
                        $p_pdmis,
                        $p_stoploss,
                        $p_phdmis,
                        $p_pcadd10,
                        $p_pcadd15,
                        $p_pcadd20,
                        $p_pcadd25,
                        $p_pcadd30,
                        $p_pcadd35,
                        $p_phvar,
                        $p_fahmm #,
                     #   $peigen_pred10,
                     #   $peigen_pred15,
                    #    $peigen_pc10,
                    #    $peigen_pc15#,
                        #@p_cnn
			
			));
		return $rs;
	  
	  
	}
	
	



	#print "Input: perl sanity.samocha.pl output mappability.bed.gz\n";

	my $ftranscript=$ARGV[0]; #"$path/Resources/Transcript/GencodeV19_splice.gz";MY $QMAP=1;
	my $fout=$ARGV[1];	
	my $QMAP=0;
	if(@ARGV>2){
		$QMAP=$ARGV[2];
	}

	print "when QMAP!=0, load mappability from $file_map' Currently qmap=$QMAP\n";
	print "load 3mer mutation rate from $file_rate\n";
	print "load anminoacid codon from $file_aminoacid\n";

	print "load annotation from $fmeta\n";

	load_info();

	main($ftranscript,$fout,$QMAP);
	print "output to $fout\n"

use strict;
use Bio::TreeIO;
my $tree_clade = $ARGV[0];
open (OUT, ">>IMPORTS.".$tree_clade.".WithCountries.txt") or die $!;


my $in = Bio::TreeIO -> new(-file => $tree_clade.'/nonbinary_full.nwk',
                            -format => 'newick');

my $tree = $in -> next_tree;
my $root = $tree->get_root_node;
$root -> id("NODE_root");
my $root_id = $root -> id;
#print $root_id."\n";

my @leafnodes = $tree->get_leaf_nodes;
my @leafs;

foreach my $leafnode(@leafnodes)  {
    my $leaf = $leafnode->id;
    # my $brlen = $leafnode -> branch_length;
    # print $leaf."\t".$brlen."\n";
    push(@leafs,$leaf);
}



my %level;

my $m = 0;
my $q = 0;
push (@{$level{$m}},$root_id);
lev ($root_id,$q);
my %node_rus;
my %node_nonrus;
#my @russian_ancestors;   #internal nodes with only-russian descendents
my %node_days_date;
my %node_russeqn;
my %node_countries;

foreach my $n_id(@leafs) {
    my $date = "noinf";
    my $days = 0;
    if ($n_id =~ m/\|20\d\d\-\d+\-\d+/) {
	$n_id =~ m/\|20\d\d\-\d+\-\d+/;
	$date = $&;
	$date =~ s/\|//;
	my @date = $date =~ m/\d+/g;
	my $year = $date[0];
	my $month = $date[1];
	my $day = $date[2];
	my @mon = $month =~ m/\d/g;
	if ($#mon == 0) {
	    $month = "0".$month;
	}
	
	if($#mon>1) {
	    $date = "noinf";
	}
	else {
	    $date = $year."-".$month."-".$day;
	    if ($year == 2020) {
		$days += 365;
	    }
	    $days += ($month-1)*30;
	}
    }
    
    elsif ($n_id =~ m/\|20\d\d\-\d+/) {
	$n_id =~ m/\|20\d\d\-\d+/;
	$date = $&;
	$date =~ s/\|//;
	my @date = $date =~ m/\d+/g;
	my $year = $date[0];
	my $month = $date[1];
	my $day = 30;
	my @mon = $month =~ m/\d/g;
	if ($#mon == 0) {
	    $month = "0".$month;
	}
	if ($month eq "02") {
	    $day = 28;
	}
	
	if($#mon>1) {
	    $date = "noinf";
	}
	else {
	    $date = $year."-".$month."\-".$day;
	    if ($year == 2020) {
		$days += 365;
	    }
	    $days += ($month-1)*30;
	}
    }
    
    if ($date eq "noinf") {
	$days = 1000;
	$date = "3000-01-01";
    }
    my $days_date = $days."_".$date;
    $node_days_date{$n_id} = $days_date;
}


until ($m == -1)  {

    my @nds = @{$level{$m}};
    foreach my $n_id(@nds)  {
	#	    print "n_id\t".$n_id."\n";
	my $rus_seq_number = 0;
	if ($n_id ~~ @leafs)  {
	
	    $n_id =~ m/\|\d+$/;
	    my $pre_seqn = $&;
	    $pre_seqn =~ m/\d+/;
	    my $seqn = $&;
	    
	    $n_id =~ m/\//;
	    my $pre_country = "$'";
	    $pre_country =~ m/\//;
	    my $country = $`;

	    my $node = $tree -> find_node(-id => $n_id);
	    if ($n_id =~ m/Russia/ | $n_id =~ m/FMBA/) {
		$node_rus{$n_id} += 1;
		$node_russeqn{$n_id} += $seqn;
	    }
	    else {
		$node_nonrus{$n_id} += 1;
		
		${$node_countries{$n_id}}{$country} += 1;
	     }   
	}
	
	else  {
	    my $days = 1000;
	    my $date = "noinf";
	    my $node = $tree->find_node(-id=>$n_id);

	    my $nonrus_days = 1000;
	    my $nonrus_date = "noinf";
		

	    my @des = $node->each_Descendent;
	    my $desn = $#des+1;
	    my $rusn = 0;
	    my @des_ids;
	    my @onlyrus;
	    my @others;
	    foreach my $des(@des) {
		my $n1 = $des->id;
		if (exists $node_russeqn{$n1}) {
		    my $rn = $node_russeqn{$n1};
		    $node_russeqn{$n_id} += $rn;
		}

		if (exists $node_countries{$n1}) {
		    my %countries =  %{$node_countries{$n1}};
		    foreach my $country(keys %countries) {
			my $count = $countries{$country};
			${$node_countries{$n_id}}{$country} += $count;
		    }
		}
		
		my $days_date = $node_days_date{$n1};
		$days_date =~ m/\_/;
		my $des_days = $`;
		my $des_date = "$'";
		if ($des_days < $days) {
		    $days = $des_days;
		    $date = $des_date;
		}
		push(@des_ids,$n1);
		
		
		my $nonrusn = 0;
		my $rusn = 0;
		
		if (exists $node_nonrus{$n1}) {
		    $nonrusn = $node_nonrus{$n1};
		    $node_nonrus{$n_id} += $nonrusn;
		}
		
		if (exists $node_rus{$n1}) {
		    $rusn = $node_rus{$n1};
		    $node_rus{$n_id} += $rusn;
		}
		my $numb = $rusn + $nonrusn;

		if ($numb == $nonrusn) {
		    push(@others,$n1);
		    if ($des_days< $nonrus_days) {
			$nonrus_days = $des_days;
			$nonrus_date = $des_date;
		    }
		}
		
		elsif ($numb == $rusn) {
		    push(@onlyrus,$n1);
		    $rusn += 1;
		}

		
	    }
	    my $days_date = $days."_".$date;
	    $node_days_date{$n_id} = $days_date;
	    if ($#onlyrus>-1 && $#others >-1) {
		my $anc = $node -> ancestor;
		my @m_des = $anc->each_Descendent;
		my $desn = 0;
		my $rnumber = 0;
		my $m_days = 1000;
		my $m_date = "no_inf";
		foreach my $m_des(@m_des) {
		    my $m_des_id = $m_des -> id;
		    if ($m_des_id ~~ @leafs) {
			$desn += 1;
			if ($m_des_id =~ m/Russia/ | $m_des_id =~ m/FMBA/) {
			    $rnumber += 1;   
			}
			else {
			  #  $nonrusn += 1;
			    my $days_date = $node_days_date{$m_des_id};
			    $days_date =~ m/\_/;
			    my $des_days = $`;
			    my $des_date = "$'";
			    if ($des_days < $m_days) {
				$m_days = $des_days;
				$m_date = $des_date;
			    }
			}
			#else {
			 
			#}
		    }
			
		    
		    elsif ($m_des ne $node) {
		    #my $m_des_id = $m_des -> id;
			my @k_des = $m_des->get_all_Descendents;
			$desn += $#k_des;
			$desn += 1;
			
			foreach my $k_des(@k_des) {
			    my $n1 = $k_des->id;
			    if ($n1 =~ m/NODE/) {
			    }
			    elsif ($n1 =~ m/Russia/ | $n1 =~ m/FMBA/) {
				
			   # }
			    #else {
				$rnumber += 1;
			    }
			   
			    #if (exists $node_nonrus{$n1} && $node_nonrus{$n1}>0) {
				#$nonrusn += 1;
				my $days_date = $node_days_date{$n1};
				$days_date =~ m/\_/;
				my $des_days = $`;
				my $des_date = "$'";
				if ($des_days < $m_days) {
				    $m_days = $des_days;
				    $m_date = $des_date;
				}
			    #}
			}
		    }
		}
		if ($rnumber == 0) {
		    my $early_days=1000;
		    my $early_date = "no_inf"; 
		    foreach my $des_id(@onlyrus) {
			my $days_date = $node_days_date{$des_id};
			$days_date =~ m/\_/;
			my $des_days = $`;
			my $des_date = "$'";
			if ($des_days < $early_days) {
			    $early_days = $des_days;
			    $early_date = $des_date;
			}
		    }
		    #my @countries = keys %{$node_countries{$n_id}};
		    my $countries;
		    my %countries = %{$node_countries{$n_id}};
		    foreach my $country (keys %countries) {
			my $number = $countries{$country};
			$countries = $countries.$country."-".$number.",";
		    }
			
		    my $rus_seq_number = $node_russeqn{$n_id};
		    print OUT $n_id."\t".$rus_seq_number."\t".$nonrus_days."\t".$nonrus_date."\t".$early_days."\t".$early_date."\t".$m_days."\t".$m_date."\t".$tree_clade."\t".$countries."\n";
		}
	    }
	}
    }
    $m=$m-1;
}



sub lev  {

    my ($par_id,$n) = @_;

    $n++;

    my $par = $tree->find_node(-id => $par_id);

  CHILD: for my $child ($par -> each_Descendent)  {

      my $ch_id = $child->id;
      push (@{$level{$n}},$ch_id);
      #     print $n."\n";
      #    foreach my $val(@{$level{$n}})  {
      #  print $val."\t";
      #     }
      #    print "\n";

      if ($n > $m)  {

	  $m = $n;
      }



      if($ch_id ~~ @leafs) {

      }

      else {
	  lev ($ch_id,$n);
      }
      next CHILD;
  }

}

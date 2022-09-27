#!/usr/bin/perl -w

use strict;

my $verbose = 0;

#seq_id  agez    ss_dist bp_seq  bp_scr  y_cont  ppt_off ppt_len ppt_scr svm_scr
#chr10:104627330-104680644:-     46      79      gcctgagag       -0.183638673693 0.608108108108  12      14      22      0.010586033
#chr10:104627330-104680644:-     46      40      tcttgaaat       -1.26859042483  0.771428571429  6       10      21      0.087537836
#chr10:104627330-104680644:-     46      24      ttatcaatc       -1.8557743933   0.894736842105  2       18      41      0.2338105
#chr10:104627330-104680644:-     46      20      caatcactt       -0.0031532448193        1.0     1       15      38      1.0466396

#chr10:107651562-107667492:+     39      87      tggtaatct       1.11110115461   0.487804878049  1       16      25      1.2842952
#chr10:107651562-107667492:+     39      82      atcttattt       -2.139714333    0.467532467532  1       11      16      0.014234349
#chr10:107651562-107667492:+     39      26      aaataatgt       0.5337489328    0.761904761905  1       7       11      1.1257521

#chr10:107667652-107671054:+     37      58      ttctgacta       1.85732996461   0.641509433962  4       8       14      1.3789546
#chr10:107667652-107671054:+     37      20      atctaattc       1.52510243521   0.733333333333  15      0       0       0.44494301
#chr10:107667652-107671054:+     37      16      aattcactg       0.479817158006  0.727272727273  11      0       0       0.33676822

#chr10:107671199-107677160:+     65      57      tgttgatgt       1.28629479991   0.692307692308  3       9       16      1.2537323
#chr10:107671199-107677160:+     65      52      atgtcattt       -2.19939817752  0.702127659574  1       6       13      0.061740983
#chr10:107671199-107677160:+     65      36      tattgactt       1.15128577414   0.774193548387  1       9       20      1.3831855

my %bps;
while(<>){
    chomp;
    next if /seq_id/;
    my ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr ) = split;
    
    # keep only positive SVM scores
    next unless $svm_scr > 0;

    # store all predictions per intron:
    push( @{$bps{$seq_id}}, [$seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr] );
}

foreach my $seq_id (keys %bps){
    
    # Is there a BP inside the AGEZ with positive BP score?
    my @filtered_bps;
    my @soft_filtered_bps;
    foreach my $bp ( @{$bps{$seq_id}} ){
	my ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr) =  @{$bp};

	if ($verbose){
	    my $s = join "\t", ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr);
	    print $s."\n";
	}
	if ($ss_dist < $agez + 9 && $bp_scr > 0){
	    push( @filtered_bps, $bp );
	}
	if ($ss_dist < $agez + 9 ){
	    push( @soft_filtered_bps, $bp );
	}
    }
    if ( scalar(@filtered_bps) == 1 ){
	my ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr) = @{$filtered_bps[0]};
	print "Selected (1)\n" if $verbose;
	my $s = join "\t", ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr);
	print $s."\n";
    }
    elsif ( scalar(@filtered_bps) > 1 ){
	# select the one with the highest svm score:
	my @sorted_bps = sort {$b->[9] <=> $a->[9]} @filtered_bps;
	my ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr) = @{$sorted_bps[0]};
	print "Selected (2)\n" if $verbose;
	my $s = join "\t", ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr);
	print $s."\n";
    }
    elsif( scalar(@soft_filtered_bps) > 0 ){
	# we drop the condition of positive bp score
	# and take the best svm_score that are close to AGEZ
	my @sorted_bps = sort {$b->[9] <=> $a->[9]} @soft_filtered_bps;
        my ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr) = @{$sorted_bps[0]};
        print "Selected (3)\n" if $verbose;
        my $s = join "\t", ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr);
        print $s."\n";
    }
    else{
	# we take all BPs and choose the one with the best svm score:
	my @sorted_bps = sort {$b->[9] <=> $a->[9]} @{$bps{$seq_id}};
        my ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr) = @{$sorted_bps[0]};
        print "Selected (4)\n" if $verbose;
        my $s = join "\t", ($seq_id, $agez, $ss_dist, $bp_seq, $bp_scr, $y_cont, $ppt_off, $ppt_len, $ppt_scr, $svm_scr);
        print $s."\n";
    }
    #my $s = join "\t", ($seq_id, scalar(@filtered_bps) );
    #print $s."\n";
}

    

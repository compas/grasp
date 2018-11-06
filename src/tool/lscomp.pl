#!/usr/bin/perl
#use strict;
#use warnings;

print "   LSCOMP.PL\n";
print "   This PERL script creates files lscomp.tex and energylabel.latex\n";
print "   \n";
print "   File lscomp.tex contains energy level data with up to \n";
print "   three LS components with a contribution > 0.02 of the \n";
print "   total wave function. \n";
print "   \n";
print "   File energylabel.latex may be used by RTABTRANS2 to produce\n";
print "   LaTeX tables of transition data.\n";
print "   \n";
print "   Input files : state1.lsj.lbl and state2.lsj.lbl\n";
print "                 state1.(c)h and state2.(c)h (optional for g_J factors)\n";
print "   Output files: lscomp.tex and energylabel.latex\n";
print "                                  Jorgen Ekman Sep. 2015\n";
print "   \n";

#$state1 = "even_5_b";
#$state2 = "odd_5_b";

print "   State 1?  ";
chomp($state1 = <>);
print "   State 2?  ";
chomp($state2 = <>);
if (length($state2) == 0) {
    $state2 = $state1;
} 


$infile1 = $state1.".lsj.lbl";
$infile2 = $state2.".lsj.lbl";
$infile3 = $state1.".ch";
$infile4 = $state2.".ch";
$infile5 = $state1.".h";
$infile6 = $state2.".h";

#Look for *.lsj.lbl files. If not found - exit program.
if (-e $infile1 && -e $infile2) {
    printf("   Necessary input file(s) exist!\n");
    printf("\n");    
}else{
    printf("   Necessary input file(s) do not exist! Program terminates!\n");
    printf("\n");    
    exit;
}

printf("   Do you want to include Lande g_J factors in the energy table? (y/n) ");
chomp($gJinclude = <>);
if($gJinclude eq "y"){
    $gJinclude = 1;
    print "   Lande g_J factors from a CI calculation? (y/n) ";
    chomp($ci = <>);
    if($ci eq "y"){
	$ci = 1;
    }else{
	$ci = 0;
    }    
}else{
    $gJinclude = 0;
}

if($ci == 1 && $gJinclude ==1){
    if (-e $infile3 && -e $infile4) {
	printf("   File(s) with g_J factors exist!\n");    
    }else{
	printf("   File(s) with g_J factors do not exist!\n");    
	$gJinclude = 0;    
    }
}elsif($ci == 0 && $gJinclude == 1){
    if (-e $infile5 && -e $infile6) {
	printf("   File(s) with g_J factors exist!\n");    
    }else{
	printf("   File(s) with g_J factors do not exist!\n");    
	$gJinclude = 0;    
    }
    
}

printf("\n");    
printf("   Do you want an extra empty column for e_obs in the energy table? (y/n) ");
chomp($eobsinclude = <>);
if($eobsinclude eq "y"){
    $eobsinclude = 1;
}else{
    $eobsinclude = 0;
}
printf("\n");
printf("   Inspect the labels of the states and \n");
printf("   determine how many positions should be skipped in \n");
printf("   the string that determines the label. For example\n");
printf("   if all the states have a common core 1s(2) in the \n");
printf("   label then 6 positions should be skipped\n");
    
printf("\n");
printf("   How many positions should be skipped? ");
chomp($nochar = <>);

$ec = 219474.6313702;   #1 a.u = 219474.6313702 cm^-1

#READ CONTENT IN FILE 1
#----------------------
open(INPUTFILE1, $infile1);
$i=0;
while(<INPUTFILE1>){
    my($line) = $_;
    chomp($line);
    $linecontent[$i] = $line;
    $size[$i] = scalar(split(' ', $linecontent[$i]));
    $states1[$i] = $state1;
    #printf("$linecontent[$i]\n");
    $i++;
}
$imax = $i;
close(INPUTFILE1);

if($infile1 ne $infile2){
    
    #READ CONTENT IN FILE 2
    #----------------------
    open(INPUTFILE2, $infile2);
    $ii=0;
    while(<INPUTFILE2>){
	my($line) = $_;
	chomp($line);
	$linecontent2[$ii] = $line;
	$size2[$ii] = scalar(split(' ', $linecontent2[$ii]));
	$states2[$ii] = $state2;
	#printf("$linecontent2[$ii]\n");
	$ii++;
    }
    $iimax = $ii;
    close(INPUTFILE2);
    
    #MERGE CONTENT OF THE FILES 1 & 2
    #--------------------------------
    for($i=0; $i<$iimax; $i++){
	$linecontent[$i+$imax] = $linecontent2[$i];
	$size[$i+$imax] = $size2[$i];
	$states1[$i+$imax] = $states2[$i];
    }
    $imax = $imax + $iimax;
}

if($gJinclude == 1) {
    
    #READ CONTENT IN FILE 3
    #----------------------
    if($ci == 1){
	open(INPUTFILE3, $infile3);
    }else{
	open(INPUTFILE3, $infile5);
    }
    $i=0;
    $j=0;
    while(<INPUTFILE3>){
	if($i > 8){
	    my($line) = $_;
	    chomp($line);
	    $linecontent3[$j] = $line;
	    #printf("$j  $linecontent3[$j]\n"); 
	    $j++;
	}
	$i++;
    }
    $iiimax = $j;
    close(INPUTFILE3);

    if($infile3 ne $infile4) {
	
	#READ CONTENT IN FILE 4
	#----------------------
	if($ci == 1){
	    open(INPUTFILE4, $infile4);
	}else{
	    open(INPUTFILE4, $infile6);
	}
	$i=0;
	$j=0;
	while(<INPUTFILE4>){
	    if($i > 8){
		my($line) = $_;
		chomp($line);
		$linecontent4[$j] = $line;
		#printf("$j  $linecontent4[$j]\n"); 
		$j++;
	    }
	    $i++;
	}
	$iiiimax = $j;
	close(INPUTFILE4);
	
	#MERGE CONTENT OF FILES 3 & 4
	#------------------------------
	for($i=0; $i<$iiiimax; $i++){
	    $linecontent3[$i+$iiimax] = $linecontent4[$i];
	}
	$iiimax = $iiimax + $iiiimax;
	
    }
    
    #EXTRACT RELEVANT DATA FROM CONTENT IN FILES 3 & 4
    #--------------------------------------------------
    #printf("\n");
    for($i=0; $i<$iiimax; $i++){
	@linesplit = split(' ', $linecontent3[$i]);
	$or[$i] = $linesplit[0];
	$sp[$i] = $linesplit[1];
	$pa[$i] = $linesplit[2];
	$gj[$i] = $linesplit[5];
	my $substring = "D+00";
	my $substring2 = "D-01";
	my $substring3 = "D-02";
	if ($gj[$i] =~ /\Q$substring\E/) {
	    $gj[$i] = substr($gj[$i], 0, 7);
	    $gj[$i] = sprintf("%.5f", $gj[$i]); 
	}
	if ($gj[$i] =~ /\Q$substring2\E/) {
	    $gj[$i] = substr($gj[$i], 0, 8);
	    $gj[$i] = $gj[$i]/10.0;
	    $gj[$i] = sprintf("%.5f", $gj[$i]); 
	}
	if ($gj[$i] =~ /\Q$substring3\E/) {
	    $gj[$i] = substr($gj[$i], 0, 8);
	    $gj[$i] = $gj[$i]/100.0;
	    $gj[$i] = sprintf("%.5f", $gj[$i]); 
	}
    }
}

#CONTINUE PROCESSING DATA
#--------------------------------------------------
printf("\n");
$j=0;
for($i=0; $i<$imax; $i++){
    @linesplit = split(' ', $linecontent[$i]);
    $size[$i] = scalar(@linesplit);
    if($size[$i] == 5){
	$states[$j] = $states1[$i];
	$order[$j]  = $linesplit[0];
	$spin[$j]   = $linesplit[1];
	$parity[$j] = $linesplit[2];
	$energy[$j] = $linesplit[3];
	$energytot[$j] = $energy[$j];
	
	if($gJinclude == 1) {
	    for($n=0; $n<$iiimax; $n++){
		if($order[$j] == $or[$n] && $spin[$j] eq $sp[$n] && $parity[$j] eq $pa[$n]){
		    $gjval[$j] = $gj[$n]; 
		    #printf("$order[$j]   $or[$n]          $spin[$j]  $sp[$n]       $parity[$j]  $pa[$n]        $gjval[$j]\n");
		}
	    }
	}

	#printf("$energy[$j]\n");
	$k=$i+1;
	$l=0;
	while($size[$k] == 3){
	    @linesplit = split(' ', $linecontent[$k]);
	    if($linesplit[1] < 0.02 && $l == 0){
		$wavemag[$j][$l]  = $linesplit[1];
		$configuration[$j][$l]  = $linesplit[2];

		$conflength = length($configuration[$j][$l]);
		$configuration2[$j][$l] = substr($configuration[$j][$l], 0, $conflength - 3);

		#printf("$order[$j] $spin[$j] $parity[$j] $energy[$j] $wavemag[$j][$l]   $configuration[$j][$l]\n");
		$l++;
	    }
	    if($linesplit[1] >= 0.02 && $l < 3){
		$wavemag[$j][$l]  = $linesplit[1];
		$configuration[$j][$l]  = $linesplit[2];

		$conflength = length($configuration[$j][$l]);
		$configuration2[$j][$l] = substr($configuration[$j][$l], 0, $conflength - 3);

		#printf("$order[$j] $spin[$j] $parity[$j] $energy[$j] $wavemag[$j][$l]   $configuration[$j][$l]\n");
		$l++;
	    }
	    $k++;
	}
	$lmax[$j] = $l;
	$j++;
	#printf("\n");
    }
    $jmax = $j;
}

#SORT DATA
#---------
@energysort = sort { $a <=> $b } @energy;
#printf("@energysort\n");
#printf("$energysort[0]\n");
for($j=0; $j<$jmax; $j++){
    for($k=0; $k<$jmax; $k++){
	if($energy[$k] == $energysort[$j]){
	    $energysortr7[$j] = sprintf "%.7f", $energysort[$j];	    

	    $eenergysort[$j] = $ec*($energysort[$j] - $energysort[0]);
	    $eenergysortr[$j] = sprintf "%.0f", $eenergysort[$j];
	    $eenergysortr2[$j] = sprintf "%.2f", $eenergysort[$j];
	    $eenergysortrsep[$j] = thousandsep($eenergysortr[$j]);
	    $statessort[$j] = $states[$k];
	    $ordersort[$j]  = $order[$k];
	    $spinsort[$j]   = $spin[$k];

	    if($gJinclude == 1) {
		$gjvalsort[$j]   = $gjval[$k];
	    }

            $spinsort2[$j] = $spinsort[$j];

	    $paritysort[$j] = $parity[$k];
	    $lmaxsort[$j] = $lmax[$k];
	    for($l=0; $l<$lmax[$k]; $l++){
		$wavemagsort[$j][$l] = $wavemag[$k][$l];
		$wavemagsortr[$j][$l] = sprintf "%.2f", $wavemagsort[$j][$l];
		$configurationsort[$j][$l]  = $configuration[$k][$l];
		$configurationsort2[$j][$l]  = $configuration2[$k][$l];
		#printf("$ordersort[$j] $spinsort[$j] $paritysort[$j] $eenergysortr[$j] $eenergysort[$j] $wavemagsortr[$j][$l] $wavemagsort[$j][$l]   $configurationsort[$j][$l]\n");
	    }
	    #printf("$k $j $energy[$k] $energysorted[$j]\n");
	}
    }
}


#JE testing
$nom = 0;
$testconfje = level2($configurationsort[$nom][1], $spinsort[$nom], $paritysort[$nom], 1, $nochar);
#printf("-----------------------\n");
#printf("$configurationsort[$nom][1]\n");
#printf("$testconfje\n");
#printf("-----------------------\n");

#PRODUCE LATEX STYLE CONFIGURATIONS
#----------------------------------
for($i=0; $i<$jmax; $i++){
    for($l=0; $l<$lmaxsort[$i]; $l++){  #energy sorted
	if($l==0){
	    $testconfsort[$i][$l] = level2($configurationsort[$i][$l], $spinsort[$i], $paritysort[$i], 0, $nochar);
	    $teststring[$i] = $wavemagsortr[$i][$l];
	    $llss[$i] = lsj($configurationsort[$i][$l]);
            $llss[$i] = "\$".$llss[$i]."_".$spinsort2[$i];
	    if($paritysort[$i] eq "-"){        #add negative parity symbol "o" to LSJ term
		$llss[$i] .= "^{\\circ}\$";
	    }else{
		$llss[$i] .= "\$";
	    }
	}else{
	    $testconfsort[$i][$l] = level2($configurationsort[$i][$l], $spinsort[$i], $paritysort[$i], 1, $nochar);
	    $llss2 = lsj($configurationsort[$i][$l]);
	    $teststring[$i] = $teststring[$i]." + ".$wavemagsortr[$i][$l]."~".$testconfsort[$i][$l];
	}
    }
}

for($i=0; $i<$jmax; $i++){
    $test = 0;
    for($j=$i+1; $j<$jmax; $j++){
	if(($i != $j) && ($configurationsort[$i][0] eq $configurationsort[$j][0]) && ($spinsort[$i] eq $spinsort[$j])){
	    if($test == 0) {
		$testconfsort[$i][0] = $testconfsort[$i][0]."\$_a\$";
		$testconfsort[$j][0] = $testconfsort[$j][0]."\$_b\$";
	    }
	    if($test == 1) {$testconfsort[$j][0] = $testconfsort[$j][0]."\$_c\$";}
	    if($test == 2) {$testconfsort[$j][0] = $testconfsort[$j][0]."\$_d\$";}
	    if($test == 3) {$testconfsort[$j][0] = $testconfsort[$j][0]."\$_e\$";}
	    if($test == 4) {$testconfsort[$j][0] = $testconfsort[$j][0]."\$_f\$";}
	    if($test == 5) {$testconfsort[$j][0] = $testconfsort[$j][0]."\$_g\$";}	    	    
	    $test++;
	}
    }
}

#PRODUCE ENERGYLABEL FILE CALLED energylabel
#-------------------------------------------
open (MYOUTPUTFILE, '>energylabel.latex');
for($j=0; $j<9; $j++){
    printf MYOUTPUTFILE "\n";
}
for($j=0; $j<3; $j++){
    printf MYOUTPUTFILE "--------------\n";
}
for($j=0; $j<$jmax; $j++){
    printf MYOUTPUTFILE "%3s %2s %3s %2s %15s %11s  %-10s          %-50s\n" , $j+1, $ordersort[$j], $spinsort[$j], $paritysort[$j], $energysortr7[$j], $eenergysortr2[$j], $statessort[$j],$testconfsort[$j][0];
}
close(MYOUTPUTFILE);



#PRODUCE LATEX TABLE WITH LEVEL INFORMATION
#------------------------------------------
$header_gJ_eobs = "No. & State & \$LS\$-composition & \$E(CI) \$ & \$E(OBS) \$  & \$g_J \$   \\\\ \n";
$header_gJ = "No. & State & \$LS\$-composition & \$E(CI) \$ & \$g_J \$   \\\\ \n";
$header_eobs = "No. & State & \$LS\$-composition & \$E(CI) \$ & \$E(OBS) \$  \\\\ \n";
$header = "No. & State & \$LS\$-composition & \$E(CI) \$ \\\\ \n";    

open (MYOUTPUTFILE2, '>lscomp.tex');
print MYOUTPUTFILE2 "\\documentclass[12pt]{article}\n";
print MYOUTPUTFILE2 "\\usepackage{longtable}\n";
print MYOUTPUTFILE2 "\\usepackage[cm]{fullpage}\n";
print MYOUTPUTFILE2 "\\thispagestyle{empty}\n";
print MYOUTPUTFILE2 "\\begin{document}\n";
print MYOUTPUTFILE2 "\{\\scriptsize\n";
if ($gJinclude == 1 && $eobsinclude == 1) {print MYOUTPUTFILE2 "\\begin{longtable}{\@\{\}rllrrr}\n";}
if ($gJinclude == 1 && $eobsinclude == 0) {print MYOUTPUTFILE2 "\\begin{longtable}{\@\{\}rllrr}\n";}
if ($gJinclude == 0 && $eobsinclude == 1) {print MYOUTPUTFILE2 "\\begin{longtable}{\@\{\}rllrr}\n";}
if ($gJinclude == 0 && $eobsinclude == 0) {print MYOUTPUTFILE2 "\\begin{longtable}{\@\{\}rllr}\n";}            
print MYOUTPUTFILE2 "\\caption\{Energies.....\}\\\\  \n";
print MYOUTPUTFILE2 "\\hline\n";
if ($gJinclude == 1 && $eobsinclude == 1) {print MYOUTPUTFILE2 "$header_gJ_eobs\n";}
if ($gJinclude == 1 && $eobsinclude == 0) {print MYOUTPUTFILE2 "$header_gJ\n";}
if ($gJinclude == 0 && $eobsinclude == 1) {print MYOUTPUTFILE2 "$header_eobs\n";}
if ($gJinclude == 0 && $eobsinclude == 0) {print MYOUTPUTFILE2 "$header\n";}            
print MYOUTPUTFILE2 "\\hline\n";
print MYOUTPUTFILE2 "\\endfirsthead\n";
print MYOUTPUTFILE2 "\\caption\{Continued.\}\\\\  \n";
print MYOUTPUTFILE2 "\\hline\n";
if ($gJinclude == 1 && $eobsinclude == 1) {print MYOUTPUTFILE2 "$header_gJ_eobs\n";}
if ($gJinclude == 1 && $eobsinclude == 0) {print MYOUTPUTFILE2 "$header_gJ\n";}
if ($gJinclude == 0 && $eobsinclude == 1) {print MYOUTPUTFILE2 "$header_eobs\n";}
if ($gJinclude == 0 && $eobsinclude == 0) {print MYOUTPUTFILE2 "$header\n";}                
print MYOUTPUTFILE2 "\\hline\n";
print MYOUTPUTFILE2 "\\endhead\n";
print MYOUTPUTFILE2 "\\hline\n";
print MYOUTPUTFILE2 "\\endfoot\n";

for ($i=0; $i<$jmax; $i++)
{
    if ($gJinclude == 1 && $eobsinclude == 1) {printf MYOUTPUTFILE2 "%-3s & %-50s & %-90s & %-12s & %-2s & %-7s\\\\ \n", $i+1, $testconfsort[$i][0], $teststring[$i], $eenergysortrsep[$i], " ", $gjvalsort[$i];}
    if ($gJinclude == 1 && $eobsinclude == 0) {printf MYOUTPUTFILE2 "%-3s & %-50s & %-90s & %-12s & %-7s\\\\ \n", $i+1, $testconfsort[$i][0], $teststring[$i], $eenergysortrsep[$i], $gjvalsort[$i];}
    if ($gJinclude == 0 && $eobsinclude == 1) {printf MYOUTPUTFILE2 "%-3s & %-50s & %-90s & %-12s & %-2s \\\\ \n", $i+1, $testconfsort[$i][0], $teststring[$i], $eenergysortrsep[$i], " ";}
    if ($gJinclude == 0 && $eobsinclude == 0) {printf MYOUTPUTFILE2 "%-3s & %-50s & %-90s & %-12s \\\\ \n", $i+1, $testconfsort[$i][0], $teststring[$i], $eenergysortrsep[$i]}            
}

print MYOUTPUTFILE2 "\\hline \n";
print MYOUTPUTFILE2 "\\end{longtable}\n";
print MYOUTPUTFILE2 "\}\n";

print MYOUTPUTFILE2 "\\end{document}\n";
close(MYOUTPUTFILE2);

print "   Files lscomp.tex and energylabel.latex written to disc.  \n";
print "\n";

#FUNCTION LSJ
#--------------
sub lsj {
    ($confstring) = @_;
    $conflength = length($confstring);
 
    $term = substr($confstring, $conflength - 2, 2);
    $sterm = substr($term, 0, 1);                                     #extract 2S+1
    $lterm = substr($term, 1, 1);                                     #extract L
    $lsj2 = "^{".$sterm."}".$lterm;                                   #form LS term

    #$lsj2 = "\$".$lsj2."\$";

    return($lsj2);
}

#FUNCTION LEVEL 2 (NEW)
#----------------------
sub level2 {
    ($confstring, $jvalue, $parity2, $levelflag, $rmnochar) = @_;
    $conflength = length($confstring);
    $conf = substr($confstring, 0, $conflength - 3);                    # remove LS term at the end (ex. _3D)
    $conf = substr($conf, $rmnochar, $conflength - 3);
    @confsplit = split('\.', $conf);                                    # split configuration by dots "."
    $size = scalar(@confsplit);                                         # number of dot separeted strings
    #printf("$conf\n");
    for ($j=$confflag; $j<$size; $j++){                                 # loop over dot seprated strings
	#printf("\n");	
	#printf("$confsplit[$j]\n");
	#printf("--------------------\n");	
        #if($confsplit[$j] =~ /\Q(\E/ && $confsplit[$j] !~ /\Q_\E/){     # if more than nl electron
        if($confsplit[$j] =~ /\Q(\E/){                                   # if more than nl electron	    
	    #printf("check\n");	
		$confpartlength = length($confsplit[$j]);
		if($confpartlength < 8){
		    $nl = substr($confsplit[$j], 0, 2);                 # extract nl
		    $checkpow = substr($confsplit[$j], 4, 1); 
		    if($checkpow eq ")"){                               # check if char 4 is ")"  example 5s(2)
			$pow = substr($confsplit[$j], 3, 1);            # if so extract power = 1 char
		    }else{
			$pow = substr($confsplit[$j], 3, 2);            # if not extract power = 2 chars
		    }
		    $nlpow = $nl."^{".$pow."}";                         # contruct nl^pow latex style
		    #printf("$nlpow\n");
		}else{
		    $nl = substr($confsplit[$j], 0, 2);
		    $checkpow = substr($confsplit[$j], 4, 1); 
		    if($checkpow eq ")"){                               # check if char 4 is ")"  example 5s(2)
			$pow = substr($confsplit[$j], 3, 1);            # if so extract power = 1 char
		    }else{
			$pow = substr($confsplit[$j], 3, 2);            # if not extract power = 2 chars
		    }
		    if($confsplit[$j] =~ /\Q_\E/){
			#printf("end\n");
			$intterm3 = substr($confsplit[$j], $confpartlength - 2, 1);
			$intterm4 = substr($confsplit[$j], $confpartlength - 1, 1);
			$intterm5 = substr($confsplit[$j], $confpartlength - 6, 1);
			$intterm6 = substr($confsplit[$j], $confpartlength - 5, 1);
			$intterm7 = substr($confsplit[$j], $confpartlength - 4, 1);		    
			$nlpow = $nl."^{".$pow."}"."(^{$intterm5}_{$intterm7}$intterm6)"."~^{$intterm3}$intterm4";       # contruct nl^pow latex style
		    }else{
			#printf("not end\n");
			$intterm5 = substr($confsplit[$j], $confpartlength - 3, 1);
			$intterm6 = substr($confsplit[$j], $confpartlength - 2, 1);
			$intterm7 = substr($confsplit[$j], $confpartlength - 1, 1);		    
			$nlpow = $nl."^{".$pow."}"."(^{$intterm5}_{$intterm7}$intterm6)";       # contruct nl^pow latex style
			#$nlpow = $nl."(^{".$intterm3."}".$intterm4.")";
			#printf("$nlpow\n");
		    }
		}
        }else{
	    #printf("check 2\n");	                                                      #if exactly 1 nl electron
	    if($confsplit[$j] =~ /\Q_\E/){
		$nl = substr($confsplit[$j], 0, 2);
		$intterm1 = substr($confsplit[$j], 3, 1);
		$intterm2 = substr($confsplit[$j], 4, 1);
		$nlpow = $nl."~^{".$intterm1."}".$intterm2;
	    }else{
		$nlpow = substr($confsplit[$j], 0, 2);
	    }
        }
        if($j == $confflag){
            $conf2 = $nlpow;
        }else{
            $conf2 = $conf2."\\,".$nlpow;
        }
    }
    $term = substr($confstring, $conflength - 2, 2);
    $sterm = substr($term, 0, 1);                                     #extract 2S+1
    $lterm = substr($term, 1, 1);                                     #extract L
    $lsj = "^{".$sterm."}".$lterm."\_{".$jvalue."}";                  #form LSJ term
    $lsj2 = "^{".$sterm."}".$lterm;                                   #form LS term
    if($parity2 eq "-"){                                              #add negative parity symbol "o" to LSJ term
        $lsj .= "^{\\circ}";
        $lsj2 .= "^{\\circ}";
    }
    if($levelflag == 0){
	$conf2 = $conf2."~".$lsj;
    }else{
	$conf2 = $conf2."~".$lsj2;
    }
    $conf2 = "\$".$conf2."\$";

    return($conf2);
}

#FUNCTION THOUSANDSEP
#--------------------
sub thousandsep {
    ($unsepenergy) = @_;
    $lenergy = length($unsepenergy);
    if($lenergy < 4){
	$sepenergy = $unsepenergy;
    }elsif($lenergy > 3 && $lenergy < 7){
	$lastpart = substr($unsepenergy, $lenergy - 3, 3);
	$firstpart = substr($unsepenergy, 0, $lenergy - 3);
	$sepenergy = $firstpart."~".$lastpart;
    }elsif($lenergy > 6){
	$lastpart = substr($unsepenergy, $lenergy - 3, 3);
	$middlepart = substr($unsepenergy, $lenergy - 6, 3);
	$firstpart = substr($unsepenergy, 0, $lenergy - 6);
	$sepenergy = $firstpart."~".$middlepart."~".$lastpart;
    }
    return($sepenergy);
}

#!/usr/local/bin/perl
# For calculating energy and variance after the single-step Lanczos 
  #input!!
  $Sample=5; # # of sample

#===============Read input file====================
  for($i=1;$i<=$Sample;$i++){
    #$cnt=0;
    $file=sprintf("zvo_aft_Lz_ls_00%d.dat",$i);
    #print "$file \n";
    open(INPUTFILE,$file);
    while($name=<INPUTFILE>){
      chomp $name;
      #DELETE EMPTY
      $_=$name; 
      s/^\s+//;
      $name=$_; 
      @tmp = split /\s+/, $name;
      #printf "$tmp[0]\n";
      $r_H[$i]      = $tmp[0];
      $r_H2_1[$i]   = $tmp[1];
      $r_H2_2[$i]   = $tmp[2];
      $r_H3[$i]     = $tmp[3];
      $r_H4[$i]     = $tmp[4];
    }
    close(INPUTFILE);
  }
#== Determination of suitable alpha =========
  for($i=1;$i<=$Sample;$i++){
    $a   = $r_H[$i];
    $b   = $r_H2_1[$i];
    $c   = $r_H2_2[$i];
    $d   = $r_H3[$i];
    $tmp_AA  = $b*($b+$c)-2*$a*$d;
    $tmp_BB  = -$a*$b+$d;
    $tmp_CC  = $b*(($b+$c)**2)-($a**2)*$b*($b+2*$c)+4*($a**3)*$d-2*$a*(2*$b+$c)*$d+$d**2;
    $tmp_xp  = ($tmp_BB+sqrt($tmp_CC))/$tmp_AA;
    $tmp_xm  = ($tmp_BB-sqrt($tmp_CC))/$tmp_AA;
    $alpha_p[$i] = $tmp_xp ;
    $alpha_m[$i] = $tmp_xm ;
    printf("$i: %lf %lf %lf: %lf %lf\n",$a,$b,$c,$tmp_xp,$tmp_xm);
  }
#=== calculating energy and variance ========================
  $fname="Result_Lz.dat";
  open(FILE,">$fname");

  $Ave_alpha_p = 0.0;
  $Ave_alpha_m = 0.0;
  $Err_alpha_p = 0.0;
  $Err_alpha_m = 0.0;
  $Ave_Lz_V_p  = 0.0;
  $Err_Lz_V_p  = 0.0;

  for($i=1;$i<=$Sample;$i++){
    $alpha        = $alpha_p[$i];
    $Ave_alpha_p += $alpha_p[$i]; 
    $Err_alpha_p += $alpha_p[$i]**2; 
    $tmp_1        = $r_H[$i]+$alpha*($r_H2_1[$i]+$r_H2_2[$i])+$alpha*$alpha*$r_H3[$i];
    $tmp_2        = 1.0+2*$alpha*$r_H[$i]+$alpha*$alpha*$r_H2_1[$i];
    $tmp_3        = $r_H2_1[$i]+2*$alpha*$r_H3[$i]+$alpha*$alpha*$r_H4[$i];
    $tmp_V        =  (($tmp_3/$tmp_2)-($tmp_1/$tmp_2)**2)/($tmp_1/$tmp_2)**2;
    $Ave_Lz_E_p  += ($tmp_1/$tmp_2);
    $Err_Lz_E_p  += ($tmp_1/$tmp_2)**2;
    $Ave_Lz_V_p  += $tmp_V;
    $Err_Lz_V_p  += ($tmp_V)**2;

    $alpha        = $alpha_m[$i];
    $Ave_alpha_m += $alpha_m[$i]; 
    $Err_alpha_m += $alpha_m[$i]**2; 
    $tmp_1        = $r_H[$i]+$alpha*($r_H2_1[$i]+$r_H2_2[$i])+$alpha*$alpha*$r_H3[$i];
    $tmp_2        = 1.0+2*$alpha*$r_H[$i]+$alpha*$alpha*$r_H2_1[$i];
    $tmp_3        = $r_H2_1[$i]+2*$alpha*$r_H3[$i]+$alpha*$alpha*$r_H4[$i];
    $tmp_V        =  (($tmp_3/$tmp_2)-($tmp_1/$tmp_2)**2)/($tmp_1/$tmp_2)**2;
    $Ave_Lz_E_m  += ($tmp_1/$tmp_2);
    $Err_Lz_E_m  += ($tmp_1/$tmp_2)**2;
    $Ave_Lz_V_m  += $tmp_V;
    $Err_Lz_V_m  += ($tmp_V)**2;
  }
  $H_p     =  $Ave_Lz_E_p/(1.0*$Sample);
  $H_p_err =  sqrt($Err_Lz_E_p/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$H_p**2);
  $V_p     =  $Ave_Lz_V_p/(1.0*$Sample);
  $V_p_err =  sqrt($Err_Lz_V_p/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$V_p**2);
  $a_p     =  $Ave_alpha_p/(1.0*$Sample);
  $a_p_err =  sqrt($Err_alpha_p/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$a_p**2);

  $H_m     =  $Ave_Lz_E_m/(1.0*$Sample);
  $H_m_err =  sqrt($Err_Lz_E_m/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$H_m**2);
  $V_m     =  $Ave_Lz_V_m/(1.0*$Sample);
  $V_m_err =  sqrt($Err_Lz_V_m/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$V_m**2);
  $a_m     =  $Ave_alpha_m/(1.0*$Sample);
  $a_m_err =  sqrt($Err_alpha_m/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$a_m**2);

  printf("%lf %lf %lf %lf %lf %lf    \n",$a_p,$a_p_err,$H_p,$H_p_err,$V_p,$V_p_err);
  printf("%lf %lf %lf %lf %lf %lf    \n",$a_m,$a_m_err,$H_m,$H_m_err,$V_m,$V_m_err);
  printf FILE ("%lf %lf %lf %lf %lf %lf    \n",$a_p,$a_p_err,$H_p,$H_p_err,$V_p,$V_p_err);
  printf FILE ("%lf %lf %lf %lf %lf %lf    \n",$a_m,$a_m_err,$H_m,$H_m_err,$V_m,$V_m_err);

  close(FILE); 
 
#===============for check====================
  $fname="Lz_energy.dat";
  open(FILE,">$fname");
  $fname_2="Lz_variance.dat";
  open(FILE_2,">$fname_2");

#===============changing alpha |psi>=(1+alpha*H)|phi> ===================
# H_b,V_b -> better
# H,V     -> simple

  $min_alpha=0.0;
  $max_alpha=0.0;
  $delta_alpha=0.0;
  $period = 20.0;
  for($i=1;$i<=$Sample;$i++){
    if($alpha_p[$i]>$alpha_m[$i]){
      $min_alpha   += $alpha_m[$i];
      $max_alpha   += $alpha_p[$i];
      $delta_alpha += ($alpha_p[$i]-$alpha_m[$i])/$period;
    }else{
      $min_alpha   += $alpha_p[$i];
      $max_alpha   += $alpha_m[$i];
      $delta_alpha += ($alpha_m[$i]-$alpha_p[$i])/$period;
    }
  }
  $delta_alpha = $delta_alpha/(1.0*$Sample);
  $min_alpha   = $min_alpha/(1.0*$Sample)-10*$delta_alpha;
  $max_alpha   = $max_alpha/(1.0*$Sample)+10*$delta_alpha;
  printf("$min_alpha $max_alpha $delta_alpha\n");
  for($alpha= $min_alpha  ;$alpha<= $max_alpha  ;$alpha+=$delta_alpha){
    #printf("$alpha \n");
    $Ave_Lz_E = 0.0;
    $Err_Lz_E = 0.0;
    $Ave_Lz_E_b = 0.0;
    $Err_Lz_E_b = 0.0;

    $Ave_Lz_V = 0.0;
    $Err_Lz_V = 0.0;
    $Ave_Lz_V_b = 0.0;
    $Err_Lz_V_b = 0.0;
    for($i=1;$i<=$Sample;$i++){
      $a   = $r_H[$i];
      $b   = $r_H2_1[$i];
      $c   = $r_H2_2[$i];
      $d   = $r_H3[$i];
      $tmp_AA  = $b*($b+$c)-2*$a*$d;
      $tmp_BB  = -$a*$b+$d;
      $tmp_CC  = $b*(($b+$c)**2)-($a**2)*$b*($b+2*$c)+4*($a**3)*$d-2*$a*(2*$b+$c)*$d+$d**2;
      $tmp_xp  = ($tmp_BB+sqrt($tmp_CC))/$tmp_AA;
      $tmp_xm  = ($tmp_BB-sqrt($tmp_CC))/$tmp_AA;
      $alpha_p[$i] = $tmp_xp ;
      $alpha_m[$i] = $tmp_xm ;
      #printf("$i: %lf %lf %lf: %lf %lf\n",$a,$b,$c,$tmp_xp,$tmp_xm);
      $tmp_1   = $r_H[$i]+2*$alpha*$r_H2_1[$i]+$alpha*$alpha*$r_H3[$i];
      $tmp_2   = 1.0+2*$alpha*$r_H[$i]+$alpha*$alpha*$r_H2_1[$i];
      $tmp_3   = $r_H2_1[$i]+2*$alpha*$r_H3[$i]+$alpha*$alpha*$r_H4[$i];
      $tmp_V   =  (($tmp_3/$tmp_2)-($tmp_1/$tmp_2)**2)/($tmp_1/$tmp_2)**2;
      $Ave_Lz_E += ($tmp_1/$tmp_2);
      $Err_Lz_E += ($tmp_1/$tmp_2)**2; 
      $Ave_Lz_V += ($tmp_V);
      $Err_Lz_V += ($tmp_V)**2; 

      $tmp_1   = $r_H[$i]+$alpha*($r_H2_1[$i]+$r_H2_2[$i])+$alpha*$alpha*$r_H3[$i];
      $tmp_2   = 1.0+2*$alpha*$r_H[$i]+$alpha*$alpha*$r_H2_1[$i];
      $tmp_3   = $r_H2_1[$i]+2*$alpha*$r_H3[$i]+$alpha*$alpha*$r_H4[$i];
      $tmp_V   =  (($tmp_3/$tmp_2)-($tmp_1/$tmp_2)**2)/($tmp_1/$tmp_2)**2;
      $Ave_Lz_E_b += ($tmp_1/$tmp_2);
      $Err_Lz_E_b += ($tmp_1/$tmp_2)**2; 
      $Ave_Lz_V_b += $tmp_V;
      $Err_Lz_V_b += ($tmp_V)**2; 
    }
    $H     =  $Ave_Lz_E/(1.0*$Sample);
    $H_err =  sqrt($Err_Lz_E/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$H**2);
    $V     =  $Ave_Lz_V/(1.0*$Sample);
    $V_err =  sqrt($Err_Lz_V/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$V**2);

    $H_b     =  $Ave_Lz_E_b/(1.0*$Sample);
    $H_b_err =  sqrt($Err_Lz_E_b/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$H_b**2);
    $V_b     =  $Ave_Lz_V_b/(1.0*$Sample);
    $V_b_err =  sqrt($Err_Lz_V_b/(1.0*($Sample-1.0))-$Sample/($Sample-1.0)*$V_b**2);

    #printf("%lf %lf %lf %lf %lf   \n",$alpha,$H,$H_err,$H_b,$H_b_err);
    #printf("%lf %lf %lf %lf %lf   \n",$alpha,$V,$V_err,$V_b,$V_b_err);
    printf FILE ("%lf %lf %lf %lf %lf   \n",$alpha,$H_b,$H_b_err,$H,$H_err);
    printf FILE_2 ("%lf %lf %lf %lf %lf   \n",$alpha,$V_b,$V_b_err,$V,$V_err);
  }
  close(FILE); 
  close(FILE_2); 


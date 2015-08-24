#!/usr/bin/perl -w
  $PI=3.14159265358979;
  #input!!
  &input;
  #input!!
  $orb_num=$tmp_orb;
  $L_x=$tmp_Lx;
  $L_y=$tmp_Ly;
  $L=$L_x*$L_y;
  $Ns=$L;
  $All_N=$Ns*$orb_num;
  $Particle=1;
  printf("CHECK L_x=$L_x L_y=$L_y \n");

  $tmp=2*$All_N;
  $fname=sprintf("zinitial.def");
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NInitial $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========initial Green functions ====== \n";
  printf FILE "========misawa======== \n";

  $Qx        = 1.0;
  $Qy        = 1.0;
  $Intensity = 0.3;
  for($all_i=0;$all_i<$All_N;$all_i++){
    $orb_i  = $all_i%$orb_num;  
    $site_i = ($all_i-$orb_i)/$orb_num;  
    $int_x  = $site_i%$L_x;
    $int_y  = ($site_i-$int_x)/$L_x;
    $mg     = $Intensity*cos($Qx*$PI*$int_x+$Qy*$PI*$int_y);
    $tmp    = 0.5*($Particle+$mg);   
    printf FILE (" %4d %4d  0 %f\n",$all_i,$all_i,$tmp); 
  }
  for($all_i=0;$all_i<$All_N;$all_i++){
    $orb_i  = $all_i%$orb_num;  
    $site_i = ($all_i-$orb_i)/$orb_num;  
    $int_x  = $site_i%$L_x;
    $int_y  = ($site_i-$int_x)/$L_x;
    $mg     = $Intensity*cos($Qx*$PI*$int_x+$Qy*$PI*$int_y);
    $tmp    = 0.5*($Particle-$mg);   
    printf FILE (" %4d %4d  1 %f\n",$all_i,$all_i,$tmp); 
  }

  close(FILE);
 sub input{
  #input START 
  $Lx_cnt=0;
  $Ly_cnt=0;
  $orb_cnt=0;
  $file=sprintf("input.txt");
  open(INPUTFILE,$file);
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    @tmp = split /\s+/, $name;
    #printf "$tmp[0] $tmp[1] \n";
    if($tmp[0] eq 'Lx'){
      #printf "AA $tmp[0] $tmp[1] \n";
      $tmp_Lx = $tmp[1];
      $Lx_cnt=1;
    } 
    if($tmp[0] eq 'Ly'){
      #printf "AA $tmp[0] $tmp[1] \n";
      $tmp_Ly = $tmp[1];
      $Ly_cnt=1;
    } 
    if($tmp[0] eq 'orb_num'){
      #printf "AA $tmp[0] $tmp[1] \n";
      $tmp_orb = $tmp[1];
      $orb_cnt=1;
    } 
  }
  if($Lx_cnt==0 || $Ly_cnt==0||$orb_cnt==0){
    printf "FAITAL ERROR IN input.txt !!!!!!!!!!!!! \n";
  }
  #input FINISH
  close(INPUTFILE);
 }

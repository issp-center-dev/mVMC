#!/usr/bin/perl
  #print "start !! \n";
  $PI=3.14159265358979;
  #input!!
  &input;
  $orb_num=$tmp_orb;
  #input!!
  $L_x=$tmp_Lx;
  $L_y=$tmp_Ly;
  $L=$L_x*$L_y;
  $Ns=$L;
  $All_N=$Ns*$orb_num;
  printf("CHECK $All_N L_x=$L_x L_y=$L_y  orb=$orb_num \n");

  $cnt=0;
  $file=sprintf("zvo_UHF_cisajs.dat");
  #print "$file \n";
  open(INPUTFILE,$file);
  while($name=<INPUTFILE>){
     chomp $name;
     #DELETE EMPTY
     $_=$name; 
     s/^\s+//;
     $name=$_; 
     #DELETE EMPTY FINISH
     @foo = split /\s+/, $name;
     #printf "$cnt $foo[0] $foo[1] $foo[2] $foo[3] \n";
     $all_i                       = $foo[0]; 
     $all_j                       = $foo[1]; 
     $spin                        = $foo[2]; 
     $Green[$all_i][$all_j][$spin] = $foo[3];
     $cnt+=1;
  }
  close(INPUTFILE);
  #printf "$cnt\n";

 for($orb_cnt=0;$orb_cnt<$orb_num;$orb_cnt++){
   $fname=sprintf("Result_UHF_local_$orb_cnt.dat");
   open(FILE,">$fname");
   for($i_y=0;$i_y<$L_y;$i_y++){
     for($i_x=0;$i_x<$L_x;$i_x++){
       $all_i = $orb_cnt+($i_x+$i_y*$L_x)*$orb_num;
       $charge = $Green[$all_i][$all_i][0]+$Green[$all_i][$all_i][1];
       $spin   = $Green[$all_i][$all_i][0]-$Green[$all_i][$all_i][1];
       $TSpin[$orb_cnt][$i_x][$i_y]   = $spin; 
       $TCharge[$orb_cnt][$i_x][$i_y] = $charge; 
       printf FILE ("%4d %4d %f %f \n",$i_x,$i_y,$charge,$spin);
     }
     printf FILE ("\n");
   }
   close(FILE);
 } 

 for($orb_cnt=0;$orb_cnt<$orb_num;$orb_cnt++){
   $fname=sprintf("Result_UHF_Q_$orb_cnt.dat");
   open(FILE,">$fname");
   for($kx=0;$kx<=$L_x;$kx+=1){
     for($ky=0;$ky<=$L_y;$ky+=1){
       $All_Spin     = 0.0;
       $All_Charge   = 0.0;
       $Total_Charge = 0.0;
       for($i_x=0;$i_x<$L_x;$i_x++){
         for($i_y=0;$i_y<$L_y;$i_y++){
           $diff_x        = 2*$PI*($i_x)/$L_x;
           $diff_y        = 2*$PI*($i_y)/$L_y;
           $All_Spin     += $TSpin[$orb_cnt][$i_x][$i_y]*cos($diff_x*$kx+$diff_y*$ky);
           $All_Charge   += $TCharge[$orb_cnt][$i_x][$i_y]*cos($diff_x*$kx+$diff_y*$ky);
           $Total_Charge += $TCharge[$orb_cnt][$i_x][$i_y];
         }
       }
       $All_Spin     = $All_Spin/($L_x*$L_y);
       if($kx%$L_x==0 && $ky%$L_y==0 ){
         $All_Charge   = 0.0;
         $Total_Charge   = $Total_Charge/($L_x*$L_y);
         printf FILE ("%4d %4d %f %f %f\n",$kx,$ky,$All_Charge,$All_Spin,$Total_Charge);
       }else{
         $All_Charge   = $All_Charge/($L_x*$L_y);
         printf FILE ("%4d %4d %f %f \n",$kx,$ky,$All_Charge,$All_Spin);
       } 
     } 
     printf FILE ("\n");
   } 
   close(FILE);
 }
 


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
 }

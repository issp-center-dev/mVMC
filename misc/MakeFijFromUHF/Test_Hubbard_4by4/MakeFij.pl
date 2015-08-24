#!/usr/bin/perl
  #input!!
  &input;
  #input!!
  $orb_num=$tmp_orb;
  $L_x=$tmp_Lx;
  $L_y=$tmp_Ly;
  $L=$L_x*$L_y;
  $Ns=$L;
  $All_N=$Ns*$orb_num;
  #readdef !!!
  &read_def;
  #doublecount !!!

  printf("CHECK: dc_fij = $dc_fij \n");
  printf("CHECK: L_x=$L_x L_y=$L_y orb=$orb_num \n");
  printf("CHECK: N_G = $N_G, N_J = $N_J, N_DH2 = $N_DH2, N_DH4 = $N_DH4, N_fij = $N_fij\n");
  $N_all    = 2+$N_G+$N_J+$N_DH2+$N_DH4+$N_fij;

  printf("\n");

  $cnt = 0;
  $file=sprintf("zvo_fij.dat");
  #print "$file \n";
  open(INPUTFILE,$file);
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    #DELETE EMPTY FINISH
    @tmp = split /\s+/, $name;
    $all_i                               = $tmp[0]; 
    $all_j                               = $tmp[1]; 

    $Input_Fij[$all_i][$all_j]           = $tmp[2];
    $cnt                                += 1;  
  }
  close(INPUTFILE);

  $N_wo_G   = 2;
  $N_wo_J   = 2+$N_G;
  $N_wo_DH  = 2+$N_G+$N_J;
  $N_wo_fij = 2+$N_G+$N_J+$N_DH2+$N_DH4;
  $N_all    = 2+$N_G+$N_J+$N_DH2+$N_DH4+$N_fij;
  for($cnt=0;$cnt<2*2;$cnt+=1){
     $var[$cnt] = 0;
  }
  #Gutz
  for($cnt=2*$N_wo_G;$cnt<2*(2+$N_G);$cnt+=1){
    if($cnt%2==0){
     $cnt_tmp    = (($cnt-2*$N_wo_G)/2); 
     #$Gutz[$cnt_tmp] = $tmp_3[$cnt];
     $var[$cnt] = 0.0;#$Small_Gutz[$cnt_tmp];
    }else{
      $var[$cnt] = 0;
    }
  }
  #Jast
  for($cnt=2*$N_wo_J;$cnt<2*(2+$N_G+$N_J);$cnt+=1){
    if($cnt%2==0){
     $cnt_tmp    = (($cnt-2*$N_wo_J)/2); 
     #$Jast[$cnt_tmp] = $tmp_3[$cnt];
     $var[$cnt] = 0.0;#$Small_Jast[$cnt_tmp];
    }else{
      $var[$cnt] = 0;
    }
  }
  # DH
  #printf("AAAA:$N_wo_DH B:$N_wo_fij C:$N_DH2 D:$N_DH4\n");
  for($cnt=2*$N_wo_DH;$cnt<2*($N_wo_fij);$cnt+=1){
    #printf ("$cnt \n");
    if($cnt%2==0){
     $cnt_tmp    = (($cnt-2*$N_wo_DH)/2); 
     #$Jast[$cnt_tmp] = $tmp_3[$cnt];
     $var[$cnt] = 0.0;#$Small_Jast[$cnt_tmp];
    }else{
      $var[$cnt] = 0;
    }
  }


  #printf("$S_Small_fij[-1][-1][4][4][3] \n");
  for($cnt=2*$N_wo_fij;$cnt<2*$N_all;$cnt+=1){
    if($cnt%2==0){
      $cnt_tmp    = (($cnt-2*$N_wo_fij)/2); 
      $all_i      = $inv_fij_i[0][$cnt_tmp];
      $all_j      = $inv_fij_j[0][$cnt_tmp];

      #$orb_i  = $all_i % $orb_num;
      #$site_i = ($all_i-$orb_i)/$orb_num;
      #$x_i    = $site_i % $L_x;    
      #$y_i    = ($site_i-$x_i) / $L_x;    

      #$orb_j  = $all_j % $orb_num;
      #$site_j = ($all_j-$orb_j)/$orb_num;
      #$x_j    = $site_j % $L_x;    
      #$y_j    = ($site_j-$x_j) / $L_x;    

      $var[$cnt] = $Input_Fij[$all_i][$all_j];#+0.001*(1-2*rand());
      #printf("$all_i $all_j;  $var[$cnt] \n");
    }else{
      $var[$cnt]   = 0.0;
    }
  }

  $fname="VMC_zinitial.def";
  open(FILE,">$fname");
  for($cnt=0;$cnt<$N_all*2;$cnt+=1){
    printf FILE (" $var[$cnt] ");
  }
  printf FILE (" \n ");
  close(FILE);

  $fname="result.dat";
  open(FILE,">$fname");
  for($cnt=0;$cnt<$N_all*2;$cnt+=1){
    printf FILE (" $cnt $var[$cnt] \n");
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
  close(INPUTFILE);
  #input FINISH
 }

 sub read_def{
  #input START 
  $file=sprintf("zgutzwilleridx.def");
  open(INPUTFILE,$file);
  $cnt = 0;
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    @tmp = split /\s+/, $name;
    #printf "$tmp[0] $tmp[1] \n";
    if($cnt == 1){
      $N_G = $tmp[1];
      last;
    } 
    $cnt += 1;
  }
  close(INPUTFILE);
  #input FINISH

  #input START 
  $file=sprintf("zjastrowidx.def");
  open(INPUTFILE,$file);
  $cnt = 0;
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    @tmp = split /\s+/, $name;
    #printf "$tmp[0] $tmp[1] \n";
    if($cnt == 1){
      $N_J = $tmp[1];
      last;
    } 
    $cnt += 1;
  }
  close(INPUTFILE);
  #input FINISH

  #input START 
  $file=sprintf("zdoublonholon2siteidx.def");
  open(INPUTFILE,$file);
  $cnt = 0;
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    @tmp = split /\s+/, $name;
    #printf "$tmp[0] $tmp[1] \n";
    if($cnt == 1){
      $N_DH2 = 6*$tmp[1];
      last;
    } 
    $cnt += 1;
  }
  close(INPUTFILE);
  #input FINISH

  #input START 
  $file=sprintf("zdoublonholon4siteidx.def");
  open(INPUTFILE,$file);
  $cnt = 0;
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    @tmp = split /\s+/, $name;
    #printf "$tmp[0] $tmp[1] \n";
    if($cnt == 1){
      $N_DH4 = 10*$tmp[1];
      last;
    } 
    $cnt += 1;
  }
  close(INPUTFILE);
  #input FINISH


  #input START 
  $file=sprintf("zorbitalidx.def");
  open(INPUTFILE,$file);
  $cnt = 0;
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    @tmp = split /\s+/, $name;
    #printf "$tmp[0] $tmp[1] \n";
    if($cnt == 1){
      $N_fij = $tmp[1];
      last;
    } 
    $cnt += 1;
  }
  close(INPUTFILE);
  #Make inv
  $dc_fij  = $All_N*$All_N/$N_fij;
  for($cnt = 0;$cnt<$N_fij;$cnt++){
    $inv_fij_cnt[$cnt] = 0;
  }
  #input START 
  $file=sprintf("zorbitalidx.def");
  open(INPUTFILE,$file);
  $cnt = 0;
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    @tmp = split /\s+/, $name;
    #printf "$tmp[0] $tmp[1] \n";
    if($cnt > 4 && $cnt < $All_N**2+5){
      $all_i                         = $tmp[0]; 
      $all_j                         = $tmp[1]; 
      $num                           = $tmp[2]; 
      #$num_fij[$all_i][$all_j]       = $num;

      $tmp_cnt                       = $inv_fij_cnt[$num];  
      $inv_fij_cnt[$num]            += 1;      
      
      $inv_fij_i[$tmp_cnt][$num]     = $all_i;
      $inv_fij_j[$tmp_cnt][$num]     = $all_j;
    } 
    if($cnt >= $All_N**2+5){
      last;
    } 
    $cnt += 1;
  }
  close(INPUTFILE);
  #input FINISH

  for($cnt = 0;$cnt<$N_fij;$cnt++){
    $tmp          = $inv_fij_cnt[$cnt];
    if($tmp != $dc_fij){
      printf("fij FAITAL ERROR \n");
    }
  }
} 

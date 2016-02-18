#!/usr/local/bin/perl
  use Term::ANSIColor qw(:constants);
  #$Term::ANSIColor::AUTORESER = 1;
  
  #printf("$ARGV[0] \n");

  if($ARGV[0]==144){
    $rack_name[0] = "r3";
    $rack_name[1] = "r4";
    $rack_name[2] = "r7";
    $rack_name[3] = "r8";

    $max_node[0] = 288;
    $max_node[1] = 288;
    $max_node[2] = 288;
    $max_node[3] = 144;

    $tmp_cnt_max = 4;
    $nodename    = "F$ARGV[0]cpu";
  }elsif($ARGV[0]==72){
    $rack_name[0] = "r6";
    $max_node[0]  = 144;
    $tmp_cnt_max  = 1;
    $nodename     = "F$ARGV[0]acc";
  }elsif($ARGV[0]==36){
    $rack_name[0] = "r2";
    $max_node[0]  = 288;
    $tmp_cnt_max  = 1;
    $nodename     = "F$ARGV[0]cpu";
  }elsif($ARGV[0]==18){
    $rack_name[0] = "r5";
    $max_node[0]  = 108;
    $tmp_cnt_max  = 1;
    $nodename     = "F$ARGV[0]acc";
  }elsif($ARGV[0]==4){
    $rack_name[0] = "r1";
    $max_node[0]  = 216;
    $tmp_cnt_max  = 1;
    $nodename     = "F$ARGV[0]cpu";
  }else{
    printf BOLD RED;
    printf("INPUT ERROR \n");
    printf RESET;
    exit;
  }

  $home = $ENV{"HOME"};
  $td="$home/Scripts";
  #printf "$PWD \n";
#===============def input====================
  $cnt=0;
  $file=sprintf("$td/tmp_f$ARGV[0]");
  open (INPUTFILE, $file);
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name;
    s/^\s+//;
    $name=$_;
    @tmp = split /\s+/, $name;
    
    $job[$cnt][0] = $tmp[0]; 
    $job[$cnt][1] = $tmp[1]; 
    $job[$cnt][2] = $tmp[2]; 
    $job[$cnt][3] = $tmp[3]; 
    $job[$cnt][4] = $tmp[4]; 
    $job[$cnt][5] = $tmp[5]; 

    @tmp_2 = split /\./, $tmp[0];
    #printf "$tmp[0] $tmp[5]\n";
    #printf "$tmp_2[0] $tmp_2[1]\n";
    $jobid[$cnt] = $tmp_2[0];
    #system("qstat -f $tmp_2[0] ");
    $cnt += 1;
  }
  close(INPUTFILE);

  $cnt_max = $cnt;
  
  for($tmp_cnt=0;$tmp_cnt<$tmp_cnt_max;$tmp_cnt++){
     $rack_used{$rack_name[$tmp_cnt]} = 0;
  }

  $run_cnt = 0;
  for($cnt=0;$cnt<$cnt_max;$cnt++){
    if($job[$cnt][4] eq "R"){
      system("finger $job[$cnt][2] > $td/info_finger");
      $file=sprintf("$td/info_finger");
      open (INPUTFILE, $file);
      $tmp_cnt = 0;
      while($name=<INPUTFILE>){
        if($tmp_cnt==1){
          last;
        }
        chomp $name;
        #DELETE EMPTY
        $_=$name;
        s/^\s+//;
        $name=$_;
        @tmp = split /\:+/, $name;
        #printf("$tmp[2] \n");  
        $tmp_cnt += 1;
      }
      close(INPUTFILE);

      $run_jobid[$run_cnt]  = $jobid[$cnt] ;
      $run_userid[$run_cnt] = $job[$cnt][2] ;
      $run_type[$run_cnt]   = $job[$cnt][5] ;
      $run_name[$run_cnt]   = $tmp[2];
      #printf("jobid=$jobid[$cnt] user id=$job[$cnt][2]: name=$tmp[2]\n");

      #system("qstat -f $jobid[$cnt] |grep -e nodect -e  exec_host > info_$jobid[$cnt]"); 
      #$file=sprintf("info_$jobid[$cnt]");
      system("qstat -f $jobid[$cnt] |grep -e nodect -e  exec_host > $td/tmp_info"); 
      system("qstat -f $jobid[$cnt] |grep resources_used.walltime >> $td/tmp_info"); 
      $file=sprintf("$td/tmp_info");
      open (INPUTFILE, $file);
      $tmp_cnt = 0;
      while($name=<INPUTFILE>){
        chomp $name;
        #DELETE EMPTY
        $_=$name;
        s/^\s+//;
        $name=$_;
        @tmp = split /\s+/, $name;
        if($tmp_cnt == 0){
          @tmp_2 = split /i+/, $tmp[2];
          #printf("$tmp_2[0] $tmp_2[1] \n");  
          $rack = $tmp_2[0];
          #@tmp_3 = split /\r/, $tmp[0];
          #printf("$tmp_3[0] \n");  
        }elsif($tmp_cnt==1){
          $rack_used{"$rack"} += $tmp[2];
          $run_used[$run_cnt]  = $tmp[2];
          $run_rack[$run_cnt]  = $rack;
        }elsif($tmp_cnt==2){
          $run_time[$run_cnt]  = $tmp[2];
        }
        $tmp_cnt += 1;
      }
      close(INPUTFILE);
      $run_cnt+=1;
    }
  }
  $run_cnt_max =$run_cnt;
  printf(" \n");
  printf("$nodename analysis: \n");
  for($run_cnt=0;$run_cnt<$run_cnt_max;$run_cnt++){
     printf("jobid = $run_jobid[$run_cnt] userid = $run_userid[$run_cnt] ($run_name[$run_cnt]) :nodes = $run_used[$run_cnt]: rackname = $run_rack[$run_cnt]: time = $run_time[$run_cnt] ($run_type[$run_cnt]) \n");
  }
  printf(" \n");

  printf("$nodename nodes analysis: \n");
  for($tmp_cnt=0;$tmp_cnt<$tmp_cnt_max;$tmp_cnt++){
    $rest = $max_node[$tmp_cnt]-$rack_used{$rack_name[$tmp_cnt]};
    printf("rack name = $rack_name[$tmp_cnt] (max nodes = $max_node[$tmp_cnt]): used = $rack_used{$rack_name[$tmp_cnt]} rest = $rest \n");
  }
  printf(" \n");

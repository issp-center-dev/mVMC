qstat|grep -e F$1 -e L$1 -e B$1 > $HOME/Scripts/tmp_f$1
perl -w $HOME/Scripts/AllNodes.pl $1
date >> $HOME/Scripts/rec_f$1.txt
perl -w $HOME/Scripts/AllNodes.pl $1 >> $HOME/Scripts/rec_f$1.txt


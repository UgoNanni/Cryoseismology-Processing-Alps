#!/bin/bash

# ==========================================================

# ----------------------------------------------------------

# DATA
path2dat=/media/ugonanni/2E3D-8697/DATA/ARGENTIERE/ARG_Borehole_Ugo
path2mseed=$path2dat/MSEED
path2datsac=$path2dat/SACt
path2sac03=$path2datsac/03files
path2sacraw=$path2datsac/RAW
path2saccorr=$path2datsac/CORR
path2sacday=$path2datsac/DAY

path2PZ=$path2dat/Information

home=/home/ugonanni/Share/PhD/Processing/scripts_BASH

# ----------------------------------------------------------
NAME_NET=ARG
NAME_STAT=B01

declare -a arrayy=("2000" "2017" "2018")
year=${#arrayy[@]}

declare -a arrayc=("Z")
component=${#arrayc[@]}

ls $path2mseed   > $path2mseed/mseed_list
cd $path2sacraw
# ----------------------------------------------------------
for (( yy=1; yy<${year}+1; yy++ ))
do  YEAR=${arrayy[$yy-1]}
    
    # DAY LOOP
    for (( jj=1; jj<${component}+1; jj++ ))
    do  COOR=${arrayc[$jj-1]}
	
	# COMPONENT LOOP
	for DAY in `cat $path2mseed/mseed_list`
	do  	    
	    file=$NAME_NET.$NAME_STAT..1$COOR.D.$YEAR.$DAY.????00.SAC
	    echo $file
sac <<EOF
r $file
merge overlap average gap interp 
w $path2sacday/$NAME_NET.$NAME_STAT..1$COOR.D.$YEAR.$DAY.DAY.SAC
quit

EOF
	done
    done
done

	    

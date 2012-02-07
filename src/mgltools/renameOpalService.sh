
SERVERNAME=$1
INPUTFILE="VinaVSpublic_0.1_net.py"
SERVICENAME=`echo $SERVERNAME | sed 's/\./_/g'`

sed "s/AutodockVina_Screening_kryptonite_nbcr_net/autodockvina_screening_1_1_2_$SERVICENAME/g" $INPUTFILE > CorrectHost_net.py
sed -i "s/kryptonite\.nbcr\.net/$SERVERNAME/g" CorrectHost_net.py 




# LC 
# script that can be used to change the server name 
# pointed by the opal services in the Vision Network
# 

#the new rocks server offering the server
SERVERNAME=$1
SERVICENAME=`echo $SERVERNAME | sed 's/\./_/g'`
#basepath of the mgltools
BASEPATH="/opt/mgltools/MGLToolsPckgs"

#cadd workflow
INPUTFILE="$BASEPATH/CADD/workflows/virtualScreening/VinaVSpublic_0.1_net.py"
sed -i -e "s/AutodockVina_Screening_kryptonite_nbcr_net/autodockvina_screening_1_1_2_$SERVICENAME/g" \
       -e "s/kryptonite\.nbcr\.net/$SERVERNAME/g" \
       -e '61,64d' $INPUTFILE 

INPUTFILE="$BASEPATH/CADD/workflows/MDanalysis/GROMOSClustering_0.1_net.py"
sed -i -e "s/GROMOS_ClusterFiles_kryptonite_nbcr_net/gromos_cluster_4_5_5_$SERVICENAME/g" \
       -e "s/kryptonite\.nbcr\.net/$SERVERNAME/g" $INPUTFILE

#modifying the adt macro node
INPUTFILE="$BASEPATH/AutoDockTools/VisionInterface/Adt/Macro/Vina.py"
sed -i -e "s/AutodockVina_Screening_kryptonite_nbcr_net/autodockvina_screening_1_1_2_$SERVICENAME/g" \
       -e "s/kryptonite\.nbcr\.net/$SERVERNAME/g" $INPUTFILE 

INPUTFILE="$BASEPATH/AutoDockTools/VisionInterface/Adt/Macro/GromosCluster.py"
sed -i -e "s/GROMOS_ClusterFiles_kryptonite_nbcr_net/gromos_cluster_4_5_5_$SERVICENAME/g" \
       -e "s/kryptonite\.nbcr\.net/$SERVERNAME/g" $INPUTFILE

INPUTFILE="$BASEPATH/AutoDockTools/VisionInterface/Adt/Macro/PrepareReceptor.py"
sed -i -e "s/Pdb2pqrOpalService_ws_nbcr_net/pdb2pqr_1_8_$SERVICENAME/g" \
       -e "s/PrepareReceptorOpalService_ws_nbcr_net/prepare_receptor_1_5_6_$SERVICENAME/g" \
       -e "s/ws\.nbcr\.net/$SERVERNAME/g" $INPUTFILE


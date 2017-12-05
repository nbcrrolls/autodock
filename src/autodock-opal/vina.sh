#!/bin/bash
#
# run vina in a "vina" SGE queue 

# fixed settings 
VINA=/opt/mgltools/bin/vina
LIBRARYBASEDIR=/share/opal/libraries
SCREEN_ROOT=/opt/cadd/bin

. /etc/profile.d/sge-binaries.sh

echo Inputs arguments: $@

function usage
{
  echo "  ./vina_screening.sh [options] "
  echo 
  echo "      --config file       Used to specify the file with the config option used by vina"
  echo "      --user email        email address to use for sending notificaiton at the end of the simulation"
  echo "      --flex file         flexible part of the receptor"
  echo "      --receptor file     receptor pdbqt file"
  echo "      --ligand file       ligand pdbqt file"
}

while [ "$1" != "" ]; do
    case $1 in
	--config )              shift 
	                        config=$1
                                ;;
	--user )                shift 
	                        email=$1
                                ;;
        --flex )                shift
                                flex=$1
                                ;;
        --ligand )              shift
                                ligand=$1
                                ;;
        --receptor )            shift
                                receptor=$1
                                ;;
        --out )                 shift
                                out=$1
                                ;;
        * )                     pair=`echo $1 $2 | awk -F"--" '{print $2}'`
                                echo $pair | awk '{print $1 " = " $2}'   >> tmp.config
                                shift                                 
    esac
    shift
done

if [ ! "$receptor" ]; then
    usage
    echo
    echo "ERROR: must provide receptor file"
    echo ""
    exit 1;
fi

if [ ! "$ligand" ]; then
    usage
    echo
    echo "ERROR: must provide ligand file"
    echo ""
    exit 1;
fi

if [ "$flex" ]; then
    FLAGS="--flex $flex"
fi

if [ "$out" ]; then
    OUT="--out $out"
fi

# test for config file
if [ "$config"  == "" ] ; then
    config="config"
    touch $config
fi
cp $config $config.org

# add any other form parameters to config
if test -e tmp.config; then
  cat $config >> tmp.config
  cp tmp.config $config
fi

#  -----   Job Submission   -----
s=vina.sub
echo "#!/bin/bash" > $s
echo "#$ -clear" >> $s
echo "#$ -cwd" >> $s
echo "#$ -S /bin/bash" >> $s
echo "#$ -q vina" >> $s
echo "#$ -o vina.out" >> $s
echo "#$ -e vina.err" >> $s
echo "" >> $s
echo $VINA --config $config --receptor $receptor --ligand $ligand --cpu 1 $FLAGS $OUT >& vina.log >> $s

chmod +x $s
qsub -sync y $s &
wait

exit 0

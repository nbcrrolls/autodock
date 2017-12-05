#!/bin/bash

# autodock src files  are available on google drive 
# on nbcr.ucsd@gmail.com # from nbcrrolls/autodock repo

if [ ! -d $NBCRDEVEL ]; then
    echo "Need to install roll 'nbcr' first"
    exit 1
fi

# set google url for dowloading sources. botostrap-values is in updated nbcr roll
# make sure have latest rocks-devel RPM with gget.sh and get_sources.sh updated
#. $NBCRDEVEL/bootstrap-values.sh
export SURL="https://googledrive.com/host/0BwytB2hwlJ6Vbmx3TFNHY2I5OUk"

# download sources
. $ROLLSROOT/etc/bootstrap-functions.sh

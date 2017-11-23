#!/bin/bash

# To set the enviornment, run the following:
# source setenv.sh [-q]

# If -q flag or --quick is set, the workspace and makefiles will not be created.
# If -r flag or --rucio is set, Rucio will be set (you will be promted with a passphrase to enter.  

# Defaults:
QUICK='no'
RUCIO='no'
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -q|--quick)
    QUICK="yes"
    shift # past argument
    #shift # past value - only flags so far
    ;;
    -r|--rucio)
    RUCIO="yes"
    shift # past argument
    #shift # past value - only flags so far
    ;;
	
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo "Setting enviornment with the following parameters:"
echo quick = "${QUICK}"
echo rucio = "${RUCIO}"

# Setup Atlas stuff
setupATLAS
asetup AthAnalysisBase,2.4.30,here
asetup AtlasProduction,20.7.8.4,here


# If not quick - also create the workarea and the makefiles.
# This is not necessary unless dependencies have changed or something of that sort.
if [[ "$QUICK" == "no" ]]; then
	setupWorkArea.py
	cd WorkArea/cmt
	cmt bro cmt config
	cd ../../
fi

# Setup rucio if flag is set
if [[ "$RUCIO" == "yes" ]]; then
	lsetup rucio
	voms-proxy-init -voms atlas
fi

lsetup panda

echo 'Enviornment set up successfully!'


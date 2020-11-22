#!/bin/bash

_usage() {
    echo 'Metis'
    echo '   options: --int64 : use 64 bit integer'
    echo '            --real64  : use 64 bit real'
    echo '            --help'    
    
}

SCRIPT=$(dirname "$0")/env_${TwoPiDevice}.sh
source $SCRIPT

SCRIPT=$(dirname "$0")/subs/ostype.sh
source $SCRIPT

GIT=$(command -v git)
SRCDIR=${TwoPiRoot}/src
CMAKE=$(command -v cmake)
MAKE=$(command -v make)

DO_INT64=0
DO_REAL64=0

while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    --int64)
    DO_INT64=1
    shift # past argument    
    ;;
    --real64)
    DO_REAL64=1
    shift # past argument    
    ;;
    --help)
    _usage
    exit 1
    ;;
    *)
    echo "Unknown option " $key
    exit 2  #  error_code=2
    ;;
esac
done

cd ${SRCDIR}/metis-5.1.0

if [ "$DO_INT64" -eq "1" ]; then
    echo "!!!! 64 bit int mode"
    sed -i 's/#define IDXTYPEWIDTH 32/#define IDXTYPEWIDTH 64/g' include/metis.h
else
    sed -i 's/#define IDXTYPEWIDTH 64/#define IDXTYPEWIDTH 32/g' include/metis.h
fi

if [ "$DO_REAL64" -eq "1" ]; then
    echo "!!!! 64 bit real mode"
    sed -i 's/#define REALTYPEWIDTH 32/#define REALTYPEWIDTH 64/g' include/metis.h
else
    sed -i 's/#define REALTYPEWIDTH 64/#define REALTYPEWIDTH 32/g' include/metis.h
fi

$MAKE config shared=1 prefix=$TwoPiRoot cc=${CC}
$MAKE $MAKEOPT
$MAKE install

OSTYPE=$(ostype)
if [ "${OSTYPE}" == "Mac" ]; then
   install_name_tool -id ${TwoPiRoot}/lib/libmetis.dylib ${TwoPiRoot}/lib/libmetis.dylib
fi




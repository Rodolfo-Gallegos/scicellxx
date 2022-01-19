#! /bin/sh

#====================================================================
# A few helper functions
#====================================================================

# Usage function
usage()
{
cat << EOF
usage: $0 [OPTIONS]

This script stops any container running SciCell++

OPTIONS:
   -h      Show this message
   -v      Verbose
EOF
}

#====================================================================
# Variables (and default values)
#====================================================================
# Verbose option
verbose=TRUE

#====================================================================
# Variables (and default values)
#====================================================================
project_name=scicellxx
container_name=$project_name

#====================================================================
# Parse arguments
# ====================================================================
while getopts “hv” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         v)
             verbose=TRUE
             ;;
         ?)
             echo ""
             echo "============================================================= "
             echo ""
             echo "[ERROR] No recognized parameter"
             echo ""
             echo "============================================================= "
             echo ""
             usage
             exit 1
             ;;
     esac
done

#####################################################################
                   ## STOP THE DOCKER CONTAINER ##
#####################################################################

if ! docker rm --force -v $container_name ; then
    echo ""
    echo "[ERROR] - Removing the docker container"
    echo ""
    exit 1
fi

echo ""
echo "============================================================= "
echo ""
echo "[DONE] - Finishing SciCell++ in docker container"
echo ""
echo "============================================================= "
echo ""

# Success
exit 0

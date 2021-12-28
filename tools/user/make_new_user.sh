#! /bin/sh

#====================================================================
# A few helper functions
#====================================================================

# Usage function
usage()
{
cat << EOF
usage: $0 [OPTIONS]

This script creates the folder structure and files to register a new
user of SciCell++

OPTIONS:
   -h      Show this message
   -u      User name
EOF
}

# An small function 'borrowed' from the oomph-lib installation
# script...
OptionPrompt() 
{ 
 printf "%s " "$1" 
}

# Another small function 'borrowed' from the oomph-lib installation
# script...
OptionRead()
{
 read Opt
 if test "$Opt" = "" ; then
  Opt=$1
 fi
 echo $Opt
}

#====================================================================
# Variables (and default values)
#====================================================================
private_dir=private

# The name of the library
lib_name=SciCell++
# Indicates whether a user name was given, otherwise the script
# prompts for a username
user_name_given=FALSE
user_name=default

#====================================================================
# Parse arguments
# ====================================================================
# A letter followed by a ':' indicates that it requires an argument
# value. Example: 'u:' indicates that 'u' is followed by an argument
# value. 'h' do not require arguments.
while getopts “hu:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         u)
             user_name=$OPTARG
             user_name_given=TRUE
             ;;
         ?)
             echo ""
             echo "============================================================= "
             echo ""
             echo "ERROR: No recognized parameter"
             echo ""
             echo "============================================================= "
             echo ""
             usage
             exit 1
             ;;
     esac
done

#====================================================================
# The building script
#====================================================================

echo ""
echo "============================================================= "
echo ""
echo "This script creates a new user folder for "$lib_name
echo ""
echo "============================================================= "
echo ""
echo ""

#====================================================================
# User name
#====================================================================

# Check whether a user name was given
if test "$user_name_given" = "FALSE" ; then
    
    echo "Type your user name:"
    OptionPrompt "[Do not use spaces or special characters!]"
    user_name=`OptionRead`

    echo ""
    echo "============================================================= "
    echo ""
    
fi

#====================================================================
# Check whether the given user is already in the users list
#====================================================================
if (test -d  $private_dir/$user_name); then
    echo ""
    echo "$user_name"
    echo ""
    echo "============================================================= "
    echo "[ERROR] The user name already exist!"
    echo "   Choose another user name and try again"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

#====================================================================
# Your user name
#====================================================================

echo ""
echo "============================================================= "
echo ""
echo "I am going to create a new folder for you:"
echo ""
echo "$private_dir/$user_name"
echo ""
echo "============================================================= "
echo ""

if ! mkdir -p $private_dir/$user_name ; then
    echo ""
    echo "============================================================= "
    echo "[ERROR] The new user folder could not be created!"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

#====================================================================
# Add the new user folder into the CMakeLists.txt file
#====================================================================

echo ""
echo "============================================================= "
echo ""
echo "Add the new user to the folders tree"
echo ""
echo "============================================================= "
echo ""

echo "ADD_SUBDIRECTORY("$user_name")" >> $private_dir/CMakeLists.txt

# Done
exit 0

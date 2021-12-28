#! /bin/sh

#====================================================================
# A few helper functions
#====================================================================

# Usage function
usage()
{
cat << EOF
usage: $0 [OPTIONS]

This script creates a new project structure for a given SciCell++ user

OPTIONS:
   -h      Show this message
   -u      User name
   -p      Project name
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

# Indicates whether a project name was given, otherwise the script
# prompts for a project name
project_name_given=FALSE
project_name=default

#====================================================================
# Parse arguments
# ====================================================================
# A letter followed by a ':' indicates that it requires an argument
# value. Example: 'u:' indicates that 'u' is followed by an argument
# value. 'h' do not require arguments.
while getopts “hu:p:” OPTION
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
         p)
             project_name=$OPTARG
             project_name_given=TRUE
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
echo "This script creates a new project for a given "$lib_name" user"
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
if ! (test -d  $private_dir/$user_name); then
    echo ""
    echo "$user_name"
    echo ""
    echo "============================================================= "
    echo "[ERROR] The user name does not exist!"
    echo "   Make sure you correctly typed your user name and try again"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

#====================================================================
# Project name
#====================================================================

# Check whether a user name was given
if test "$project_name_given" = "FALSE" ; then
    
    echo "Type your project name:"
    OptionPrompt "[Do not use spaces or special characters!]"
    project_name=`OptionRead`

    echo ""
    echo "============================================================= "
    echo ""
    
fi

#====================================================================
# Check whether the project name already exist
#====================================================================
if (test -d  $private_dir/$user_name/$project_name); then
    echo ""
    echo "$project_name"
    echo ""
    echo "============================================================= "
    echo "[ERROR] The project name already exist!"
    echo "   Choose another project name and try again"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

#====================================================================
# Your user and project name
#====================================================================

echo ""
echo "============================================================= "
echo ""
echo "I am going to create a new folder for user $user_name:"
echo ""
echo "$private_dir/$user_name/$project_name"
echo ""
echo "============================================================= "
echo ""

if ! mkdir -p $private_dir/$user_name/$project_name ; then
    echo ""
    echo "============================================================= "
    echo "[ERROR] The new project folder could not be created!"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

if ! mkdir -p $private_dir/$user_name/$project_name/RESLT ; then
    echo ""
    echo "============================================================= "
    echo "[ERROR] The RESLT folder for the new project could not be created!"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

#====================================================================
# Copy template files into the new project folder
#====================================================================

echo ""
echo "============================================================= "
echo ""
echo "Copy template files into the new project folder"
echo ""
echo "============================================================= "
echo ""

project_code_filename=$user_name_$project_name.cpp
SRC_user_project_tag=SRC_$user_name_$project_name
user_project_name=$user_name_$project_name
LIB_user_project_tag=LIB_$user_name_$project_name

if ! cp ./tools/templates/demo_template.cpp $private_dir/$user_name/$project_name/$project_code_filename ; then
    echo ""
    echo "============================================================= "
    echo "[ERROR] The template code [.cpp] project could not be copied!"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

if ! cp ./tools/templates/CMakeLists.txt.private_template $private_dir/$user_name/$project_name/CMakeLists.txt ; then
    echo ""
    echo "============================================================= "
    echo "[ERROR] The template CMakeLists.txt.private_template file could not be copied!"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

if ! sed -i 's/SRC_demo_john/$SRC_user_project_tag/g' $private_dir/$user_name/$project_name/CMakeLists.txt ; then
    echo ""
    echo "============================================================= "
    echo "[ERROR] The CMakeLists.txt could not be modified!"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

if ! sed -i 's/demo_john.cpp/$project_code_filename/g' $private_dir/$user_name/$project_name/CMakeLists.txt ; then
    echo ""
    echo "============================================================= "
    echo "[ERROR] The CMakeLists.txt could not be modified!"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

if ! sed -i 's/demo_john/$user_project_name/g' $private_dir/$user_name/$project_name/CMakeLists.txt ; then
    echo ""
    echo "============================================================= "
    echo "[ERROR] The CMakeLists.txt could not be modified!"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

if ! sed -i 's/LIB_demo_john/$LIB_user_project_tag/g' $private_dir/$user_name/$project_name/CMakeLists.txt ; then
    echo ""
    echo "============================================================= "
    echo "[ERROR] The CMakeLists.txt could not be modified!"
    echo "============================================================= "
    echo ""
    exit 1 # Error flag
fi

echo ""
echo "============================================================= "
echo "Make sure to include any additional modules to your project"
echo "prior compilation"
echo "============================================================= "
echo ""

echo ""
echo "[DONE]"
echo ""

# Done
exit 0

#! /bin/sh

#--------------------------------------------------------------
# Script to remove licence (all lines starting 
# with //LIC//) from all SciCell++ *.h *.cpp files.
#--------------------------------------------------------------


#--------------------------------------------------------------
# Find *.h and *.cpp files in SciCell++ distribution (exclude
# external_src)
#--------------------------------------------------------------
scicellxx_h_and_cpp_files=`find demos private src \( -name '*.h' -o -name '*.cpp' \) -exec ls  {} \; `


#--------------------------------------------------------------
# Loop over all of these files
#--------------------------------------------------------------
for file in $scicellxx_h_and_cpp_files; do
    echo "Delicencing:" $file
    grep -v "//LIC//" $file > $file.delicenced
    #gawk -f tools/development/licence/remove_licence.awk $file > $file.delicenced
    mv -f $file.delicenced $file
done

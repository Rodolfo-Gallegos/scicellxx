#! /bin/sh

#--------------------------------------------------------------
# Script to apply licence prefix to all SciCell++ *.h and *.cpp
# files.
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
    echo "Licencing:" $file
    cat tools/development/licence/cpp_and_h_licence_block.txt $file > $file.licenced 
    gawk -f tools/development/licence/remove_licence.awk $file.licenced > $file.delicenced
    if (test `diff $file $file.delicenced | wc -w` != 0); then
        echo "Original and licenced/delicenced file don't match!"
        echo "Please run bin/remove_licence.sh first."
        echo "Filename: " $file
        exit 1
    else
        echo "OK"
        mv -f $file.licenced $file
        rm  $file.delicenced
#        svn propset svn:keywords "LastChangedDate LastChangedRevision" $file
    fi
#    echo " "
done

#! /bin/sh

#====================================================================
# A few helper functions
#====================================================================

# A little function 'borrowed' from the oomph-lib installation
# script...
OptionPrompt() 
{ 
 printf "%s " "$1" 
}

# Another little function 'borrowed' from the oomph-lib installation
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
# Variables
#====================================================================
build_dir=build
lib_dir=lib
lib_name=developing_library

#====================================================================
# The building script
#====================================================================

echo " "
echo "============================================================= "
echo "            "$lib_name" installation script" 
echo "============================================================= "
echo " "

echo ""
echo "============================================================= "
echo ""
echo "I am going to build "$lib_name
echo "I will use the ./"$build_dir" directory as the temp build folder"
echo "The complete built library will be stored in the ./"$lib_dir" folder"
echo ""
echo "============================================================= "
echo ""

#====================================================================
# Checking if lib_dir exists
#====================================================================
if (test -d $lib_dir); then 
	continue
else
	mkdir $lib_dir
fi

#====================================================================
# Going to the build directory
#====================================================================
if (test -d  $build_dir); then 
	cd $build_dir
	echo "Cleaning up ..."
	rm -r *
	echo "Done"
else
	mkdir $build_dir
	cd $build_dir
fi
echo ""
echo ""

#====================================================================
# Library version
#====================================================================

echo "Which library type do you want to build?"
OptionPrompt "The a)STATIC or the b)SHARED type library? [a/b -- default: a]"
static_or_shared=`OptionRead`
    if test "$static_or_shared" = "b" -o "$static_or_shared" = "B" ; then 
		lib_type=SHARED
		lib_ext=.so
	else
		lib_type=STATIC
		lib_ext=.a
	fi
	
echo ""
echo "============================================================= "
echo ""
	
echo "Which library version do you want to build?"
OptionPrompt "The a)DEBUG or the b)RELEASE versionlibrary ? [a/b -- default: a]"
debug_or_release=`OptionRead`
    if test "$debug_or_release" = "b" -o "$debug_or_release" = "B" ; then 
		lib_version=RELEASE
	else
		lib_version=DEBUG
	fi
	
echo ""
echo ""
echo "Building the " $lib_type "/" $lib_version " version of the library ..."

#====================================================================
# Calling CMake
#====================================================================
# Go one folder up since we did a cd into ./build
cmake ../ -Dlib_type=$lib_type -DCMAKE_BUILD_TYPE=$lib_version
make clean
make
cp lib$lib_name$lib_ext ../$lib_dir
echo "Done"

#====================================================================
# Finishing up !!!
#====================================================================

cd ..

echo ""
echo "============================================================= "
echo ""
echo "Finishing library built process (autogen.sh has finished!)"
echo "If you can't spot any error messages above this, the" 
echo "multirate library should now be ready to use... " 
echo " "
echo "Please contact the developers if you encountered any"
echo "building problem!"
echo ""
echo "============================================================= "
echo ""

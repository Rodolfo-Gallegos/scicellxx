# Add directories with source files (in order of dependence, the ones
# with no dependency first then the others)
ADD_SUBDIRECTORY(numerical_recipes)
ADD_SUBDIRECTORY(argparse)
# Once all subdirectories have compiled their stuff then create a VARIABLE to hold their names
SET(LIB_external_src_lib numerical_recipes_lib argparse_lib PARENT_SCOPE)

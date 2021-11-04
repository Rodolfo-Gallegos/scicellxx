Add your project to the ``demos`` folder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you add a new feature to SciCell++ we encourage you to
:doc:`create_a_tutorial` and a demo showing these new features. Here
we detail the process to include your project as part of the demos of
SciCell++. We divide this process in two parts, the first one guides
you to create your folder and your validation files, the second part
shows you how to configure the SciCell++ to build and execute your
demo. In both sections we suppose that your demo is called
``demo_sophy``.

**Create your demo and validation folder for your demo**

The initial steps to include your demo as part of SciCell++ involve
create a folder in the SciCell++ demos folder structure and to
generate the validation files.

1. Execute your project and save its output into a file. We encorage
   you to execute it using single and double precision so that we have
   two different outputs. The files that you generate should be named:
   
   * ``validate_demo_sophy.dat`` for the single precision generated
     output.
   * ``validate_double_demo_sophy.dat`` for the double precision
     generated output.

2. Create a new folder into the ``demos`` folder structure. Use a name
   that captures the intent of your project.

   .. code-block:: shell

      mkdir <your-folder-name>

3. Add the following line at the end of the ``CMakeLists.txt`` file
   that lives at the same level of the folder that you created:
   
   .. code-block:: shell
      
      ADD_SUBDIRECTORY(your-folder-name)

4. Step into your demo folder and create a folder called
   ``validate``.

5. Copy the two output files (or copy all of them if you have more
   than two) generated at step 1 into the ``validate`` folder.

**Configure SciCell++ to build and execute your demo**

Once you have created your folder and copied the validation files
there you are ready to configure SciCell++ to build and execute your
demo.

1. Copy the source code for your project into your demo folder, in
   this case we suppose that the source code for your project is
   the file ``demo_sophy.cpp``.

2. Copy the ``CMakeLists.txt.demo_template`` from the ``/tools/``
   folder into your demo folder. Rename this file as
   ``CMakeLists.txt``.

3. Change the content of the ``CMakeLists.txt`` file as follow:
   
   * Change all the instances of the tag ``SRC_demo_john`` for your
     own tag to identify your source code. For example:
     ``SRC_demo_sophy``.

   * Change all the instances of ``demo_john.cpp`` for the name of
     your source code file. For example: ``demo_sophy.cpp``.

   * Change all the instances of ``demo_john``, this will be the name
     of your executable and the name you need to type at the terminal
     to compile your project. For example:``demo_sophy``.
     
   * Change all the instances of the tag ``LIB_demo_john`` for your
     own tag to identify libraries required for your code. For
     example: ``LIB_demo_sophy``.

   * Include the modules you need. In the template we only include the
     ``general_lib`` and the ``problem_lib`` modules. Check the
     :doc:`modules` document for the full list of module and their
     details.
    
4. In the same file perform the following changes in the ``Test
   section``.
   
   * Change all the instances of ``TEST_demo_john_run`` by the name of
     your demo. For example: ``TEST_demo_sophy_run``.

     .. important:: Make sure to keep the ``TEST`` and ``_run`` prefix
                    and postfix, respectively.
  
   * Change all the instances of ``demo_john`` with the name of your
     demo. For example: ``demo_sophy``.

   * Change all the instances of ``VALIDATE_FILENAME_demo_john`` with
     the name of your tag for the validation file. For example:
     ``VALIDATE_FILENAME_demo_sophy``.

   * Change the name of the validation file
     ``validate_double_demo_john.dat`` by yours. Recall that this file
     should store the output of your project executed using double
     precision. For example: ``validate_double_demo_sophy.dat``.

   * Change the name of the validation file ``validate_demo_john.dat``
     by yours. Recall that this file should store the output of your
     project executed using single double precision. For example:
     ``validate_demo_sophy.dat``.
  
   * Change all instances of ``TEST_demo_john_check_output`` with the
     name of your demo. For example: ``TEST_demo_sophy_check_output``.

   .. important:: Make sure to keep the ``TEST`` and ``_output``
                  prefix and postfix, respectively.

5. Make sure that the computations of your demo are stored in an
   output file. If the file that you generate is called differently
   than ``output_test.dat`` then modify any instance of that name in
   the ``CMakeLists.txt`` file.

6. Go to the root folder of SciCell++ and execute the ``./autogen.sh``
   script and enable the execution of the demos. If you find errors
   please make sure you correctly changed all the tags indicated in
   the previous steps. Your project should be automatically built,
   executed and validated.


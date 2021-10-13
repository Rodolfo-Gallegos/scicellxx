A quick starting-up guide
=========================

In this section we provide you with the basic tool to start working
with SciCell++, we present the main folder structure and how to create
your very first project. There is also a section for including your
project as a demo for further reference of new collaborators.

Running demos
-------------

SciCell++ is released with a set of demos that show you some of its
main features. We recommend you to explore the demos section of the
documentation and the demos folder. Whenever you want to run a demo
just go to the demo folder which you are interested, create a folder
called ``RESLT`` if it is not already there and type ``./bin/``
followed by the name of the demo.

* **Example:** Lets say you want to run the Lotka-Volterra demo in the
  folder ``/demos/lotka_volterra/``, once you are in that folder
  create the ``RESLT`` folder where the output is stored (all the
  demos are configured to store its output in a folder with that name,
  if the folder does not exist then the output is not generated) and
  run the demo:

  .. code-block:: shell

     mkdir RESLT
     ./bin/demo_lotka_volterra

  Once the demo has started you should see output messages on the
  terminal with general information about the results of the
  computations. You can check the produced results in the ``RESLT``
  folder.

.. note:: Observe that some demos are equipped with Python or GNUPlot
          script to visualise the results. Try to run them as ``python
          <name-of-the-python-script.py>`` or ``gnuplot
          <name-of-the-gnu-script.gp>``.

Input arguments
^^^^^^^^^^^^^^^

Some demos require input arguments to run, if you try to run one of
those and pass nothing you will get a message indicating what you need
to pass. You can also check what input arguments a demo needs by
passing the ``--help`` or ``-h`` options at running time.

Create your ``private`` folder
------------------------------

Every user has its own private folder, use this folder to store all of
your work, in-development demos and any of your new developed features
for SciCell++. One of the first things that you should do in order to
start developing new features for SciCell++ is to create your private
folder, to do so follow theses instructions:

1. Open a terminal and go to the ``private`` folder of SciCell++ and
   typet the following (make sure to substitute ``john_cool`` by
   your name):

   .. code-block:: shell

      cd private
      mkdir john_cool
      cd john_cool

2. Update the ``CMakeLists.txt`` file in the private folder by adding
   your folder name at the end of the file as follow (make sure to
   substitute ``john_cool`` by your name):

   .. code-block:: shell

      ADD_SUBDIRECTORY(john_cool)

3. Run the ``autogen.sh`` script at the root folder of SciCell++ and
   make sure no problems are found. If there are any problem
   double-check that you added your folder inside the ``private``
   folder of SciCell++ and that you are modifying the correct
   ``CMakeLists.txt`` file.

Creating your own project
-------------------------

The easiest way to start a new project is to use a demo as a
template. For this example we are going to copy the demo driver
``demo_basic_interpolation.cpp`` from the folder
``demos/interpolation/basic_interpolation``.

1. Open a terminal and go to your private folder.

2. Type the following to copy the demo driver into your private folder:

   .. code-block:: shell

      cp ../../demos/interpolation/basic_interpolation/demo_basic_interpolation.cpp demo_john.cpp

3. Copy the ``CMakeLists.txt.private_template`` file from the
``tools`` folder into your private directory and change its name to
``CMakeLists.txt``

   .. code-block:: shell

      cp ../../tools/CMakeLists.txt.private_template CMakeLists.txt

4. Change the content of the ``CMakeLists.txt`` file as follow:

  * Change all the instances of the tag ``SRC_demo_john`` for your own
    tag to identify your source code. For example: ``SRC_project_sophy``.

  * Change all the instances of ``demo_john.cpp`` for the name of your
    source code file. For example: ``project_sophy.cpp``.

  * Change all the instances of ``demo_john``, this will be the name
    of your executable and the name you need to type at the terminal
    to compile your project. For example:``project_sophy``.
    
  * Change all the instances of the tag ``LIB_demo_john`` for your own
    tag to identify libraries required for your code. For example:
    ``LIB_project_sophy``.

  * Include the modules you need. In the template we only include the
    ``general_lib`` and the ``problem_lib`` modules. Check the
    :doc:`modules` document for the full list of module and their
    details.

5. Go to the root folder of SciCell++ and execute the ``./autogen.sh``
   script. If you find errors please make sure you correctly changed
   all the tags indicated in the previous step. Once building has
   finished without errors you can build your own project.

Building and executing your project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Open a terminal and follow these instructions:

1. Go to the ``build`` folder in the root SciCell++ folder and type

   .. code-block:: shell
   
      make demo_sophy
      
   The building output should be displayed at your screen. Once no
   errors have been reported you may run your code.

2. Go to your ``private`` folder, create a ``RESLT`` folder if you
   have no one, and type:

   .. code-block:: shell

      ./bin/demo_sophy
                   
3. You should see the output of your project at the terminal.

.. important:: As you noticed, the generation and execution of your
               project is performed in two different folders:

               * the ``build`` folder (building)
               * your ``private`` folder (execution)

               We use this two-folders strategy to avoid cluttering
               the folder structure of SciCell++ with files
               automatically generated by CMake. By following this
               strategy we keep a clean folder structure for SciCell++
               and group all files generated by CMake in the ``build``
               folder. This help us to keep track for changes easily
               since we can exclude the whole ``build`` folder from
               the git repository.

               **Just keep in mind the following:**

               * Whenever you want to build your project you need to do so in the ``build`` folder, inthere just type ``make`` followed by the name of your project.

               * Whenever you want to execute your project go to your ``private`` folder and type ``./bin/the-name-of-your-project``.
            
Add your project to the ``demos`` folder
----------------------------------------

If you add a new feature to SciCell++ we encourage you to
:doc:`create_a_tutorial` and a demo showing these new features. Here
we detail the process to include your project as part of the demos of
SciCell++. We divide this process in two parts, the first one guides
you to create your folder and your validation files, the second part
shows you how to configure the SciCell++ to build and execute your
demo. In both sections we suppose that your demo is called
``demo_sophy``.

Create your demo and validation folder for your demo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Configure SciCell++ to build and execute your demo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

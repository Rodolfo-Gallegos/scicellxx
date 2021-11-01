Create your ``private`` folder
==============================

Every SciCell++ user has its own private folder. You should use this
folder to store all of your work, this should include in-development
demos, and any non-released features you are developing for
SciCell++.

.. note:: This workflow should be executed only once by every new
          SciCell++ user, however, if you are collaborating with other
          SciCell++ users then you may need a private folder for your
          collaborative work.

If you require to move your code into another folder then have a look
at the :ref:`folder structure
<additional-features-folder-structure_initial_steps.rst>` section of
SciCell++ to find the folder that better fits your needs.


Workflow
--------

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

3. Run the ``autogen.sh`` script at the main folder of SciCell++ and
   make sure no problems are found. If there are any problem
   double-check that you added your folder into the ``private`` folder
   of SciCell++ and that you modified the correct ``CMakeLists.txt``
   file.

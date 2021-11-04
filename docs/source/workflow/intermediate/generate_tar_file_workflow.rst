Generate a ``.tar.gz`` file to distribute SciCell++
===================================================

The easiest way to distributed SciCell++ is by means of the official
GitHub repository, however, if you need to move your current copy of
SciCell++ to a computer with no Internet access (ex. an isolated
cluster of computers or a SuperComputer) this is an easy way to do
so. Follow the steps in this section to create a ``.tar.gz`` package
file with your current version of SciCell++.

Workflow
--------

1. Save all of your work, including source files, images and
   documentation.

2. Make sure that your current version has neither errors nor broken
   demos. Verify this by running the ``./autogen.sh`` script at the
   main folder of SciCell++.

3. At the main folder of SciCell++ type the following:

   .. code-block:: shell

                   ./make_clean_distro.sh

   The full folder containing SciCell++ will be copied into a
   temporary location, all the control version information generated
   by Git will be removed. You will be prompted to remove all files
   with the extension ``.dat, .png, .tar.gz, .fig, .bin, .rar, .vtu,
   .ubx, .gp, .m`` (only those in the ``demos`` folder will be
   keep). The process of creating a compressed file will start.

4. Once finished a file named ``SciCell++.tar.gz`` will be created at
   the main folder of SciCell++.


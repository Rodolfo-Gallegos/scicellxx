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

                   ./tools/user/make_clean_distro.sh


   This script will perform a set of instructions to generate a
   ``.tar.gz`` file package. You will be prompted to whether remove or
   not all files with extensions ``.dat, .png, .jpeg, .jpg, .tar.gz,
   .fig, .bin, .rar, .vtu, .ubx, .gp, .m`` (only those in the
   ``demos`` and ``private`` folders will be ignored for
   deletion). The process of creating a compressed file will start.

4. Once finished a file named ``SciCell++.tar.gz`` will be created at
   the main folder of SciCell++.


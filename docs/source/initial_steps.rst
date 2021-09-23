Initial steps
=============

This document shows you how to :ref:`install
<installation-label_initial_steps.rst>` and :ref:`configure
<configuration-label_initial_steps.rst>` SciCell++ on a Windows or
Linux system. The default installation strategy is based on containers
so that all software dependencies are preinstalled and ready to use by
SciCell++. The :ref:`configuration section
<configuration-label_initial_steps.rst>` shows you how to compile and
enable special features for SciCell++. :ref:`Additional features
<additional-features-label_initial_steps.rst>` are included to ease
your journey with SciCell++.

.. _installation-label_initial_steps.rst:

Installation
------------

The default installations use a container strategy, thus that all
software dependencies are preinstalled and ready to use for SciCell++.

:doc:`installation/windows_installation`
     For Windows systems users.

:doc:`installation/linux_installation`
     For Linux systems users. We also provide instructions for a
     non-container based installation in case you have enough spare
     time.

     :ref:`Advanced installation
     <advanced_installation-label_linux_installation.rst>`: Use this
     installation only if you are familiar with Unix based systems,
     within this installation you have full control over the versions
     of the third-part packages used by SciCell++.
         
.. _configuration-label_initial_steps.rst:

Configuration
-------------

This section guides you through the configuration process of
SciCell++. The configuration is performed with help of the
``autogen.sh`` script which lives in the main SciCell++ folder.

1. Open a terminal and go to the ``scicellxx`` folder.
2. Execute the automatic generator script by typing:

   .. code-block:: shell

                   ./autogen.sh

   .. important::

      This commands executes a full compilation of SciCell++ and runs
      all the demos and tests to make sure you are working with an
      stable copy. If you want a full list of available parameters for
      this script then add the ``-h`` parameter or review the
      :ref:`options section for autogen.sh
      <autogen.sh-options-label_initial_steps.rst>`.

   .. important::

      If you are using a docker container to run SciCell++ then do not
      forget to pass the ``-c ./configs/container`` option to the
      ``autogen.sh`` script.
      
   A summary of the compilation and testing process is shown once they
   have finished. If no errors were reported then SciCell++ is ready
   to go. We recommend you to have a look at the :doc:`tutorials` and
   :doc:`demos` as follow up.

.. _autogen.sh-options-label_initial_steps.rst:
        
Options for the ``autogen.sh`` script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to fully configure the compilation process use the ``-i``
parameter. Check the options for the ``autogen.sh`` section.
sectionsYou will be able to specify the number of processors to
compile SciCell++, also the number of processor to run the demos, use
predefined configuration files for access to third-party libraries and
many more. For a full list of available options use the ``-h``
parameter.

.. _additional-features-label_initial_steps.rst:

Additional features
-------------------

In this section we present some additional features that may help you
to generate the full documentation of SciCell++ from source code, and
to move SciCell++ to a computer with no Internet access.

Generate ``doxygen`` documentation for SciCell++
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This allows you to create class diagrams and browseable documentation
directly from the source code of SciCell++.

**Requirements**

* `Doxygen <https://www.doxygen.nl/index.html>`_ and `Latex
  <https://www.latex-project.org/>`_ to generate documentation from
  source code.

  Check :ref:`this section <doxygen-installation-label_initial_steps.rst>` for doxygen installation.
  
**Steps**
  
1. Open a command line and go to the upper level folder of the
   project, probably called ``scicellxx``.

2. In the command line type the following:
  
   .. code-block:: shell

                   ./make_doc.sh

   Voila! The documentation will be automatically generated into the
   ``docs/doxy_doc/html`` folder.

3. Open the file ``index.html`` within your favorite web-browser to
   read the documentation.

Generate a ``.tar.gz`` file to distribute SciCell++
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to distributed SciCell++ is by means of the official
GitHub repository, however, if you need to move your current copy of
SciCell++ to a computer with no Internet access (ex. an isolated
cluster of computers or a SuperComputer) this is an easy way to do
so. Follow the steps in this section to create a ``.tar.gz`` package
file with your current version of SciCell++.

**Requirements**

* Save all of your work
* Make sure that your current version has neither errors nor broken
  demos. You can verify this by running the ``./autogen.sh`` script at
  the root directory of SciCell++.

**Steps**

1. Go to the upper level folder of the project, probably called
   ``scicellxx``.

2. Open a command line and type

   .. code-block:: shell

                   ./make_clean_distro.sh

   The full folder containing SciCell++ will be copied into a
   temporary location, all the control version information generated
   by Git will be removed. You will be prompted to remove all files
   with the extension ``.dat, .png, .tar.gz, .fig, .bin, .rar, .vtu,
   .ubx, .gp, .m`` (only those in the ``demos`` folder will be
   keep). The process of creating a compressed file will start.

3. Once finished a file named ``SciCell++.tar.gz`` will be created in
   the root folder of SciCell++.

Workflow
--------

The main differences on the workflow for Windows and Linux users are
on the graphic interfaces. We provide details only for the graphic
interfaces for those steps that may be required.

:doc:`workflow/windows_workflow`
     For Windows systems users.

:doc:`workflow/linux_workflow`
     For Linux systems users.

   

Create your ``private`` folder
==============================

Every SciCell++ user has its own private folder. You should use this
folder to store all of your work, this should include in-development
demos, and any non-released features you are developing for
SciCell++.

.. note:: This workflow should be executed only once by every new
          SciCell++ user. However, if you are collaborating with other
          SciCell++ users or you require an additional private folder
          then you may need to execute this workflow again.

.. _folder-structure_create_your_private_folder_workflow.rst:
   
Folder structure
----------------

Prior to create your own private folder we encourage you to explore
the folder structure of SciCell++. In this section we briefly mention
what each folder is about:
    
* ``build``, this folder is automatically generated when compiling
  SciCell++, all compilation files are stored inhere. You do not need
  to deal with the files within this folder, just leave them alone.
  
* ``configs``, store configuration files where each file corresponds
  to an specialised configuration of the framework. For example, you
  can indicate to use ``Armadillo``, ``VTK``, double precision
  arithmetic, panic mode, etc. Have a look at the :ref:`options for
  configuration files
  <options_for_the_configuration_file-label_initial_steps.rst>`. The
  current configuraton is stored in the ``current`` file. If you want
  to use an specialised configuration of the framework you should copy
  the configuration into the ``current`` file or choose it as the
  `configuration file` when running the ``autogen.sh`` script (use the
  ``-c`` option followed by the configuration filename). Try any other
  of the configurations in this folder to reduce the compilation time
  or to improve the performance of SciCell++.
  
* ``demos``, stores a large set of demos that you may want to use as
  templates or starting points for your project. These demos provide a
  good insight on the features available in SciCell++. This folder
  also helps for testing the framework and reporting any issues found
  when new features are implemented. Whenever you want to update your
  contributions to the SciCell++ repository make sure all of the demos
  compile, run and pass the tests. Check the corresponding
  :doc:`demos` documentation for details.

* ``docs``, the source files for this documentation. If you add a demo
  to SciCell++ you will be requested to write the documentation to
  your demo within this folder.
       
* ``external_src``, this folder stores any external software packages
  used within the framework to provide extra features. You should not
  modify this folder unless you are providing new functionalities that
  depend on external software packages. If you are using a container
  to run SciCell++ then most of the software within this folder is not
  used.
  
* ``private``, stores private files for each user or
  collaborator. Each one should have its own private folder here, this
  should be used as the development folder for each one. We encorage
  you to fully document your projects so that it can later be included
  in the ``demos`` folder to shown any specialised features of the
  framework that you contributed with.
  
* ``src``, this is SciCell++'s soul, here lives all the source
  code. Prior to including files in this folder you should test them
  in your ``private`` folder. Any addition into this folders requires
  the aprovement of the main developers team since you would be mainly
  extending SciCell++'s capabilities.
  
* ``tools``, a set of tools used for the library, stores scripts used
  by the framework at compilation time or to generate *clean
  distributions* of the framework.

Once you have a glance of what each folder is about on SciCell++
organization lets create your own private folder.

Workflow
--------

1. Open a terminal and on the main folder of SciCell++ execute the
   following script:

   .. code-block:: shell

      ./tools/user/make_new_user

   The script will prompt you for a user name, this should not include
   whitespaces or any special character.

2. Run the ``autogen.sh`` script at the main folder of SciCell++ and
   check for any errors.

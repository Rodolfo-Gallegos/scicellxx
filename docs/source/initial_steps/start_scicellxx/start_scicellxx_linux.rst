Starting SciCell++ on Linux
===========================

Overview
--------

1. :ref:`Run Docker container with SciCell++
   <run-docker-container-label_start_scicellxx_linux.rst>`
2. :ref:`Troubleshooting
   <troubleshooting-label_start_scicellxx_linux.rst>`
   
.. _run-docker-container-label_start_scicellxx_linux.rst:
   
Run Docker container with SciCell++
-----------------------------------

1. Open a terminal and go to the folder where you ran the `git clone`
   command, then type the following:

   .. code-block:: shell

      sudo ./scicellxx/tools/common/run_scicellxx_on_docker.sh

      
 If you spot no errors then follow the instructions on section
 :ref:`Configuration <configuration-label_initial_steps.rst>`.

.. _troubleshooting-label_start_scicellxx_linux.rst:
   
Troubleshooting
---------------

I am getting an error when executing the script run_scicellxx_on_docker.sh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get an error when executing the ``run_scicellxx_on_docker.sh``
script stating that the name ``/scicellxx`` is already in use by
another container then you need to ``DELETE`` that container prior to
start another one with the same name.

 - Open a terminal and type the folowing:

   .. code-block:: shell

      sudo ./scicellxx/tools/common/stop_scicellxx_on_docker.sh
                 

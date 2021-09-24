Starting SciCell++ on Windows
=============================

Overview
--------

1. :ref:`Run Docker container with SciCell++
   <run-docker-container-label_start_scicellxx_windows.rst>`
2. :ref:`Troubleshooting
   <troubleshooting-label_start_scicellxx_windows.rst>`

.. _run-docker-container-label_start_scicellxx_windows.rst:
   
Run Docker container with SciCell++
-----------------------------------
   
1. Run the Docker Desktop application. If you installed it with the
   default options then it should be already running on the
   backgroud. Open the interface by double clicking the docker icon at
   the botton-right menu of your task bar and check no errors are
   reported.

2. Open a ``Windows PowerShell`` (there is no need to do so with
   administrative rights) and type the following.

   .. code-block:: shell
                   
                   docker run --name=scicellxx -v C:\Users\tachi\Documents\GitHub\scicellxx:/home/scicellxx -w /home/scicellxx/ -it scicellxx/scicellxx-base-all:0.1
   
   .. warning:: 
   
      Make sure to change
      ``C:\Users\tachi\Documents\GitHub\scicellxx`` by the path where
      you cloned the SciCell++ repository in your local machine.

   You should have a similar output as that shown in the image. Wait
   for completion.

   .. image:: figures/01.png
      :width: 700
   
3. Once finished, you should have a prompt as that shown in the
   image. That means SciCell++ is ready to run.

   .. image:: figures/02.png
      :width: 700

   You could also check the docker interface that should show a
   running image with the name ``scicellxx`` as shown below:

   .. image:: figures/03.png
      :width: 600

4. Continue with the :ref:`configuration step
   <configuration-label_initial_steps.rst>` at the initial steps
   document.

.. _troubleshooting-label_start_scicellxx_windows.rst:
   
Troubleshooting
---------------

I am getting an error when running the docker run command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get an error when running the ``docker run`` command stating
that the name ``scicellxx`` is already in use by another container you
need to ``DELETE`` the container in your docker interface. Open the
docker interface and in the ``Containers/Apps`` section find the
``scicellxx`` container and click on the ``Trash can`` icon to delete
it.

   .. image:: figures/troubleshoot/01.png
      :width: 600
                 

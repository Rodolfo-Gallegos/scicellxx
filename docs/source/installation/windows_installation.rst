Windows installation
====================

This section provides instructions for the installation of SciCell++
on a Windows system. We tested these instructions on Windows 10 but we
(hopefully) expect them to work on recent versions too.

Overview
--------

1. :ref:`Enable virtualisation on Windows
   <enable-virtualisation-label_windows_installation.rst>`
2. :ref:`Install Docker Desktop
   <install-docker-desktop-label_windows_installation.rst>`
3. :ref:`Install GitHub Desktop
   <install-github-desktop-label_windows_installation.rst>`
4. :ref:`Test your installation
   <test-installation-label_windows_installation.rst>`

.. _enable-virtualisation-label_windows_installation.rst:
   
Enable virtualisation on Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following instructions are based on this `YouTube video
<https://youtu.be/6cVBG9BHibo>`_ and the `official webpage
<https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_ for
WSL2 installation for Windows 10.

1. Look for ``bios`` on the windows search tool and select the
   ``Change advanced startup options`` option.
   
   .. image:: figures_windows_installation/01_virtualisation/01.png
      :width: 300

2. On the Advanced startup section click on the ``Restart now`` button.

   .. image:: figures_windows_installation/01_virtualisation/02.png
      :width: 300

3. Click on ``Troubleshoot``, then on ``Advanced options``, followed
   by ``UEFI Firmware Settings`` and finally click on the ``Restart``
   button.

4. Once your computer has launched you may see a screen similar to the
   one below (the specific screen depends on your vendor's
   machine). Look for an ``Advanced Menu`` and make sure that
   ``Virtualization Technologies`` are enabled. **Save your changes**
   and restart your computer. You may have a similar option if you are
   using a different processor brand.

   .. image:: figures_windows_installation/01_virtualisation/03.jpg
      :width: 500

5. Once your computer has restarted look for ``powershell`` on the
   search bar, right-click on the ``Windows PowerShell`` option and
   select ``Run as administrator``.

   .. image:: figures_windows_installation/01_virtualisation/04.png
      :width: 300
              
6. On the command line type (or copy-paste) the following and wait for completion.

   .. code-block:: shell

                   dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart

7. Check your Windows version by typing the ``winver`` command in the
   ``Run`` dialog, press ``Windows Key+R`` to open the ``Run`` dialog.

   .. image:: figures_windows_installation/01_virtualisation/05.png
      :width: 400

8. In the ``About Windows`` dialog check you fullfill the following
   requirements (as indicated in `Step 2 on this webpage
   <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_):

   * For x64 systems: Version 1903 or higher, with Build 18362 or higher.
   * For ARM64 systems: Version 2004 or higher, with Build 19041 or higher.

   If that is not the case then use the ``Windows Update settings``
   tool to install system updates.

   .. image:: figures_windows_installation/01_virtualisation/06.png
      :width: 300

9. Once again open a ``Windows PowerShell`` with administrative
   rights, type (or copy-paste) the following and wait for completion.

   .. code-block:: shell
                   
                   dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart

10. Download and install the ``WSL2 Linux kernel update package for
    x64 machines`` as indicated on `Step 4 on this page
    <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_

   .. image:: figures_windows_installation/01_virtualisation/07.png
      :width: 400

11. Once more open a ``Windows PowerShell`` with administrative
    rights, type (or copy-paste) the following and wait for
    completion.

   .. code-block:: shell
                   
                   wsl --set-default-version 2

12. Install a Linux distribution as indicated on `Step 6 on this page
    <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_. We
    recommend to install the latest Ubuntu available distribution
    (20.04 LTS or 18.04 LTS). **Do not forget to launch and set a
    password for your newly installed linux distribution.**

    .. image:: figures_windows_installation/01_virtualisation/08.png
      :width: 400

.. _install-docker-desktop-label_windows_installation.rst:
   
Install Docker Desktop
^^^^^^^^^^^^^^^^^^^^^^

1. Download `Docker Desktop
   <https://www.docker.com/products/docker-desktop>`_ for windows (at
   the writing of this document lastest version was 3.5.2).

2. Install Docker Desktop with the default options.

   .. image:: figures_windows_installation/02_docker_installation/01.png
      :width: 400
   
3. Once the installation process finish you need to restart your
   computer. Click on the ``Close and restart`` button.

4. (Optional) Open docker, go to ``Settings>General`` and make sure
   the ``Use the WSL2 based engine`` check box is ticked.
   
.. _install-github-desktop-label_windows_installation.rst:
   
Install GitHub Desktop
^^^^^^^^^^^^^^^^^^^^^^

1. Download `GitHub Desktop <https://desktop.github.com/>`_ (you will
   need lo sign up on `GitHub <https://github.com/>`_).

2. Install GitHub Desktop and select the ``Sign in to GitHub.com``
   option.

3. In the browser use your GitHub credentials to login. If prompted,
   select the ``open on GitHub desktop`` option.

4. On the ``Configure Git`` dialog select the ``Use my GitHub account
   name and email address`` option and click on ``Finish``.
         
5. Select the ``Clone a repository from the Internet...`` option.

6. Look for the ``scicellxx`` repository and select it. Use the
   default location to clone the repository or choose one in your
   local drive (make sure to remember this location since you will
   need it to use SciCell++).

7. Click on the ``Clone`` button and wait for completion.

8. Create a new branch on the Github Desktop application. Go to the
   menu ``Branch`` and select ``New branch...``. This will open a
   dialog where you specify the new branch name, use your name in
   lowercase as the branch name. For example `john_cool`.

   .. image:: figures_windows_installation/03_github_desktop/01.png
      :width: 400
 
   .. note::

         Whenever you start to work with SciCell++ you should ensure
         that you are working on your own branch. In case you are on a
         different branch you can switch to your branch (or any other)
         by selecting it on the popup menu (`current branch`).

   .. note::
      
         Any commits to SciCell++ must be done to your own branch, so
         make sure the ``Commit to ..`` button spells your branch
         name.
         
.. _test-installation-label_windows_installation.rst:
   
Starting SciCell++
^^^^^^^^^^^^^^^^^^

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

   .. image:: figures_windows_installation/04_start_scicellxx/01.png
      :width: 700
   
3. Once finished, you should have a prompt as that shown in the
   image. That means SciCell++ is ready to run.

   .. image:: figures_windows_installation/04_start_scicellxx/02.png
      :width: 700

   You could also check the docker interface that should show a
   running image with the name ``scicellxx`` as shown below:

   .. image:: figures_windows_installation/04_start_scicellxx/03.png
      :width: 600

4. Continue with the :ref:`configuration step
   <configuration-label_initial_steps.rst>` at the initial steps
   document.

**What if I am getting an error when running the ``docker`` command**

    
.. warning:: Not finished section.
             

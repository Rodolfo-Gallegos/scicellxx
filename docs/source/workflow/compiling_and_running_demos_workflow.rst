Compiling and running demos
===========================

SciCell++ is released with a set of demos that show you some of its
main features. We recommend you to explore the demos section of the
documentation and the demos folder to find a demo that be of your
interest. Here we show you how to run a basic demo, however, these
workflow should work for most of the SciCell++ demos. Carefully review
the documentation associated for the specific demo you are interested
for any additional requirements.

.. note:: Whenever you want to run a demo you need to create a
          ``RESLT`` folder, this is where the demo output will be
          located. If the folder does not existe when you run the demo
          then no output will be generated. If the folder does exists
          then delete or move its content prior to running the demo to
          avoid overwritting.

Workflow          
--------

Suppose you want to run the Lotka-Volterra demo in the folder
``/demos/odes/lotka_volterra/``.

Compiling
^^^^^^^^^

Running
^^^^^^^

You need to move into that folder, do
it as follow:

  .. code-block:: shell

     cd demos/odes/lotka_volterra

once in the folder create the ``RESLT`` folder to store the output of
the demo.

  .. code-block:: shell

     mkdir RESLT

Run the demo by typping its name after the ``./bin/`` string as follow:

  .. code-block:: shell

     ./bin/demo_lotka_volterra

You should see output messages on the terminal with general
information about the results of the computations. Once finished check
the results in the ``RESLT`` folder.

.. note:: Some demos are equipped with ``Python`` or ``GNUPlot``
          script to visualise the results. Try to run them as ``python
          <name-of-the-python-script.py>`` or ``gnuplot
          <name-of-the-gnu-script.gp>``.

Input arguments
---------------

Some demos require input arguments to run, if you try to run one of
those and pass nothing you will get a message with the list of
arguments that you need to pass. You can also check what input
arguments a demo needs by passing the ``--help`` or ``-h`` options
when executing the demo. Example:

  .. code-block:: shell

     ./bin/demo_lotka_volterra --help

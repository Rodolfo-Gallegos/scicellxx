Running demos
=============

SciCell++ is released with a set of demos that show you some of its
main features. We recommend you to explore the demos section of the
documentation and the demos folder.

Whenever you want to run a demo just go to the demo folder which you
are interested, create a folder called ``RESLT`` if it is not already
there and type ``./bin/`` followed by the name of the demo.

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
---------------

Some demos require input arguments to run, if you try to run one of
those and pass nothing you will get a message with th elist of
arguments that you need to pass it. You can also check what input
arguments a demo needs by passing the ``--help`` or ``-h`` options
when executing the demo.

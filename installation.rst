Ostrich parallel processing
---------------------------

If you are on a modern computer, you can make use of its parallel
processing capabilities to save time. You can use any mpi solution, the
steps outlined here will focus on ``openmpi``, though.

System dependencies
~~~~~~~~~~~~~~~~~~~

You will need:

.. code:: bash

   sudo apt-get openmpi-bin
   sudo apt-get install libopenmpi-dev

Compiling Ostrich
~~~~~~~~~~~~~~~~~

Get the Ostrich source code from
https://www.civil.uwaterloo.ca/envmodelling/Ostrich.html . The download
link is next to the green icon of an ostrich. Next, execute the
following commands to unzip and compile the source code.

.. code:: bash

   unzip Ostrich_v17.12.19.zip Source/* -d $HOME/build
   cd $HOME/build/Source/
   make OMPI

This will create OstrichMPI inside the Source folder. These build
instructions slightly differ from the ones on the new
`Ostrich <https://usbr.github.io/ostrich/pages/development/solution/building.html>`__
webpage.

Running Ostrich in MPI mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy OstrichMPI to your project’s root folder. The executable in run
with ``mpirun ./OstrichMPI``. If you prefer to use the same executable
for multiple projects, you can create a shell script ``Ostrich_MPI.sh``
that links to the Ostrich executable:

.. code:: bash

   #!/bin/bash

   # match assignment to location of OSTRICH installation
   OSTRICH_MPI=./OstrichMPI

   mpirun $OSTRICH_MPI

Passing the flag ``-n [number of processes to run`` allows to set the
number of processes. Omitting this flag will let ``mpirun`` decide by
itself. For further information on openmpi, see the
`README <https://github.com/open-mpi/ompi/blob/v4.1.x/README>`__ file in
GitHub (for versions <v5.0) or `Read the
Docs <https://docs.open-mpi.org/en/v5.0.x/index.html>`__ (since v5.0)

Please make sure that you have selected a parallel processing compatible
``ProgramType`` in Ostrich’s config file ``ostIn.txt``.  [1]_

Ostrich Configuration
---------------------

Project Setup
=============

A generic directory tree for a model could be set up as follows: -
{model_name} - model > Contains the RAVEN model files. -
{model_name}.rv\* - data_obs - output - OstrichMPI - {model_name}.rvp -
{model_name}.rvp.tpl - save_best.sh - Ostrich_MPI.sh - ostIn.txt

Sphinx
======

See `Sphinx Tutorial <https://sphinx-tutorial.readthedocs.io/start/>`__
and
`Automodule <https://stackoverflow.com/questions/67065530/how-to-add-automodule-to-sphinx-toctree>`__.
Also: `Thorough
example <https://thomas-cokelaer.info/tutorials/sphinx/docstring_python.html>`__
and `Module
Docstrings <https://coderslegacy.com/python/module-docstrings/>`__

.. [1]
   See Table 1: Catalog of Algorithms Implemented in OSTRICH in the
   Ostrich_Manual_17_12_19.pdf for a list of compatible algorithms.

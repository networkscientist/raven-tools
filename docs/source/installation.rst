.. role:: bash(code)
  :language: bash
  :class: highlight

Installation
============

Version Control
---------------
The entire code is hosted in this repo on GitHub. To clone it, simply use:

.. code-block:: bash

    git clone https://github.com/networkscientist/raven-tools.git

Poetry
------
Package management is done through Poetry. To install Poetry, execute in the terminal [#ref_poetry1]_:

.. code-block:: bash


    curl -sSL https://install.python-poetry.org | python3 -

Afterwards, simply :bash:`cd` into the *raven-tools* folder and run

.. code-block:: bash

    poetry install

to install the dependencies. This will create a new Poetry environment, which you can then activate from within the same folder with

.. code-block:: bash

    poetry shell

Alternatively, it is possible to run a script file using this environment without activating it by issuing

.. code-block:: bash

    poetry run python your_script.py

.. rubric:: Footnotes

.. [#ref_poetry1] See `Poetry documentation <https://python-poetry.org/docs/#installation>`_, from where the Poetry installation commands have been taken.
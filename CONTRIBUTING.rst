.. image:: https://raw.githubusercontent.com/Spectre5/section-properties/dev/docs/source/_static/logo.png
  :width: 100 %
  :alt: sectionproperties
  :align: left

This document outlines how to use ``sectionproperties`` locally and get setup for making a pull request.
The ``runner`` script is written for ``bash``, but should work with ``git-bash`` on Windows as well.
However, make will also be needed on Windows to buid the docs.

Guide
-----

Clone the code, then execute the runner script which will use `poetry` to setup a virtual envivonrment, format the code, check it, and test it.
If your code changes cleanly pass the ``runner`` script without any arguments, then the code will likely pass any CI checks when the PR is opened.

.. code-block:: shell

   $ ./runner

Only run the tests with pytest:

.. code-block:: shell

   $ ./runner clean

Clean tempoary files:

.. code-block:: shell

   $ ./runner clean

Clean tempoary files and remove the virtual environment:

.. code-block:: shell

   $ ./runner clean -a

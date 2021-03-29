.. image:: https://raw.githubusercontent.com/Spectre5/section-properties/dev/docs/source/_static/logo.png
  :width: 100 %
  :alt: sectionproperties
  :align: left

This document outlines how to use ``sectionproperties`` locally and get setup for making a pull request.
The ``runner`` script is written for ``bash``, but should work with ``git-bash`` on Windows as well.

Setup
-----

This library makes use of ``poetry`` and it must be accessible on your path.
It be installed with ``pipx`` as shown below.

.. code-block:: shell

   $ python3 -m pip install --user --upgrade pipx
   $ pipx install poetry
   $ pipx upgrade poetry
   $ pipx ensurepath

``runner`` Script
-----------------

This is a convenience script to perform some common actions.
Executing the runner script without any arguments will use `poetry` to do the following:

- Create and setup a local virtual environment
- Install pre-commit hooks
- Format the code (``black`` and ``isort``)
- Lint the code (``flake8``, ``pylint``, ``rstcheck``)
- Build the docs (``Sphinx``)
- Test the code (``pytest`` with ``coverage``)
- Build and check the sdist and wheel distribution (``poetry`` and ``twine``)

Be sure that your code passes the ``runner`` script without any arguments before submitting a pull request as these checks are also enforced there.

.. code-block:: shell

   $ ./runner

See the help for various tasks to run.

.. code-block:: shell

   $ ./runner help

Benchmarking
------------

There are some benchmarking tests in the pytest suite.
These benchmarks are not run by default with the ``runner`` script nor are they run for the CI processes, since they are significantly slower than normal tests.
The benchmark tests are intended to be used by developers to compare the execution speeds when looking for improvements.
Benchmarks can be run by executing:

.. code-block:: shell

   $ ./runner benchmark

Only relative comparisons are meaningful since the computer hardware and other factors can impact the timing.
Thus, when comparing benchmarks, be sure you use the same computer and setup.
When running benchmarks repeatedly to compare changes in the code, it may be useful to run a subset of benchmarks (or only 1) so that it completes more quickly.


Testing Arguments
-----------------

Testing is performed with ``pytest`` and when running ``./runner test`` or ``./runner benchmark``, any additional arguments are passed on to ``pytest``.
This can be useful for various purposes, such as running only a subset (or only 1) test case.
A few examples are shown below, please consult the ``pytest`` documentation for many more examples and details.

Skip all of the tests marked with the "slow" marker:

.. code-block:: shell

   $ ./runner test -m "not slow"

Only run tests the match the expression "Rectangle":

.. code-block:: shell

   $ ./runner test -k "Rectangle"

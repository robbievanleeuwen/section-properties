.. image:: https://raw.githubusercontent.com/Spectre5/section-properties/dev/docs/source/_static/logo.png
  :width: 100 %
  :alt: sectionproperties
  :align: left

This document outlines how to use ``sectionproperties`` locally and get setup for making a pull request.
The ``runner`` script is written for ``bash``, but should work with ``git-bash`` on Windows as well.
However, you'll need to install ``make`` to build the docs on Windows.
On Windows, an error when building the docs is skipped since ``make`` may not be installed.

Setup
-----

This library makes use of ``poetry`` and it must be accessible on your path.
I recommend that it be installed with ``pipx`` as shown below.

.. code-block:: shell

   $ python3 -m pip install --user --upgrade pipx
   $ pipx install poetry
   $ pipx upgrade poetry
   $ pipx ensurepath

``runner`` Script
-----------------

Executing the runner script which will use `poetry` to:

- Create and setup a local virtual envivonrment
- Install pre-commit hooks
- Format the code (``isort`` and ``black``)
- Lint the code (``rstcheck``, ``flake8``, ``pylint``)
- Build the docs
- Test the code
- Build the sdist and wheel distributions

Be sure that your code passes the ``runner`` script without any arguments before submitting a pull request.

.. code-block:: shell

   $ ./runner

See the help for various tasks to run.

.. code-block:: shell

   $ ./runner help

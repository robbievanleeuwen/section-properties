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

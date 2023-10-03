.. _label-contributing:

Contributor Guide
=================

Thank you for your interest in improving this project. This project is
open-source under the `MIT
license <https://opensource.org/licenses/MIT>`__ and welcomes
contributions in the form of bug reports, feature requests, and pull
requests.

Here is a list of important resources for contributors:

-  `Source
   Code <https://github.com/robbievanleeuwen/section-properties>`__
-  `Documentation <https://sectionproperties.readthedocs.io/>`__
-  `Issue
   Tracker <https://github.com/robbievanleeuwen/section-properties/issues>`__
-  :doc:`Code of Conduct <codeofconduct>`

How to report a bug
-------------------

Report bugs on the `Issue
Tracker <https://github.com/robbievanleeuwen/section-properties/issues>`__.

When filing an issue, make sure to answer these questions:

-  Which operating system and Python version are you using?
-  Which version of this project are you using?
-  What did you do?
-  What did you expect to see?
-  What did you see instead?

The best way to get your bug fixed is to provide a test case, and/or
steps to reproduce the issue.

How to request a feature
------------------------

Features that improve ``sectionproperties`` can be suggested on the
`Issue
Tracker <https://github.com/robbievanleeuwen/section-properties/issues>`__.
It's a good idea to first submit the proposal as a feature request prior
to submitting a pull request as this allows for the best coordination of
efforts by preventing the duplication of work, and allows for feedback
on your ideas.

How to set up your development environment
------------------------------------------

You need Python 3.9, 3.10 or 3.11, and the following tools:

-  `Poetry <https://python-poetry.org/>`__
-  `Nox <https://nox.thea.codes/>`__
-  `nox-poetry <https://nox-poetry.readthedocs.io/>`__

Recommended dependency installation method:

#. Install `pipx <https://pypa.github.io/pipx/installation/>`_:

   .. code:: shell

      python3 -m pip install --user pipx
      python3 -m pipx ensurepath


#. Install `Poetry <https://python-poetry.org/>`__:

   .. code:: shell

      pipx install poetry
      poetry --version

#. Install `Nox <https://nox.thea.codes/>`__ and
   `nox-poetry <https://nox-poetry.readthedocs.io/>`__:

   .. code:: shell

      pipx install nox
      pipx inject nox nox-poetry

#. If you do not have ``pandoc`` installed, it will be required to build the docs. The
   `installation method <https://pandoc.org/installing.html>`_ depends on which OS you
   are running.

Now that you have all the dependencies up and running, you can install
``sectionproperties`` with development requirements:

.. code:: shell

   poetry install

Install with the ``rhino`` and ``cad`` extras:

.. code:: shell

   poetry install --all-extras

You can now run an interactive Python session, or the command-line interface:

.. code:: shell

   poetry run python
   poetry run sectionproperties

How to test the project
-----------------------

Run the full test suite:

.. code:: shell

   nox

List the available Nox sessions:

.. code:: shell

   nox --list-sessions

You can also run a specific Nox session. For example, invoke the unit test suite like
this:

.. code:: shell

   nox --session=tests

Unit tests are located in the *tests* directory, and are written using
the `pytest <https://pytest.readthedocs.io/>`__ testing framework.

How to submit changes
---------------------

Open a `pull
request <https://github.com/robbievanleeuwen/section-properties/pulls>`__
to submit changes to this project.

Your pull request needs to meet the following guidelines for acceptance:

-  The Nox test suite must pass without errors and warnings.
-  Include unit tests. This project aims for a high code coverage.
-  If your changes add functionality, update the documentation
   accordingly.

To run linting and code formatting checks before committing your change,
you can install pre-commit as a Git hook by running the following
command:

.. code:: shell

   nox --session=pre-commit -- install

It is recommended to open an issue before starting work on anything.
This will allow a chance to talk it over with the owners and validate
your approach.

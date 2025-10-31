Getting Started
===============

Installation
------------

.. code-block:: bash

   cd python
   poetry install --with docs

Local Build
-----------

Generate static HTML locally before publishing:

.. code-block:: bash

   cd python
   poetry run sphinx-build -b html ../docs ../docs/_build/html

The output is written to ``docs/_build/html``.

Live Preview
------------

Install the optional live-reload server and preview the docs as you edit them:

.. code-block:: bash

   cd python
   poetry run sphinx-autobuild ../docs ../docs/_build/html

Open ``http://127.0.0.1:8000`` in your browser to view the live site.

Publishing
----------

Pushing to the ``main`` branch triggers the GitHub Pages workflow, which installs the
project, builds the Sphinx documentation, and deploys it automatically.

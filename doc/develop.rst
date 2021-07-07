Development
===========

- Installation_
- `Conda environment files`_
- Testing_

Installation
------------

Install ``dbspro`` for development using

..  code-block:: bash

    pip install -e .[dev]

This installs ``dbspro`` in editable mode along with ``pytest`` for testing and 
``flake8`` for stylechecks.

Conda environment files
-----------------------

Strongly inspired by from: https://github.com/FrickTobias/BLR/doc/develop.rst

There are two types of files that describe Conda environments.

- The file ``environment.yml`` contains abstract dependencies such as ``pysam``.
  This file is managed manually and needs to be updated whenever there are new
  dependencies or when the required version for a dependency changes.

- The ``environment.linux.lock.yml`` and ``environment.osx.lock.yml`` files
  (lock files) contain a fully specified description of the entire environment,
  with locked-down versions.  These files are used to create the test
  environment.

Use the script ``misc/condalock.sh`` to update the lock files whenever you make
changes to ``environment.yml``.

Testing
-------

To run tests `pytest` need to be installed. 

Before running also make sure that the environement variable ``PYTHONHASHSEED`` 
is set to ``1``. This is required to make UMI-tools (used in ``dbspro.cli.splitcluster``) 
consistent. This can be done using:

..  code-block:: bash

    export PYTHONHASHSEED=1 

To run test use the following line from the base directory:

..  code-block:: bash

    pytest -v tests/

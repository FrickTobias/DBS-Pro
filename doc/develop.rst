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

There are three kinds of environement files in the repository:

- ``environement.yml``: Flexible base YAML for dependencies
- ``environement.linux-64.lock``: Lock file for fully reproducible Linux environement
- ``environement.osx-64.lock``: Lock file for fully reproducible MacOS environement

Whenever the file ``environment.yaml`` is updated new lock-files need to be generated. 
This is done using `conda-lock <https://pypi.org/project/conda-lock/>`_ which can be install 
using:

..  code-block:: bash

    pip install conda-lock

To update the lock files run:

..  code-block:: bash

    conda-lock -f environment.yml -p linux-64 -p osx-64 --filename-template "environment.{platform}.lock"

    
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

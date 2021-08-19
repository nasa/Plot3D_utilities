Installation
============

.. note::
    We do not recommend installation as root user on your system python.
    Please setup an `Anaconda/Miniconda <https://conda.io/docs/user-guide/install/index.html/>`_ environment or create a `Docker image <https://www.docker.com/>`_.

Please follow the steps below for a successful installation.

Installation via Pip
-------------------------

#. Install the relevant packages:

.. code-block:: none

    pip install plot3d


Installation from Source
-------------------------

Download the git repository and navigate to the python folder where you see setup.py. Run the following command `python setup.py install` and that should install the library in your python environment. 


.. note:: 

    For CentOS you need to run this command prior to installation *scl enable devtoolset-7 bash*

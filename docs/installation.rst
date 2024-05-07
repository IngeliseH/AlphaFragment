Installation
============

Requirements
------------
- Python 3.6 or newer
- pandas>=1.3.5
- requests>=2.31.0
- matplotlib>=3.5.3

Installation from Source
------------------------

To install AlphaFragment, you will need to clone the repository and install it manually. Here are the steps:

1. Clone the repository:

.. code-block:: bash

    git clone https://github.com/IngeliseH/AlphaFragment.git
    cd AlphaFragment

2. Install the package:

.. code-block:: bash

    pip install .

This will install AlphaFragment along with its required dependencies. Ensure that you have `git` installed on your system to clone the repository.

Alternatively, if you have direct access to the package's archive, you can install it using:

.. code-block:: bash

    pip install /path/to/AlphaFragment.tar.gz

Replace `/path/to/AlphaFragment.tar.gz` with the actual path to the source archive.

Verifying Installation
----------------------

After installation, you can verify that AlphaFragment has been installed correctly by checking its version:

.. code-block:: bash

    python -c "import alphafragment; print(alphafragment.__version__)"

This should print the installed version of AlphaFragment.

Troubleshooting
---------------

If you encounter any issues during installation, ensure that you have the latest version of `pip` and `setuptools`:

.. code-block:: bash

    pip install --upgrade pip setuptools

Also, check that you meet all the Python version and dependency requirements.
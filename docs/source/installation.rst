Installation and updating
=========================

Requirements
------------

PEAKachu was developed using Python 3.5 and the user is advised to run PEAKachu
with this or a higher version. The following third party software is required
to run PEAKachu with all features:

- The Python packages `setuptools <https://pythonhosted.org/setuptools>`_,
  `pip <http://www.pip-installer.org>`_,
  `biopython <http://biopython.org>`_,
  `pysam <http://pysam.readthedocs.org/en/latest/api.html>`_,
  `pandas <http://pandas.pydata.org>`_,
  `matplotlib <http://matplotlib.org/>`_,
  `bcbio-gff <https://github.com/chapmanb/bcbb/tree/master/gff>`_,
  `statsmodels <http://statsmodels.sourceforge.net/>`_ and
  `rpy2 <http://rpy2.bitbucket.org/>`_.

- The R bioconductor packages
  `edgeR <https://bioconductor.org/packages/release/bioc/html/edgeR.html>`_ and
  `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_.
- The tool
  `blockbuster <http://hoffmann.bioinf.uni-leipzig.de/LIFE/blockbuster.html>`_.


Installing on a clean Ubuntu system
-----------------------------------

The following installation procedure was tested on a newly generated Docker container
based on an Ubuntu 15.10 image.

1. Install required Debian/Ubuntu packages (root privileges required)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before starting it is a good idea to update the package list::

  apt-get update

Now you can install the packages::

  apt-get install wget python3-dev build-essential zlib1g-dev libfreetype6-dev pkg-config \
                  gfortran libopenblas-dev liblapack-dev r-base r-base-dev r-cran-xml

Some comments:

- Some packages not available in the Docker image are already present in the
  installation version of Ubuntu 15.10.
- ``python3-dev`` and ``build-essential`` are required for ``pip3``
  installation and package compilation
- ``gfortran``, ``libopenblas-dev`` and ``liblapack-dev`` are required for ``scipy``
- ``zlib1g-dev`` is required for ``pysam``
- ``libfreetype6-dev`` and ``pkg-config`` are required for ``matplotlib``
- ``r-base`` is required for ``rpy2``
- ``r-base-dev`` and ``libxml2-dev`` are required for R package dependencies

2. Install R packages
~~~~~~~~~~~~~~~~~~~~~
Start ``R``::

    R

and install the R packages inside of the interactive command line
interface. You might be asked to confirm the installation path::

  source("http://bioconductor.org/biocLite.R")
  biocLite(c("edgeR", "DESeq2"))

Leave ``R``::

    quit(save = "no")


3. Install blockbuster
~~~~~~~~~~~~~~~~~~~~~~
Download, extract and compile the source code::

  wget http://www.bioinf.uni-leipzig.de/~david/LIFE/blockbuster.tar.gz
  tar -xvf blockbuster.tar.gz
  cd blockbuster && make && cd ..

Copy the executable to a directory that is part of your ``PATH``, e.g. to
``/usr/local/bin`` (with root privileges)::

    cp blockbuster/blockbuster.x /usr/local/bin/

... or the bin folder of your home directory::

    mkdir ~/bin
    cp blockbuster/blockbuster.x ~/bin/

4. Install Python 3 packages via the `Python Package Index (PyPI) <https://pypi.python.org/pypi>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Install pip3::

  wget https://bootstrap.pypa.io/get-pip.py
  python3 get-pip.py

Install PyPI packages::

  pip3 install Cython numpy
  pip3 install scipy
  pip3 install pysam pandas biopython matplotlib bcbio-gff statsmodels numexpr rpy2

5. Install PEAKachu
~~~~~~~~~~~~~~~~~~~
::

  pip3 install PEAKachu

Updating PEAKachu
--------------------

Once you have installed PEAKachu as described above you can easily
update it to the newest version by running::

  pip3 install PEAKachu -U

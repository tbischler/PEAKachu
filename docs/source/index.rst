PEAKachu - peak calling based on RIP/CLIP-seq data
**************************************************
Table of content
================

.. toctree::
   :maxdepth: 1

   installation
   license
   versions

PEAKachu in a nutshell
=========================
PEAKachu is a peak calling tool for the identification of RNA tagets of specific
proteins based on RIP/CLIP-seq data. The input consists of mapped read files in
BAM format for ideally several replicates of experiment and control libraries.
While experiment libraries for CLIP-seq data consist of samples where the
protein was tagged and cross-linked to the target RNA, the control libraries
can either result from non-tagged cross-linked or tagged non-cross-linked
samples.
For RIP-seq experiments the experiment libraries use tagged proteins while the
controls just consist of the wild-type.
In both cases the RNA is pulled-down using a tag-specific antibody and
subsequently subjected to library preparation and deep sequencing.

PEAKachu implements two alternative approaches which both rely on the
comparison of read coverage in experiment and control libraries along a
reference genome. The window approach calculates the coverage using a sliding
window of predefined size where neighboring significantly enriched windows
locations are merged into peaks. The predefined_peak approach is based on
potential peaks derived from read clusters identified by blockbuster and
subsequent comparison to the control libraries via DESeq2.


Download
========

PEAKachu can be download from `its PyPI page
<https://pypi.python.org/pypi/PEAKachu/>`_. Please read the
`installation instructions <installation.html>`_.


Source code
===========

The source code of PEAKachu can be found at https://github.com/tbischler/PEAKachu.

Cite
====


Contact
=======

For question and requests feel free to contact Thorsten Bischler
<thorsten.bischler@uni-wuerzburg.de>

Peak calling tool for CLIP-seq data

Installation
------------

PEAKachu is available via `bioconda <https://bioconda.github.io/>`__ and
can be installed using

::

   conda create -n peakachu peakachu python=3

Alternativaly, PEAKachu can be installed from source

::

   $ make readme_rst
   $ make package
   $ pip3 install --user dist/PEAKachu-0.2.0.tar.gz

License
-------

`ICSL <https://en.wikipedia.org/wiki/ISC_license>`__ (Internet Systems
Consortium license ~ simplified BSD license) - see LICENSE

Development
-----------

-  The git braching model is very close to the one proposed
   `here <http://nvie.com/posts/a-successful-git-branching-model/>`__.
   There two main branches:

   -  master
   -  dev(elopment)

   And there are further supporting branches:

   -  feature branches - branched off and back to the dev branch
   -  release branches - branched off from dev and merged back into dev
      and master
   -  hotfix branches - branched off from master and merged back into
      dev and master

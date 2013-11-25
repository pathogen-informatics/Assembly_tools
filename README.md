Assembly_tools
==============

Scripts relating to genome assemblies

Installation
------------

Prerequisites:

  * [Fastaq] [Fastaq] >= v0.1
  * [Bowtie2] [Bowtie2] (only required for some scripts)

Once the prerequisites are installed, run the tests (these do not check if bowtie2 is installed and in your path):

    python3 setup.py test

If all was OK, then install:

    python3 setup.py install


[Fastaq]: https://github.com/sanger-pathogens/Fastaq
[Bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

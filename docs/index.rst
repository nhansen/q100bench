.. Q100bench documentation master file, created by
   sphinx-quickstart on Fri Aug  4 13:33:33 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Q100bench
==========

.. toctree::
   :hidden:
   :maxdepth: 2
   :name: mastertoc

   Install <install>

`Q100bench <http://github.com/nhansen/q100bench>`_ is a software package for reporting
differences between an assembly or read dataset and a benchmark assembly, such as the
Q100 assembly of the sample HG002.

Install
========

Software dependencies
---------------------

Q100bench requires bedtools (https://bedtools.readthedocs.io/en/latest/), and R with Bioconductor (https://www.bioconductor.org/).

Until q100bench is available on PyPi and bioconda, the easiest way to use it is to install it locally. First clone this github repository:

::

  git clone https://github.com/nhansen/q100bench.git
  cd q100bench

Create a virtual environment for the project:

::

  python3 -m venv venv
  source venv/bin/activate

Finally use python's pip installer to install and test a development copy for yourself to run:

::

  python3 -m pip install -e .
  pytest



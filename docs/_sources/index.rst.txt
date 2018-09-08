.. test documentation master file, created by
   sphinx-quickstart on Sat Jan  7 22:10:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: ../bib/ref.txt
   
What is mVMC?
=============

High-accuracy analyses of theoretical models for quantum many-body
systems are expected to play important role for clarifying the nature of
novel quantum phases such as high-temperature superconductivities and
quantum spin liquids. Furthermore, recent theoretical progress enables
us to obtain low-energy effective models for real materials with
non-empirical way [ImadaMiyake_ ]. To clarify the
electronic structures of real materials and control the physical
properties, it is important to perform the high-accuracy analyses for
the low-energy effective models. One of the most reliable theoretical
tools for treating the strongly correlated electron systems is the exact
diagonalization method. However, applicable range of system size is
strongly limited in the exact diagonalization method. Variational Monte
Carlo method [Gros_ ] is one of promising way to perform
the high-accuracy calculations for the larger system sizes beyond exact
diagonalization method. Although the strong limitation of the
variational wave function is the origin of the poor accuracy of the
variational Monte Carlo method, recent development of the theoretical
method and computers relaxes the limitation of the variational wave
functions and enable us to perform the highly-accurate
calculations [Tahara2008_ , Misawa2014_ , Morita2015_ ].

mVMC (many-variable Variational Monte Carlo method) is a software for
performing the highly-accurate variational Monte Carlo calculations with
the simple and flexible user interface. mVMC also supports the
large-scale parallelization. For the conventional models in strongly
correlated electron systems such as the Hubbard model, the Heisenberg
model, and the Kondo-lattice model, users can perform the calculation by
preparing the one input files whose length is shorter than ten lines. By
using the same input file, users can perform the exact diagonalization
through HPhi [HPhi_ ]. Thus, it is easy to examine the accuracy of the
variational calculation for small system sizes and to perform the
calculations for large system sizes that can not be treated by the exact
diagonalization. A broad spectrum of users including experimental
scientists is cordially welcome.

Overview of mVMC
----------------

By using mVMC, the following calculation can be done:

-  The variational wave function which gives the minimum value of the
   expected value of energy in the range of the degree of freedoms of
   variational parameters is numerically generated. The calculation
   limited to the partial space divided by quantum numbers is also
   possible.

-  The expected values of the physical quantities such as correlation
   functions can be calculated by using the generated variational wave
   functions.

The calculation flow in mVMC is shown as follows:

#. Read input files (\*.def)

#. Optimize variational parameters :math:`\vec{\alpha}` to minimize
   :math:`\langle {\cal H} \rangle`

#. Calculate one body and two body Green functions

#. Output variational parameters and expected values

In calculation, the simple parallelization for the generation of the
real space arrangement :math:`|x\rangle` and the collecting samples and
the calculation result of expected energies is implemented. Following
the procedure for each cluster computers, the parallelized calculation
using MPI is automatically done by indicating the parallel number.
However, mVMC cannot execute under the environment where the MPI job is
forbidden such as the front-end of system B at ISSP. In mVMC, we use
PFAPACK [PFAPACK_ ] to compute the Pfaffian matrix.

License
-------

The distribution of the program package and the source codes for mVMC
follows GNU General Public License version 3 (GPL v3).

We would appreciate if you cite the following article in your research with mVMC:  
mVMC - Open-source software for many-variable variational Monte Carlo method, Takahiro Misawa, Satoshi Morita, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Yuichi Motoyama, Kota Ido, Takahiro Ohgoe, Masatoshi Imada, Takeo Kato, arXiv:1711.11418. [mVMC_ ]

Copyright
---------

    *©2016- The University of Tokyo.* *All rights reserved.*

This software is developed under the support of “\ *Project for
advancement of software usability in materials science* " by The
Institute for Solid State Physics, The University of Tokyo.

Contributors
------------

This software is developed by following contributors.

-  ver.1.0.3 (released at 2018/7/23)

-  ver.1.0.2 (released at 2017/8/25)

-  ver.1.0.1 (released at 2017/6/8)

-  ver.1.0.0 (released at 2017/5/23)

-  ver.0.2.0 (released at 2017/3/16)

-  ver.0.1.1 (released at 2016/12/16)

-  ver.0.1.0 (released at 2016/10/26)

   -  Developers

      -  Takahiro Misawa
         (The Institute for Solid State Physics, The University of
         Tokyo)

      -  Satoshi Morita
         (The Institute for Solid State Physics, The University of
         Tokyo)

      -  Takahiro Ogoe
         (Department of Applied Physics, The University of Tokyo)

      -  Kota Ido
         (Department of Applied Physics, The University of Tokyo)

      -  Masatoshi Imada
         (Department of Applied Physics, The University of Tokyo)

      -  Yuichi Motoyama
         (The Institute for Solid State Physics, The University of
         Tokyo)

      -  Mitsuaki Kawamura
         (The Institute for Solid State Physics, The University of
         Tokyo)

      -  Kazusyohi Yoshimi
         (The Institute for Solid State Physics, The University of
         Tokyo)

   -  Project coordinator

      -  Takeo Kato
         (The Institute for Solid State Physics, The University of
         Tokyo)

Operating environment
---------------------

mVMC is tested in the following platform:

-  The supercomputer system-B “sekirei”S in ISSP

-  K computer

-  OpenMPI + Intel Compiler + MKL

-  MPICH + Intel Compiler + MKL

-  MPICH + GNU Compiler + MKL

   
   
Manual for mVMC
===============

.. toctree::
   :maxdepth: 2

   start.rst
   tutorial.rst
   standard.rst
   expert.rst
   output.rst
   algorithm.rst
   appendix.rst
   fourier/index
   wannier/index
   acknowledge.rst

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


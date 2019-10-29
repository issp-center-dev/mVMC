How to use mVMC?
================

Prerequisite
------------

mVMC requires the following packages:

-  C compiler (intel, Fujitsu, GNU, etc. )

-  MPI library

-  Option: ScaLAPACK library (intel MKL, Fujitsu, ATLAS, etc.)

.. note::
   
   **Settings of intel compiler**

   When you use the intel compiler, you can use easily scripts attached to
   the compiler. In the case of the bash in 64 bit OS, write the following
   in your ``~/.bashrc``:

   .. code-block:: bash

      $ source /opt/intel/bin/compilervars.sh intel64

   or

   .. code-block:: bash

      $ source /opt/intel/bin/iccvars.sh intel64
      $ source /opt/intel/mkl/bin/mklvars.sh

   Please read manuals of your compiler/library for more information.

Installation
------------

You can download mVMC in the following place.
https://github.com/issp-center-dev/mVMC/releases

You can obtain the mVMC directory by typing

.. code-block:: bash

   $ tar xzvf mVMC-xxx.tar.gz

There are two kind of procedures to install mVMC.

Using ``mVMCconfig.sh``
~~~~~~~~~~~~~~~~~~~~~~~

Please run ``mVMCconfig.sh`` script in the mVMC directory as follow (for
ISSP system-B "sekirei"):

.. code-block:: bash

   $ bash mVMCconfig.sh sekirei

Then environmental configuration file ``make.sys`` is generated in
``src/`` directory. The command-line argument of ``mVMCconfig.sh`` is as
follows:

-  ``sekirei`` : ISSP system-B "sekirei"

-  ``kei`` : K computer and ISSP system-C "maki"

-  ``intel-openmpi`` : Intel compiler + OpenMPI

-  ``intel-mpich`` : Intel compiler + MPICH2

-  ``intel-intelmpi`` : Intel compiler + IntelMPI

-  ``gcc-openmpi`` : GCC + OpenMPI

-  ``gcc-mpich-mkl`` : GCC + MPICH + MKL

``make.sys`` is as follows (for ISSP-system-B "sekirei"):

.. code-block:: makefile

   CC = mpicc
   F90 = mpif90
   CFLAGS = -O3 -no-prec-div -xHost -qopenmp -Wno-unknown-pragmas
   FFLAGS = -O3 -implicitnone -xHost
   LIBS = -L $(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
           -lmkl_intel_thread -lmkl_core -lmkl_blacs_sgimpt_lp64 -lpthread -lm
   SFMTFLAGS = -no-ansi-alias -DHAVE_SSE2

We explain macros of this file as:

-  ``CC`` : C compiler (``mpicc``, ``mpifccpx``)

-  ``F90`` : fortran compiler (``ifort``, ``frtpx``)

-  ``Libs`` : Linker option

-  ``CFLAGS`` : C compile option

-  ``FFLAGS`` : fortran compile option

Then you are ready to compile mVMC. Please type

.. code-block:: bash

   $ make mvmc

and obtain ``vmc.out`` and ``vmcdry.out`` in ``src/`` directory; you
should add this directory to the ``$PATH``.

You can make a PATH to mVMC as follows:
``$ export PATH=${PATH}:``\ *mVMC\_top\_directory*\ ``/src/``
If you keep this PATH, you should write above in ``~/.bashrc`` (for
``bash`` as a login shell)

Using ``cmake``
~~~~~~~~~~~~~~~

We can compile mVMC as

.. code-block:: bash

   cd $HOME/build/mvmc
   cmake -DCONFIG=gcc $PathTomvmc
   make

Here, we set a path to mVMC as ``$PathTomvmc`` and to a build directory
as `` $HOME/build/mvmc``. After compilation, ``src`` directory is
constructed below a ``$HOME/build/mvmc`` directory and we obtain an
executable ``vmc.out`` in ``src/`` directory.

In the above example, we compile mVMC by using a gcc compiler. We can
select a compiler by using the following options:

-  ``sekirei`` : ISSP system-B "sekirei"

-  ``fujitsu`` : Fujitsu compiler

-  ``intel`` : Intel compiler + Linux PC

-  ``gcc`` : GCC compiler + Linux PC.

An example of compiling mVMC by using the Intel compiler is shown as
follows:

.. code-block:: bash

   mkdir ./build
   cd ./build
   cmake -DCONFIG=intel ../
   make

After compilation, ``src/`` directory is created below the ``build/`` directory
and an execute ``vmc.out`` in the ``src/`` directory. We can select
ScaLAPACK instead of LAPACK for vmc calculation by adding the following
as the cmake option

::

    -DUSE_SCALAPACK=ON -DSCALAPACK_LIBRARIES="xxx".

Here, xxx is the libraries to use ScaLAPACK. Please note that we must
delete the ``build/`` directory and repeat the above operations when we
change the compiler.

.. note::

   Before using cmake for sekirei, you must type

   .. code-block:: bash

      $ source /home/issp/materiapps/tool/env.sh

   When we type the following in sekirei,

   .. code-block:: bash

      $ cmake -DCONFIG=sekirei ../ -DUSE_SCALAPACK=ON ,

   ``-DSCALAPACK_LIBRARIES`` is automatically set as

   ::

      -DSCALAPACK_LIBRARIES="\${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 
      -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
      -lmkl_blacs_sgimpt_lp64".

   When the path to libraries for ScaLAPACK is different in your
   circumstance, please set ``-DSCALAPACK_LIBRARIES`` as the correct path.

Directory structure
-------------------
 When mVMC-xxx.tar.gz is unzipped, the following directory structure is composed.

::
   
   | |──COPYING
   | |──mVMCconfig.sh
   | |──doc/
   | │      |──bib/
   | │      │      |──elsart-num\_mod.bst
   | │      │     |──userguide.bib
   | │      |──figs/
   | │      │      |──\*.pdf
   | │      │     |──\*.xbb
   | │      |──fourier/
   | │      │      |──en/
   | │      │      |──figs/
   | │      │     |──ja/
   | │      |──jp/
   | │      │     |──\*.tex
   | │     |──en/
   | │            |──\*.tex
   | |──sample/
   | │     |──Standard/
   | │                  |──Hubbard/
   | │                  │      |─square/
   | │                  │      │      |──StdFace.def
   | │                  │      │     |──reference/
   | │                  │      │                |──\*\*.dat
   | │                  │     |─triangular/
   | │                  │           |──\ :math:`\cdots`
   | │                  |──Kondo/
   | │                  │     |─chain/
   | │                  │           |──\ :math:`\cdots`
   | │                 |──Spin/
   | │                              |─HeisenbergChain/
   | │                              │     |──\ :math:`\cdots`
   | │                              |─HeisenbergSquare/
   | │                              │     |──\ :math:`\cdots`
   | │                             |─Kagome/
   | │                                    |──\ :math:`\cdots`
   | |──src/
   | │          |──mVMC/
   | │          │      |─ \*\*.c
   | │          │     |──include/
   | │          │             |──\*\*.h
   | │          |──ComplexUHF/
   | │          │      |─ \*\*.c
   | │          │     |──include/
   | │          │             |──\*\*.h
   | │          |──StdFace/
   | │          │       |──\*\*.c
   | │          │      |──\*\*.h
   | │          |──pfapack/
   | │          │       |──makefile\_pfapack
   | │          │      |──\*\*.f
   | │          |──sfmt/
   | │                  |──makefkie\_sfmt
   | │                  |──\*\*.c
   | │                 |──\*\*.h
   | |──tool/
   |           |──\*\*.f90
   |           |──makefile\_tool

Basic usage
-----------

mVMC works as whether the following two modes:

-  Expert mode

   mVMC supports the arbitrary fermion-/spin-lattice system; we can
   specify the hopping, etc. at each site independently. Although this
   makes us able to specify flexibly the target this requires many
   input-files, and the setup of the calculation is complicated.

-  Standard mode

   For some typical models (such as the Heisenberg model on the square
   lattice), we can start calculation with a few parameters (for
   example, the size of the simulation cell, the common coupling
   parameter). In this case, the input-files for Expert mode are
   automatically generated. Although the number of available systems is
   smaller than that number of Expert mode, the setup of the calculation
   is easier than in Expert mode.

We can calculate by using these modes as follows:

#. Prepare a minimal input file

   You can choose a model (the Heisenberg model, the Hubbard model,
   etc.) and a lattice (the square lattice, the triangular lattice,
   etc.) from ones provided; you can specify some parameters (such as
   the first/second nearest neighbor hopping integrals, the on-site
   Coulomb integral, etc.) for them. The input file format is described
   in :ref:`HowToStandard`.

#. Run

   Run a executable ``vmc.out`` in terminal by specifying the name of
   input file written in previous step (option ``-s`` is required).

   .. code-block:: bash

      $ mpiexec -np number-of-processes Path/vmc.out -s Input-file-name

   When you use a queuing system in workstations or super computers,
   sometimes the number of processes is specified as an argument for the
   job-submitting command. If you need more information, please refer
   manuals for your system.

#. Watch calculation logs

   Log files are outputted in the ``output/`` directory which is
   automatically made in the directory for a calculation scenario. The
   details of output files are shown in :ref:`OutputFile`.

#. Results

   If the calculation is finished normally, the result files are
   outputted in the ``output/`` directory. The details of output files
   are shown in :ref:`OutputFile`.

#. Prepare and run Expert mode

   In the above case, the calculation starts as soon as input files for
   Expert mode are generated. If we only generate files without starting
   the calculation, we can use a executable ``vmcdry.out`` as follows
   (MPI is not used in this step):

   .. code-block:: bash

      $ Path/vmcdry.out Input-file-name

   Then, we can edit generated files by hand and run a executable
   ``vmc.out`` with ``namelist.def`` as an argument (option ``-e`` is
   required) as follows:

   .. code-block:: bash

      $ mpiexec -np number-of-processes Path/vmc.out -e namelist.def

.. note::
      
   **The number of threads for OpenMP**

   If you specify the number of OpenMP threads for mVMC, you should set it
   as follows (in case of 16 threads) before the running:

   .. code-block:: bash

      $ export OMP_NUM_THREADS=16

Printing version ID
-------------------

By using ``-v`` option as follows, you can check which version of mVMC
you are using.

.. code-block:: bash

   $ PATH/vmcdry.out -v
   $ PATH/vmc.out -v

Tutorial
========

List of sample files
--------------------

There are following tutorials in ``samples/Standard/``.

-  The Hubbard model on the two dimensional square lattice

   (``samples/Standard/Hubbard/square/``)

-  The Hubbard model on the two dimensional triangular lattice

   (``samples/Standard/Hubbard/triangular/``)

-  The one dimensional Kondo chain

   (``samples/Standard/Kondo/chain/``)

-  The one dimensional antiferromagnetic Heisenberg chain

   (``samples/Standard/Spin/HeisenbergChain/HeisenbergChain/``)

-  The antiferromagnetic Heisenberg model on the two dimensional square
   lattice

   (``samples/Standard/Spin/HeisenbergSquare/``)

-  The antiferromagnetic Heisenberg model on the two dimensional Kagome
   lattice

   (``samples/Standard/Spin/Kagome/``)

We can perform these tutorials in the same way. In the following, the
tutorial of the one dimensional antiferromagnetic Heisenberg chain is
shown.

Heisenberg model
----------------

This tutorial should be performed in

::

    sample/Standard/Spin/HeisenbergChain/.

This directory contains the following things:

The input file: StdFace.def

reference outputs: reference/

In this case, we treat the one dimensional antiferromagnetic Heisenberg
chain which has a nearest neighbor spin coupling.

.. math::

   \begin{aligned}
   {\hat H} = J \sum_{i=1}^{L} {\hat {\boldsymbol S}}_i \cdot {\hat {\boldsymbol S}}_{i+1}
   \end{aligned}

The input file is as follows:

::

    L = 16
    Lsub=4
    model = "Spin"
    lattice = "chain lattice"
    J = 1.0
    2Sz = 0
    NMPtrans=1

In this tutorial, J and the number of sites are set to 1 (arbitrary
unit) and 16 respectively.

Running
^^^^^^^

We execute the following command.

.. code-block:: bash

   $ mpiexec -np number-of-processes Path/vmcdry.out -s StdFace.def

The MPI command depends on your system (such as ``mpiexec``, ``mpirun``,
``mpijob``, ``poe``, etc.). Then, the calculation starts and the
following standard message is outputted in the terminal.

::

    ######  Standard Intarface Mode STARTS  ######

      Open Standard-Mode Inputfile StdFace.def 

      KEYWORD : l                    | VALUE : 16 
      KEYWORD : lsub                 | VALUE : 4 
      KEYWORD : model                | VALUE : spin 
      KEYWORD : lattice              | VALUE : chain 
      KEYWORD : j                    | VALUE : 1.0 
      KEYWORD : nmptrans             | VALUE : 1 

    #######  Parameter Summary  #######

      @ Lattice Size & Shape

                    L = 16 
                 Lsub = 4         
                    L = 16        
                    W = 1         
               phase0 = 0.00000     ######  DEFAULT VALUE IS USED  ######

      @ Hamiltonian 

                   2S = 1           ######  DEFAULT VALUE IS USED  ######
                    h = 0.00000     ######  DEFAULT VALUE IS USED  ######
                Gamma = 0.00000     ######  DEFAULT VALUE IS USED  ######
                    D = 0.00000     ######  DEFAULT VALUE IS USED  ######
                  J0x = 1.00000   
                  J0y = 1.00000   
                  J0z = 1.00000   

      @ Numerical conditions

                 Lsub = 4         
                 Wsub = 1         
          ioutputmode = 1           ######  DEFAULT VALUE IS USED  ######

    ######  Print Expert input files  ######

        qptransidx.def is written.
             filehead = zvo         ######  DEFAULT VALUE IS USED  ######
             filehead = zqp         ######  DEFAULT VALUE IS USED  ######
          NVMCCalMode = 0           ######  DEFAULT VALUE IS USED  ######
         NLanczosMode = 0           ######  DEFAULT VALUE IS USED  ######
        NDataIdxStart = 1           ######  DEFAULT VALUE IS USED  ######
          NDataQtySmp = 1           ######  DEFAULT VALUE IS USED  ######
          NSPGaussLeg = 8           ######  DEFAULT VALUE IS USED  ######
             NMPTrans = 1         
        NSROptItrStep = 1000        ######  DEFAULT VALUE IS USED  ######
         NSROptItrSmp = 100         ######  DEFAULT VALUE IS USED  ######
           NVMCWarmUp = 10          ######  DEFAULT VALUE IS USED  ######
         NVMCInterval = 1           ######  DEFAULT VALUE IS USED  ######
           NVMCSample = 1000        ######  DEFAULT VALUE IS USED  ######
        NExUpdatePath = 2         
              RndSeed = 123456789   ######  DEFAULT VALUE IS USED  ######
           NSplitSize = 1           ######  DEFAULT VALUE IS USED  ######
               NStore = 0           ######  DEFAULT VALUE IS USED  ######
         DSROptRedCut = 0.00100     ######  DEFAULT VALUE IS USED  ######
         DSROptStaDel = 0.02000     ######  DEFAULT VALUE IS USED  ######
         DSROptStepDt = 0.02000     ######  DEFAULT VALUE IS USED  ######
              NSPStot = 0           ######  DEFAULT VALUE IS USED  ######
          ComplexType = 0           ######  DEFAULT VALUE IS USED  ######
        locspn.def is written.
        trans.def is written.
        interall.def is written.
        jastrowidx.def is written.
        coulombintra.def is written.
        coulombinter.def is written.
        hund.def is written.
        exchange.def is written.
        orbitalidx.def is written.
        gutzwilleridx.def is written.
        namelist.def is written.
        modpara.def is written.
        greenone.def is written.
        greentwo.def is written.

    ######  Input files are generated.  ######
    -----------
    Start: Read *def files.
      Read File namelist.def .
      Read File 'modpara.def' for ModPara.
      Read File 'locspn.def' for LocSpin.
      Read File 'trans.def' for Trans.
      Read File 'coulombintra.def' for CoulombIntra.
      Read File 'coulombinter.def' for CoulombInter.
      Read File 'hund.def' for Hund.
      Read File 'exchange.def' for Exchange.
      Read File 'gutzwilleridx.def' for Gutzwiller.
      Read File 'jastrowidx.def' for Jastrow.
      Read File 'orbitalidx.def' for Orbital.
      Read File 'qptransidx.def' for TransSym.
      Read File 'greenone.def' for OneBodyG.
      Read File 'greentwo.def' for TwoBodyG.
    End  : Read *def files.
    Start: Read parameters from *def files.
    End  : Read parameters from *def files.
    Start: Set memories.
    End  : Set memories.
    Start: Initialize parameters.
    End  : Initialize parameters.
    Start: Initialize variables for quantum projection.
    End  : Initialize variables for quantum projection.
    Start: Optimize VMC parameters.
    End  : Optimize VMC parameters.
    -----------

In the beginning of this run, files describing the detail of considered
Hamiltonian

-  ``locspin.def``

-  ``trans.def``

-  ``coulombinter.def``

-  ``coulombintra.def``

-  ``exchange.def``

-  ``hund.def``

-  ``namelist.def``

-  ``modpara.def``

and files for setting variational parameters

-  ``gutzwilleridx.def``

-  ``jastrowidx.def``

-  ``orbitalidx.def``

-  ``qptransidx.def``

and files specifying elements of correlation functions that will be
calculated

-  ``greenone.def``

-  ``greentwo.def``

are generated. The details of these files are shown in :ref:`HowToExpert`.

During the calculation, the following files are outputted in
``output`` directory:

::

    zvo_SRinfo.dat
    zvo_out_001.dat
    zvo_time_001.dat
    zvo_var_001.dat
    zvo_CalcTimer.dat

In ``zvo_out_001.dat``, the following quantities are outputted at each
bins

.. math::

   \langle H \rangle, \langle H^2 \rangle,
   \frac{\langle H^2 \rangle- \langle H \rangle^2 }{\langle H \rangle^2} \nonumber.

By seeing these informations, the conversion of the calculation can be
judged. By using gnuplot, we can check the evolution of
:math:`\langle H \rangle` as follows:

.. code-block:: gnuplot

   gnuplot> plot "zvo_out_001.dat" u 1

The details of these outputted files are shown in :ref:`OutputFile`.

Output results
^^^^^^^^^^^^^^

After finishing calculation normally, the files for the energy, the
deviation, the optimized variational parameters and the time of
execution for each calculation steps are outputted in ``output/``
directory. In the following, the outputted files are shown

::

    gutzwiller_opt.dat
    jastrow_opt.dat
    orbital_opt.dat
    zqp_opt.dat
    ClacTimer.dat

The details of these outputted files are shown in :ref:`OutputFile`.

Calculation of Green functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After changing the value of ``NVMCCalMode`` from 0 to 1 in
``modpara.def`` file, we execute the following command. When we add
``"zqp_opt.dat"`` after ``"namelist.dat"`` as a command-line argument as
follows, the calculation of Green functions is done by using the
optimized variational parameters.

.. code-block:: bash
                
   $ Path/vmc.out -e namelist.def output/zqp_opt.dat

After the calculation finishes, the following files are outputted in
``output/`` directory.

::

    zvo_cisajs_001.dat
    zvo_cisajscktalt_001.dat

The details of these outputted files are shown in :ref:`OutputFile`.

Input files for Expert mode
---------------------------

In mVMC, the calculation is done by reading input files categorized by
the following six parts.

(1) List:
    Specify the kinds and names of input files.

(2) Basic parameters:
    Specify the basic parameters.

(3) Set Hamiltonian:
    Specify the Hamiltonian.

(4) Set condition of variational parameters :
    Specify the variational parameters to be optimized.

(5) Initial variational parameters:
    Specify the initial values of the variational parameters.

(6) Output:
    Specify the components of one-body and two-body Green’s functions to
    be outputted.

The calculation for complex models can be done by directly making above
input files. The details for each files are shown in :ref:`HowToExpert`.

Fourier transformation of correlation functions
-----------------------------------------------

This package has a utility which performs the Fourier transformation of
the correlation function and plots that function.
For more details, please see :ref:`fourier`.

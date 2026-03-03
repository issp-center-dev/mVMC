.. include:: ../../bib/ref.txt

.. _HowToExpert:
             
Input files for Expert mode
===========================

In this section, detailed input files (\*.def) are explained. Input
files are categorized by the following six parts. The files that are
listed in parentheses correspond to the file made by vmcdry.out.

(1) List:

    No keyword (namelist.def): This file is a list of input file names
    with keywords. Each keywords is fixed, but file names are free to
    be determined.

(2) Basic parameters:
      
    **ModPara** (modpara.def): Set the parameters for basic parameters
    such as site number, electron number, Lanczos step *etc*.  
    **LocSpin (locspn.def)**: Set the location of local spin.

(3) Hamiltonian:
      
    Hamiltonian for mVMC is denoted by

    .. math::

       \begin{aligned}
       {\cal H}&={\cal H}_T+{\cal H}_U+{\cal H}_V+{\cal H}_H+{\cal H}_E+{\cal H}_P+{\cal H}_I,\\
       {\cal H}_T&=-\sum_{i, j}\sum_{\sigma_1, \sigma_2}t_{ij\sigma_1\sigma_2} c_{i\sigma_1}^{\dagger}c_{j\sigma_2},\\
       {\cal H}_U&=\sum_{i} U_i n_ {i \uparrow}n_{i \downarrow},\\
       {\cal H}_V&=\sum_{i,j} V_{ij}n_ {i}n_{j},\\
       {\cal H}_H&=-\sum_{i,j}J_{ij}^{\rm Hund} (n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}),\\
       {\cal H}_E&=\sum_{i,j}J_{ij}^{\rm Ex} (c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{j \downarrow}^{\dagger}c_{i  \downarrow}+c_ {i \downarrow}^{\dagger}c_{j\downarrow}c_{j \uparrow}^{\dagger}c_{i  \uparrow}),\\
       {\cal H}_P&=\sum_{i,j}J_{ij}^{\rm Pair} c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{i \downarrow}^{\dagger}c_{j  \downarrow},\\
       {\cal H}_I&=\sum_{i,j,k,l}\sum_{\sigma_1,\sigma_2, \sigma_3, \sigma_4}
       I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}, \end{aligned}

    as the format of interactions for electron system. Here, we
    define the charge density operator with spin :math:`\sigma` at
    site :math:`i` as
    :math:`n_{i \sigma}=c_{i\sigma}^{\dagger}c_{i\sigma}` and the total
    charge density operator at site :math:`i` as
    :math:`n_i=n_{i\uparrow}+n_{i\downarrow}`. Each parameters are
    specified by the following files, respectively;

    **Trans (trans.def)**: :math:`t_{ij\sigma_1\sigma_2}` in
    :math:`{\cal H}_T`,
          
    **CoulombIntra (coulombintra.def)**: :math:`U_i` in
    :math:`{\cal H}_U`,
          
    **CoulombInter (coulombinter.def)**: :math:`V_{ij}` in
    :math:`{\cal H}_V`,
          
    **Hund (hund.def)**: :math:`J_{ij}^{\rm Hund}` in
    :math:`{\cal H}_H`,
          
    **Exchange (exchange.def)**: :math:`J_{ij}^{\rm Ex}` in
    :math:`{\cal H}_E`,
          
    **PairHop**: :math:`J_{ij}^{\rm Pair}` in :math:`{\cal H}_P`,

    **InterAll**: :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}` in
    :math:`{\cal H}_I`.

(4) Variational parameters to be optimized:
      
    The variational parameters to be optimized are specified by using
    this categorized files. In mVMC, the variational wave function is
    given as

    .. math::

       \begin{aligned}
       |\psi \rangle &= {\cal N}_{General RBM} {\cal P}_G{\cal P}_J{\cal P}_{d-h}^{(2)}{\cal P}_{d-h}^{(4)}{\cal L}^S{\cal L}^K{\cal L}^P |\phi_{\rm pair} \rangle,\\
       {\cal P}_G&=\exp\left[ \sum_i g_i n_{i\uparrow} n_{i\downarrow} \right],\\
       {\cal P}_J&=\exp\left[\frac{1}{2} \sum_{i\neq j} v_{ij} (n_i-1)(n_j-1)\right],\\
       {\cal P}_{d-h}^{(2)}&= \exp \left[ \sum_t \sum_{n=0}^2 (\alpha_{2nt}^d \sum_{i}\xi_{i2nt}^d+\alpha_{2nt}^h \sum_{i}\xi_{i2nt}^h)\right],\\
       {\cal P}_{d-h}^{(4)}&= \exp \left[ \sum_t \sum_{n=0}^4 (\alpha_{4nt}^d \sum_{i}\xi_{i4nt}^d+\alpha_{4nt}^h \sum_{i}\xi_{i4nt}^h)\right],\\
       {\cal N}_{\rm General RBM}&= \exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] \prod_k^{N_h} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right],\\
       {\cal L}_S&=\frac{2S+1}{8 \pi^2}\int d\Omega P_s(\cos \beta) \hat{R}(\Omega),\\
       {\cal L}_K&=\frac{1}{N_s}\sum_{{\boldsymbol R}}e^{i {\boldsymbol K} \cdot{\boldsymbol R} } \hat{T}_{\boldsymbol R},\\
       {\cal L}_P&=\sum_{\alpha}p_{\alpha} \hat{G}_{\alpha},\end{aligned}

    where :math:`\Omega=(\alpha, \beta, \gamma)` is the Euler angle,
    :math:`\hat{R}(\Omega)` is the rotational operator, :math:`P_S(x)`
    is the :math:`S`-th polynomial, :math:`{\boldsymbol K}` is the momentum
    operator of the whole system and :math:`\hat{T}_{\boldsymbol R}` is the
    translational operators corresponding to the translational vector
    :math:`{\boldsymbol R}`, :math:`\hat{G}_{\alpha}` is the point group
    operator, and :math:`p_\alpha` is the parity operator,
    respectively. The details of :math:`{\cal P}_{d-h}^{(2)}` and
    :math:`{\cal P}_{d-h}^{(4)}` are shown in
    [Tahara2008_ ]. The one body part of the
    wavefunction is represented as the pair function of the real
    space:

    .. math::

       |\phi_{\rm pair} \rangle = \left[\sum_{i, j=1}^{N_s}
       \sum_{\sigma_1, \sigma_2}f_{i\sigma_1j \sigma_2}
       c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger} \right]^{N/2}|0 \rangle,

    where :math:`N` is the number of electrons and :math:`N_s` is the
    number of sites. The setting for optimizing variational parameters
    or not is given by the following files (the parameters for
    :math:`{\cal L}_S` are specified in the **ModPara** file).
             
    **Gutzwiller (gutzwilleridx.def)**: Set the target parameters
    :math:`g_i` in :math:`{\cal P}_G` to be optimized.
          
    **Jastrow (jastrowidx.def)**: Set the target parameters
    :math:`v_{ij}` in :math:`{\cal P}_J` to be optimized.
          
    **DH2**: Set the target 2-site doublon-holon correlation factor
    :math:`\alpha_{2nt}^{d(h)}` in :math:`{\cal P}_{d-h}^{(2)}` to be
    optimized.
          
    **DH4**: Set the target 4-site doublon-holon correlation factor
    :math:`\alpha_{4nt}^{d(h)}` in :math:`{\cal P}_{d-h}^{(4)}` to be
    optimized.

    **GeneralRBM_PhysLayer**: Set the target variational parameters :math:`a_{i\sigma}` in the RBM correlation factor :math:`{\cal N}_{\rm General RBM}`.

    **GeneralRBM_HiddenLayer**: Set the target variational parameters :math:`h_{k}` in the RBM correlation factor :math:`{\cal N}_{\rm General RBM}`.

    **GeneralRBM_PhysHidden**: Set the target variational parameters :math:`W_{i\sigma k}` in the RBM correlation factor :math:`{\cal N}_{\rm General RBM}`.

    **Orbital/OrbitalAntiParallel (orbitalidx.def)**: Set the pair
    orbital with anti-parallel spins :math:`f_{i\uparrow j\downarrow}`
    in :math:`|\phi_{\rm pair} \rangle` to be optimized.
    
    **OrbitalParallel**: Set the pair orbital with anti-parallel spins
    :math:`f_{i\sigma j\sigma}` in :math:`|\phi_{\rm pair} \rangle` to
    be optimized.  

    **OrbitalGeneral**: Set the pair orbital with anti-parallel spins
    :math:`f_{i\sigma j\sigma'}` in :math:`|\phi_{\rm pair} \rangle`
    to be optimized.

    **TransSym (qptransidx.def)**: Set the the momentum projection
    operators :math:`{\cal L}_K` and the lattice translational
    projection operators :math:`{\cal L}_P`.

(5) Initial variational parameters:
      
    Set the initial values of the variational parameters. When the
    keyword is not setting, the corresponding parameters are given by
    random values as default values.

    **InGutzwiller**: Set the initial values of :math:`g_i` in
    :math:`{\cal P}_G`.
          
    **InJastrow**: Set the initial values of :math:`v_{ij}` in
    :math:`{\cal P}_J`.
          
    **InDH2**: Set the initial values of :math:`\alpha_{2nt}^{d(h)}`
    in :math:`{\cal P}_{d-h}^{(2)}`.
    
    **InDH4**: Set the initial values of :math:`\alpha_{4nt}^{d(h)}`
    in :math:`{\cal P}_{d-h}^{(4)}`.

    **InGeneralRBM_PhysLayer**: Set the initial values of :math:`a_{i\sigma}` in the RBM correlation factor :math:`{\cal N}_{\rm General RBM}`.

    **InGeneralRBM_HiddenLayer**: Set the initial values of :math:`h_{k}` in the RBM correlation factor :math:`{\cal N}_{\rm General RBM}`.

    **InGeneralRBM_PhysHidden**: Set the initial values of :math:`W_{i\sigma k}` in the RBM correlation factor :math:`{\cal N}_{\rm General RBM}`.
    
    **InOrbital /InOrbitalAntiParallel**: Set the initial values of
    :math:`f_{i\uparrow j\downarrow}` in
    :math:`|\phi_{\rm pair} \rangle`.
      
    **InOrbitalParallel**: Set the initial values of
    :math:`f_{i\sigma j\sigma}` in :math:`|\phi_{\rm pair} \rangle`.
       
    **InOrbitalGeneral**: Set the initial values of
    :math:`f_{i\sigma j\sigma'}` in :math:`|\phi_{\rm pair} \rangle`.

(6) Output:
      
    **OneBodyG (greenone.def)**: Set the components of one-body green
    functions to output.

    **TwoBodyG (greentwo.def)**: Set the components of two-body green
    functions to output.

    **Twist (twist.def)**: Set the components of twist operators to
    output.

(7) Others:

    **Lattice (lattice.def)**: Set the position of each site.

.. _InputFileList:
    
List file for Input files (namelist.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file determines input filenames correlated with keywords. An
example of the file format is shown as follows.

::

    ModPara  modpara.def
    LocSpin  zlocspn.def
    Trans    ztransfer.def
    InterAll zinterall.def
    Orbital orbitalidx.def
    OneBodyG zcisajs.def
    TwoBodyG	zcisajscktaltdc.def

File format
^^^^^^^^^^^

[string01] [string02]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type

   **Description :** Select a word from keywords.

-  [ string02 ]

   **Type :** string-type

   **Description :** An input filename which is correlated with
   keywords.

User rules
^^^^^^^^^^

-  After setting keywords at [string 01], half-width state is needed for
   writing a filename. You can set the filename freely.

-  Keywords for input files are shown in Table [Table:Defs].

-  Essential keywords are "CalcMod", "ModPara" , "LocSpin", "Orbital"
   and "TransSym".

-  Keywords can be set in random order.

-  If keywords or filenames are incorrect, the program is terminated.

-  When the head of line is ``#``, the line is skipped.

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Keywords                 
     - Details for corresponding files
   * - ModPara :math:`^*`       
     - Parameters for calculation.
   * - LocSpin :math:`^*`       
     - Configurations of the local spins for Hamiltonian.
   * - Trans                    
     - Transfer and chemical potential for Hamiltonian.
   * - InterAll                 
     - Two-body interactions for Hamiltonian.
   * - CoulombIntra             
     - CoulombIntra interactions.
   * - CoulombInter             
     - CoulombInter interactions.
   * - Hund                     
     - Hund couplings.
   * - PairHop                  
     - Pair hopping couplings.
   * - Exchange                 
     - Exchange couplings.
   * - Gutzwiller               
     - Gutzwiller factors.
   * - Jastrow                  
     - Charge Jastrow factors.
   * - DH2                      
     - 2-site doublon-holon correlation factors.
   * - DH4                      
     - 4-site doublon-holon correlation factors.
   * - GeneralRBM_PhysLayer     
     - A part of the general RBM correlation factor including variational parameters in a physical layer :math:`\alpha_{i\sigma}`.
   * - GeneralRBM_HiddenLayer   
     - A part of the general RBM correlation factor including variational parameters in a hidden layer :math:`h_{k}`.
   * - GeneralRBM_PhysHidden    
     - A part of the general RBM correlation factor including variational parameters which connect physical and hidden layers :math:`W_{i\sigma k}`.
   * - Orbital :math:`^*`       
     - Pair orbital factors with anti-parallel spins :math:`f_{i\uparrow j\downarrow}`.
   * - OrbitalAntiParallel      
     - Pair orbital factors with anti-parallel spins :math:`f_{i\uparrow j\downarrow}`.
   * - OrbitalParallel          
     - Pair orbital factors with parallel spins :math:`f_{i\sigma j\sigma}`.
   * - OrbitalGeneral           
     - Pair orbital factors :math:`f_{i\sigma_1 j\sigma_2}`.
   * - TransSym :math:`^*`      
     - Translational and lattice symmetry operation.
   * - InGutzwiller             
     - Initial values of Gutzwiller factors.
   * - InJastrow                
     - Initial values of charge Jastrow factors.
   * - InDH2                    
     - Initial values of 2-site doublon-holon correlation factors.
   * - InDH4                    
     - Initial values of 4-site doublon-holon correlation factors.
   * - InGeneralRBM_PhysLayer   
     - Initial values of variational parameters of the general RBM correlation factors in a physical layer :math:`\alpha_{i\sigma}`.
   * - InGeneralRBM_HiddenLayer 
     - Initial values of variational parameters of the general RBM correlation factors in a hidden layer :math:`h_{k}`.
   * - InGeneralRBM_PhysHidden  
     - Initial values of variational parameters of the general RBM correlation factors which connect physical and hidden layers :math:`W_{i\sigma k}`.
   * - InOrbital                
     - Initial values of pair orbital factors :math:`f_{i\uparrow j\downarrow}`.
   * - InOrbitalAntiParallel    
     - Initial values of pair orbital factors :math:`f_{i\uparrow j\downarrow}`.
   * - InOrbitalParallel        
     - Initial values of pair orbital factors :math:`f_{i\sigma j\sigma}`.
   * - InOrbitalGeneral         
     - Initial values of pair orbital factors :math:`f_{i\sigma_1 j\sigma_2}`.
   * - OneBodyG                 
     - Output components for Green functions :math:`\langle c_{i\sigma}^{\dagger}c_{j\sigma}\rangle`
   * - TwoBodyG                 
     - Output components for Correlation functions :math:`\langle c_{i\sigma}^{\dagger}c_{j\sigma}c_{k\tau}^{\dagger}c_{l\tau}\rangle`
   * - Twist                    
     - Output components for Twist operators :math:`\langle \exp ( i 2\pi \sum_{i\sigma} \sum_{\mu=x,y,z} c^{(\alpha)\mu }_{i\sigma } \mu_{i} n_{i\sigma} ) \rangle`
   * - Lattice                  
     - Position of each site.

ModPara file (modpara.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~

This file determines parameters for calculation. An example of the
file format is shown as follows.

::

    --------------------
    Model_Parameters   0
    --------------------
    VMC_Cal_Parameters
    --------------------
    CDataFileHead  zvo
    CParaFileHead  zqp
    --------------------
    NVMCCalMode    0
    NLanczosMode   0
    --------------------
    NDataIdxStart  1
    NDataQtySmp    1
    --------------------
    Nsite          16
    Nelectron      8
    NSPGaussLeg    1
    NSPStot        0
    NMPTrans       1
    NSROptItrStep  1200
    NSROptItrSmp   100
    DSROptRedCut   0.001
    DSROptStaDel   0.02
    DSROptStepDt   0.02
    NVMCWarmUp     10
    NVMCInterval   1
    NVMCSample     1000
    NExUpdatePath  0
    RndSeed        11272
    NSplitSize     1
    NStore         1  

File format
^^^^^^^^^^^

-  Lines 1 - 5: Header

-  Line 6: [string01] [string02]

-  Line 7: [string03] [string04]

-  Line 8: Header

-  Lines 9 - : [string05] [int01] (or [double01])

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** Set a keyword for header of output files.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** Set a header of output files. For example, the
   output file of one-body green’s functions are named as
   **xxx\_cisajs.dat**, where **xxx** is [ string02 ].

-  [ string03 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** Set a keyword for header of output files for
   variational parameters.

-  [ string04 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** Set a header of output files for variational
   parameters. For example, the output file of optimized variational
   parameters are named as **xxx\_opt.dat**, where **xxx** is
   [ string04 ].

-  [ string05 ]

   **Type :** string-type

   **Description :** Select a word from keywords.

-  [ int01 ] ([double01])

   **Type :** int (double)-type (blank parameter not allowed)

   **Description :** A parameter which is correlated with a keyword.

User rules
^^^^^^^^^^

-  From Line 9: After setting keywords at [string 01], a half-width
   blank is needed for setting a parameter.

-  From Line 9: When the first character of the line is "-", the line is
   not read and skipped. 

Keywords and parameters
^^^^^^^^^^^^^^^^^^^^^^^

-  ``NVMCCalMode``

   **Type :** int-type (default value: 0)

   **Description :** [0] Optimization of variational parameters, [1]
   Calculation of one body and two body Green’s functions.

-  ``NLanczosMode``

   **Type :** int-type (default value: 0)

   **Description :** [0] Not using single Lanczos step, [1] Calculating
   energy by using Single Lanczos Step, [2] Calculating one body and two
   body Green’s functions by using Single Lanczos Step (Condition: The
   options 1 and 2 can be selected when ``NVMCCalMode`` = 1).

-  ``NDataIdxStart``

   **Type :** int-type (default value: 0)

   **Description :** An integer for numbering of output files. For
   ``NVMCCalMode`` = 0 , ``NDataIdxStart`` is added at the end of the
   output files. For ``NVMCCalMode`` = 1, the files are outputted with
   the number from ``NDataIdxStart`` to
   ``NDataIdxStart`` + ``NDataQtySmp``-1.

-  ``NDataQtySmp``

   **Type :** int-type (default value: 1)

   **Description :** The set number for outputted files (only used for
   ``NVMCCalMode`` = 1).

-  ``Nsite``

   **Type :** int-type (Positive integer)

   **Description :** The number of sites.

-  ``Nelectron``

   **Type :** int-type (Positive integer)

   **Description :** The number of electron pairs (the electron number
   is given by 2 ``Nelectron``).

-  ``Ncond``

   **Type :** int-type (greater than 0)

   **Description :** The number of conduction electrons.

-  ``2Sz``

   **Type :** int-type

   **Description :** The value of :math:`2S_z`. Since the electrons form
   pair, :math:`2S_z` must be even number.

-  ``NSPGaussLeg``

   **Type :** int-type (Positive integer, default value: 8)

   **Description :** The mesh number for the Gauss-legendre quadrature
   about :math:`\beta` integration (:math:`S_y` rotation) for the spin
   quantum-number projection in actual numerical calculation.

-  ``NSPStot``

   **Type :** int-type ( greater than 0, default value: 0)

   **Description :** The spin quantum-number.

-  ``NMPTrans``

   **Type :** int-type (default value: 1)

   **Description :** The absolute value gives the number of the momentum
   and lattice translational quantum-number projection. When the value
   is negative, the mode of anti-periodic condition turns on. The
   quantum-number projection is used from the top to ``NMPTrans`` with
   the specified weight indicated in ``TransSym`` file. In the case of
   not applying the projection, this value must be equal to 1.

-  ``NSROptItrStep``

   **Type :** int-type (Positive integer, default value: 1000)

   **Description :** The whole step number to optimize variational
   parameters by SR method. Only used for ``NVMCCalMode`` =0.

-  ``NSROptItrSmp``

   **Type :** int-type (Positive integer, default value:
   ``NSROptItrStep``/10)

   **Description :** In the ``NSROptItrStep`` step, the average values
   of the each variational parameters at the ``NSROptItrStep`` step are
   adopted as the optimized values. Only used for ``NVMCCalMode`` =0.

-  ``DSROptRedCut``

   **Type :** double-type (default value: 0.001)

   **Description :** The stabilized factor for the SR method by
   truncation of redundant directions corresponding to
   :math:`\varepsilon_{\rm wf}` in the ref.
   [Tahara2008_ ].

-  ``DSROptStaDel``

   **Type :** double-type (default value: 0.02)

   **Description :** The stabilized factor for the SR method by
   modifying diagonal elements in the overwrap matrix corresponding to
   :math:`\varepsilon` in the ref. [Tahara2008_ ].

-  ``DSROptStepDt``

   **Type :** double-type

   **Description :** The time step using in the SR method.

-  ``NSROptCGMaxIter``

   **Type :** int-type (default value: 0)

   **Description :** The maximum number of CG steps for the SR method.
   If this is zero or negative, CG steps will be run as many as the size
   of :math:`S` matrix at maximum. Only used for ``NSRCG``!=0.

-  ``DSROptCGTol``

   **Type :** double-type (default value: 1.0e-10)

   **Description :** The convergence condition of a CG step in the SR
   method. CG method runs until the root mean square of the residues
   becomes below this value. Only used for ``NSRCG``!=0.

-  ``NVMCWarmUp``

   **Type :** int-type (Positive integer, default value: 10)

   **Description :** Idling number for the Malkov chain Montecarlo
   Methods.

-  ``NVMCInterval``

   **Type :** int-type (Positive integer, default value: 1)

   **Description :** The interval step between samples. The local update
   will be performed ``Nsite`` × ``NVMCInterval`` times.

-  ``NVMCSample``

   **Type :** int-type (Positive integer, default value: 1000)

   **Description :** The sample numbers to calculate the expected
   values.

-  ``NExUpdatePath``

   **Type :** int-type (Positive integer)

   **Description :** The option for local update about exchange terms.
   0: not update, 1: update for electron system. For Spin system, the
   value must be 2.

-  ``RndSeed``

   **Type :** int-type

   **Description :** The initial seed of generating random number. For
   MPI parallelization, the initial seeds are given by ``RndSeed`` +my
   rank+1 at each ranks.

-  ``NSplitSize``

   **Type :** int-type (Positive integer, default value: 1)

   **Description :** The number of processes of MPI parallelization.

-  ``NStore``

   **Type :** int-type (0 or 1, default value: 1)

   **Description :** The option of applying matrix-matrix product to
   calculate expected values :math:`\langle O_k O_l \rangle` (0: off, 1:
   on). This speeds up calculation but increases the amount of memory
   usage from :math:`O(N_\text{p}^2)` to
   :math:`O(N_\text{p}^2) + O(N_\text{p}N_\text{MCS})`, where
   :math:`N_\text{p}` is the number of the variational parameters and
   :math:`N_\text{MCS}` is the number of Monte Carlo sampling.

-  ``NSRCG``

   **Type :** int-type (0 or 1, default value: 0)

   **Description :** The option of solving :math:`Sx=g` in the SR method
   without constructing :math:`S`
   matrix [NeuscammanUmrigarChan_ ]. (0: off, 1: on).
   This reduces the amount of memory usage from
   :math:`O(N_\text{p}^2) + O(N_\text{p}N_\text{MCS})` to
   :math:`O(N_\text{p}) + O(N_\text{p}N_\text{MCS})` when
   :math:`N_\text{p} > N_\text{MCS}`.

-  ``useDiagScale``

   **Type :** int-type (0 or 1, default value: 0)

   **Description :** The option of using the point Jacobi method (scaling by diagonal elements of :math:`S` matrix) when solving the linear equation :math:`Sx=g` in the SR method by CG method (0: off, 1: on, ``NSRCG`` must be 1).

-  ``NneuronGeneral``

   **Type :** int-type (default value: 0)

   **Description :** The number of neurons :math:`N_{\rm General RBM}` in the hidden layer of RBM.

LocSpin file (locspn.def)
~~~~~~~~~~~~~~~~~~~~~~~~~

This file determines sites with localized spins. An example of the
file format is shown as follows.

::

    ================================ 
    NlocalSpin     6  
    ================================ 
    ========i_0LocSpn_1IteElc ====== 
    ================================ 
        0      1
        1      0
        2      1
        3      0
        4      1
        5      0
        6      1
        7      0
        8      1
        9      0
       10      1
       11      0

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 -: [int02] [int03]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of localized spins. You
   can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of localized spins.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 ] :math:`<` ``Nsite``).

-  [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer for selecting an electron state whether
   localized spin or itinerant electron states (0: Itinerant electron
   state, 1: localized spin state with :math:`S=1/2`).

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of localized spins indicated by
   [ int03 ].

-  A program is terminated, when [ int02 ] is
   different from the total number of sites.

-  A program is terminated under the condition
   [ int02 ] :math:`<0` or
   ``Nsite`` :math:`\leq` [ int02 ].

Trans file (trans.def)
~~~~~~~~~~~~~~~~~~~~~~

The Hamiltonian for general one-body interactions

.. math::

   \begin{aligned}
   {\cal H}_{T} =-\sum_{ij\sigma_1\sigma_2}
   t_{ij\sigma_1\sigma_2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2},\end{aligned}

is added to the whole Hamiltonian by setting the parameters
:math:`t_{ij\sigma_1\sigma2}`. An example of the file format is shown
as follows.

::

    ======================== 
    NTransfer      24  
    ======================== 
    ========i_j_s_tijs====== 
    ======================== 
        0     0     2     0   1.000000  0.000000
        2     0     0     0   1.000000  0.000000
        0     1     2     1   1.000000  0.000000
        2     1     0     1   1.000000  0.000000
        2     0     4     0   1.000000  0.000000
        4     0     2     0   1.000000  0.000000
        2     1     4     1   1.000000  0.000000
        4     1     2     1   1.000000  0.000000
        4     0     6     0   1.000000  0.000000
        6     0     4     0   1.000000  0.000000
        4     1     6     1   1.000000  0.000000
        6     1     4     1   1.000000  0.000000
        6     0     8     0   1.000000  0.000000
        8     0     6     0   1.000000  0.000000
    …

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3-5: Header

-  Lines 6-: [int02]  [int03]  [int04]  [int05]  [double01]  [double02]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of transfer integrals.
   You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of transfer
   integrals.

-  [ int02 ], [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 ],  [ int04 ] :math:`<` ``Nsite``).

-  [ int03 ], [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a spin index,

   0: up-spin,

   1: down-spin.

-  [ double01 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** A value for a real part of
   :math:`t_{ij\sigma_1\sigma_2}`.

-  [ double02 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** A value for an imaginary part of
   :math:`t_{ij\sigma_1\sigma_2}`.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  Blank line is not allowed.

-  A program is terminated, when [ int01 ] is
   different from the total number of transfer integrals defined in this
   file.

-  A program is terminated, when
   [ int02 ]-[ int05 ] are out of
   range from the defined values.

-  Since Hamiltonian must be Hermitian, the following relation must be
   satisfied,
   :math:`t_{ij\sigma_1\sigma_2}=t_{ji\sigma_2\sigma_1}^{\dagger}`.

InterAll file
~~~~~~~~~~~~~

The Hamiltonian for general two-body interactions

.. math::

   \begin{aligned}
   {\cal H}_I =\sum_{i,j,k,l}\sum_{\sigma_1,\sigma_2, \sigma_3, \sigma_4}
   I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}
   c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}.
   \end{aligned}

is added to the whole Hamiltonian by setting the parameters
:math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`. An example of file
format is shown as follows.

::

    ====================== 
    NInterAll      36  
    ====================== 
    ========zInterAll===== 
    ====================== 
    0    0    0    1    1    1    1    0   0.50  0.0
    0    1    0    0    1    0    1    1   0.50  0.0
    0    0    0    0    1    0    1    0   0.25  0.0
    0    0    0    0    1    1    1    1  -0.25  0.0
    0    1    0    1    1    0    1    0  -0.25  0.0
    0    1    0    1    1    1    1    1   0.25  0.0
    2    0    2    1    3    1    3    0   0.50  0.0
    2    1    2    0    3    0    3    1   0.50  0.0
    2    0    2    0    3    0    3    0   0.25  0.0
    2    0    2    0    3    1    3    1  -0.25  0.0
    2    1    2    1    3    0    3    0  -0.25  0.0
    2    1    2    1    3    1    3    1   0.25  0.0
    4    0    4    1    5    1    5    0   0.50  0.0
    4    1    4    0    5    0    5    1   0.50  0.0
    4    0    4    0    5    0    5    0   0.25  0.0
    4    0    4    0    5    1    5    1  -0.25  0.0
    4    1    4    1    5    0    5    0  -0.25  0.0
    4    1    4    1    5    1    5    1   0.25  0.0
    ...

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 -:
   [int02] [int03] [int04] [int05] [int06] [int07] [int08] [int09] [double01] [double02]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of generalized two body
   interactions. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of generalized two
   body interactions.

-  [ int02 ], [ int04 ],
   [ int06 ], [ int08 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 :math:`],  [` int04 ], [ int06 ], [ int08 ] :math:`<` ``Nsite``).

-  [ int03 ], [ int05 ],
   [ int07 ], [ int09 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a spin index,

   0: up-spin,

   1: down-spin.

-  [ double01 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** A value for a real part of
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`.

-  [ double02 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** A value for an imaginary part of
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  Since Hamiltonian must be Hermitian, the following relation must be
   satisfied,
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}=I_{lkji\sigma_4\sigma_3\sigma_2\sigma_1}^{\dagger}`.

-  A program is terminated, when [ int01 ] is
   different from the total number of generalized two body interactions
   defined in this file.

-  A program is terminated, when
   [ int02 ]-[ int09 ] are out of
   range from the defined values.

CoulombIntra file (coulombintra.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Hamiltonian for the coulombintra interactions

.. math::

   {\cal H}_U =\sum_{i}U_i n_ {i \uparrow}n_{i \downarrow}

is added to the whole Hamiltonian by setting :math:`U_i`. An example of
the file format is shown as follows.

::

    ====================== 
    NCoulombIntra 6  
    ====================== 
    ========i_0LocSpn_1IteElc ====== 
    ====================== 
       0  4.000000
       1  4.000000
       2  4.000000
       3  4.000000
       4  4.000000
       5  4.000000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 -: [int02] [double01]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of on-site interactions.
   You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of on-site
   interactions.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 ] :math:`<` ``Nsite``).

-  [ double01 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** A value for :math:`U_i`.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of on-site interactions defined in
   this file.

-  A program is terminated, when [ int02 ] is out of
   range from the defined values.

CoulombInter file (coulombinter.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Hamiltonian for the coulombinter interactions

.. math::

   {\cal H}_V = \sum_{i,j}V_{ij} n_ {i} n_{j}

is added to the whole Hamiltonian by setting :math:`V_{ij}`. An example
of the file format is shown as follows.

::

    ====================== 
    NCoulombInter 6  
    ====================== 
    ========CoulombInter ====== 
    ====================== 
       0     1  1.0000
       1     2  1.0000
       2     3  1.0000
       3     4  1.0000
       4     5  1.0000
       5     0  1.0000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 -: [int02] [int03] [double01]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of off-site
   interactions. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of off-site
   interactions.

-  [ int02 ], [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 ], [ int03 ] :math:`<` ``Nsite``).

-  [ double01 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** A value for :math:`V_{ij}`.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of off-site interactions defined in
   this file.

-  A program is terminated, when either [ int02 ] or
   [ int03 ] are out of range from the defined values.

Hund file (hund.def)
~~~~~~~~~~~~~~~~~~~~

The Hamiltonian for Hund couplings

.. math::

   {\cal H}_H =-\sum_{i,j}J_{ij}^{\rm Hund} (n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow})

is added to the whole Hamiltonian by setting the parameters
:math:`J_{ij}^{\rm Hund}`. An example of the file format is shown as
follows.

::

    ====================== 
    NHund 6  
    ====================== 
    ========Hund ====== 
    ====================== 
       0     1 -0.250000
       1     2 -0.250000
       2     3 -0.250000
       3     4 -0.250000
       4     5 -0.250000
       5     0 -0.250000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 -: [int02] [int03] [double01]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of Hund couplings. You
   can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of Hund couplings.

-  [ int02 ], [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 ], [ int03 ] :math:`<` ``Nsite``).

-  [ double01 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** A value for :math:`J_{ij}^{\rm Hund}`.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of Hund couplings defined in this
   file.

-  A program is terminated, when either [ int02 ] or
   [ int03 ] are out of range from the defined values.

PairHop file
~~~~~~~~~~~~

The Hamiltonian for PairHop couplings

.. math::

   {\cal H}_P=\sum_{i,j}J_{ij}^{\rm Pair} (c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{i \downarrow}^{\dagger}c_{j  \downarrow}
   +c_{j \downarrow}^{\dagger}c_{i  \downarrow}c_ {j \uparrow}^{\dagger}c_{i\uparrow})

is added to the whole Hamiltonian by setting the parameters
:math:`J_{ij}^{\rm Pair}`. An example of the file format is shown as
follows.

::

    ====================== 
    NPairhop 6  
    ====================== 
    ========Pairhop ====== 
    ====================== 
       0     1  0.50000
       1     2  0.50000
       2     3  0.50000
       3     4  0.50000
       4     5  0.50000
       5     0  0.50000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 -: [int02] [int03] [double01]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of PairHop couplings.
   You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of PairHop
   couplings.

-  [ int02 ], [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 ], [ int03 ] :math:`<` ``Nsite``).

-  [ double01 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** A value for :math:`J_{ij}^{\rm Pair}`.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of PairHop couplings defined in this
   file.

-  A program is terminated, when either [ int02 ] or
   [ int03 ] are out of range from the defined values.

Exchange file (exchange.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Hamiltonian for exchange couplings

.. math::

   {\cal H}_E  =\sum_{i,j}J_{ij}^{\rm Ex}
   (c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{j \downarrow}^{\dagger}c_{i  \downarrow}
   +c_ {i \downarrow}^{\dagger}c_{j\downarrow}c_{j \uparrow}^{\dagger}c_{i  \uparrow})

is added to the whole Hamiltonian by setting :math:`J_{ij}^{\rm Ex}`.
Note that, for spin systems, the definition of the Exchange is different from that in HPhi.
An example of the file format is shown as follows.

::

    ====================== 
    NExchange 6  
    ====================== 
    ========Exchange ====== 
    ====================== 
       0     1  0.50000
       1     2  0.50000
       2     3  0.50000
       3     4  0.50000
       4     5  0.50000
       5     0  0.50000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3-5: Header

-  Lines 6-: [int02] [int03] [double01]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of Exchange couplings.
   You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of Exchange
   couplings.

-  [ int02 ], [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 ], [ int03 ] :math:`<` ``Nsite``).

-  [ double01 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** A value for :math:`J_{ij}^{\rm Ex}`.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of Exchange couplings defined in this
   file.

-  A program is terminated, when either [ int02 ] or
   [ int03 ] are out of range from the defined values.

Gutzwiller file (gutzwiller.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file sets the calculation conditions of Gutzwiller factors

.. math::

   {\cal P}_G=\exp\left[ \sum_i g_i n_{i\uparrow} n_{i\downarrow} \right].

A site number :math:`i` and the variational parameters :math:`g_i` are
specified. An example of the file format is shown as follows.

::

    ======================
    NGutzwillerIdx 2  
    ComplexType 0
    ====================== 
    ====================== 
       0     0
       1     0
       2     0
       3     1
    (continue...)
      12     1
      13     0
      14     0
      15     0
       0     1
       1     0

File format
^^^^^^^^^^^

In the following, we define the whole number of sites as :math:`N_s` and
variational parameters as :math:`N_g`, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+ :math:`N_s`): [int03] [int04]

-  Lines (6+ :math:`N_s`) - (5+ :math:`N_s` + :math:`N_g`):
   [int05] [int06]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters :math:`g_i`. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters :math:`g_i`.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters :math:`g_i`. You can freely give a name of
   the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters :math:`g_i` (0: double, 1: complex).

-  [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int03 ] :math:`<` ``Nsite``).

-  [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   :math:`g_i` (0 :math:`\leq` [ int04 ] :math:`<` [ int01]).

-  [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters
   (0 :math:`\leq` [ int05 ] :math:`<` [ int01]).

-  [ int06 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int05] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when components of variational parameters
   are double counted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int06 ] are out of range from the defined values.

Jastrow file (jastrow.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~

This file sets the calculation conditions of Jastrow factors

.. math::

   {\cal P}_J=\exp\left[\frac{1}{2} \sum_{i\neq j} v_{ij} (n_i-1)(n_j-1)\right]

Site numbers :math:`i` :math:`j`, and the variational parameters
:math:`v_{ij}` are specified. An example of the file format is shown as
follows.

::

    ======================
    NJastrowIdx 5  
    ComplexType 0
    ====================== 
    ======================
       0     1     0 
       0     2     1 
       0     3     0 
     (continue...)
       0    1 
       1    1 
       2    1 
       3    1 
       4    1 

File format
^^^^^^^^^^^

In the following, we define the total number of sites as :math:`N_s` and
variational parameters as :math:`N_j`, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+ :math:`N_s\times (N_s-1))`: [int03] [int04] [int05]

-  Lines (6+ :math:`N_s\times (N_s-1)`) -
   (5+ :math:`N_s\times (N_s-1)` + :math:`N_j`): [int06] [int07]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters :math:`v_{ij}`. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters :math:`v_{ij}`.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters :math:`v_{ij}`. You can freely give a name
   of the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters :math:`v_{ij}` (0: double, 1: complex).

-  [ int03 ], [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int03 ], [ int04 ] :math:`<` ``Nsite``).

-  [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   :math:`v_{ij}` (0 :math:`\leq` [ int05 ] :math:`<` [ int01]).

-  [ int06 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters
   (0 :math:`\leq` [ int06 ] :math:`<` [ int01]).

-  [ int07 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int06] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int07 ] are out of range from the defined values.

DH2 file
~~~~~~~~

This file sets the calculation conditions of 2-site doublon-holon
correlation factors

.. math::

   {\cal P}_{d-h}^{(2)}= \exp \left[ \sum_t \sum_{n=0}^2
   (\alpha_{2nt}^d \sum_{i}\xi_{i2nt}^d+\alpha_{2nt}^h \sum_{i}\xi_{i2nt}^h)\right].

A site number :math:`i`, the two sites around :math:`i` site and the
variational parameters :math:`\alpha_{2nt}^{d(h)}` which have :math:`t`
kinds at each sites are specified. The details of the parameters
:math:`\alpha_{2nt}^{d(h)}` and the operator :math:`\xi_{i2nt}^{d(h)}`
are shown in ref. [Tahara2008_ ]. An example of the file
format is shown as follows.

::

    ====================================
    NDoublonHolon2siteIdx 2  
    ComplexType 0
    ====================================
    ====================================
       0     5   15    0
       0    13    7    1
     (continue...)
      15     8    2    1
       0     1 
    (continue...)
      11     1 

File format
^^^^^^^^^^^

In the following, we define the total number of sites as :math:`N_s` and
variational parameters as :math:`N_{\rm dh2}`, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+ :math:`N_s\times N_{\rm dh2})`:
   [int03] [int04] [int05] [int06]

-  Lines (6+ :math:`N_s\times N_{\rm dh2}`) -
   (5+ :math:`(N_s+6) \times N_{\rm dh2})`: [int07] [int08]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters. You can freely give a name of the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters (0: double, 1: complex).

-  [ int03 ], [ int04 ],
   [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int03 ], [ int04 ], [ int05 ] :math:`<` ``Nsite``).

-  [ int06 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   (0 :math:`\leq` [ int06 ] :math:`<` [ int01]).

-  [ int07 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters.
   The value is :math:`(2n+s)\times` [int01]:math:`+t`, where
   :math:`n`, :math:`s` and :math:`t` are given by the following
   relation:

   -  :math:`n`: The number of doublon (holon) around the center site
      (0, 1, 2),

   -  :math:`s`: When the center is doublon (holon), s=0 (1),

   -  :math:`t`: The kind of variational parameters (0, :math:`\cdots`
      [int1]-1).

-  [ int08 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int07] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int08 ] are out of range from the defined values.

DH4 file
~~~~~~~~

This file sets the calculation conditions of 4-site doublon-holon
correlation factors

.. math::

   {\cal P}_{d-h}^{(4)}= \exp \left[ \sum_t \sum_{n=0}^4
   (\alpha_{4nt}^d \sum_{i}\xi_{i4nt}^d+\alpha_{4nt}^h \sum_{i}\xi_{i4nt}^h)\right]

A site number :math:`i`, the four sites around :math:`i` site and the
variational parameters :math:`\alpha_{4nt}^{d(h)}` which have :math:`t`
kinds at each sites are specified. The details of the parameters
:math:`\alpha_{4nt}^{d(h)}` and the operator :math:`\xi_{i4nt}^{d(h)}`
are shown in ref. [Tahara2008_ ]. An example of the file
format is shown as follows.

::

    ====================================
    NDoublonHolon4siteIdx 1  
    ComplexType 0
    ====================================
    ====================================
       0     1    3    4   12    0
       1     2    0    5   13    0
     (continue...)
      15    12   14    3   11    0
       0     1 
    (continue...)
       9     1 

File format
^^^^^^^^^^^

In the following, we define the total number of sites as :math:`N_s` and
variational parameters as :math:`N_{\rm dh4}`, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+ :math:`N_s\times N_{\rm dh4}`):
   [int03] [int04] [int05] [int06] [int07] [int08]

-  Lines (6+ :math:`N_s\times N_{\rm dh4}`) -
   (5+ :math:`(N_s+10) \times N_{\rm dh4}`): [int09] [int10]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters. You can freely give a name of the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters (0: double, 1: complex).

-  [ int03 ], [ int04 ],
   [ int05 ], [ int06 ],
   [ int07 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int03 :math:`], \cdots,  [` int07 ] :math:`<` ``Nsite``).

-  [ int08 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   (0 :math:`\leq` [ int08 ] :math:`<` [ int01]).

-  [ int09 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters.
   The value is :math:`(2n+s)\times` [int01]:math:`+t`, where
   :math:`n`, :math:`s` and :math:`t` are given by the following
   relation:

   -  :math:`n`: The number of doublon (holon) around the center site
      (0, 1, 2, 3, 4),

   -  :math:`s`: When the center is doublon (holon), s=0 (1),

   -  :math:`t`: The kind of variational parameters (0, :math:`\cdots`
      [int1]-1).

-  [ int10 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int09] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when components of variational parameters
   are double counted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int10 ] are out of range from the defined values.


GeneralRBM_PhysLayer file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the general RBM correlation factors, 

.. math::

   {\cal N}_{\rm General RBM}= \exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] \prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right]


this file sets the calculation conditions of  :math:`\exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] ` .
A site index :math:`i`, a spin index :math:`\sigma` and the variational parameters
:math:`a_{i \sigma}` are specified. An example of the file format is shown as
follows.

::

    --------------------
    NRBM_PhysLayerIdx	1
    ComplexType	1
    i s RBM_PhysLayer_Idx
    --------------------
     0	0	 0
     0	1	 0
     1	0	 0
     1	1	 0
     (continue...)
     0    1 

File format
^^^^^^^^^^^

In the following, we define the total number of sites as :math:`N_s` and
variational parameters as :math:`N_{v}`, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+ :math:`2N_s`): [int03] [int04] [int05] 

-  Lines (6+ :math:`2N_s` ) -
   (5+  :math:`2N_s` + :math:`N_v`): [int06] [int07]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters :math:`a_{i\sigma}`. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters :math:`a_{i\sigma}`.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters :math:`a_{i\sigma}`. You can freely give a name
   of the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters :math:`a_{i\sigma}` (0: double, 1: complex).

-  [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int03 ] :math:`<` ``Nsite``).

-  [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a spin index
   ([int04] = 0 or 1).

-  [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   :math:`v_{ij}` (0 :math:`\leq` [ int05 ] :math:`<` [ int01]).

-  [ int06 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters
   (0 :math:`\leq` [ int06 ] :math:`<` [ int01]).

-  [ int07 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int06] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  This function is a beta version. Please use it carefully. The file format and the implementation may be changed in the future release.
  
-  This function is available only for the VMC calculation with ``ComplexType=1`` and ``Orbital``. This function is not available for the Power Lanczos calculations.

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int07 ] are out of range from the defined values.


GeneralRBM_HiddenLayer file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the general RBM correlation factors, 

.. math::

   {\cal N}_{\rm General RBM}= \exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] \prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right]


this file sets the calculation conditions of :math:`b_k` in :math:`\prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right] `.
A hidden neuron index :math:`k` and the variational parameters
:math:`h_{k}` are specified. An example of the file format is shown as
follows.

::

    --------------------
    NRBM_HiddenLayerIdx	2
    ComplexType	1
    k RBM_HiddenLayer_Idx
    --------------------
     0	0
     1	0
     2	0
     3	0
     (continue...)
     0    1 
     1    1 

File format
^^^^^^^^^^^

In the following, we define the total number of neurons as :math:`N_{\rm neuronGeneral}` and
variational parameters as :math:`N_{v}`, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+:math:`N_{\rm neuronGeneral}`): [int03] [int04] 

-  Lines (6+ :math:`N_{\rm neuronGeneral}` ) -
   (5+ :math:`N_{\rm neuronGeneral}` + :math:`N_v`): [int05] [int06]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters :math:`h_{k}`. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters :math:`h_{k}`.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters :math:`h_{k}`. You can freely give a name
   of the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters :math:`h_{k}` (0: double, 1: complex).

-  [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a neuron index
   (0 :math:`\leq` [ int03 ] :math:`<` ``NneuronGeneral``).

-  [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   :math:`h_{k}` (0 :math:`\leq` [ int05 ] :math:`<` [ int01]).

-  [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters
   (0 :math:`\leq` [ int06 ] :math:`<` [ int01]).

-  [ int06 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int06] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  This function is a beta version. Please use it carefully. The file format and the implementation may be changed in the future release.
  
-  This function is available only for the VMC calculation with ``ComplexType=1`` and ``Orbital``. This function is not available for the Power Lanczos calculations.

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int06 ] are out of range from the defined values.


GeneralRBM_PhysHidden file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the general RBM correlation factors, 

.. math::

   {\cal N}_{\rm General RBM}= \exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] \prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right]


this file sets the calculation conditions of :math:`W_{i\sigma k}` in :math:`\prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right] `.
A site index :math:`i`, a spin index :math:`\sigma`, a hidden neuron index :math:`k` and the variational parameters
:math:`W_{i\sigma k}` are specified. An example of the file format is shown as
follows.

::

    --------------------
    NRBM_HiddenLayerIdx	32
    ComplexType	1
    i s k RBM_PhysHidden_Idx
    --------------------
     0	0   0   0
     0	1   0   1
     1	0   0   2
     1	1   0   3
     2	0   0   4
     2	1   0   5
     (continue...)
     0    1 
     1    1
     (continue...) 

File format
^^^^^^^^^^^

In the following, we define the total number of sites as :math:`N_s`, the total number of neurons as :math:`N_{\rm neuronGeneral}` and
variational parameters as :math:`N_{v}`.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+:math:`2 N_s N_{\rm neuronGeneral}`): [int03] [int04] [int05] [int06]

-  Lines (6+ :math:`2 N_s N_{\rm neuronGeneral}` ) -
   (5+ :math:`2 N_s N_{\rm neuronGeneral}` + :math:`N_v`): [int07] [int08]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters :math:`W_{i \sigma k}`. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters :math:`W_{i \sigma k}`.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters :math:`W_{i \sigma k}`. You can freely give a name
   of the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters :math:`W_{i \sigma k}` (0: double, 1: complex).

-  [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int03 ] :math:`<` ``Nsite``).

-  [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a spin index
   ([ int03 ] = 0 or 1  ).

-  [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a neuron index
   (0 :math:`\leq` [ int03 ] :math:`<` ``NneuronGeneral``).

-  [ int06 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   :math:`W_{i \sigma k}` (0 :math:`\leq` [ int05 ] :math:`<` [ int01]).

-  [ int07 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters
   (0 :math:`\leq` [ int06 ] :math:`<` [ int01]).

-  [ int08 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int06] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  This function is a beta version. Please use it carefully. The file format and the implementation may be changed in the future release.
  
-  This function is available only for the VMC calculation with ``ComplexType=1`` and ``Orbital``. This function is not available for the Power Lanczos calculations.

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int06 ] are out of range from the defined values.

Orbital/OrbitalAntiParallel file (orbitalidx.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file sets the calculation conditions of pair orbitals

.. math::

   |\phi_{\rm pair} \rangle = \left[\sum_{i, j=1}^{N_s}
   f_{ij}c_{i\uparrow}^{\dagger}c_{j\downarrow}^{\dagger} \right]^{N/2}|0 \rangle.


Site numbers :math:`i, j` and the variational parameters :math:`f_{ij}`
are indicated. An example of the file format is shown as follows.

::

    ====================================
    NOrbitalIdx 64   
    ComplexType 0
    ====================================
    ====================================
       0     0     0 
       0     1     1 
       0     2     2 
       0     3     3 
     (continue...)
      15     9    62 
      15    10    63 
       0    1 
       1    1 
    (continue...)
      62    1 
      63    1 

File format
^^^^^^^^^^^

In the following, we define the total number of sites as :math:`N_s` and
variational parameters as :math:`N_{\rm o}`, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+ :math:`N_s^2`): [int03] [int04] [int05] [int06]

-  Lines (6+ :math:`N_s^2` )- (5+ :math:`N_s^2+N_{\rm o}`):
   [int06] [int07]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters. You can freely give a name of the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters (0: double, 1: complex).

-  [ int03 ], [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int03 ], [ int04 ] :math:`<` ``Nsite``).

-  [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   (0 :math:`\leq` [ int05 ] :math:`<` [ int01]).

-  [ int06 ]

   **Type :** int-type

   **Description :** When the mode of the anti-periodic condition turns
   on (the mode turns on when the value of ``NMPTrans`` in ``ModPara``
   file is negative), the sign of :math:`f_{ij}` is specified by setting
   [ int06 :math:`]=\pm1`. This term can be omitted when the
   mode of the anti-periodic condition is off.

-  [ int07 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters
   (0 :math:`\leq` [ int06 ] :math:`<` [ int01]).

-  [ int08 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int06] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int09 ] are out of range from the defined values.

OrbitalParallel file
~~~~~~~~~~~~~~~~~~~~

This file sets the calculation conditions of pair orbitals

.. math::

   |\phi_{\rm pair} \rangle = \left[\sum_{i, j=1}^{N_s} \sum_{\sigma}
   f_{i\sigma j\sigma}c_{i\sigma}^{\dagger}c_{j\sigma}^{\dagger} \right]^{N/2}|0 \rangle.

Site numbers :math:`i, j`, the spin index :math:`\sigma` and the
variational parameters :math:`f_{i\sigma j\sigma}` are indicated. The
indexes of :math:`f_{i\sigma j\sigma}` must satisfy the condition
:math:`i < j `, where :math:`\sigma = 0` or :math:`1` and process will
terminate when the condition is broken. An example of the file format is
shown as follows.

::

    ====================================
    NOrbitalIdx 120  
    ComplexType 0
    ====================================
    ====================================
       0     1     0
       0     2     1 
       0     3     2 
     (continue...)
      15    13    118
      15    14    119
       0    1 
       1    1 
    (continue...)
      118    1 
      119    1 

File format
^^^^^^^^^^^

In the following, we define the total number of sites as :math:`N_s` and
variational parameters as :math:`N_{\rm o}`, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+ :math:`N_s*(N_s-1)/2`): [int03] [int04] [int05] [int06]

-  Lines (6+ :math:`N_s*(N_s-1)/2` )-
   (5+ :math:`N_s*(N_s-1)/2+N_{\rm o}`): [int06] [int07]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters. You can freely give a name of the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters (0: double, 1: complex).

-  [ int03 ], [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int03 ], [ int04 ] :math:`<` ``Nsite``).

-  [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   (0 :math:`\leq` [ int05 ] :math:`<` [ int01]).

-  [ int06 ]

   **Type :** int-type

   **Description :** When the mode of the anti-periodic condition turns
   on (the mode turns on when the value of ``NMPTrans`` in ``ModPara``
   file is negative), the sign of :math:`f_{ij}` is specified by setting
   [ int06 :math:`]=\pm1`. This term can be omitted when the
   mode of the anti-periodic condition is off.

-  [ int07 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters
   (0 :math:`\leq` [ int06 ] :math:`<` [ int01]).

-  [ int08 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int06] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int09 ] are out of range from the defined values.

OrbitalGeneral file
~~~~~~~~~~~~~~~~~~~

This file sets the calculation conditions of pair orbitals

.. math::

   |\phi_{\rm pair} \rangle = \left[\sum_{i, j=1}^{N_s} \sum_{\sigma_1, \sigma_2}
   f_{i\sigma_1j\sigma_2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger} \right]^{N/2}|0 \rangle.

Site numbers :math:`i, j`, spin indexes :math:`\sigma_1, \sigma_2` and
the variational parameters :math:`f_{i\sigma_1j\sigma_2}` are indicated.
The indexes of :math:`f_{i\sigma_1j\sigma_2}` must satisfy the condition
:math:`i+\sigma_1 N_s < j+\sigma_2 N_s`, where :math:`\sigma_i = 0` or
:math:`1` and process will terminate when the condition is broken. An
example of the file format is shown as follows.

12.5cm

::

    ====================================
    NOrbitalIdx 255
    ComplexType 0
    ====================================
    ====================================
       0  0  0  1  0
       0  0  1  1  1 
     (continue...)
      14  0  15  1  253
      15  0  15  1  254
       0    1 
       1    1 
    (continue...)
      253   1 
      254   1 

File format
^^^^^^^^^^^

In the following, we define the total number of sites as :math:`N_s` and
variational parameters as :math:`N_{\rm o}`, respectively. A total
number of variational parameters :math:`N_p` is given by :math:`N_s^2`,
:math:`2 N_s^2 -N_s` in :math:`Sz` conserved and :math:`S_z`
unconserved system, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Line 3: [string02] [int02]

-  Lines 4 - 5: Header

-  Lines 6 - (5+ :math:`N_p`):
   [int03] [int04] [int05] [int06] [int07] [int08]

-  Lines (6+ :math:`N_p`) - (5+ :math:`N_p+N_{\rm o}`): [int09] [int10]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters.

-  [ string02 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for indicating the double or complex type
   of variational parameters. You can freely give a name of the keyword.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer indicates the double or complex type of
   variational parameters (0: double, 1: complex).

-  [ int03 ], [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int03 ], [ int04 ] :math:`<` ``Nsite``).

-  [ int04 ], [ int06 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a spin index (:math:`0:\uparrow`
   spin, :math:`1:\downarrow` spin).

-  [ int07 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   (0 :math:`\leq` [ int07 ] :math:`<` [ int01]).

-  [ int08 ]

   **Type :** int-type

   **Description :** When the mode of the anti-periodic condition turns
   on (the mode turns on when the value of ``NMPTrans`` in ``ModPara``
   file is negative), the sign of :math:`f_{i\sigma_1j\sigma_2}` is
   specified by setting [ int08 :math:`]=\pm1`. This term can
   be omitted when the mode of the anti-periodic condition is off.

-  [ int09 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of variational parameters
   (0 :math:`\leq` [ int09 ] :math:`<` [ int01]).

-  [ int10 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer to select the target of variational
   parameters indicated at [int09] to be optimized or not (0: not
   optimize, 1: optimize).

User rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int10 ] are out of range from the defined values.

TransSym file (qptransidx.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file sets the weight and corresponding site numbers of momentum
projection
:math:`{\cal L}_K=\frac{1}{N_s}\sum_{{\boldsymbol R}}e^{i {\boldsymbol K} \cdot{\boldsymbol R} } \hat{T}_{\boldsymbol R}`
and lattice translational projection
:math:`{\cal L}_P=\sum_{\alpha}p_{\alpha} \hat{G}_{\alpha}`. The
patterns of projection are indicated by :math:`(\alpha, {\boldsymbol R})`. We
note that the weight must be equal to :math:`1.0` when the projection is
not done. An example of the file format is shown as follows.

12.5cm

::

    ====================================
    NQPTrans 4  
    ====================================
    == TrIdx_TrWeight_and_TrIdx_i_xi  ==
    ====================================
       0  1.000000  0.000000
       1  1.000000  0.000000
       2  1.000000  0.000000
       3  1.000000  0.000000
       0     0    0
     (continue...)
       3    12    1
       3    13    2 

File format
^^^^^^^^^^^

In the following, we define the total number of sites as :math:`N_s` and
projection patterns as :math:`N_{\rm TS}`, respectively.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 - (5+ :math:`N_{\rm TS})`: [int02] [double01] [double02]

-  Lines (6+ :math:`N_{\rm TS}`) - (5+ :math:`(N_s+1) \times N_{\rm TS}`):
   [int03] [int04] [int05] [int06]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of projection patterns.
   You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of projection
   patterns.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving the projection pattern
   :math:`(\alpha, {\boldsymbol R})`
   (0 :math:`\leq` [ int02 ] :math:`<` [ int01]).

-  [ double01 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** The weight
   Real part of :math:`p_{\alpha}\cos ({\boldsymbol K}\cdot {\boldsymbol R})` of the projection
   pattern :math:`(\alpha, {\boldsymbol R})`.

-  [ double02 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** The weight
   Imaginary part of :math:`p_{\alpha}\cos ({\boldsymbol K}\cdot {\boldsymbol R})` of the projection
   pattern :math:`(\alpha, {\boldsymbol R})`.


-  [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving kinds of the projection pattern
   :math:`(\alpha, {\boldsymbol R})`
   (0 :math:`\leq` [ int03 ] :math:`<` [ int01]).

-  [ int04 ], [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int04 ], [ int05 ] :math:`<` ``Nsite``).
   The site number [ int05 ] is given by applying the
   translation and point group transformation indicated by
   [ int03 ] to the site [ int04 ].

-  [ int06 ]

   **Type :** int-type

   **Description :** When the mode of the anti-periodic condition turns
   on (the mode turns on when the value of ``NMPTrans`` in ``ModPara``
   file is negative), the sign of the translational operator is
   specified by setting [ int06 ] :math:`=\pm1`. This term can
   be omitted when the mode of the anti-periodic condition is off.

User rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of projection patterns defined in
   this file.

-  A program is terminated, when [ int02 ] -
   [ int06 ] are out of range from the defined values.

.. _InputParam:
   
Files to set initial values of variational parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file sets the initial values of variational parameters. The kinds
of variational parameters are specified by setting the following
keywords in ``List`` file (namelist.def):
``InGutzwiller``, ``InJastrow``, ``InDH2``, ``InDH4``, ``InGeneralRBM_PhysLayer``, ``InGeneralRBM_HiddenLayer``, ``InGeneralRBM_PhysHidden``, ``InOrbital``,
``InOrbitalAntiParallel``, ``InOrbitalParallel``,
``InOrbitalGeneral``.
The file format is common and an example of the ``InJastrow`` file is
shown as follows.

::

    ======================
    NJastrowIdx  28
    ====================== 
    == i_j_JastrowIdx  ===
    ====================== 
    0 -8.909963465082626488e-02  0.000000000000000000e+00
    1  5.521681211878626955e-02  0.000000000000000000e+00
    (continue...)
    27 -9.017586139930480749e-02  0.000000000000000000e+00

File format
^^^^^^^^^^^

In the following, we define the total number of variational parameters
as :math:`N_v`.

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 - (5+ :math:`N_v`): [int03] [double01] [double02]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of variational
   parameters. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of variational
   parameters.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer setting kinds of variational parameters
   (0 :math:`\leq` [ int02 ] :math:`<` [ int01]).

-  [ double01 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** The real part of the variational parameter
   indicated by [int01].

-  [ double02 ]

   **Type :** double-type

   **Description :** The imaginary part of the variational parameter
   indicated by [int01].

User rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of variational parameters defined in
   this file.

OneBodyG file (greenone.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file determines the target components to calculate and output
one-body Green’s function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`. An example
of file format is shown as follows.

::

    ===============================
    NCisAjs         24
    ===============================
    ======== Green functions ======
    ===============================
        0     0     0     0
        0     1     0     1
        1     0     1     0
        1     1     1     1
        2     0     2     0
        2     1     2     1
        3     0     3     0
        3     1     3     1
        4     0     4     0
        4     1     4     1
        5     0     5     0
        5     1     5     1
        6     0     6     0
        6     1     6     1
        7     0     7     0
        7     1     7     1
        8     0     8     0
        8     1     8     1
        9     0     9     0
        9     1     9     1
       10     0    10     0
       10     1    10     1
       11     0    11     0
       11     1    11     1

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 -: [int02]  [int03]  [int04]  [int05]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of one-body Green’s
   functions. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of one-body Green’s
   functions.

-  [ int02 ], [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 ], [ int04 ] :math:`<` ``Nsite``).

-  [ int03 ], [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a spin index,

   0: up-spin,

   1: down-spin.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of one-body Green’s functions defined
   in this file.

-  A program is terminated, when
   [ int02 ]-[ int05 ] are out of
   range from the defined values.

TwoBodyG file (greentwo.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file determines the target components to calculate and output
two-body Green’s function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`.
For Spin, the condition :math:`i=j` and :math:`k=l` must be satisfied.
An example of file format is shown as follows.

::

    =============================================
    NCisAjsCktAltDC        576
    =============================================
    ======== Green functions for Sq AND Nq ======
    =============================================
        0     0     0     0     0     0     0     0
        0     0     0     0     0     1     0     1
        0     0     0     0     1     0     1     0
        0     0     0     0     1     1     1     1
        0     0     0     0     2     0     2     0
        0     0     0     0     2     1     2     1
        0     0     0     0     3     0     3     0
        0     0     0     0     3     1     3     1
        0     0     0     0     4     0     4     0
        0     0     0     0     4     1     4     1
        0     0     0     0     5     0     5     0
        0     0     0     0     5     1     5     1
        0     0     0     0     6     0     6     0
        0     0     0     0     6     1     6     1
        0     0     0     0     7     0     7     0
        0     0     0     0     7     1     7     1
        0     0     0     0     8     0     8     0
        0     0     0     0     8     1     8     1
        0     0     0     0     9     0     9     0
        0     0     0     0     9     1     9     1
        0     0     0     0    10     0    10     0
        0     0     0     0    10     1    10     1
        0     0     0     0    11     0    11     0
        0     0     0     0    11     1    11     1
        0     1     0     1     0     0     0     0
        ...

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 -:
   [int02]  [int03]  [int04]  [int05]  [int06]  [int07]  [int08]  [int09]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of two-body Green’s
   functions. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving total number of two-body Green’s
   functions.

-  [ int02 ],
   [ int04 ], [ int06 ],
   [ int08 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving a site index
   (0 :math:`\leq` [ int02 ], [ int04 ], [ int06 ], [ int08 ] :math:`<` ``Nsite``).

-  [ int03 ],
   [ int05 ], [ int07 ],
   [ int09 ]

   **Type :** int-type (blank parameter not allowed)

   | **Description :** An integer giving a spin index,
   | 0: up-spin,
   | 1: down-spin.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated, when [ int01 ] is
   different from the total number of two-body Green’s functions defined
   in this file.

-  A program is terminated, when
   [ int02 ]-[ int09 ] are out of
   range from the defined values.


Twist file (twist.def)
~~~~~~~~~~~~~~~~~~~~~~

This file sets the target components to calculate and output the
Twist operator
:math:`P^{(\alpha)} = \langle \exp ( i 2\pi \sum_{i\sigma} \sum_{\mu=x,y,z} c^{(\alpha)\mu }_{i\sigma } \mu_{i} n_{i\sigma} ) \rangle`.
Position operator :math:`\mu_i` is defined in Lattice definition file (lattice.def).
An example of file format is shown as follows. 

::

    --------------------
    NTwist	2
    --------------------
    idx_site_s_x_y_z
    --------------------
    0	0	0	0.0	0.0	0.0
    0	1	0	0.0	0.0	0.0
    0	2	0	0.0	0.0	0.0
    0	3	0	0.0	0.0	0.0
    0	4	0	0.0	0.0	0.0
    0	5	0	0.0	0.0	0.0
    0	0	1	0.3333333333	0.0	0.0
    0	1	1	0.3333333333	0.0	0.0
    0	2	1	0.3333333333	0.0	0.0
    0	3	1	0.3333333333	0.0	0.0
    0	4	1	0.3333333333	0.0	0.0
    0	5	1	0.3333333333	0.0	0.0
    1	0	0	0.3333333333	0.0	0.0
    1	1	0	0.3333333333	0.0	0.0
        …

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3 - 5: Header

-  Lines 6 -:
   [int02]  [int03]  [int04]  [double01]  [double02]  [double03]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword for total number of Twist operators. You can freely give a name of the keyword.

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving the total number of Twist operators.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :**
   An integer giving a Twist operator index :math:`\alpha`. You can specify it as an integer from 0 to [ int01 ]-1.

-  [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :**
   An integer giving a site index (0 :math:`\leq` [ int03 ] :math:`<` ``Nsite``).

-  [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :**
   An integer giving a spin index,
   0: up-spin,
   1: down-spin.

-  [ double01 ], [ double02 ],
   [ double03 ]

   **Type :** double-type (blank parameter not allowed)

   **Description :** :math:`c^{(\alpha)\mu}_{i\sigma}` for :math:`\mu = x, y, z`, :math:`c^{(\alpha)\mu}_{i\sigma}`. You need to specify the :math:`x` direction component as [ double01 ], the :math:`y` direction component as [ double02 ], and the :math:`z` direction component as [ double03 ].

Use rules
^^^^^^^^^^

-  Headers cannot be omitted.

-  All the combinations of sites and spins must be specified to designate each Twist operator.

Lattice file (lattice.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~

This file defines the position operator :math:`\mu_i` corresponding to site number :math:`i` and the orbital number in the unit cell.
An example of file format is shown as follows.

::

    --------------------
    NLattice  4 4 4 2
    --------------------
    i_x_y_z_orb
    --------------------
    0 0 0 0 0
    1 0 0 0 1
    2 1 0 0 0
    3 1 0 0 1
    4 2 0 0 0
    5 2 0 0 1
        …

File format
^^^^^^^^^^^^

- Line 1: Header

- Line 2: [string01] [int01] [int02]  [int03]  [int04]

- Lines 3 - 5: Header

- Lines 6 -: [int05]  [int06]  [int07]  [int08]  [int09]

Parameters
^^^^^^^^^^

-  [ string01 ]

   **Type :** string-type (blank parameter not allowed)

   **Description :** A keyword name (optional).

-  [ int01 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving the maximum value in the :math:`x` direction for the position operator :math:`\mu`.

-  [ int02 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving the maximum value in the :math:`y` direction for the position operator :math:`\mu`.

-  [ int03 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving the maximum value in the :math:`z` direction for the position operator :math:`\mu`.

-  [ int04 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving the number of orbitals in the unit cell.

-  [ int05 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :**
   An integer giving a site index (0 :math:`\leq` [ int05 ] :math:`<` ``Nsite``).

-  [ int06 ], [ int07 ],
   [ int08 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving the position operator :math:`\mu_i` for site number [ int05 ] in each coordinate direction :math:`\mu = x, y, z`. You can specify the :math:`x, y, z` direction components as 0 or greater and less than [ int01 ], [ int02 ], [ int03 ].

-  [ int09 ]

   **Type :** int-type (blank parameter not allowed)

   **Description :** An integer giving the orbital index in the unit cell corresponding to site number [ int05 ]. You can specify it as an integer from 0 to [ int04 ]-1.

Usage rules
^^^^^^^^^^^

The rules for using this file are as follows:

-  Headers cannot be omitted.

-  All the site indices must be specified.

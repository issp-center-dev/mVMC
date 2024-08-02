.. include:: ../../bib/ref.txt

.. _OutputFile:
             
Output files
============

The list of output files are shown in Table [Table:Output], where \*\*\*
and xxx are the header indicated by
``CParaFileHead``,``CDataFileHead`` in ``ModPara`` file, respectively.
yyy is a number given by
``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp``,
where both ``NDataIdxStart`` and ``NDataQtySmp`` are defined in
``ModPara`` file and zzz is a number given by ``NDataIdxStart`` in
``ModPara``.

+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| Name                                     | Details for corresponding files                                                                                        |
+==========================================+========================================================================================================================+
| \*\*\*\_opt.dat                          | All optimized parameters.                                                                                              |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| \*\*\*\_gutzwiller\_opt.dat              | Optimized gutzwiller factors.                                                                                          |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| \*\*\*\_jastrow\_opt.dat                 | Optimized jastrow factors.                                                                                             |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| \*\*\*\_doublonHolon2site\_opt.dat       | Optimized 2-site doublon-holon correlation factors.                                                                    |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| \*\*\*\_doublonHolon4site\_opt.dat       | Optimized 4-site doublon-holon correlation factors.                                                                    |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| \*\*\*\_generalRBM\_physlayer\_opt.dat   | Optimized variational parameters of the general RBM correlation factors only in a physical layer.                      |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| \*\*\*\_generalRBM\_hiddenlayer\_opt.dat | Optimized variational parameters of the general RBM correlation factors only in a hidden layer.                        |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
|  \*\*\*\_orbital\_opt.dat                | Optimized variational parameters of the general RBM correlation factors which connect physical and hidden layers.      |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| \*\*\*\_orbital\_opt.dat                 | Optimized pair orbital factors.                                                                                        |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| xxx\_out\_yyy.dat                        | Energy and deviation.                                                                                                  |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| xxx\_var\_yyy.dat                        | Progress information for optimizing variational parameters.                                                            |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| xxx\_CalcTimer.dat                       | Computation time for each processes.                                                                                   |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| xxx\_time\_zzz.dat                       | Progress information for Monte Carlo samplings.                                                                         |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| xxx\_cisajs\_yyy.dat                     | One body Green’s functions.                                                                                            |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| xxx\_cisajscktalt\_yyy.dat               | Correlation functions.                                                                                                 |
+------------------------------------------+------------------------------------------------------------------------------------------------------------------------+

Output file for variational parameters (\*\*\*\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The average and deviation values of variational parameters and the
energy optimized by the SR method are outputted in the following order:

.. math::

   \langle H \rangle, \langle H^2 \rangle, g_i, v_{ij},
   \alpha_{2nt}^{d(h)}, \alpha_{4nt}^{d(h)}, f_{ij} \nonumber.

The type of average values is a complex number, while that of the
deviation is a real number. Since the initial values of all variational
parameters are specified at the beginning of the calculation, the
calculation of physical quantities is done by using this file file.
Here, \*\*\* is the header indicated by ``CParaFileHead`` in ``ModPara``
file.

Output files for variational parameters at each steps (xxx\_var\_yyy.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The average and deviation values of variational parameters and the
energy optimized by the SR method are postscripted by the order same as
$$$\_opt.dat file (the deviation is always outputted as 0). Here, xxx is
the header indicated by ``CDataFileHead`` in ``ModPara`` file and yyy is
a number given by
``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp``,
where both ``NDataIdxStart`` and ``NDataQtySmp`` are defined in
``ModPara`` file.

Output file for gutzwiller factors (\*\*\*\_gutzwiller\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optimized Gutzwiller factors by SR method are outputted. The file
format is same as the ``InGutzwiller`` file defined in Sec.
:ref:`InputParam`.

Output file for jastrow factors (\*\*\*\_jastrow\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optimized Jastrow factors by SR method are outputted. The file
format is same as the ``InJastrow`` file defined in Sec.
:ref:`InputParam`.

Output file for doublonHolon 2-site factors (\*\*\*\_doublonHolon2site\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optimized 2-site doublon-holon correlation factors by SR method are
outputted. The file format is same as the ``InDH2`` file defined in Sec.
:ref:`InputParam`.

Output file for doublonHolon 4-site factors (\*\*\*\_doublonHolon4site\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optimized 4-site doublon-holon correlation factors by SR method are
outputted. The file format is same as the ``InDH4`` file defined in Sec.
:ref:`InputParam`.

Output files for general RBM factors(\*\*\*\_generalRBM\_physlayer\_opt.dat, \*\*\*\_generalRBM\_hiddenlayer\_opt.dat, \*\*\*\_generalRBM\_physhidden\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optimized general RBM correlation factors by SR method are outputted. 
The file format is same as the ``InGeneralRBM_PhysLayer``, ``InGeneralRBM_HiddenLayer`` and ``InGeneralRBM_PhysHidden`` files defined in Sec.
:ref:`InputParam`.

Output file for pair orbitals (\*\*\*\_orbital\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optimized pair orbitals by SR method are outputted. The file format
is same as the ``InOrbital`` file defined in Sec. :ref:`InputParam`.

xxx\_out\_yyy.dat
~~~~~~~~~~~~~~~~~

The calculation information at each bins are outputted in the order:

.. math::

   \langle H \rangle, \langle H^2 \rangle,
   \frac{\langle H^2 \rangle- \langle H \rangle^2 }{\langle H \rangle^2},
   \langle S^z \rangle, \langle (S^z)^2 \rangle
   \nonumber.

The type of :math:`\langle H \rangle` is a complex number, and that of
the others is a real number. Here, xxx is the header indicated by
``CDataFileHead`` in ``ModPara`` file and yyy is a number given by
``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp``,
where both ``NDataIdxStart`` and ``NDataQtySmp`` are defined in
``ModPara`` file. An example of outputted file is shown as follows.

::

    1.151983765704212992e+01  8.124622418360909482e-01  \
    1.619082955438887268e+02  2.019905203939084959e-01 
    1.288482613817423150e+01  5.006903733262847433e-01  
    1.972000325276957824e+02  1.824505193695792893e-01
    1.308897206011880421e+01  5.701244886956570168e-01  \
    2.072610167083121837e+02  2.029162857569105916e-01
    ...

xxx\_CalcTimer.dat 
~~~~~~~~~~~~~~~~~~~

After finishing calculation, the processing time is outputted in the
order of the name, the number assigned by the process and the seconds at
each processes. An example of outputted file is shown as follows.

::

    All                         [0]     15.90724
    Initialization              [1]      0.04357
      read options             [10]      0.00012
      ReadDefFile              [11]      0.00082
      SetMemory                [12]      0.00002
      InitParameter            [13]      0.03026
    VMCParaOpt                  [2]     15.86367
      VMCMakeSample             [3]     12.85650
    ...

xxx\_time\_zzz.dat 
~~~~~~~~~~~~~~~~~~~

The calculation information at each bins are outputted in the order of
the sampling number, the acceptance ratio for hopping and exchange term
(acc\_hopp, acc\_ex), trial numbers to update for hopping and exchange
term (n\_hopp, n\_ex) and the time stamp. Here, xxx is the header
indicated by ``CDataFileHead`` in ``ModPara`` file and zzz is a number
given by ``NDataIdxStart`` in ``ModPara``. An example of outputted file
is shown as follows.

::

    00000  acc_hop acc_ex  n_hop    n_ex     : Mon Jul 25 14:03:29 2016
    00001  0.59688 0.00000 320      0        : Mon Jul 25 14:03:30 2016
    00002  0.47727 0.00000 176      0        : Mon Jul 25 14:03:30 2016
    00003  0.50000 0.00000 176      0        : Mon Jul 25 14:03:30 2016
    00004  0.49432 0.00000 176      0        : Mon Jul 25 14:03:30 2016
    00005  0.57386 0.00000 176      0        : Mon Jul 25 14:03:30 2016
    00006  0.55114 0.00000 176      0        : Mon Jul 25 14:03:30 2016    
    …

xxx\_cisajs\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~

This file is the outputted files for one-body Green’s function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`. The target
components are set in the input file with the keyword "OneBodyG". Here,
xxx is the header indicated by ``CDataFileHead`` in ``ModPara`` file and
yyy is a number given by
``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp``,
where both ``NDataIdxStart`` and ``NDataQtySmp`` are defined in
``ModPara`` file.

An example of the file format is as follows.

::

        0    0    0    0 0.4452776740 0.0000000000
        0    1    0    1 0.4452776740 0.0000000000
        1    0    1    0 0.5000000000 0.0000000000
        1    1    1    1 0.5000000000 0.0000000000
        2    0    2    0 0.4452776740 0.0000000000
        2    1    2    1 0.4452776740 0.0000000000
        3    0    3    0 0.5000000000 0.0000000000
        3    1    3    1 0.5000000000 0.0000000000
        ...

File format
^^^^^^^^^^^

-  [int01]  [int02]  [int03]  [int04]  [double01]  [double02]

Parameters
^^^^^^^^^^

-  [int01], [int03]

   **Type :** Int

   **Description :** The integer of the site number.
   [int01] and [int03] show the
   :math:`i` and :math:`j` site numbers, respectively.

-  [int02], [int04]

   **Type :** Int

   | **Description :** The integer of the spin index:
   | 0: Up-spin
   | 1: Down-spin.
   | [int02] and [int04] show
     :math:`\sigma_1` and :math:`\sigma_2`, respectively.

-  [double01], [double02]

   **Type :** Double

   | **Description :** The value of
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`.
   | [double01] and [double02]
     show the real and imaginary part of
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`,
     respectively.

xxx\_cisajscktalt\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file is the outputted files for the two-body Green’s function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`.
The target components are set in the input file with the keyword
"TwoBodyG". Here, xxx is the header indicated by ``CDataFileHead`` in
``ModPara`` file and yyy is a number given by
``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp``,
where both ``NDataIdxStart`` and ``NDataQtySmp`` are defined in
``ModPara`` file. An example of the file format is as follows.

::

        0    0    0    0    0    0    0    0 0.4452776740 0.0000000000
        0    0    0    0    0    1    0    1 0.1843355815 0.0000000000
        0    0    0    0    1    0    1    0 0.1812412105 0.0000000000
        0    0    0    0    1    1    1    1 0.2640364635 0.0000000000
        0    0    0    0    2    0    2    0 0.0279690007 0.0000000000
        0    0    0    0    2    1    2    1 0.2009271524 0.0000000000
        0    0    0    0    3    0    3    0 0.2512810778 0.0000000000
        0    0    0    0    3    1    3    1 0.1939965962 0.0000000000
        ...

File format
^^^^^^^^^^^

-  [int01]  [int02]  [int03]  [int04]  [int05]  [int06]  [int07]  [int08]  [double01]  [double02].

Parameters
^^^^^^^^^^

-  [int01],
   [int03],[int05],
   [int07]

   **Type :** Int

   **Description :** The integer of the site number.
   [int01], [int03],
   [int05], and [int07] show the
   :math:`i`, :math:`j`, :math:`k`, and :math:`l` site numbers,
   respectively.

-  [int02],
   [int04],[int06],
   [int08]

   **Type :** Int

   **Description :** The integer of the spin index:

   0: Up-spin

   1: Down-spin.

   [int02], [int04],
   [int06], and [int08] show
   :math:`\sigma_1`, :math:`\sigma_2`, :math:`\sigma_3`, and
   :math:`\sigma_4`, respectively.

-  [double01], [double02]

   **Type :** Double

   **Description :** The value of
   :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`.
   [double01] and [double02]
   show the real and imaginary part of
   :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`,
   respectively.

xxx\_ls\_out\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~~

This file is the outputted files for :math:`\langle H \rangle`,
:math:`(\langle H^2\rangle - \langle H \rangle^2)/\langle H \rangle^2`, and the optimized parameter :math:`\alpha`
obtained by Power Lanczos method. This file is outputted when
``NVMCCalMode`` = 1, ``NLanczosmode`` = 1 or 2 are set in ``ModPara``
file. Here, xxx is the header indicated by ``CDataFileHead`` in
``ModPara`` file and yyy is a number given by
``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp``,
where both ``NDataIdxStart`` and ``NDataQtySmp`` are defined in
``ModPara`` file.

xxx\_ls\_cisajs\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~~~~~

This file is the outputted files for one-body Green’s function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle` obtained by
Power Lanczos method. The file format is same as the
xxx\_cisajs\_yyy.dat file. This file is outputted when ``NVMCCalMode`` =
1, ``NLanczosmode`` = 2 are set in ``ModPara`` file.

xxx\_ls\_cisajscktalt\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file is the outputted files for the two-body Green’s function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`
obtained by Power Lanczos method. The file format is same as the
xxx\_cisajscktalt\_yyy.dat file. This file is outputted when
``NVMCCalMode`` = 1, ``NLanczosmode`` = 2 are set in ``ModPara`` file.


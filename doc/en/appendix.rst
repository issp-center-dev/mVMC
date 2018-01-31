Program for the unrestricted Hartree-Fock approximation
=======================================================

In mVMC package, there is a program to calculate the initial values of
the pair orbital parameters :math:`f_{ij}` by using the unrestricted
Hartree-Fock (UHF) approximation (relation between Pfaffian Slater
determinant and single Slater determinant is explained in sec.
:ref:`PuffAndSlater`). It is noted that the target system of this program
is the itinerant electron system.

Overview
--------

In UHF approximation, two-body interaction terms are approximated as
one-body interaction terms by taking into account of the fluctuation,
:math:`\delta A \equiv A-\langle A \rangle`, up to the first order. As
an example, we consider the inter-site coulomb interactions

.. math::

   {\cal H}_V = \sum_{i,j}V_{ij} n_ {i}n_{j},

where we define :math:`i\equiv (i, \sigma)`,
:math:`j\equiv (j, \sigma')` for simplicity. Then, the interaction terms
can be approximated as

.. math::

   \begin{aligned}
   n_ {i}n_{j} &= (\langle n_{i} \rangle +\delta n_i) (\langle n_{j} \rangle +\delta n_j)
   - \left[ \langle c_{i}^{\dagger}c_j \rangle +\delta (c_{i}^{\dagger}c_j ) \right]
     \left[ \langle c_{j}^{\dagger}c_i \rangle +\delta (c_{j}^{\dagger}c_i )\right]
   \nonumber\\
   &\sim
   \langle n_{i} \rangle n_j+\langle n_{j} \rangle  n_i - \langle c_{i}^{\dagger}c_j \rangle  c_{j}^{\dagger}c_i
   - \langle c_{j}^{\dagger}c_i \rangle c_{i}^{\dagger}c_j 
   - \langle n_{i} \rangle \langle n_j \rangle +  \langle c_{j}^{\dagger}c_i \rangle \langle c_{i}^{\dagger}c_j \rangle.
   \end{aligned}

Also for other types of interaction, the problem can be attributed to a
one-body problem by using a similar approximation. Actual calculation is
performed iteratively until that self-consistent solution for the mean
values of the above observables are obtained.

Source code
~~~~~~~~~~~

A set of source codes are included in the directory
``src/ComplexUHF/src``.

How to compile
~~~~~~~~~~~~~~

To compile source codes, move to the directory just below the main
directory of mVMC, and execute

.. code-block:: bash

    $ make mvmc

in a similar way as the compile of mVMC. After compiling, an executable
file ``UHF`` is generated in ``src/ComplexUHF/src``.

Input files
~~~~~~~~~~~

A file for assigning input files (namelsit.def)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following files are needed to use the program of UHF. The format
of ``namelist.def`` is the same as defined in :ref:`InputFileList`.

-  ``ModPara``

-  ``LocSpin``

-  ``Trans``

-  ``CoulombIntra``

-  ``CoulombInter``

-  ``Hund``

-  ``PairHop``

-  ``Exchange``

-  ``Orbital`` / ``OrbitalAntiParallel``

-  ``OrbitalParallel``

-  ``OrbitalGeneral``

-  ``Initial``

Although the format of these files are the same as those for mVMC
basically, the following items are different:

-  Parameters assigned in ``ModPara`` file.

-  Addition of ``Initial`` file.

We explain details of the format of these files as follows.

Parameters assigned in ModPara file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters needed in the program of UHF are as follows:

-  ``Nsite``

-  ``Ne``

-  ``Mix``

-  ``EPS``

-  ``IterationMax``

The parameters, ``Nsite`` and ``Ne``, are common as mVMC. The other
three parameters are specific to UHF:

-  ``Mix``
   
   Linear mixing is assigned by double-type. When mix=1, a new Green’s
   function is fully updated without using a old one.

-  ``EPS``
      
   A condition for convergence is assigned by int-type. When a
   residual error between a new Green’s function and a previous one is
   less than :math:`10^{-\verb|eps|}`, the iteration of calculation is
   stopped.

-  ``IterationMax``
      
   A maximum number of the loop is assigned by int-type.

If there are the other parameters for mVMC in this file , warning is
output to the standard output (the calculation is not stopped).

Initial file
^^^^^^^^^^^^

Initial values of Green’s function
:math:`G_{ij\sigma_1\sigma_2}\equiv \langle c_{i\sigma_1}^\dagger c_{j\sigma_2}\rangle`
are given. The format is the same as ``Trans`` file, and instead of
:math:`t_{ij\sigma_1\sigma_2}`, values of :math:`G_{ij\sigma_1\sigma_2}`
are described. Green’s function is set as zero if values are not given.

Usage
-----

Calculation of UHF is performed by the same way as mVMC, i.e., by
executing the command

::

    $ UHF namelist.def

The routine of the calculation is as follows.

#. Reading files

#. Construction of a Hamiltonian

#. Self-consistent calculation of Green’s function

#. Output of :math:`f_{ij}` and other files

Examples of output after calculation are as follows.

-  zvo\_result.dat:

   The energy and the particle number are output.

   ::

        energy -15.2265348135
        num    36.0000000000

-  zvo\_check.dat:

   The step number of the iteration, the mean of the
   absolute value of the residual error in Green’s function, the energy
   in convergence process, and the particle number are output in order.

   ::

        0  0.004925645652 -544.963484605164 36.000000
        1  0.002481594941 -278.304285708488 36.000000
        2  0.001274395448 -147.247026925130 36.000000
        3  0.000681060599 -82.973664527606 36.000000
       ...

-  zvo\_UHF\_cisajs.dat:

   Convergent one-body Green’s function
   :math:`G_{ij\sigma_1\sigma_2}\equiv\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`
   is output.
   For all the components,
   :math:`i, \sigma_1, j, \sigma_2, {\rm Re}\left[G_{ij\sigma_1\sigma_2}\right], {\rm Im}\left[G_{ij\sigma_1\sigma_2}\right]`
   are output in order.

   ::

           0    0    0    0 0.5037555283 0.0000000000
           0    0    0    1 0.4610257618 0.0003115503
           0    1    0    0 0.4610257618 -0.0003115503
           0    1    0    1 0.4962444717 0.0000000000
        ...

-  zvo\_eigen.dat:
   
   Convergent eigenvalues of the Hamiltonian are
   output in ascending order.

   ::

        1  -2.9425069199
        2  -2.9425069198
        3  -1.5005359205 
        ...

-  zvo\_gap.dat:
   
   For the total electron number :math:`N_{\rm tot}`, the
   energy difference :math:`\Delta E= E(N_{\rm tot}+1)-E(N_{\rm tot})`
   is output.

   ::

         5.2208232631

-  zvo\_orbital\_opt.dat:

   :math:`f_{ij}` generated from the Slater
   determinant. The file with the same format as
   ``InOrbital``,\ ``InOrbitalAntiParallel``,
   ``InOrbitalParallel``,\ ``InOrbitalAntiGeneral`` file is output. By
   referring ``Orbital``, ``OrbitalAntiParallel``,
   ``OrbitalParallel``,\ ``OrbitalAntiGeneral`` file, :math:`f_{ij}` is
   calculated (for the same type of parameters, the averaged value is
   calculated).

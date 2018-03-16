Overview
========

In this document, we introduce how we compute downfolded models
with mVMC or :math:`{\mathcal H}\Phi` in conjunction to
`RESPACK <https://sites.google.com/view/kazuma7k6r>`_.

Prerequisite
------------

We compute the Kohn-Sham orbitals with
`QuantumESPRESSO <http://www.quantum-espresso.org/>`_
or
`xTAPP <http://xtapp.cp.is.s.u-tokyo.ac.jp/>`_,
and obtain the Wannier function, the dielectric function,
the effective interaction with RESPACK,
and simulate quantum lattice models with
mVMC or :math:`{\mathcal H}\Phi`.
Therefore, these programs must be available in our machine.

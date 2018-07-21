.. include:: ../bib/ref.txt

Algorithm
=========

Variational Monte Calro Method
------------------------------

The variational Monte Carlo (VMC) method is a method for calculating
approximate wave functions of a ground state and low-lying excited
states by optimizing variational parameters included in a trial wave
function. In calculating expectation values of physical quantities for
the trial wave functions, the Markov chain Monte Carlo method is applied
for efficient important sampling.

In the mVMC package, we choose a spatial configuration for electrons as
a complete set of bases in sampling:

.. math::

   | x\rangle =  \prod_{n=1}^{N/2} c_{r_{n\uparrow}}^{\dagger} \prod_{n=1}^{N/2} c_{r_{n\downarrow}}^{\dagger} |0 \rangle,

where :math:`r_{n\sigma}` is a position of :math:`n`\ -th electron with
:math:`\sigma (=\uparrow \rm{or} \downarrow)` spin, and
:math:`c_{r_{n\sigma}}^{\dagger}` is a creation operator of electrons. By
using this basis set, the expectation value of an operator :math:`A` is
expressed as

.. math::

   \langle A \rangle =\frac{\langle \psi| A| \psi \rangle}{\langle \psi | \psi \rangle} 
   =\sum_x \frac{\langle \psi| A | x\rangle \langle x| \psi \rangle}{\langle \psi |\psi \rangle}.

If we define a weight of the Markov chain Monte Carlo method as

.. math::

   \rho(x)=\frac{|\langle x| \psi \rangle|^2}{\langle \psi | \psi \rangle} \ge 0, \quad \sum_{x} \rho(x)=1,

we can rewrite :math:`\langle A \rangle` in the following form:

.. math::

   \langle A \rangle =\sum_x \rho(x) \frac{\langle \psi| A | x\rangle }{\langle \psi |x \rangle}.

By using this form, the Markov chain Monte Carlo method is performed
for sampling with respect to :math:`x`. The local Green’s function
:math:`G_{ij\sigma\sigma'}(x)`, which is defined as

.. math::

   G_{ij\sigma\sigma'}(x)=\frac{\langle \psi | c_{i\sigma}^{\dagger} c_{j\sigma'} | \psi \rangle}{\langle \psi | x \rangle},

is also evaluated by the same sampling method by taking
:math:`A = c_{i\sigma}^{\dagger} c_{j\sigma'}`. We adopt the Mersenne
twister method as a random number generator for
sampling [Mutsuo2008_ ].

Bogoliubov representation
-------------------------

In the VMC calculation for spin systems, we use the Bogoliubov
representation. In the input files defining the one-body term
(``transfer``) and the two-body term (``InterAll``), and the output
files for correlation functions, the indices must be assigned by the
Bogoliubov representation, in which the spin operators are generally
expressed by creation/annihilation operators of fermions as

.. math::

   \begin{aligned}
   S_{i z} &= \sum_{\sigma = -S}^{S} \sigma c_{i \sigma}^\dagger c_{i \sigma},
   \\
   S_{i}^+ &= \sum_{\sigma = -S}^{S-1} 
   \sqrt{S(S+1) - \sigma(\sigma+1)} 
   c_{i \sigma+1}^\dagger c_{i \sigma},
   \\
   S_{i}^- &= \sum_{\sigma = -S}^{S-1} 
   \sqrt{S(S+1) - \sigma(\sigma+1)} 
   c_{i \sigma}^\dagger c_{i \sigma+1}.
   \end{aligned}

Since the present package support only :math:`S=1/2` spin systems, the
Bogoliubov representation obtained by substituting :math:`S=1/2` into
the above equations is used.

.. _PuffAndSlater:

Properties of the Pfaffian-Slater determinant
---------------------------------------------

In this section, we explain some properties of the Pfaffian-Slater
determinant. We derive the general relation between a Pfaffian-Slater
determinant and a single Slater determinant in :ref:`Antiparallel Pfaffian <PfaffianAP>`
and :ref:`General Pfaffian <PfaffianP>` . We also discuss meaning of the singular value
decomposition of coefficients :math:`f_{ij}` in
:ref:`SVD <PfaffianSingular>`.

.. _PfaffianAP:

Relation between :math:`f_{ij}` and :math:`\Phi_{in\sigma}` (the case of the anti-parallel pairing)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the many-variable variational Monte Carlo (mVMC) method, the one-body
part of the trial wave function is expressed by the Pfaffian Slater
determinant defined as

.. math::

   |\phi_{\rm Pf}\rangle=\Big(\sum_{i,j=1}^{N_{s}}f_{ij}
   c_{i\uparrow}^{\dagger}c_{j\downarrow}^{\dagger}\Big)^{N_{\rm e}/2}|0\rangle,

where :math:`N_{s}` is number of sites, :math:`N_{e}` is number of
total particles, and :math:`f_{ij}` are variational parameters. For
simplicity, we assume that :math:`f_{ij}` are a real number. The single
Slater determinant is defined as

.. math::

   \begin{aligned}
   |\phi_{\rm SL}\rangle&=\Big(\prod_{n=1}^{N_{e}/2}\psi_{n\uparrow}^{\dagger}\Big)
   \Big(\prod_{m=1}^{N_{e}/2}\psi_{m\downarrow}^{\dagger}\Big)|0\rangle, \\
   \psi_{n\sigma}^{\dagger}&=\sum_{i=1}^{N_{s}}\Phi_{in\sigma}c^{\dagger}_{i\sigma},
   \end{aligned}

Here, :math:`\Phi_{in\sigma}` is an orthonormal basis, i.e., satisfies

.. math:: \sum_{i=1}^{N_{s}}\Phi_{in\sigma}\Phi_{im\sigma}=\delta_{nm},

where :math:`\delta_{nm}` is the Kronecker’s delta. From this
orthogonality, we can prove the relation

.. math::

   \begin{aligned}
   [\psi^{\dagger}_{n\sigma},\psi_{m\sigma}]_{+}&=\delta_{nm},\\
   G_{ij\sigma}=\langle c_{i\sigma}^{\dagger}c_{j\sigma}\rangle 
   &=\frac{\langle \phi_{\rm SL}| c_{i\sigma}^{\dagger}c_{j\sigma} | \phi_{\rm SL}\rangle}{\langle \phi_{\rm SL}|\phi_{\rm SL}\rangle } \\
   &=\sum_{n} \Phi_{in\sigma} \Phi_{jn\sigma}.
   \end{aligned}

Next, let us prove the relation between :math:`f_{ij}` and
:math:`\Phi_{in\sigma}` by modifying :math:`|\phi_{\rm SL}\rangle`. By
the commutation relation for :math:`\psi^{\dagger}_{n\sigma}`,
:math:`|\phi_{\rm SL}\rangle` is rewritten as

.. math::

   \begin{aligned}
   |\phi_{\rm SL}\rangle \propto \prod_{n=1}^{N_{e}/2}\Big(\psi_{n\uparrow}^{\dagger}\psi_{\mu(n)\downarrow}^{\dagger}\Big)|0\rangle,
   \end{aligned}

where :math:`\mu(n)` represents permutation of a sequence of natural
numbers, :math:`n= 1, 2, \cdots, N_{e}/2`. For simplicity, let us take
identity permutation (:math:`\mu(n) = n`). By defining
:math:`K_{n}^{\dagger}=\psi_{n\uparrow}^{\dagger}\psi_{n\downarrow}^{\dagger}`,
and by using the relation
:math:`K_{n}^{\dagger}K_{m}^{\dagger}=K_{m}^{\dagger}K_{n}^{\dagger}`,
we can derive the relation

.. math::

   \begin{aligned}
   |\phi_{\rm SL}\rangle &\propto \prod_{n=1}^{N_{e}/2}\Big(\psi_{n\uparrow}^{\dagger}\psi_{n\downarrow}^{\dagger}\Big)|0\rangle
   =\prod_{n=1}^{N_{e}/2} K_{n}^{\dagger}|0\rangle \\
   &\propto\Big(\sum_{n=1}^{\frac{N_{e}}{2}}K_{n}^{\dagger}\Big)^{\frac{N_{e}}{2}} |0\rangle
   =\Big(\sum_{i,j=1}^{N_{s}}\Big[\sum_{n=1}^{\frac{N_{e}}{2}}\Phi_{in\uparrow}\Phi_{jn\downarrow}\Big]
   c_{i\uparrow}^{\dagger}c_{j\downarrow}^{\dagger}\Big)|0\rangle.
   \end{aligned}

This result indicates that :math:`f_{ij}` is expressed by the
coefficients of the single Slater determinant as

.. math::

   \begin{aligned}
   f_{ij}=\sum_{n=1}^{\frac{N_{e}}{2}}\Phi_{in\uparrow}\Phi_{jn\downarrow}.
   \end{aligned}

We note that this is one of a number of possible expressions of
:math:`f_{ij}` derived from one single Slater determinant. Since
:math:`f_{ij}` depends not only on the choice of the pairing degrees of
freedom (i.e., the choice of :math:`\mu(n)`) but also on the choice of
the gauge degrees of freedom (i.e., the sign of
:math:`\Phi_{in\sigma}`), the parameter :math:`f_{ij}` has huge
redundancy.

.. _PfaffianP:

Relation between :math:`F_{IJ}` and :math:`\Phi_{In\sigma}` (the case of the general pairing)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We extend the relation between the Pfaffian-Slater wave function and the
single Slater wave function into the general pairing case including the
spin-parallel pairing. We define the Pfaffian-Slater wave function and
the single Slater wave function as

.. math::

   \begin{aligned}
   |\phi_{\rm Pf}\rangle&=\Big(\sum_{I,J=1}^{2N_{s}}F_{IJ}c_{I}^{\dagger}c_{J}^{\dagger}\Big)^{N_{\rm e}/2}|0\rangle, \\
   |\phi_{\rm SL}\rangle&=\Big(\prod_{n=1}^{N_{e}}\psi_{n}^{\dagger}\Big)|0\rangle,~~\psi_{n}^{\dagger}=\sum_{I=1}^{2N_{s}}\Phi_{In}c^{\dagger}_{I},
   \end{aligned}

respectively, where :math:`I`, :math:`J` denote the site index
including the spin degrees of freedom. By the similar argument as the
anti-parallel pairing case, we can derive the following relation:

.. math::

   \begin{aligned}
   F_{IJ}=\sum_{n=1}^{\frac{N_{e}}{2}}\Big(\Phi_{In}\Phi_{Jn+1}-\Phi_{Jn}\Phi_{In+1}\Big).
   \end{aligned}

Because this relation hold for the case of anti-parallel pairing, we
employ this relation in mVMC ver 1.0 and later.

.. _PfaffianSingular:

Singular value decomposition of :math:`f_{ij}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We define matrices :math:`F`, :math:`\Phi_{\uparrow}`,
:math:`\Phi_{\downarrow}`, and :math:`\Sigma` as

.. math::

   \begin{aligned}
   &(F)_{ij}=f_{ij},~~~ 
   (\Phi_{\uparrow})_{in}=\Phi_{in\uparrow},~~~ 
   (\Phi_{\downarrow})_{in}=\Phi_{in\downarrow}, \\
   &\Sigma={\rm diag}[\underbrace{1,\cdots,1}_{N_e/2},0,0,0].
   \end{aligned}

When :math:`f_{ij}` (i.e., the matrix :math:`F`) is related with a
single Slater determinant of the wave function, we can show that the
singular value decomposition of :math:`F` becomes

.. math::

   \begin{aligned}
   F=\Phi_{\uparrow}\Sigma\Phi_{\downarrow}^{t}.
   \end{aligned}

This result indicates that when the number of nonzero singular values
is :math:`N_{e}/2`, and when all the nonzero singular values of
:math:`F` are one in the singular value decomposition of :math:`F`, the
Pfaffian-Slater wave function parametrized by :math:`f_{ij}` coincides
with a single Slater determinant (i.e. a solution of the mean-field
approximation). In other words, the numbers of the nonzero singular
values and their difference from one offer a quantitative criterion how
the Pfaffian-Slater determinant deviates from the single Slate
determinant.

Power Lanczos method
--------------------

In this section, we show how to determine :math:`\alpha` in the
power-Lanczos method. We also explain the calculation of physical
quantities after the single-step Lanczos method.

Determination of :math:`\alpha`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we briefly explain the sampling procedure of the variational
Monte Carlo (VMC) method. Physical properties :math:`\hat{A}` are
calculated as follows:

.. math::

   \begin{aligned}
   &\langle \hat{A}\rangle = \frac{\langle \phi| \hat{A}|\phi \rangle}{\langle \phi| \phi \rangle} = \sum_{x} \rho(x) F(x, {\hat{A}}),\\
   & \rho(x)=\frac{|\langle \phi|x\rangle|^2}{\langle \phi | \phi \rangle}, ~~~~F(x,  {\hat{A}}) =  \frac{\langle x| \hat{A}|\phi \rangle}{\langle x| \phi \rangle}.
   \end{aligned}

There are two ways to calculate the product of the operators
:math:`\hat{A}\hat{B}`.

.. math::

   \begin{aligned}
   &\langle \hat{A} \hat{B}\rangle = \sum_{x} \rho(x) F(x, {\hat{A}\hat{B}}),\\
   &\langle \hat{A} \hat{B}\rangle = \sum_{x} \rho(x) F^{\dagger}(x, {\hat{A})F(x, \hat{B}}).
   \end{aligned}

As we explain later, in general, the latter way is numerical stable
one. For example, we consider the expectation value of the variance,
which is defined as
:math:`\sigma^2=\langle (\hat{H}-\langle \hat{H}\rangle)^2\rangle`.
There are two ways to calculate the variance.

.. math::

   \begin{aligned}
   \sigma^2&=\sum_{x} \rho(x) F(x,  (\hat{H}-\langle \hat{H}\rangle)^2) = \sum_{x} \rho(x) F(x,  \hat{H}^2) - \left[ \sum_{x} \rho(x) F(x,  \hat{H})\right]^2 ,\\
   \sigma^2&=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}-\langle \hat{H}\rangle)F(x,  \hat{H}-\langle \hat{H}\rangle) \nonumber \\
   &= \sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{H})- \left[ \sum_{x} \rho(x) F(x,  \hat{H})\right]^2
   \end{aligned}

From its definition, the latter way gives the positive definitive
variance even for the finite sampling while the former way does not
guarantee the positive definitiveness of the variance. Here, we consider
the expectation values of energy and variance for the (single-step)
power Lanczos wave function
:math:`|\phi\rangle =(1+\alpha \hat{H}) |\psi \rangle`. The energy is
calculated as

.. math::

   \begin{aligned}
   E_{LS}(\alpha) =\frac{\langle \phi| \hat{H} |\phi\rangle}{\langle \phi|\phi\rangle}=\frac{h_1 + \alpha(h_{2(20)} + h_{2(11)}) + \alpha^2 h_{3(12)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}},
   \end{aligned}

where we define :math:`h_1`, :math:`h_{2(11)},~h_{2(20)},` and
:math:`h_{3(12)}` as

.. math::

   \begin{aligned}
   &h_1 =\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}),\\
   &h_{2(11)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{H}),\\
   &h_{2(20)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}^2),\\
   &h_{3(12)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H})F(x,  \hat{H}^2).
   \end{aligned}

From the condition
:math:`\frac{\partial E_{LS}(\alpha)}{\partial \alpha}=0`, i.e., by
solving the quadratic equations, we can determine the :math:`\alpha`.
The variance is calculate in the similar way.

Calculation of physical quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By using the optimized parameter :math:`\alpha`, we can calculate the
expected value of the operator :math:`\hat{A}` as

.. math::

   \begin{aligned}
   A_{LS}(\alpha) =\frac{\langle \phi| \hat{A} |\phi\rangle}{\langle \phi|\phi\rangle}=\frac{A_0 + \alpha(A_{1(10)} + A_{1(01)}) + \alpha^2 A_{2(11)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}},
   \end{aligned}

where we define :math:`A_0`, :math:`A_{1(10)},~A_{1(01)},` and
:math:`A_{2(11)}` as

.. math::

   \begin{aligned}
   &A_0 =\sum_{x} \rho(x) F(x,  \hat{A}),\\
   &A_{1(10)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{A}),\\
   &A_{1(01)}=\sum_{x} \rho(x) F(x, \hat{A}\hat{H}),\\
   &A_{2(11)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H})F(x,  \hat{A}\hat{H}).
   \end{aligned}



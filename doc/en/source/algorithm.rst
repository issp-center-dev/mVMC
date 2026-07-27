.. include:: ../../bib/ref.txt

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

   | x\rangle =  \prod_{n=1}^{N_e/2} c_{r_{n\uparrow}}^{\dagger} \prod_{n=1}^{N_e/2} c_{r_{n\downarrow}}^{\dagger} |0 \rangle,

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

Relation between :math:`F_{IJ}` and :math:`\Phi_{In}` (the case of the general pairing)
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
   F_{IJ}=\sum_{n=1}^{\frac{N_{e}}{2}}\Big(\Phi_{I,2n-1}\Phi_{J,2n}-\Phi_{J,2n-1}\Phi_{I,2n}\Big).
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

The power-Lanczos method improves an optimized VMC wave function
:math:`|\psi\rangle` by using the trial wave function

.. math::

   |\phi_p\rangle
   =P_p(\hat H)|\psi\rangle
   =\sum_{k=0}^{p}c_k\hat H^k|\psi\rangle.

The coefficients :math:`c_k` minimize the energy in the Krylov subspace
spanned by
:math:`\{|\psi\rangle,\hat H|\psi\rangle,\ldots,\hat H^p|\psi\rangle\}`.
The first- and second-step trial wave functions used by mVMC are

.. math::

   \begin{aligned}
   |\phi_1\rangle &= (1+\alpha\hat H)|\psi\rangle,\\
   |\phi_2\rangle &= (1+\alpha\hat H+\beta\hat H^2)|\psi\rangle.
   \end{aligned}

The first step requires ``NVMCCalMode=1``, ``NLanczosMode>=1``, and
``NLanczosStep=1``. The second step requires ``NVMCCalMode=1``,
``NLanczosMode=1``, and ``NLanczosStep=2``.

The following subsections first explain the first-step optimization of
:math:`\alpha`. They then separate the variational principle that
determines :math:`\alpha,\beta` in the second step from the numerical
stabilization used by the implementation. The final subsection describes
the calculation of observables after the first power-Lanczos step.

First power-Lanczos step and determination of :math:`\alpha`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we briefly explain the sampling procedure of the variational
Monte Carlo (VMC) method. Here, :math:`x` denotes a configuration sampled
from the base VMC wave function :math:`|\psi\rangle`. The expectation
value of an observable :math:`\hat A` is

.. math::

   \begin{aligned}
   &\langle \hat{A}\rangle_\psi
   = \frac{\langle \psi|\hat A|\psi\rangle}{\langle \psi|\psi\rangle}
   = \sum_x \rho_\psi(x)F(x,\hat A),\\
   &\rho_\psi(x)=\frac{|\langle x|\psi\rangle|^2}{\langle\psi|\psi\rangle},
   \qquad
   F(x,\hat A)=\frac{\langle x|\hat A|\psi\rangle}
                     {\langle x|\psi\rangle}.
   \end{aligned}

The quantity :math:`F(x,\hat A)` is the local estimator of :math:`\hat A`
at configuration :math:`x`. The expectation value of a product
:math:`\hat A\hat B` has the following two representations:

.. math::

   \begin{aligned}
   &\langle \hat A\hat B\rangle_\psi
   = \sum_x\rho_\psi(x)F(x,\hat A\hat B),\\
   &\langle \hat A\hat B\rangle_\psi
   = \sum_x\rho_\psi(x)F^\dagger(x,\hat A)F(x,\hat B).
   \end{aligned}

In the former representation, the product :math:`\hat A\hat B` is
evaluated as a single local estimator. In the latter, the complete set
:math:`\sum_{x} |x\rangle \langle x|` is inserted between
:math:`\hat A` and :math:`\hat B`, i.e., the product is split into the
bra part :math:`\hat A` and the ket part :math:`\hat B`.
The latter expression assumes that :math:`\hat A` is Hermitian. For a
general :math:`\hat A`, :math:`F^\dagger(x,\hat A)` must be replaced by
:math:`F^\dagger(x,\hat A^\dagger)`.
Both representations then give the same expectation value in the limit of an
infinite number of samples. They differ for a finite number of samples,
and the latter is numerically more stable in general. For compact
notation, :math:`\rho_\psi(x)` is abbreviated as :math:`\rho(x)` below.

For example, the variance of the energy,
:math:`\sigma^2=\langle (\hat{H}-\langle \hat{H}\rangle)^2\rangle`, can
be calculated in the two ways:

.. math::

   \begin{aligned}
   \sigma^2&=\sum_{x} \rho(x) F(x,  (\hat{H}-\langle \hat{H}\rangle)^2) = \sum_{x} \rho(x) F(x,  \hat{H}^2) - \left[ \sum_{x} \rho(x) F(x,  \hat{H})\right]^2 ,\\
   \sigma^2&=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}-\langle \hat{H}\rangle)F(x,  \hat{H}-\langle \hat{H}\rangle) \nonumber \\
   &= \sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{H})- \left[ \sum_{x} \rho(x) F(x,  \hat{H})\right]^2
   \end{aligned}

The latter way guarantees that the variance is non-negative even for a
finite number of samples, while the former way does not.

Next, we consider the energy expectation value of the first-step wave
function :math:`|\phi_1\rangle =(1+\alpha \hat{H}) |\psi \rangle`.
Expanding it with respect to :math:`\alpha` (the numerator and the
denominator are divided by :math:`\langle \psi | \psi \rangle`) gives

.. math::

   \begin{aligned}
   E_{LS}(\alpha) =\frac{\langle \phi_1| \hat{H} |\phi_1\rangle}{\langle \phi_1|\phi_1\rangle}
   =\frac{\langle \hat{H} \rangle_\psi + 2\alpha \langle \hat{H}^2 \rangle_\psi + \alpha^2 \langle \hat{H}^3 \rangle_\psi}
   {1 + 2\alpha \langle \hat{H} \rangle_\psi + \alpha^2 \langle \hat{H}^2 \rangle_\psi},
   \end{aligned}

which requires the moments of :math:`\hat{H}` up to the third order. To
estimate these moments with the splitting described above, we introduce
the following notation:

.. math::

   \begin{aligned}
   h_{n(ij)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}^i) F(x, \hat{H}^j), ~~~~ i+j=n.
   \end{aligned}

Namely, :math:`h_{n(ij)}` is an estimator of the :math:`n`-th moment
:math:`\langle \hat{H}^n \rangle_\psi`, in which :math:`\hat{H}^i` is
assigned to the bra side and :math:`\hat{H}^j` to the ket side (we set
:math:`F(x, \hat{H}^0)=1` and abbreviate :math:`h_{1(10)}` as
:math:`h_1`). The quantities used in the calculation are

.. math::

   \begin{aligned}
   &h_1 =\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}),\\
   &h_{2(11)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{H}),\\
   &h_{2(20)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}^2),\\
   &h_{3(12)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H})F(x,  \hat{H}^2),\\
   &h_{4(22)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}^2)F(x,  \hat{H}^2).
   \end{aligned}

Note that both :math:`h_{2(11)}` and :math:`h_{2(20)}` are estimators
of :math:`\langle \hat{H}^2 \rangle_\psi`, but they differ for a finite
number of samples because of the different splittings. By assigning the
powers of :math:`\hat{H}` to the bra and ket sides as evenly as
possible, only :math:`F(x, \hat{H}^2)` appears as the local estimator,
i.e., it is not necessary to evaluate the local estimators of
:math:`\hat{H}^3` and :math:`\hat{H}^4` directly.

With these estimators, the energy expectation value is evaluated as

.. math::

   \begin{aligned}
   E_{LS}(\alpha) =\frac{h_1 + \alpha(h_{2(20)} + h_{2(11)}) + \alpha^2 h_{3(12)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}}.
   \end{aligned}

The two :math:`\langle \hat{H}^2 \rangle_\psi` terms in the first order
of :math:`\alpha` in the numerator originate from the different factors
:math:`(1+\alpha \hat{H})` on the bra and ket sides, and thus one of
them is estimated by :math:`h_{2(11)}` and the other by
:math:`h_{2(20)}`.

The stationary condition
:math:`\partial E_{LS}(\alpha)/\partial \alpha=0` reduces to the
following quadratic equation for :math:`\alpha`:

.. math::

   \begin{aligned}
   \left[2 h_1 h_{3(12)} - h_{2(11)}\left(h_{2(11)}+h_{2(20)}\right)\right] \alpha^2
   + 2\left[h_{3(12)} - h_1 h_{2(11)}\right] \alpha
   + \left(h_{2(11)}+h_{2(20)}\right) - 2h_1^2 = 0.
   \end{aligned}

We evaluate :math:`E_{LS}(\alpha_{\pm})` for the two real roots
:math:`\alpha_{\pm}` of this quadratic equation and adopt the one that
gives the lower energy as the optimized :math:`\alpha` (if the
discriminant is negative and no real root exists, the Lanczos
calculation stops with an error). The variance of the energy for the
optimized :math:`|\phi_1\rangle` is also calculated as

.. math::

   \begin{aligned}
   \sigma^2_{LS}(\alpha) = \frac{\langle \phi_1| \hat{H}^2 |\phi_1\rangle}{\langle \phi_1|\phi_1\rangle} - E_{LS}(\alpha)^2,~~~~
   \frac{\langle \phi_1| \hat{H}^2 |\phi_1\rangle}{\langle \phi_1|\phi_1\rangle}
   = \frac{h_{2(11)} + 2\alpha h_{3(12)} + \alpha^2 h_{4(22)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}}.
   \end{aligned}

Second power-Lanczos step
~~~~~~~~~~~~~~~~~~~~~~~~~

Variational determination of :math:`\alpha` and :math:`\beta`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For ``NLanczosStep=2``, mVMC minimizes the variational energy

.. math::

   E_2(\alpha,\beta)
   =\frac{\langle\phi_2|\hat H|\phi_2\rangle}
          {\langle\phi_2|\phi_2\rangle}

of the trial wave function
:math:`|\phi_2\rangle=(1+\alpha\hat H+\beta\hat H^2)|\psi\rangle`
with respect to :math:`\alpha,\beta`. This is the direct extension of
determining :math:`\alpha` from the stationary condition of
:math:`E_{LS}(\alpha)` in the first step.

The minimization can be written as a generalized eigenvalue problem in the
Krylov basis

.. math::

   |v_a\rangle=\hat H^a|\psi\rangle,\qquad a=0,1,2.

Define the overlap and Hamiltonian matrices by

.. math::

   S_{ab}=\frac{\langle v_a|v_b\rangle}{\langle\psi|\psi\rangle},
   \qquad
   {\cal H}_{ab}
   =\frac{\langle v_a|\hat H|v_b\rangle}{\langle\psi|\psi\rangle}.

From the definition of :math:`|v_a\rangle`, these matrix elements are
nothing but the moments of the Hamiltonian with respect to the base wave
function:

.. math::

   S_{ab}=\langle \hat H^{a+b}\rangle_\psi,
   \qquad
   {\cal H}_{ab}=\langle \hat H^{a+b+1}\rangle_\psi.

The second step therefore requires the moments up to
:math:`\langle \hat H^5\rangle_\psi` (up to
:math:`\langle \hat H^6\rangle_\psi` when the variance described below
is also evaluated). For :math:`\boldsymbol{c}=(c_0,c_1,c_2)^{\mathsf T}`,

.. math::

   |\phi_2\rangle=\sum_{a=0}^{2}c_a|v_a\rangle,\qquad
   E_2(\boldsymbol{c})
   =\frac{\boldsymbol{c}^\dagger{\cal H}\boldsymbol{c}}
          {\boldsymbol{c}^\dagger S\boldsymbol{c}}.

The stationary condition is therefore

.. math::

   {\cal H}\boldsymbol{c}=E_2 S\boldsymbol{c}.

The eigenvector associated with the lowest eigenvalue gives the optimal
coefficients. Its overall scale is arbitrary; setting
:math:`\alpha=c_1/c_0` and :math:`\beta=c_2/c_0` recovers the
:math:`(1+\alpha\hat H+\beta\hat H^2)|\psi\rangle` notation above.
This generalized eigenvalue problem is the principle of the second
power-Lanczos step. Note that fixing :math:`c_2=0` and optimizing only
:math:`c_1/c_0` reduces this generalized eigenvalue problem to the
quadratic equation of the first step. In other words, the second step
applies the same variational principle with the dimension of the Krylov
subspace extended from two to three.

Numerically stable implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following operations do not change the variational principle. They
evaluate the same matrix elements from Monte Carlo samples and solve the
generalized eigenvalue problem more stably.

First, mVMC applies the local-estimator definition above to
:math:`\hat A=\hat H^a` and samples the four local powers

.. math::

   F_a(x)\equiv F(x,\hat H^a)
   =\frac{\langle x|\hat H^a|\psi\rangle}
          {\langle x|\psi\rangle},
   \qquad a=0,1,2,3.

For every sampled configuration :math:`x`, mVMC retains these four
quantities. Here, :math:`F_0(x)=1`, :math:`F_1(x)` is the usual local
energy, and :math:`F_2(x)` and :math:`F_3(x)` are the local estimators of
:math:`\hat H^2` and :math:`\hat H^3`, respectively. They are not trial
wave functions; they are sample-local quantities used to construct the
overlap matrix :math:`S`, Hamiltonian matrix :math:`{\cal H}`, and the
squared-Hamiltonian matrix used to calculate the variance.

After the sample loop, the globally weighted energy of the base VMC wave
function, :math:`E_{\rm shift}=\langle\hat H\rangle_\psi`, is computed,
and :math:`X=\hat H-E_{\rm shift}` is introduced. Because

.. math::

   {\rm span}\{|\psi\rangle,\hat H|\psi\rangle,\hat H^2|\psi\rangle\}
   =
   {\rm span}\{|\psi\rangle,X|\psi\rangle,X^2|\psi\rangle\},

this basis transformation changes neither the variational space nor its
optimum. It is introduced only for numerical stability. Before any
moment products are accumulated, each sample is transformed to

.. math::

   \widetilde F_a(x)
   =\sum_{j=0}^{a}{a\choose j}(-E_{\rm shift})^{a-j}F_j(x).

Thus :math:`\widetilde F_a` is the local estimator of :math:`X^a`, and
mVMC forms the centered moment matrix

.. math::

   \widetilde M_{ab}
   =\sum_x \rho(x)\widetilde F_a(x)^\dagger\widetilde F_b(x).

The implementation retains all combinations with
:math:`a,b=0,\ldots,3`. As in :math:`h_{n(ij)}` for the first step,
each :math:`\widetilde M_{ab}` estimates the :math:`(a+b)`-th moment
:math:`\langle X^{a+b}\rangle_\psi` by splitting the powers between the
bra :math:`\widetilde F_a^\dagger` and the ket
:math:`\widetilde F_b`. Because neither side has an exponent larger
than three, :math:`\widetilde M_{33}` represents the moment through
:math:`X^6` while requiring local powers only through :math:`F_3`.
Centering before accumulation avoids the precision loss that would occur
when raw moments such as :math:`\langle H^6\rangle` are formed for a
large extensive energy.  It requires storage proportional to four local
powers per Monte Carlo sample.
Centering primarily removes the offset from the mean energy.  When
:math:`|\langle H\rangle|` is small compared with the width of the energy
distribution, the dynamic range of the high-order local powers may remain
large; centering alone therefore does not guarantee numerical accuracy for
every system.

Because :math:`\widetilde M_{ab}` and
:math:`\widetilde M_{ba}^\dagger` are accumulated from the same
samples, their contributions are equal sample by sample; finite sampling
itself does not create a difference between them. Small anti-Hermitian
components can nevertheless arise from floating-point arithmetic and
parallel reductions, so mVMC projects the moment matrix onto its
Hermitian part. In contrast, finite-sample statistical error affects
different splittings :math:`\widetilde M_{ab}` and
:math:`\widetilde M_{cd}` with :math:`a+b=c+d`; this difference is
reported as ``hankel_residual``. mVMC then constructs the overlap,
Hamiltonian, and squared-Hamiltonian matrices in the centered Krylov basis:

.. math::

   \widetilde S_{ab}=\widetilde M_{ab},\qquad
   \widetilde{\cal H}_{ab}
   =\frac{\widetilde M_{a,b+1}+\widetilde M_{a+1,b}}{2},\qquad
   \widetilde G_{ab}=\widetilde M_{a+1,b+1},

where :math:`a,b=0,1,2`. Averaging the two moment estimates in
:math:`\widetilde{\cal H}` also preserves Hermiticity at finite sample
size.

To reduce differences among the norms of the basis vectors, mVMC applies
the diagonal scaling

.. math::

   D_{aa}=\frac{1}{\sqrt{\widetilde S_{aa}}}

and solves

.. math::

   (D\widetilde{\cal H}D)\boldsymbol{u}
   =E_X(D\widetilde S D)\boldsymbol{u}.

Only if the overlap factorization fails, a small regularization is added
to the diagonal of :math:`D\widetilde S D` before retrying. Centering,
Hermitian projection, diagonal scaling, and this fallback regularization
are numerical stabilization steps separate from the variational
principle.

The lowest eigenvalue gives the energy of the original Hamiltonian as
:math:`E=E_X+E_{\rm shift}`. If
:math:`\boldsymbol{d}=(d_0,d_1,d_2)^{\mathsf T}` denotes the coefficients
in the centered basis, the unshifted coefficients are

.. math::

   \begin{aligned}
   c_0 &= d_0-E_{\rm shift}d_1+E_{\rm shift}^2d_2,\\
   c_1 &= d_1-2E_{\rm shift}d_2,\\
   c_2 &= d_2.
   \end{aligned}

mVMC reports :math:`c_0,c_1,c_2` and reports the
:math:`\alpha,\beta` defined above as :math:`\alpha_1,\alpha_2`,
respectively.
With :math:`\boldsymbol{d}^\dagger\widetilde S\boldsymbol{d}=1`, the
variance is evaluated in the centered basis as

.. math::

   \sigma_E^2
   =\boldsymbol{d}^\dagger\widetilde G\boldsymbol{d}
    -(\boldsymbol{d}^\dagger\widetilde{\cal H}\boldsymbol{d})^2.

Computational cost
^^^^^^^^^^^^^^^^^^

The new cost introduced by the second step is the evaluation of the
local powers :math:`F_2(x)` and :math:`F_3(x)` for every sample, among
which
:math:`F_3(x)=\langle x|\hat H^3|\psi\rangle/\langle x|\psi\rangle`
dominates.

In the current second-step implementation, the Hamiltonian is limited to
spin-conserving ``Transfer`` terms and number-operator-type diagonal
interactions. The
evaluation of :math:`F_3` runs nested outer and inner loops over the
off-diagonal ``Transfer`` terms. Since the product of two ``Transfer``
terms combines into the two-body operator

.. math::

   CACA \equiv
   c_{i\sigma}^\dagger c_{j\sigma} c_{k\tau}^\dagger c_{l\tau},

the double loop generates :math:`O(N_{\mathrm{Transfer}}^2)` ``CACA``
operators, and for each of them the matrix element containing the
remaining :math:`\hat H`,
:math:`\langle \psi|\hat H\,CACA|x\rangle/\langle \psi|x\rangle`, must
be evaluated.

With the identity quantum projection (``NQPFull=1``), this matrix
element factorizes: ``CACA`` is applied to the sampled configuration
:math:`x` to form a temporary configuration :math:`x'`, and the local
estimator of :math:`\hat H` at :math:`x'` is multiplied by the overlap
ratio :math:`\langle \psi|CACA|x\rangle/\langle \psi|x\rangle`. The
double loop generates :math:`O(N_{\mathrm{Transfer}}^2)` ``CACA``
operators, but the Hamiltonian local estimator for each ``CACA`` has
another ``Transfer`` loop. Therefore, the worst-case complexity of this
path is also :math:`O(N_{\mathrm{Transfer}}^3)`. Its prefactor is
smaller than that of the nontrivial-projection path because it can use
the ordinary local estimator and a rank-two update and avoids the
``GreenFuncN`` work proportional to ``NQPFull``.

For a nontrivial quantum projection (``NQPFull>1``), this factorization
does not preserve the sum over quantum-projection components. The second
power-Lanczos implementation therefore expands
:math:`\hat H=\sum_\mu h_\mu\hat O_\mu` and treats
:math:`\hat O_\mu\,CACA` as one many-body operator for every term.
``GreenFuncN`` evaluates the corresponding many-body Green function from
the original configuration :math:`x`, including the quantum-projection
components, and the results are accumulated with weights :math:`h_\mu`.
This direct evaluation introduces one additional internal loop over
``Transfer`` terms for each :math:`\hat O_\mu\,CACA` many-body Green
function, so the dominant operator count of this path grows as
:math:`O(N_{\mathrm{Transfer}}^3)`, and the ``GreenFuncN`` work is also
proportional to ``NQPFull``. Wall-clock cost therefore depends strongly
on the projection as well as the wave function and execution
environment.

Calculation of physical quantities after the first power-Lanczos step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For ``NLanczosMode=2``, by using the optimized parameter
:math:`\alpha`, the physical quantities (one-body and two-body Green
functions) are calculated for the first-step wave function
:math:`|\phi_1\rangle` (Green functions for the second step are not
supported). The numerator of the expectation value of an operator
:math:`\hat{A}` is expanded with respect to :math:`\alpha` as

.. math::

   \begin{aligned}
   \langle \phi_1| \hat{A} |\phi_1\rangle = \langle \psi| \hat{A} |\psi\rangle
   + \alpha \langle \psi| \hat{H}\hat{A} |\psi\rangle
   + \alpha \langle \psi| \hat{A}\hat{H} |\psi\rangle
   + \alpha^2 \langle \psi| \hat{H}\hat{A}\hat{H} |\psi\rangle.
   \end{aligned}

As in the case of the energy, each term is estimated by splitting it
into the bra part :math:`\hat{H}^i` and the ket part
:math:`\hat{A}\hat{H}^j`:

.. math::

   \begin{aligned}
   &A_0 =\sum_{x} \rho(x) F(x,  \hat{A}),\\
   &A_{1(10)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{A}),\\
   &A_{1(01)}=\sum_{x} \rho(x) F(x, \hat{A}\hat{H}),\\
   &A_{2(11)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H})F(x,  \hat{A}\hat{H}).
   \end{aligned}

Here, as in the case of the energy, the subscript :math:`n(ij)` denotes
the total power :math:`n=i+j` of :math:`\hat{H}` and its assignment to
the bra side (:math:`\hat{H}^i`) and the ket side
(:math:`\hat{A}\hat{H}^j`). With these estimators, the expectation
value is calculated as

.. math::

   \begin{aligned}
   A_{LS}(\alpha) =\frac{\langle \phi_1| \hat{A} |\phi_1\rangle}{\langle \phi_1|\phi_1\rangle}=\frac{A_0 + \alpha(A_{1(10)} + A_{1(01)}) + \alpha^2 A_{2(11)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}},
   \end{aligned}

where the denominator is the same normalization factor
:math:`\langle \phi_1 | \phi_1 \rangle` as that for the energy. In the
program, the one-body and two-body Green functions are calculated based
on this expression.

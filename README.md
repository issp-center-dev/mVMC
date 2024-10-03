# mVMC

A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method

### What is mVMC ?

mVMC (many-variable Variational Monte Carlo method)
is a software for performing the highly-accurate 
variational Monte Carlo calculations
with the simple and flexible user interface.
mVMC also supports the large-scale parallelization.
For the conventional models in strongly correlated electron systems such as the Hubbard model, the Heisenberg model, and the Kondo-lattice model,
users can perform the calculation by preparing the one input files whose length is shorter than ten lines.  
By using the same input file, users can perform the exact diagonalization through [HPhi](https://github.com/QLMS/HPhi/releases).
Thus, it is easy to examine the accuracy of the variational calculation for small system sizes
and to perform the calculations 
for large system sizes that can not be treated 
by the exact diagonalization.
A broad spectrum of users including experimental scientists is cordially welcome.


### Methods

many-variable variational Monte Carlo method


### Target models

Hubbard model, Heisenberg model, Kondo lattice model, multi-orbital Hubbard model

### Available physical quantities

ground-state energy, spin/charge structure factors, superconducting correlations, etc.


## Requirement

- C compiler (intel, Fujitsu, GNU, etc. ) 
- ScaLAPACK library (intel MKL, Fujitsu, ATLAS, etc.) 
- MPI library

## Install

You can install mVMC and also get a manual for mVMC from a [release note](https://github.com/issp-center-dev/mVMC/releases).


## Licence

GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)). 

The mVMC package is developed based on the [mVMC-mini](https://github.com/fiber-miniapp/mVMC-mini) program. The license of mVMC-mini is "The BSD 3-Clause License".

We would appreciate if you cite the following article in your research with mVMC:  
mVMC - Open-source software for many-variable variational Monte Carlo method, Takahiro Misawa, Satoshi Morita, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Yuichi Motoyama, Kota Ido, Takahiro Ohgoe, Masatoshi Imada, Takeo Kato, 
[Computer Physics Communications, 235, 447-462 (2019)](https://www.sciencedirect.com/science/article/pii/S0010465518303102)


## Tutorials

Lecture notes and sample scripts used in Hands-on
are available at [mVMC-tutorial](https://github.com/issp-center-dev/mVMC-tutorial)

## Authors

Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, RuQing G. Xu, Yuichi Motoyama, Mitsuaki Kawamura, Kazuyoshi Yoshimi, Takeo Kato, Masatoshi Imada.

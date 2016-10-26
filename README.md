# mVMC

A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method

====
mVMC(many-variable variational Monte Carlo method)
is a software for performing the high-accurate 
variational Monte Carlo calculations
with the simple and flexible user interface.
mVMC also supports the massively parallel calculations.
For the conventional model in strongly correlated electron systems
such as the Hubbard model, the Heisenberg model, and the Kondo-lattice model,
users can perform the calculation by preparing the one input files whose length is
shorter the ten lines.
By using the same input file,
user can perform the exact diagonalization
through [HPhi](https://github.com/QLMS/HPhi/releases).
Thus, it is easy to perform the calculations 
for large system sizes that can not be treated 
by the exact diagonalization after examining the accuracy of the
variational calculation for small system sizes.
A broad spectrum of users including experimental scientists and 
scientists of quantum chemistry is cordially welcome.



### Methods
many-variable variational Monte Carlo method


### Target models
Hubbard model, Heisenberg model, Kondo lattice model, multi-orbital Hubbard model

### Available physical quantities
specific heat, susceptibility, ground state energy, structure factors


## Requirement
- C compiler (intel, Fujitsu, GNU, etc. ) 
- ScaLAPACK library (intel MKL, Fujitsu, ATLAS, etc.) 
- MPI library

## Install

You can install mVMC and also get a manual for mVMC from a [release note](https://github.com/issp-center-dev/mVMC/releases).


## Licence

The distribution of the program package and the source codes for mVMC follow GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)). 

The mVMC package is developed based on the [mVMC-mini](https://github.com/fiber-miniapp/mVMC-mini) program which follows "The BSD 3-Clause License".

## Author
Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido,  Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.
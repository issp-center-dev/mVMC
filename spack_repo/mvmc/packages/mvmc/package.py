# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *


class Mvmc(CMakePackage):
    """mVMC: many-variable Variational Monte Carlo method for quantum lattice models."""

    homepage = "https://www.pasums.issp.u-tokyo.ac.jp/mvmc/en/"
    git = "https://github.com/issp-center-dev/mVMC.git"

    maintainers("issp-center-dev")

    license("GPL-3.0", checked_by="freifrauvonbleifrei")

    version("1.3.0", tag="v1.3.0", submodules=True)
    version("1.2.0", tag="v1.2.0", submodules=True)
    version("develop", branch="develop", submodules=True)

    variant("scalapack", default=False, description="Enable ScaLAPACK support")
    variant(
        "pfaffian_blocked",
        default=False,
        description="Use blocked Pfaffian",
    )
    variant(
        "gemmt",
        default=True,
        description="Use GEMMT",
    )
    variant("shared", default=True, description="Build shared libraries")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")

    depends_on("cmake@3.5:", type="build")
    depends_on("mpi")
    depends_on("lapack")
    depends_on("blas")
    depends_on("scalapack", when="+scalapack")

    def cmake_args(self):
        args = [
            self.define_from_variant("USE_SCALAPACK", "scalapack"),
            self.define_from_variant("PFAFFIAN_BLOCKED", "pfaffian_blocked"),
            self.define_from_variant("USE_GEMMT", "gemmt"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define("GIT_SUBMODULE_UPDATE", False),
            self.define("Testing", self.run_tests),
        ]
        return args

Release notes
=============

Unreleased
----------

Power-Lanczos estimator selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Power-Lanczos inputs now select the corrected estimator by default. This
includes existing inputs that set ``NLanczosMode`` and omit the new
``NLanczosEstimatorMode`` keyword. ``NLanczosMode=1`` runs the full-support,
scale-normalized corrected energy/variance path with independent coefficient
and final chains and writes one compact stabilization JSON. This is numerical
stabilization evidence, not a release-statistics certification. Corrected
``NLanczosMode=2`` remains outside scope and is rejected before allocation;
additional-observable production enable is unchanged.

Corrected launchers must provide exact source, input, binary, environment, and
seed identity variables documented in :ref:`HowToExpert`. The fixed default
policy uses a 4096+4096 scale pilot, 4096 warmup and 16384 saved configurations
for each scored stage, 32 blocks, and no sample-level trace persistence.

To reproduce the previous biased base-support estimator temporarily, add the
following line to the ModPara file::

   NLanczosEstimatorMode 1

This legacy route emits a warning. Its output is diagnostic-only and must not
be treated as a corrected release result. All other new P6 estimator controls
must remain 0. See :ref:`HowToExpert` for the full keyword contract.

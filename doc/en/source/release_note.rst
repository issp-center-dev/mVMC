Release notes
=============

Unreleased
----------

Power-Lanczos estimator selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Power-Lanczos inputs now select the corrected estimator by default. This
includes existing inputs that set ``NLanczosMode`` to 1 or 2 and omit the new
``NLanczosEstimatorMode`` keyword. The corrected production route remains
disabled during the staged rollout, so such inputs fail fast with exit status
20 and an availability diagnostic. For ``NLanczosMode=2``, this happens before
the P6 observable census and prevents a census or resource-limit error from
masking the unavailable corrected route.

To reproduce the previous biased base-support estimator temporarily, add the
following line to the ModPara file::

   NLanczosEstimatorMode 1

This legacy route emits a warning. Its output is diagnostic-only and must not
be treated as a corrected release result. All other new P6 estimator controls
must remain 0. See :ref:`HowToExpert` for the full keyword contract.

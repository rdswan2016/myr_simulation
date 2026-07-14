"""Standalone mass-transfer / kinetic-fitting engine for pH-tracked titration data.

Modules:
  geometry         -- vessel/zone volume construction (2, 3, 4 compartments)
  chemistry        -- Henderson-Hasselbalch thermodynamics + titrant feed model
  transport_model   -- N-compartment convective transport ODE + least-squares fitting
  model_selection   -- automatic compartment-count selection (nested F-test, AIC/BIC,
                        identifiability, mass-balance ceiling) + plain-language rationale
  io_utils         -- flexible CSV ingestion
"""

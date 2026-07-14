# Mass-Transfer & Kinetic Fitting Tool

A standalone, local web app that fits 2-, 3-, and 4-compartment convective
mass-transfer models (coupled to Henderson-Hasselbalch chemistry) to pH-tracked
acid/base titration data, and automatically recommends a compartment count with
a plain-language rationale. Runs entirely on your machine, served over
`https://localhost` with a self-signed certificate. No external services, no
account, no network calls once dependencies are installed.

## Repository layout

```
webapp/
  app.py                    Streamlit UI (3 tabs: checklist, upload/configure, results)
  generate_cert.py          creates a self-signed TLS cert/key pair (pure Python)
  requirements.txt
  .streamlit/config.toml    HTTPS wiring + brand theme
  engine/
    geometry.py              vessel/zone volume construction (2/3/4 compartments)
    chemistry.py              Henderson-Hasselbalch + titrant feed model
    transport_model.py        N-compartment ODE engine + least-squares fitting
    model_selection.py        automatic compartment-count selection + rationale
    io_utils.py                flexible CSV ingestion
```

## Setup (one time)

```bash
cd Z:\Surya\Code\Python\myr_simulation\webapp
python -m venv .venv
.venv\Scripts\activate          # Windows.  Use: source .venv/bin/activate  on macOS/Linux
pip install -r requirements.txt
python generate_cert.py
```

`generate_cert.py` writes `certs/cert.pem` and `certs/key.pem`. These are
self-signed and machine-specific -- they are git-ignored; every user/machine
regenerates their own.

## Run

```bash
streamlit run app.py
```

Then open **https://localhost:8501**. Your browser will show a certificate
warning on first visit (expected for a self-signed, local-only certificate) --
proceed past it (e.g. "Advanced -> Proceed to localhost"). The HTTPS
certificate/key paths and the brand theme are both picked up automatically
from `.streamlit/config.toml` in this directory, so no extra CLI flags are
needed as long as you run `streamlit run app.py` **from inside `webapp/`**.

To stop the server: `Ctrl+C` in the terminal running it.

## Regenerating the certificate

Certificates are valid ~825 days. To force a fresh one (e.g. after the old one
expires, or on a new machine):

```bash
del certs\cert.pem certs\key.pem     # Windows.  Use: rm certs/*.pem  on macOS/Linux
python generate_cert.py
```

## How the Henderson-Hasselbalch loop couples to the transport ODE engine

See the full docstring at the top of `engine/transport_model.py` for the exact
per-timestep coupling; summary:

```
state per zone i:  (n_base_i, n_acid_i)   <- moles, NOT concentrations

per ODE evaluation:
    C_base_i, C_acid_i = n_base_i / V_i(t), n_acid_i / V_i(t)     # concentration
    for each junction (i, i+1):
        flux = Q[i] * (C_x_i - C_x_{i+1})    for x in {base, acid}   # Fickian exchange
        move flux from zone i to zone i+1
    zone 0 only: titrant feed converts n_base_0 <-> n_acid_0 directly (no separate
                 reaction-rate term -- neutralization assumed instantaneous relative
                 to circulation, matching the validated notebook models)

after solve_ivp integrates transport only:
    pH_probe(t) = pKa + log10( n_base_0(t) / n_acid_0(t) )   # <- H-H loop closes HERE, once,
                                                                #    only on the probe zone's output
```

Chemistry never appears inside the ODE right-hand side itself -- it only
converts the feed into a mole transfer at the source zone, and converts the
final probe-zone mole trajectory into a pH curve for comparison against the
uploaded data. Everything else in the ODE is pure conservative transport.

## Compartment-count selection logic (engine/model_selection.py)

Fits all of 2/3/4 compartments, then:

1. **Nested extra-sum-of-squares F-test** (Qᵢ₊₁ vs Qᵢ) -- is the added
   compartment's improvement in fit statistically significant, not just
   numerically nonzero (which more parameters always produce)?
2. **Identifiability check** -- did the added flow rate converge to an interior
   value, or did it get pushed to the edge of its allowed range (a sign the
   data cannot actually pin it down)?
3. **Mass-balance ceiling check** -- independent of compartment count: does the
   data's peak/plateau pH exceed the closed-system equilibrium implied by the
   total titrant/analyte charge? If so, no compartment count can close that
   gap, and the tool says so explicitly instead of recommending more zones.

AIC/BIC are computed and displayed for reference, but the automatic
recommendation is driven by the three checks above -- with unweighted
least-squares fits, AIC/BIC alone can favor an over-parameterized model whose
extra parameter isn't actually identifiable, which is exactly the failure mode
this tool is designed to catch and flag rather than silently act on.

## Validation

This engine was cross-checked against a hand-built, independently-validated
3-compartment fit for a real 160 RPM Tris/HCl titration run (see
`5L_tank_transfer_mass_3_compartments.ipynb` in the parent directory): the
general N-compartment engine here reproduces that notebook's fitted flow rates
to within numerical tolerance for 2, 3, and 4 compartments alike, and its
automatic recommendation (2 compartments; mass-balance ceiling exceeded, so no
transport-only model explains the full trace) matches the conclusion reached
manually across that notebook's extended analysis.

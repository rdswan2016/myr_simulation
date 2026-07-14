"""
Standalone web application: automated multi-compartment mass-transfer parameter
estimation and kinetic fitting from pH-tracked titration data.

Run (from this directory, after `python generate_cert.py`):
    streamlit run app.py
(the .streamlit/config.toml in this same directory wires up HTTPS automatically)
"""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from engine.chemistry import ChemistryParams, FeedParams
from engine.geometry import VesselGeometry, probe_zone_check
from engine.io_utils import load_ph_csv
from engine.model_selection import run_model_selection
from engine.transport_model import simulate_ph

# ---------------------------------------------------------------------------
# Brand palette (fixed roles per spec)
# ---------------------------------------------------------------------------
NAVY = "#003049"     # primary structures & text
RED = "#D62828"      # alerts / critical indicators
ORANGE = "#F77F00"   # secondary accents / dividers
YELLOW = "#FCBF49"   # tables / container backgrounds
DATA_COLOR = RED
FIT_COLORS = {2: NAVY, 3: ORANGE, 4: "#7A5C00"}  # 4-comp: darkened yellow for line-contrast

st.set_page_config(page_title="Mass-Transfer Fitting Tool", page_icon="🧪", layout="wide")

st.markdown(
    f"""
    <style>
    h1, h2, h3, h4, .stMarkdown p {{ color: {NAVY}; }}
    .block-container {{ padding-top: 2rem; }}
    .checklist-box {{
        background-color: {YELLOW}33;
        border-left: 6px solid {ORANGE};
        border-radius: 6px;
        padding: 1rem 1.2rem;
        margin-bottom: 1rem;
    }}
    .alert-box {{
        background-color: {RED}18;
        border-left: 6px solid {RED};
        border-radius: 6px;
        padding: 1rem 1.2rem;
        margin-bottom: 1rem;
        color: {NAVY};
    }}
    .ok-box {{
        background-color: #00830018;
        border-left: 6px solid #008300;
        border-radius: 6px;
        padding: 1rem 1.2rem;
        margin-bottom: 1rem;
        color: {NAVY};
    }}
    .stButton>button {{
        background-color: {ORANGE};
        color: white;
        font-weight: 600;
        border: none;
    }}
    .stButton>button:hover {{ background-color: {RED}; color: white; }}
    div[data-testid="stMetric"] {{
        background-color: {YELLOW}33;
        border-radius: 8px;
        padding: 0.6rem;
    }}
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🧪 Multi-Compartment Mass-Transfer & Kinetic Fitting Tool")
st.caption(
    "Automated compartment-count selection and Qᵢⱼ flow-rate estimation from "
    "pH-tracked acid/base titration data in unbaffled 5 L mixing vessels."
)

tab_home, tab_data, tab_results = st.tabs(["📋 Start Here", "📤 Upload & Configure", "📊 Results"])

# ===========================================================================
# TAB 1: Data Preparation Checklist
# ===========================================================================
with tab_home:
    st.markdown("## Data Preparation Checklist")
    st.write(
        "Before running an analysis, collect the following. The tool cannot substitute "
        "for any of these -- it fits *flow rates*, not the underlying hardware/chemistry facts."
    )

    st.markdown(
        f"""
        <div class="checklist-box">
        <h4>1. CSV file structure</h4>
        <ul>
          <li><b>Time column</b> -- elapsed time since the titrant pump started (seconds or minutes; name it something containing "time"/"sec"/"min").</li>
          <li><b>pH column</b> -- probe reading at each time point (name containing "pH").</li>
          <li><em>Optional</em>: a cumulative titrant volume column, if logged directly.</li>
          <li>One row per timestamp; no merged headers or units embedded inside data cells.</li>
        </ul>
        </div>
        <div class="checklist-box">
        <h4>2. Fixed hardware variables</h4>
        <ul>
          <li>Impeller rotation speed (RPM)</li>
          <li>Tank (vessel) inner diameter</li>
          <li>Impeller diameter</li>
          <li>Impeller height off the tank floor</li>
          <li>Initial liquid fill height</li>
        </ul>
        </div>
        <div class="checklist-box">
        <h4>3. Chemical parameters</h4>
        <ul>
          <li>Stock concentration of the analyte already in the vessel (the buffering species being titrated)</li>
          <li>Stock concentration of the titrant being pumped in</li>
          <li>Whether the titrant is an <b>acid</b> (adds H⁺) or a <b>base</b> (consumes H⁺)</li>
          <li>The species' pKₐ</li>
          <li>Total titrant volume delivered and the continuous pump feed rate</li>
        </ul>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        ### What this tool does
        1. Fits **2-, 3-, and 4-compartment** convective transport models (chained zones,
           Henderson-Hasselbalch chemistry) to your uploaded pH trace, using
           `scipy.optimize.least_squares`.
        2. Runs a **nested statistical test** to check whether each extra compartment
           actually earns its keep, rather than just always picking the model with the
           lowest raw error.
        3. Flags when a fitted flow rate is **not identifiable** (converges to a boundary
           instead of a real value) and when the data's own **mass-balance ceiling** rules
           out *any* transport-only explanation for part of the curve.
        4. Recommends a compartment count with a **plain-language rationale**, not just a number.

        ### What this tool does NOT do
        It does not replace judgment about vessel hydrodynamics, does not know your specific
        impeller's true flow number (it uses a generous literature ceiling, Nₐ ≤ 1.0, as an
        upper bound only), and cannot detect non-transport causes (CO₂ outgassing, reagent
        lot variation, electrode drift) beyond flagging that they are more likely than a
        transport explanation once the mass-balance ceiling is exceeded.
        """
    )

# ===========================================================================
# TAB 2: Upload & Configure
# ===========================================================================
with tab_data:
    st.markdown("## 1. Upload your pH tracking CSV")
    uploaded = st.file_uploader("CSV file", type=["csv"])

    df = None
    if uploaded is not None:
        try:
            df = load_ph_csv(uploaded)
            st.success(f"Loaded {len(df)} rows. Detected columns mapped to time (s) and pH.")
            c1, c2 = st.columns([2, 1])
            with c1:
                fig_preview = go.Figure()
                fig_preview.add_trace(go.Scatter(x=df["t_sec"], y=df["pH"], mode="markers+lines",
                                                  marker=dict(color=DATA_COLOR, size=5),
                                                  line=dict(color=DATA_COLOR, width=1)))
                fig_preview.update_layout(
                    title="Raw data preview", xaxis_title="Time (s)", yaxis_title="pH",
                    template="plotly_white", height=320, margin=dict(l=10, r=10, t=40, b=10),
                )
                st.plotly_chart(fig_preview, width="stretch")
            with c2:
                st.dataframe(df.head(10), height=280)
            st.session_state["df"] = df
        except ValueError as e:
            st.markdown(f'<div class="alert-box"><b>File problem:</b> {e}</div>', unsafe_allow_html=True)

    st.markdown("## 2. Hardware, chemistry & feed parameters")
    with st.form("config_form"):
        st.markdown("#### Vessel & impeller geometry (cm)")
        g1, g2, g3, g4 = st.columns(4)
        D_tank_cm = g1.number_input("Tank diameter", value=20.5, min_value=1.0, step=0.1)
        D_impeller_cm = g2.number_input("Impeller diameter", value=7.5, min_value=0.5, step=0.1)
        h_impeller_cm = g3.number_input("Impeller height off floor", value=8.0, min_value=0.1, step=0.1)
        h_liquid_initial_cm = g4.number_input("Initial liquid height", value=17.2, min_value=1.0, step=0.1)

        h1, h2, h3 = st.columns(3)
        N_RPM = h1.number_input("Impeller speed (RPM)", value=160.0, min_value=1.0, step=1.0)
        probe_height_cm = h2.number_input("pH probe height (cm, optional check)", value=12.6, min_value=0.0, step=0.1)
        Nq_max = h3.number_input("Physical Nₐ ceiling (advanced)", value=1.0, min_value=0.1, max_value=5.0, step=0.1,
                                  help="Generous upper bound on any real impeller's flow number; used only to bound the optimizer, not asserted as the true value.")

        st.markdown("#### Chemistry")
        c1, c2, c3 = st.columns(3)
        pKa = c1.number_input("pKₐ", value=8.06, min_value=0.0, max_value=14.0, step=0.01)
        C_analyte_stock_M = c2.number_input("Analyte stock conc. (mol/L)", value=0.211237, min_value=0.0, format="%.6f")
        C_titrant_stock_M = c3.number_input("Titrant stock conc. (mol/L)", value=1.0, min_value=0.0, format="%.6f")
        titrant_kind = st.radio("Titrant type", ["Acid (delivers H⁺)", "Base (consumes H⁺)"], horizontal=True)

        st.markdown("#### Continuous feed")
        f1, f2 = st.columns(2)
        V_titrant_total_mL = f1.number_input("Total titrant volume (mL)", value=266.0, min_value=0.0, step=1.0)
        pump_rate_mL_min = f2.number_input("Pump feed rate (mL/min)", value=110.0, min_value=0.01, step=1.0)

        submitted = st.form_submit_button("Save configuration")

    if submitted or "config" in st.session_state:
        config = dict(
            D_tank_cm=D_tank_cm, D_impeller_cm=D_impeller_cm, h_impeller_cm=h_impeller_cm,
            h_liquid_initial_cm=h_liquid_initial_cm, N_RPM=N_RPM, probe_height_cm=probe_height_cm,
            Nq_max=Nq_max, pKa=pKa, C_analyte_stock_M=C_analyte_stock_M,
            C_titrant_stock_M=C_titrant_stock_M, titrant_delivers_acid=titrant_kind.startswith("Acid"),
            V_titrant_total_mL=V_titrant_total_mL, pump_rate_mL_min=pump_rate_mL_min,
        )
        st.session_state["config"] = config

    if "config" in st.session_state:
        geom_preview = VesselGeometry(
            D_tank_cm=st.session_state["config"]["D_tank_cm"],
            D_impeller_cm=st.session_state["config"]["D_impeller_cm"],
            h_impeller_cm=st.session_state["config"]["h_impeller_cm"],
            h_liquid_initial_cm=st.session_state["config"]["h_liquid_initial_cm"],
        )
        from engine.geometry import build_zone_layout
        try:
            layout3 = build_zone_layout(geom_preview, 3)
            ok = probe_zone_check(layout3, st.session_state["config"]["probe_height_cm"])
            if not ok:
                st.markdown(
                    '<div class="alert-box">Warning: the pH probe height does not fall inside the '
                    "top (feed) zone under this geometry. The model assumes the probe reads the feed "
                    "zone's concentration -- double-check your probe height / liquid fill.</div>",
                    unsafe_allow_html=True,
                )
        except ValueError as e:
            st.markdown(f'<div class="alert-box">Geometry problem: {e}</div>', unsafe_allow_html=True)

    st.markdown("## 3. Run analysis")
    can_run = ("df" in st.session_state) and ("config" in st.session_state)
    if not can_run:
        st.info("Upload a CSV and save a configuration above before running.")
    if st.button("▶ Run compartment analysis", disabled=not can_run):
        cfg = st.session_state["config"]
        df_run = st.session_state["df"]
        geom = VesselGeometry(D_tank_cm=cfg["D_tank_cm"], D_impeller_cm=cfg["D_impeller_cm"],
                               h_impeller_cm=cfg["h_impeller_cm"], h_liquid_initial_cm=cfg["h_liquid_initial_cm"])
        chem = ChemistryParams(pKa=cfg["pKa"], C_analyte_stock_M=cfg["C_analyte_stock_M"],
                                C_titrant_stock_M=cfg["C_titrant_stock_M"],
                                titrant_delivers_acid=cfg["titrant_delivers_acid"])
        feed = FeedParams(V_titrant_total_mL=cfg["V_titrant_total_mL"], pump_rate_mL_min=cfg["pump_rate_mL_min"])
        N_rps = cfg["N_RPM"] / 60.0
        with st.spinner("Fitting 2-, 3-, and 4-compartment models..."):
            result = run_model_selection(geom, N_rps, chem, feed,
                                          df_run["t_sec"].to_numpy(), df_run["pH"].to_numpy(),
                                          Nq_max=cfg["Nq_max"])
        st.session_state["result"] = result
        st.session_state["result_geom"] = geom
        st.session_state["result_chem"] = chem
        st.session_state["result_feed"] = feed
        st.success("Analysis complete -- see the Results tab.")

# ===========================================================================
# TAB 3: Results
# ===========================================================================
with tab_results:
    if "result" not in st.session_state:
        st.info("Run an analysis from the 'Upload & Configure' tab first.")
    else:
        res = st.session_state["result"]
        geom = st.session_state["result_geom"]
        chem = st.session_state["result_chem"]
        feed = st.session_state["result_feed"]
        df_run = st.session_state["df"]

        st.markdown("## Recommendation")
        m1, m2, m3 = st.columns(3)
        m1.metric("Recommended compartments", res.recommended_n)
        m2.metric("Mass-balance ceiling (pH)", f"{res.ceiling_pH:.3f}")
        m3.metric("Data post-feed extreme (pH)", f"{res.data_peak_pH:.3f}")

        st.markdown(res.rationale)

        st.markdown("## Fit comparison")
        rows = []
        for n in (2, 3, 4):
            fit = res.fits[n]
            rows.append({
                "Compartments": n,
                "Flow rates fit (L/min)": ", ".join(f"Q{i+1}→{i+2}={q:.4f}{' *' if b else ''}"
                                                       for i, (q, b) in enumerate(zip(fit.Q_fit, fit.at_bound))),
                "RMSE (pH)": round(fit.rmse, 5),
                "AIC": round(res.aic[n], 2),
                "BIC": round(res.bic[n], 2),
            })
        table_df = pd.DataFrame(rows)
        st.dataframe(
            table_df.style.set_properties(**{"background-color": f"{YELLOW}22", "color": NAVY}),
            width="stretch", hide_index=True,
        )
        st.caption("`*` marks a flow rate that converged to its bound (not identifiable as a point value).")

        st.markdown("## Overlay: data vs. fitted models")
        t_smooth = np.linspace(df_run["t_sec"].min(), df_run["t_sec"].max(), 600)
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=df_run["t_sec"], y=df_run["pH"], mode="markers", name="Measured data",
                                  marker=dict(color=DATA_COLOR, size=5)))
        for n in (2, 3, 4):
            layout = res.layouts[n]
            fit = res.fits[n]
            pH_smooth = simulate_ph(fit.Q_fit, layout, chem, feed, t_smooth / 60.0, df_run["pH"].to_numpy()[0])
            fig.add_trace(go.Scatter(x=t_smooth, y=pH_smooth, mode="lines", name=f"{n}-compartment fit",
                                      line=dict(color=FIT_COLORS[n], width=2,
                                                dash="solid" if n == res.recommended_n else "dot")))
        fig.add_hline(y=res.ceiling_pH, line=dict(color=RED, dash="dash"),
                      annotation_text="Mass-balance ceiling", annotation_position="top left")
        fig.update_layout(template="plotly_white", xaxis_title="Time (s)", yaxis_title="pH",
                           height=500, legend=dict(orientation="h", y=-0.2))
        st.plotly_chart(fig, width="stretch")

        with st.expander("Zone geometry used (recommended model)"):
            layout = res.layouts[res.recommended_n]
            geom_rows = [{"Zone": lbl, "Boundary (cm)": f"{b:.2f}–{t:.2f}", "Volume (L)": f"{v:.3f}"}
                         for lbl, (b, t), v in zip(layout.labels, layout.boundaries_cm, layout.fixed_volumes_L)]
            st.table(pd.DataFrame(geom_rows))

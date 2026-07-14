"""CSV ingestion helpers: flexible column detection for non-engineer users who may
not name their columns exactly "t_sec"/"pH"."""
import io

import numpy as np
import pandas as pd

TIME_HINTS = ["time", "t (s", "t_s", "detik", "sec", "seconds", "t_min", "min"]
PH_HINTS = ["ph"]
VOL_HINTS = ["vol", "titrant", "ml", "volume"]


def _guess_column(columns, hints):
    for c in columns:
        lc = str(c).strip().lower()
        for h in hints:
            if h in lc:
                return c
    return None


def load_ph_csv(file_like) -> pd.DataFrame:
    """Read an uploaded CSV and return a DataFrame with columns guaranteed to be
    named 't_sec' and 'pH' (auto-detected; falls back to the first two numeric
    columns if no name hints match). Raises ValueError with a clear message on
    failure so the UI can surface it directly to a non-engineer user."""
    df_raw = pd.read_csv(file_like)
    if df_raw.shape[1] < 2:
        raise ValueError("CSV must have at least two columns (time and pH).")

    time_col = _guess_column(df_raw.columns, TIME_HINTS)
    ph_col = _guess_column(df_raw.columns, PH_HINTS)

    if time_col is None or ph_col is None:
        numeric_cols = [c for c in df_raw.columns if pd.api.types.is_numeric_dtype(df_raw[c])]
        if len(numeric_cols) < 2:
            raise ValueError(
                "Could not identify a time column and a pH column automatically. "
                "Please name your columns something containing 'time'/'sec' and 'pH'."
            )
        if time_col is None:
            time_col = numeric_cols[0]
        if ph_col is None:
            ph_col = numeric_cols[1] if numeric_cols[1] != time_col else numeric_cols[0]

    out = df_raw[[time_col, ph_col]].copy()
    out.columns = ["t_sec", "pH"]
    out["t_sec"] = pd.to_numeric(out["t_sec"], errors="coerce")
    out["pH"] = pd.to_numeric(out["pH"], errors="coerce")
    out = out.dropna().sort_values("t_sec").reset_index(drop=True)

    if len(out) < 5:
        raise ValueError("Fewer than 5 valid (time, pH) rows after cleaning -- check the file contents.")
    if out["pH"].between(0, 14).mean() < 0.9:
        raise ValueError("Detected pH column has values mostly outside 0-14 -- check column selection.")

    vol_col = _guess_column(df_raw.columns, VOL_HINTS)
    if vol_col is not None and vol_col not in (time_col, ph_col):
        out["titrant_vol_mL"] = pd.to_numeric(df_raw[vol_col], errors="coerce")

    return out

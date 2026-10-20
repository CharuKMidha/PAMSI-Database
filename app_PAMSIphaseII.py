import streamlit as st
import pandas as pd
import numpy as np
import zipfile
from pathlib import Path
import urllib.request
import plotly.graph_objects as go

st.set_page_config(page_title="PAMSI — Private Peptide Atlas", layout="wide")
st.title("🧬 PAMSI — Peptide Atlas for Mass Spectrometry Imaging")
st.caption("MALDI Imaging Mass Spectrometry • Extracellular Matrix Proteome • Phase 2 • Private Access")

# ================== PASSWORD ==================
if "password_correct" not in st.session_state:
    st.session_state.password_correct = False

def check_password():
    password = st.text_input("Enter Collaborator Password", type="password")
    if password == "PAMSI2026_Collab":
        st.session_state.password_correct = True
        st.rerun()
    elif password:
        st.error("Incorrect password")

if not st.session_state.password_correct:
    check_password()
    st.stop()

# ================== PATHS / URLS ==================
PSM_URL = "https://github.com/CharuKMidha/PAMSI-Database/releases/download/v1.0-data/HumanBreast_Phase2_psm_best.tsv"
FRAG_URL = "https://github.com/CharuKMidha/PAMSI-Database/releases/download/v1.0-data/HumanBreast_Phase2_fragments.tsv.zip"

DATA_DIR = Path("data")
PSM_LOCAL = DATA_DIR / "HumanBreast_Phase2_psm_best.tsv"
FRAG_ZIP_LOCAL = DATA_DIR / "HumanBreast_Phase2_fragments.tsv.zip"
FRAG_TSV_LOCAL = DATA_DIR / "HumanBreast_Phase2_fragments.tsv"

def download_with_ua(url: str, dest: Path):
    """Download from GitHub Release with a User-Agent (avoids some 403s)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": "PAMSI-Streamlit/1.0"})
    with st.spinner(f"Downloading {dest.name} (first run only)…"):
        with urllib.request.urlopen(req, timeout=600) as r, open(dest, "wb") as f:
            while True:
                chunk = r.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
    return dest

def ensure_fragments_tsv():
    """Download zip once and extract the TSV to disk (not into RAM as a full DataFrame)."""
    if FRAG_TSV_LOCAL.exists() and FRAG_TSV_LOCAL.stat().st_size > 0:
        return FRAG_TSV_LOCAL
    download_with_ua(FRAG_URL, FRAG_ZIP_LOCAL)
    with zipfile.ZipFile(FRAG_ZIP_LOCAL, "r") as zf:
        members = [m for m in zf.namelist() if m.lower().endswith((".tsv", ".csv", ".txt"))]
        if not members:
            members = zf.namelist()
        if not members:
            raise ValueError("Empty ZIP from GitHub Release")
        with zf.open(members[0]) as src, open(FRAG_TSV_LOCAL, "wb") as out:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
    return FRAG_TSV_LOCAL

@st.cache_data(show_spinner="Loading PSM table…")
def load_psm():
    download_with_ua(PSM_URL, PSM_LOCAL)
    usecols = [
        "Spectrum", "Peptide", "Modified.Peptide", "Charge", "Retention",
        "Calibrated.Observed.Mass", "Calibrated.Observed.M.Z",
        "Calculated.Peptide.Mass", "Calculated.M.Z", "Calculated.Peptide.M+H",
        "Hyperscore", "Probability", "Intensity",
        "Protein", "Gene", "Protein.Description", "Precursor", "USI",
        "Tissue", "Organism", "InstrumentModel", "ExperimentID"
    ]
    df = pd.read_csv(PSM_LOCAL, sep="\t", usecols=lambda c: c in usecols)
    df["Experiment_Table"] = "Human Breast TNBC (Phase 2)"
    return df

def load_frags_for_usi(usi: str) -> pd.DataFrame:
    """Read only rows for one USI (chunked) — never load the full 3 GB file."""
    path = ensure_fragments_tsv()
    usecols = ["USI", "Ion_Label", "Ion_mz", "Ion_Intensity", "Peptide", "Precursor"]
    parts = []
    for chunk in pd.read_csv(path, sep="\t", usecols=lambda c: c in usecols, chunksize=200_000):
        hit = chunk[chunk["USI"] == usi]
        if not hit.empty:
            parts.append(hit)
    if not parts:
        return pd.DataFrame(columns=usecols)
    return pd.concat(parts, ignore_index=True)

# ---- startup: PSM only ----
try:
    df = load_psm()
    ensure_fragments_tsv()  # extract once; do not load into memory
except Exception as e:
    st.error("Failed to prepare PAMSI data from GitHub Release.")
    st.exception(e)
    st.stop()

exp_options = sorted(df["Experiment_Table"].dropna().unique().tolist())

# ================== SIDEBAR ==================
st.sidebar.header("🔎 PAMSI Search Controls")
mass_input = st.sidebar.text_area(
    "Enter mass (comma-separated or one per line)",
    value="Add m/z here in comma-separated or one per line",
    height=100,
)
ppm_tolerance = st.sidebar.number_input("Tolerance (ppm)", 0, 1000, 15, 1)
experiments = st.sidebar.multiselect("Experiments", exp_options, default=exp_options)
protein_filter = st.sidebar.text_input("Filter by Protein / Gene ID (e.g. COL1A1)")

if st.sidebar.button("🔍 Search", type="primary", use_container_width=True):
    with st.spinner("Searching…"):
        targets = []
        if mass_input.strip():
            for x in mass_input.replace("\n", ",").split(","):
                x = x.strip()
                try:
                    targets.append(float(x))
                except ValueError:
                    pass
        mask = df["Experiment_Table"].isin(experiments)
        if targets and "Calculated.Peptide.M+H" in df.columns:
            mass_col = df["Calculated.Peptide.M+H"].astype(float).values
            mass_mask = np.zeros(len(df), dtype=bool)
            for t in targets:
                w = t * ppm_tolerance / 1e6
                mass_mask |= (mass_col >= t - w) & (mass_col <= t + w)
            mask &= mass_mask
        if protein_filter.strip():
            pf = protein_filter.strip()
            mask &= (
                df["Protein"].fillna("").astype(str).str.contains(pf, case=False)
                | df["Gene"].fillna("").astype(str).str.contains(pf, case=False)
            )
        st.session_state["filtered"] = df[mask].copy()
        st.session_state["targets"] = targets
        st.session_state["ppm"] = ppm_tolerance

# ================== RESULTS + SPECTRUM ==================
if "filtered" in st.session_state:
    filtered = st.session_state["filtered"]
    targets = st.session_state.get("targets", [])
    ppm = st.session_state.get("ppm", 15)

    if len(filtered) == 0:
        st.warning("No matches.")
    else:
        st.success(f"**{len(filtered):,} peptides** found")
        display_cols = [c for c in [
            "Calculated.Peptide.M+H", "Precursor", "Calibrated.Observed.Mass",
            "Calibrated.Observed.M.Z", "Hyperscore", "Tissue", "Peptide",
            "Protein", "Gene", "Protein.Description", "USI"
        ] if c in filtered.columns]
        st.dataframe(filtered[display_cols], use_container_width=True, hide_index=True, height=350)

        st.markdown("---")
        st.subheader("🔬 Spectrum Viewer")

        usi_map = {}
        for _, row in filtered.drop_duplicates(subset=["USI"]).iterrows():
            label = f"{row.get('Precursor','')}  |  {row.get('Peptide','')}  (score: {row.get('Hyperscore','')})"
            usi_map[label] = row["USI"]

        selected_label = st.selectbox("Select a peptide", list(usi_map.keys()))
        selected_usi = usi_map[selected_label]

        with st.spinner("Loading fragments for selected peptide…"):
            frags = load_frags_for_usi(selected_usi)

        if frags.empty:
            st.warning("No fragment ions for this USI.")
        else:
            meta = filtered[filtered["USI"] == selected_usi].iloc[0]
            st.markdown(
                f"**Peptide:** `{meta.get('Peptide','')}` | "
                f"**Precursor:** `{meta.get('Precursor','')}` | "
                f"**Protein:** `{meta.get('Protein','')}`"
            )
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=frags["Ion_mz"], y=frags["Ion_Intensity"],
                mode="markers+text", text=frags["Ion_Label"],
                textposition="top center",
                hovertemplate="<b>%{text}</b><br>m/z: %{x:.4f}<br>Intensity: %{y}<extra></extra>",
            ))
            for _, row in frags.iterrows():
                fig.add_shape(type="line", x0=row["Ion_mz"], x1=row["Ion_mz"],
                              y0=0, y1=row["Ion_Intensity"],
                              line=dict(color="rgba(31,119,180,0.4)", width=1))
            fig.update_layout(
                title=f"Annotated Spectrum – {meta.get('Precursor','')}",
                xaxis_title="m/z", yaxis_title="Intensity",
                template="plotly_white", height=500, showlegend=False
            )
            st.plotly_chart(fig, use_container_width=True)
else:
    st.info("Enter masses and click **Search**.")

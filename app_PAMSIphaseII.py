import streamlit as st
import pandas as pd
import numpy as np
import zipfile
from pathlib import Path
import urllib.request
import plotly.graph_objects as go

# ================== PAGE CONFIG (must be first Streamlit call) ==================
st.set_page_config(
    page_title="PAMSI — Private Peptide Atlas",
    page_icon="🧬",
    layout="wide",
)

# ================== HEADER ==================
LOGO = Path("assets/PAMSIlogotransparent2.png")
# Wider logo column so logo matches title visual weight
col_logo, col_title = st.columns([2.2, 5])
with col_logo:
    if LOGO.exists():
        st.image(str(LOGO), width=320)  # was 140 — ~title height
    else:
        st.markdown("### 🧬 PAMSI")
with col_title:
    st.markdown(
        '<h1 style="margin-bottom:0.15rem; font-size:2.1rem; line-height:1.2;">'
        "PAMSI — Peptide Atlas for Mass Spectrometry Imaging"
        "</h1>",
        unsafe_allow_html=True,
    )
    st.caption(
        "MALDI Imaging Mass Spectrometry • Extracellular Matrix Proteome • "
        "Human + Mouse datasets • Private Access"
    )

# ================== PASSWORD PROTECTION ==================
if "password_correct" not in st.session_state:
    st.session_state.password_correct = False


def check_password():
    password = st.text_input("Enter Collaborator Password", type="password")
    if password == "PAMSI2026_Collab":  # ← CHANGE THIS TO YOUR OWN SECRET
        st.session_state.password_correct = True
        st.rerun()
    elif password:
        st.error("Incorrect password")


if not st.session_state.password_correct:
    check_password()
    st.stop()

# ================== PATHS / URLS (GitHub Release) ==================
BASE_URL = "https://github.com/CharuKMidha/PAMSI-Database/releases/download/v1.0-data"
DATA_DIR = Path("data")

DATASETS = {
    "Human TNBC": {
        "psm_url": f"{BASE_URL}/HumanBreast_psm_best.tsv",
        "frag_url": f"{BASE_URL}/HumanBreast_fragments.tsv.zip",
        "is_zipped": True,
        "psm_local": DATA_DIR / "HumanBreast_psm_best.tsv",
        "frag_tsv_local": DATA_DIR / "HumanBreast_fragments.tsv",
        "frag_zip_local": DATA_DIR / "HumanBreast_fragments.tsv.zip",
    },
    "Mouse Cochlea": {
        "psm_url": f"{BASE_URL}/Mousecochlea_psm_best.tsv",
        "frag_url": f"{BASE_URL}/Mousecochlea_fragments.tsv",
        "is_zipped": False,
        "psm_local": DATA_DIR / "Mousecochlea_psm_best.tsv",
        "frag_tsv_local": DATA_DIR / "Mousecochlea_fragments.tsv",
        "frag_zip_local": None,
    },
}


def download_with_ua(url: str, dest: Path):
    """Download from GitHub Release with User-Agent (first run only)."""
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


def ensure_fragments_tsv(dataset_name: str):
    """Ensure the fragment TSV for a given dataset is on disk (on demand)."""
    cfg = DATASETS[dataset_name]
    frag_tsv = cfg["frag_tsv_local"]
    if frag_tsv.exists() and frag_tsv.stat().st_size > 0:
        return frag_tsv
    if cfg["is_zipped"]:
        download_with_ua(cfg["frag_url"], cfg["frag_zip_local"])
        with zipfile.ZipFile(cfg["frag_zip_local"], "r") as zf:
            members = [m for m in zf.namelist() if m.lower().endswith((".tsv", ".csv", ".txt"))]
            if not members:
                members = zf.namelist()
            if not members:
                raise ValueError(f"Empty ZIP for {dataset_name}")
            with zf.open(members[0]) as src, open(frag_tsv, "wb") as out:
                while True:
                    chunk = src.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
    else:
        download_with_ua(cfg["frag_url"], frag_tsv)
    return frag_tsv


@st.cache_data(show_spinner="Loading PSM tables…")
def load_psm():
    """Load and combine both Human and Mouse PSM tables."""
    usecols = [
        "Spectrum", "Peptide", "Modified.Peptide.Unimod", "Charge", "Retention",
        "Calibrated.Observed.Mass", "Calibrated.Observed.M.Z",
        "Calculated.Peptide.Mass", "Calculated.M.Z", "Calibrated.Observed.M+H",
        "Hyperscore", "Probability", "Intensity", "Ion.Mobility",
        "Protein", "Gene", "Protein.Description", "Precursor", "USI",
        "Tissue", "Organism", "InstrumentModel", "ExperimentID",
    ]
    dfs = []
    for exp_name, cfg in DATASETS.items():
        download_with_ua(cfg["psm_url"], cfg["psm_local"])
        df_part = pd.read_csv(cfg["psm_local"], sep="\t", usecols=lambda c: c in usecols)
        df_part["Experiment_Table"] = exp_name
        dfs.append(df_part)
    return pd.concat(dfs, ignore_index=True)


def load_frags_for_usi(usi: str, experiment: str) -> pd.DataFrame:
    """Load only fragment rows for one USI from the correct dataset file."""
    path = ensure_fragments_tsv(experiment)
    usecols = ["USI", "Ion_Label", "Ion_mz", "Ion_Intensity", "Peptide", "Precursor"]
    parts = []
    for chunk in pd.read_csv(path, sep="\t", usecols=lambda c: c in usecols, chunksize=200_000):
        hit = chunk[chunk["USI"] == usi]
        if not hit.empty:
            parts.append(hit)
    if not parts:
        return pd.DataFrame(columns=usecols)
    return pd.concat(parts, ignore_index=True)


# ---- Startup: load combined PSM only (fragments load on demand) ----
try:
    df = load_psm()
except Exception as e:
    st.error("Failed to prepare PAMSI data from GitHub Release.")
    st.exception(e)
    st.stop()

exp_options = sorted(df["Experiment_Table"].dropna().unique().tolist())

# ================== SIDEBAR FILTERS ==================
st.sidebar.header("🔎 PAMSI Search Controls")
st.sidebar.subheader("Target Calibrated Observed Mass + H values")
mass_input = st.sidebar.text_area(
    "Enter mass (Add m/z here in comma-separated or one per line)",
    value="",
    height=100,
    help="Paste a list with commas or one mass per line (from Calibrated.Observed.M+H)",
)
ppm_tolerance = st.sidebar.number_input(
    "Tolerance (ppm)", min_value=0, max_value=1000, value=15, step=1
)
experiments = st.sidebar.multiselect(
    "Experiments",
    options=exp_options,
    default=exp_options,
)
protein_filter = st.sidebar.text_input("Filter by Protein / Gene ID (e.g. COL1A1, CO1A1)")

# ================== SEARCH ==================
if st.sidebar.button("🔍 Search", type="primary", use_container_width=True):
    with st.spinner("Searching PAMSI database…"):
        targets = []
        if mass_input.strip():
            raw = [x.strip() for x in mass_input.replace("\n", ",").split(",")]
            for x in raw:
                try:
                    targets.append(float(x))
                except ValueError:
                    pass
        mask = df["Experiment_Table"].isin(experiments)
        if targets and "Calibrated.Observed.M+H" in df.columns:
            mass_col = df["Calibrated.Observed.M+H"].astype(float).values
            mass_mask = np.zeros(len(df), dtype=bool)
            for t in targets:
                window = t * ppm_tolerance / 1e6
                mass_mask |= (mass_col >= t - window) & (mass_col <= t + window)
            mask &= mass_mask
        if protein_filter.strip():
            pf = protein_filter.strip()
            mask &= (
                df["Protein"].fillna("").astype(str).str.contains(pf, case=False)
                | df["Gene"].fillna("").astype(str).str.contains(pf, case=False)
            )
        filtered = df[mask].copy()
        st.session_state["filtered"] = filtered
        st.session_state["targets"] = targets
        st.session_state["ppm"] = ppm_tolerance

# ================== RESULTS + SPECTRUM VIEWER ==================
if "filtered" in st.session_state:
    filtered = st.session_state["filtered"]
    targets = st.session_state.get("targets", [])
    ppm = st.session_state.get("ppm", 15)

    if len(filtered) == 0:
        searched = [f"{t:.4f} ±{ppm} ppm" for t in targets]
        st.warning(
            f"**No matches** for the {len(targets)} mass(es).\n\n"
            f"Searched: {', '.join(searched) if searched else '(none)'}"
        )
        st.info("Tip: Increase ppm tolerance, try different masses, or clear the protein filter.")
    else:
        st.success(
            f"**{len(filtered):,} peptides** found "
            f"({len(targets)} mass(es) at ±{ppm} ppm)"
        )

        st.subheader("📊 Peptide Hits by Experiment / Tissue")
        tissue_counts = filtered["Experiment_Table"].value_counts().reset_index()
        tissue_counts.columns = ["Tissue", "Peptide Count"]

        fig_bar = go.Figure(
            data=[
                go.Bar(
                    x=tissue_counts["Tissue"],
                    y=tissue_counts["Peptide Count"],
                    marker_color="#1f77b4",
                    text=tissue_counts["Peptide Count"],
                    textposition="outside",
                    textfont=dict(size=14, color="black"),
                )
            ]
        )
        fig_bar.update_layout(
            xaxis_title="Experiment / Tissue",
            yaxis_title="Peptide Count",
            template="plotly_white",
            height=350,
            showlegend=False,
            font=dict(size=14, color="black"),
            xaxis=dict(
                title_font=dict(size=16, color="black"),
                tickfont=dict(size=14, color="black"),
            ),
            yaxis=dict(
                title_font=dict(size=16, color="black"),
                tickfont=dict(size=14, color="black"),
            ),
            title_font=dict(size=16, color="black"),
            margin=dict(t=40, b=60),
        )
        st.plotly_chart(fig_bar, use_container_width=True)

        display_cols = [
            "Calibrated.Observed.M+H", "Calibrated.Observed.Mass", "Modified.Peptide.Unimod",
            "Charge", "Hyperscore", "Ion.Mobility", "Experiment_Table", "Tissue",
            "Peptide", "Protein", "Gene", "Protein.Description", "USI",
        ]
        display_cols = [c for c in display_cols if c in filtered.columns]

        st.subheader("Precursor table")
        st.dataframe(
            filtered[display_cols],
            use_container_width=True,
            hide_index=True,
            height=350,
        )

        st.download_button(
            label="📥 Download results as CSV",
            data=filtered[display_cols].to_csv(index=False).encode(),
            file_name=f"pamsi_results_±{ppm}ppm.csv",
            mime="text/csv",
        )

        st.markdown(
            "Files were searched with **FragPipe 24.0** (MSFragger) and pre-processed with an "
            "in-house R script available on "
            "[GitHub (PAMSI-Database)](https://github.com/CharuKMidha/PAMSI-Database). "
            "Best-scoring PSMs (by Hyperscore) were retained per spectrum. Modified peptide sequences use UNIMOD annotation. "
            "Modified residues may appear as mass tags (e.g. `N[115]`) or as UNIMOD annotations "
            "(e.g. `N[UNIMOD:7]`)."
        )
        with st.expander("Common mass-tag modifications", expanded=False):
            st.markdown(
                """
| Notation | Amino acid | Unmodified residue ≈ | Modified residue ≈ | Δ mass | Meaning | UNIMOD |
|----------|------------|----------------------|--------------------|--------|---------|--------|
| **N[115]** | Asn (N) | 114.04 | 115.03 | +0.9840 | Deamidated Asn | **UNIMOD:7** |
| **Q[129]** | Gln (Q) | 128.06 | 129.04 | +0.9840 | Deamidated Gln | **UNIMOD:7** |
| **M[147]** | Met (M) | 131.04 | 147.04 | +15.9949 | Oxidized Met | **UNIMOD:35** |
| **P[113]** | Pro (P) | 97.05 | 113.05 | +15.9949 | Hydroxylated / oxidized Pro | **UNIMOD:35** |
                """
            )
            st.caption(
                "Example: `P[113]` or `P(15.9949)` = proline with oxidation/hydroxylation "
                "(`P[UNIMOD:35]`) — frequent in collagen extracellular-matrix peptides."
            )

        # ========== SPECTRUM VIEWER ==========
        st.markdown("---")
        st.subheader("🔬 Spectrum Viewer (annotated fragments)")

        usi_map = {}
        usi_to_experiment = {}
        for _, row in filtered.drop_duplicates(subset=["USI"]).iterrows():
            usi = row["USI"]
            precursor = str(row.get("Precursor", ""))
            peptide = str(row.get("Peptide", ""))
            hyperscore = row.get("Hyperscore", "")
            exp = row.get("Experiment_Table", "Human TNBC")
            label = f"{exp} | {precursor}  |  {peptide}  (score: {hyperscore})"
            usi_map[label] = usi
            usi_to_experiment[usi] = exp

        if not usi_map:
            st.info("No peptides available in the current result set.")
        else:
            selected_label = st.selectbox(
                "Select a peptide to view its spectrum",
                options=list(usi_map.keys()),
                help="Showing: Experiment | Precursor | Peptide sequence (Hyperscore)",
            )
            selected_usi = usi_map[selected_label]
            selected_exp = usi_to_experiment[selected_usi]

            with st.spinner(f"Loading fragments for selected peptide from {selected_exp}…"):
                frags = load_frags_for_usi(selected_usi, selected_exp)

            if frags.empty:
                st.warning("No fragment ions found for the selected peptide.")
            else:
                meta = filtered[filtered["USI"] == selected_usi].iloc[0]
                st.markdown(
                    f"**Experiment:** `{meta.get('Experiment_Table', '')}` &nbsp;|&nbsp; "
                    f"**Peptide:** `{meta.get('Peptide', '')}` &nbsp;|&nbsp; "
                    f"**Precursor:** `{meta.get('Precursor', '')}` &nbsp;|&nbsp; "
                    f"**Protein:** `{meta.get('Protein', '')}` &nbsp;|&nbsp; "
                    f"**Hyperscore:** {meta.get('Hyperscore', 'N/A')}"
                )

                fig = go.Figure()
                fig.add_trace(
                    go.Scatter(
                        x=frags["Ion_mz"],
                        y=frags["Ion_Intensity"],
                        mode="markers",
                        marker=dict(size=8, color="#1f77b4"),
                        hovertemplate="<b>%{text}</b><br>m/z: %{x:.4f}<br>Intensity: %{y}<extra></extra>",
                        text=frags["Ion_Label"],
                        name="Annotated ions",
                    )
                )
                for _, row in frags.iterrows():
                    fig.add_shape(
                        type="line",
                        x0=row["Ion_mz"],
                        x1=row["Ion_mz"],
                        y0=0,
                        y1=row["Ion_Intensity"],
                        line=dict(color="rgba(31,119,180,0.4)", width=1),
                    )
                for _, row in frags.iterrows():
                    fig.add_annotation(
                        x=row["Ion_mz"],
                        y=row["Ion_Intensity"],
                        text=str(row["Ion_Label"]),
                        showarrow=False,
                        textangle=-90,
                        font=dict(size=12, color="black"),
                        xanchor="left",
                        yanchor="bottom",
                        yshift=4,
                    )
                fig.update_layout(
                    title=f"Annotated Spectrum – {meta.get('Precursor', '')} ({selected_exp})",
                    xaxis_title="m/z",
                    yaxis_title="Intensity",
                    template="plotly_white",
                    height=500,
                    showlegend=False,
                    hovermode="closest",
                    font=dict(size=14, color="black"),
                    xaxis=dict(
                        title_font=dict(size=16, color="black"),
                        tickfont=dict(size=14, color="black"),
                    ),
                    yaxis=dict(
                        title_font=dict(size=16, color="black"),
                        tickfont=dict(size=14, color="black"),
                        range=[0, float(frags["Ion_Intensity"].max()) * 1.35],
                    ),
                    title_font=dict(size=16, color="black"),
                )
                st.plotly_chart(fig, use_container_width=True)

                # ---- Zoom / pan help (BEFORE fragment ion table) ----
                st.markdown("#### Spectrum viewer — zoom & pan")
                st.markdown(
                    """
- **Scroll** on the plot to zoom in or out
- **Click and drag** to pan across m/z
- **Double-click** the plot to reset zoom
- **Hover** a peak to see ion label, m/z, and intensity

Ion labels (b/y) are angled to reduce overlap.
                    """
                )

                with st.expander("Fragment ion table"):
                    st.dataframe(
                        frags[["Ion_Label", "Ion_mz", "Ion_Intensity"]].sort_values("Ion_mz"),
                        use_container_width=True,
                        hide_index=True,
                    )
else:
    st.info(
        "Adjust parameters on the left and click **🔍 Search** to query the PAMSI database "
        "and open the spectrum viewer."
    )

st.caption(
    "✅ Private mode • Human + Mouse datasets • Data from GitHub Release only • "
    "USI-linked spectrum viewer • Author Dr. Charu Kapil Midha & Prof. Peggi Angel • 2026"
)

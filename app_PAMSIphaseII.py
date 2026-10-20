import streamlit as st
import pandas as pd
import numpy as np
import zipfile
import io
import urllib.request
import plotly.graph_objects as go

st.set_page_config(page_title="PAMSI — Private Peptide Atlas", layout="wide")
st.title("🧬 PAMSI — Peptide Atlas for Mass Spectrometry Imaging")
st.caption("MALDI Imaging Mass Spectrometry • Extracellular Matrix Proteome • Human + Mouse datasets • Private Access")

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

# ================== LOAD DATA (GitHub Release only — no local data/) ==================
PSM_URL = "https://github.com/CharuKMidha/PAMSI-Database/releases/download/v1.0-data/HumanBreast_Phase2_psm_best.tsv"
FRAG_URL = "https://github.com/CharuKMidha/PAMSI-Database/releases/download/v1.0-data/HumanBreast_Phase2_fragments.tsv.zip"


@st.cache_data(show_spinner="Loading PAMSI datasets from GitHub Release…")
def load_data():
    # ----- PSM: read TSV straight from the Release URL -----
    usecols_psm = [
        "Spectrum", "Peptide", "Modified.Peptide", "Charge", "Retention",
        "Calibrated.Observed.Mass", "Calibrated.Observed.M.Z",
        "Calculated.Peptide.Mass", "Calculated.M.Z", "Calculated.Peptide.M+H",
        "Hyperscore", "Probability", "Intensity",
        "Protein", "Gene", "Protein.Description", "Precursor", "USI",
        "Tissue", "Organism", "InstrumentModel", "ExperimentID"
    ]
    df_psm = pd.read_csv(
        PSM_URL,
        sep="\t",
        usecols=lambda c: c in usecols_psm
    )
    df_psm["Experiment_Table"] = "Human Breast TNBC (Phase 2)"

    # ----- Fragments: download zip into memory, then read the TSV inside -----
    frag_usecols = ["USI", "Ion_Label", "Ion_mz", "Ion_Intensity", "Peptide", "Precursor"]

    with urllib.request.urlopen(FRAG_URL) as response:
        zip_bytes = response.read()

    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        members = [m for m in zf.namelist() if m.lower().endswith((".tsv", ".csv", ".txt"))]
        if not members:
            members = zf.namelist()
        if not members:
            raise ValueError("ZIP archive from GitHub Release is empty")
        with zf.open(members[0]) as f:
            df_frag = pd.read_csv(f, sep="\t", usecols=lambda c: c in frag_usecols)

    return df_psm, df_frag


# ---- Call the loader ----
try:
    df, df_fragments = load_data()
except Exception as e:
    st.error("Failed to load PAMSI data from GitHub Release.")
    st.exception(e)
    st.stop()

# Now it is safe to build the sidebar
exp_options = sorted(df["Experiment_Table"].dropna().unique().tolist())

# ================== SIDEBAR FILTERS ==================
st.sidebar.header("🔎 PAMSI Search Controls")

st.sidebar.subheader("Target Calculated Peptide Mass + H values")
mass_input = st.sidebar.text_area(
    "Enter mass (comma-separated or one per line)",
    value="Add m/z here in comma-separated or one per line",
    height=100,
    help="Paste a list with commas or one mass per line (from Calculated.Peptide.M+H)"
)

ppm_tolerance = st.sidebar.number_input(
    "Tolerance (ppm)", min_value=0, max_value=1000, value=15, step=1
)

experiments = st.sidebar.multiselect(
    "Experiments",
    options=exp_options,
    default=exp_options
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

        if targets and "Calculated.Peptide.M+H" in df.columns:
            mass_col = df["Calculated.Peptide.M+H"].astype(float).values
            mass_mask = np.zeros(len(df), dtype=bool)
            for t in targets:
                window = t * ppm_tolerance / 1e6
                mass_mask |= (mass_col >= t - window) & (mass_col <= t + window)
            mask &= mass_mask

        if protein_filter.strip():
            pf = protein_filter.strip()
            mask &= (
                df["Protein"].fillna("").astype(str).str.contains(pf, case=False) |
                df["Gene"].fillna("").astype(str).str.contains(pf, case=False)
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
        st.bar_chart(tissue_counts.set_index("Tissue"))

        display_cols = [
            "Calculated.Peptide.M+H", "Precursor", "Calibrated.Observed.Mass",
            "Calibrated.Observed.M.Z", "Hyperscore", "Tissue", "Peptide",
            "Protein", "Gene", "Protein.Description", "USI"
        ]
        display_cols = [c for c in display_cols if c in filtered.columns]

        st.dataframe(
            filtered[display_cols],
            use_container_width=True,
            hide_index=True,
            height=350
        )

        st.download_button(
            label="📥 Download results as CSV",
            data=filtered[display_cols].to_csv(index=False).encode(),
            file_name=f"pamsi_results_±{ppm}ppm.csv",
            mime="text/csv"
        )

        # ========== SPECTRUM VIEWER ==========
        st.markdown("---")
        st.subheader("🔬 Spectrum Viewer (annotated fragments)")

        usi_map = {}
        for _, row in filtered.drop_duplicates(subset=["USI"]).iterrows():
            usi = row["USI"]
            precursor = str(row.get("Precursor", ""))
            peptide = str(row.get("Peptide", ""))
            hyperscore = row.get("Hyperscore", "")
            label = f"{precursor}  |  {peptide}  (score: {hyperscore})"
            usi_map[label] = usi

        if not usi_map:
            st.info("No peptides available in the current result set.")
        else:
            selected_label = st.selectbox(
                "Select a peptide to view its spectrum",
                options=list(usi_map.keys()),
                help="Showing: Precursor | Peptide sequence (Hyperscore)"
            )
            selected_usi = usi_map[selected_label]

            frags = df_fragments[df_fragments["USI"] == selected_usi].copy()

            if frags.empty:
                st.warning("No fragment ions found for the selected peptide.")
            else:
                meta = filtered[filtered["USI"] == selected_usi].iloc[0]
                st.markdown(
                    f"**Peptide:** `{meta.get('Peptide', '')}` &nbsp;|&nbsp; "
                    f"**Precursor:** `{meta.get('Precursor', '')}` &nbsp;|&nbsp; "
                    f"**Protein:** `{meta.get('Protein', '')}` &nbsp;|&nbsp; "
                    f"**Hyperscore:** {meta.get('Hyperscore', 'N/A')}"
                )

                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=frags["Ion_mz"],
                    y=frags["Ion_Intensity"],
                    mode="markers+text",
                    marker=dict(size=8, color="#1f77b4"),
                    text=frags["Ion_Label"],
                    textposition="top center",
                    textfont=dict(size=9),
                    hovertemplate="<b>%{text}</b><br>m/z: %{x:.4f}<br>Intensity: %{y}<extra></extra>",
                    name="Annotated ions"
                ))

                for _, row in frags.iterrows():
                    fig.add_shape(
                        type="line",
                        x0=row["Ion_mz"], x1=row["Ion_mz"],
                        y0=0, y1=row["Ion_Intensity"],
                        line=dict(color="rgba(31,119,180,0.4)", width=1)
                    )

                fig.update_layout(
                    title=f"Annotated Spectrum – {meta.get('Precursor', '')}",
                    xaxis_title="m/z",
                    yaxis_title="Intensity",
                    template="plotly_white",
                    height=500,
                    showlegend=False,
                    hovermode="closest"
                )
                st.plotly_chart(fig, use_container_width=True)

                with st.expander("Fragment ion table"):
                    st.dataframe(
                        frags[["Ion_Label", "Ion_mz", "Ion_Intensity"]].sort_values("Ion_mz"),
                        use_container_width=True,
                        hide_index=True
                    )

else:
    st.info("Adjust parameters on the left and click **🔍 Search** to query the PAMSI database and open the spectrum viewer.")

st.caption("✅ Private mode • Phase 2 Human Breast (TNBC) • Data from GitHub Release only • USI-linked spectrum viewer • Author Charu Kapil Midha & Prof. Peggi Angel • 2026")

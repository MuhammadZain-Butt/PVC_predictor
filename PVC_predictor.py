import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO, BytesIO
import time
import subprocess
import tempfile
import os
import io
import re
import requests
import textwrap
import networkx as nx
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import glob
import zipfile

# ============================================================
# PAGE CONFIG
# ============================================================
st.set_page_config(
    page_title="ProteomePipeline · Integrated Analysis",
    layout="wide",
    page_icon="🧬"
)

# ============================================================
# CUSTOM CSS  —  refined typography + dark bioinfo theme
# ============================================================
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=DM+Sans:ital,opsz,wght@0,9..40,300;0,9..40,400;0,9..40,500;0,9..40,600;0,9..40,700;1,9..40,300&family=DM+Serif+Display&display=swap');

:root {
    --bg:        #080d18;
    --surface:   #0f1623;
    --surface2:  #162030;
    --surface3:  #1c2840;
    --accent:    #38bdf8;
    --accent2:   #818cf8;
    --accent3:   #34d399;
    --accent4:   #fbbf24;
    --danger:    #f87171;
    --text:      #cbd5e1;
    --text-hi:   #f1f5f9;
    --text-lo:   #475569;
    --border:    rgba(56,189,248,0.12);
    --border2:   rgba(56,189,248,0.06);
}

/* ── Base ─────────────────────────────────────────── */
html, body, [class*="css"] {
    font-family: 'DM Sans', sans-serif;
    background-color: var(--bg);
    color: var(--text);
    font-size: 14px;
    line-height: 1.6;
}
.stApp {
    background: radial-gradient(ellipse at 20% 10%, rgba(56,189,248,0.04) 0%, transparent 60%),
                radial-gradient(ellipse at 80% 90%, rgba(129,140,248,0.04) 0%, transparent 60%),
                var(--bg);
}

/* ── Header ───────────────────────────────────────── */
.pipe-header {
    background: linear-gradient(110deg, rgba(56,189,248,0.06) 0%, rgba(129,140,248,0.06) 100%);
    border: 1px solid var(--border);
    border-radius: 14px;
    padding: 1.8rem 2.4rem 1.6rem;
    margin-bottom: 1.8rem;
    position: relative;
    overflow: hidden;
}
.pipe-header::before {
    content: '';
    position: absolute;
    top: 0; left: 0; right: 0; height: 2px;
    background: linear-gradient(90deg, var(--accent), var(--accent2), var(--accent3), var(--accent4));
}
.pipe-title {
    font-family: 'DM Serif Display', serif;
    font-size: 2.2rem;
    font-weight: 400;
    letter-spacing: -0.5px;
    background: linear-gradient(100deg, var(--accent) 10%, var(--accent2) 90%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    margin: 0 0 0.3rem 0;
}
.pipe-sub {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.72rem;
    font-weight: 500;
    color: var(--text-lo);
    letter-spacing: 2px;
    text-transform: uppercase;
}

/* ── Step card ────────────────────────────────────── */
.step-card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 1.4rem 1.6rem 1.2rem;
    margin-bottom: 1.2rem;
}
.step-badge {
    display: inline-block;
    background: linear-gradient(135deg, var(--accent), var(--accent2));
    color: #04090f;
    font-family: 'IBM Plex Mono', monospace;
    font-weight: 600;
    font-size: 0.65rem;
    padding: 2px 10px;
    border-radius: 20px;
    letter-spacing: 1.5px;
    text-transform: uppercase;
    margin-bottom: 0.6rem;
}
.step-title {
    font-family: 'DM Sans', sans-serif;
    font-size: 1.05rem;
    font-weight: 700;
    color: var(--text-hi);
    margin: 0.2rem 0 0.4rem 0;
    letter-spacing: -0.2px;
}
.step-desc {
    font-size: 0.82rem;
    color: var(--text-lo);
    line-height: 1.55;
    margin: 0;
    font-weight: 400;
}

/* ── Flow bar ─────────────────────────────────────── */
.flow-wrap {
    display: flex;
    align-items: center;
    gap: 4px;
    margin: 0.8rem 0 1.6rem;
    flex-wrap: wrap;
}
.flow-node {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.68rem;
    font-weight: 500;
    letter-spacing: 0.5px;
    padding: 4px 12px;
    border-radius: 6px;
    border: 1px solid rgba(71,85,105,0.4);
    background: var(--surface2);
    color: var(--text-lo);
    white-space: nowrap;
}
.flow-node.done {
    border-color: rgba(52,211,153,0.35);
    background: rgba(52,211,153,0.05);
    color: var(--accent3);
}
.flow-arrow { color: var(--text-lo); font-size: 0.7rem; padding: 0 1px; }

/* ── Metric boxes ─────────────────────────────────── */
.metric-row {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(130px, 1fr));
    gap: 0.8rem;
    margin: 1rem 0;
}
.mbox {
    background: var(--surface2);
    border: 1px solid var(--border2);
    border-radius: 10px;
    padding: 0.9rem 1rem;
    text-align: center;
}
.mbox-val {
    font-family: 'DM Sans', sans-serif;
    font-size: 1.7rem;
    font-weight: 700;
    color: var(--accent);
    line-height: 1.1;
}
.mbox-lbl {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.62rem;
    color: var(--text-lo);
    text-transform: uppercase;
    letter-spacing: 1px;
    margin-top: 4px;
}

/* ── Section labels ───────────────────────────────── */
.sec-label {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.7rem;
    font-weight: 600;
    color: var(--accent);
    letter-spacing: 1.5px;
    text-transform: uppercase;
    margin: 1.2rem 0 0.6rem;
    padding-left: 2px;
    border-left: 2px solid var(--accent);
    padding-left: 8px;
}

/* ── Sidebar ──────────────────────────────────────── */
section[data-testid="stSidebar"] {
    background: var(--surface) !important;
    border-right: 1px solid var(--border) !important;
}
section[data-testid="stSidebar"] .stMarkdown p {
    font-size: 0.78rem;
    color: var(--text-lo);
}

/* ── Buttons ──────────────────────────────────────── */
.stButton > button {
    background: linear-gradient(135deg, rgba(56,189,248,0.9), rgba(129,140,248,0.9)) !important;
    color: #04090f !important;
    font-family: 'DM Sans', sans-serif !important;
    font-weight: 600 !important;
    font-size: 0.84rem !important;
    border: none !important;
    border-radius: 8px !important;
    padding: 0.55rem 1.4rem !important;
    letter-spacing: 0.3px !important;
}
.stButton > button:hover { opacity: 0.88 !important; }

/* ── Download buttons ─────────────────────────────── */
.stDownloadButton > button {
    background: transparent !important;
    color: var(--accent) !important;
    border: 1px solid rgba(56,189,248,0.4) !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.75rem !important;
    border-radius: 6px !important;
    padding: 0.4rem 1rem !important;
}
.stDownloadButton > button:hover {
    background: rgba(56,189,248,0.06) !important;
}

/* ── Tabs ─────────────────────────────────────────── */
.stTabs [data-baseweb="tab-list"] {
    background: var(--surface2) !important;
    border-radius: 10px;
    border: 1px solid var(--border2);
    padding: 3px;
    gap: 2px !important;
}
.stTabs [data-baseweb="tab"] {
    color: var(--text-lo) !important;
    font-family: 'DM Sans', sans-serif !important;
    font-weight: 500 !important;
    font-size: 0.78rem !important;
    border-radius: 7px !important;
    padding: 0.35rem 0.9rem !important;
}
.stTabs [aria-selected="true"] {
    background: linear-gradient(135deg, rgba(56,189,248,0.12), rgba(129,140,248,0.12)) !important;
    color: var(--accent) !important;
}

/* ── Expander ─────────────────────────────────────── */
.streamlit-expanderHeader {
    font-family: 'DM Sans', sans-serif !important;
    font-weight: 600 !important;
    font-size: 0.82rem !important;
    color: var(--text) !important;
    background: var(--surface2) !important;
    border: 1px solid var(--border2) !important;
    border-radius: 8px !important;
}

/* ── Inputs ───────────────────────────────────────── */
.stNumberInput > div > div > input,
.stTextInput > div > div > input,
.stTextArea > div > div > textarea {
    background: var(--surface2) !important;
    border-color: var(--border) !important;
    color: var(--text-hi) !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.82rem !important;
    border-radius: 7px !important;
}
.stSelectbox > div > div {
    background: var(--surface2) !important;
    border-color: var(--border) !important;
    color: var(--text-hi) !important;
    font-size: 0.82rem !important;
}

/* ── Alerts ───────────────────────────────────────── */
.stSuccess > div { background: rgba(52,211,153,0.08) !important; border-color: rgba(52,211,153,0.3) !important; }
.stInfo > div { background: rgba(56,189,248,0.06) !important; border-color: rgba(56,189,248,0.25) !important; }
.stWarning > div { background: rgba(251,191,36,0.08) !important; border-color: rgba(251,191,36,0.3) !important; }
.stError > div { background: rgba(248,113,113,0.08) !important; border-color: rgba(248,113,113,0.3) !important; }

/* Alert text */
.stSuccess p, .stInfo p, .stWarning p, .stError p {
    font-size: 0.82rem !important;
    font-family: 'DM Sans', sans-serif !important;
}

/* ── Dataframe ────────────────────────────────────── */
.stDataFrame { border-radius: 8px; overflow: hidden; }
.stDataFrame th {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.72rem !important;
    text-transform: uppercase !important;
    letter-spacing: 0.5px !important;
    color: var(--text-lo) !important;
}
.stDataFrame td {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.78rem !important;
    color: var(--text) !important;
}

/* ── File uploader ────────────────────────────────── */
[data-testid="stFileUploader"] {
    background: var(--surface2) !important;
    border: 1px dashed rgba(56,189,248,0.2) !important;
    border-radius: 10px !important;
}
[data-testid="stFileUploader"] p {
    font-size: 0.8rem !important;
    color: var(--text-lo) !important;
}

/* ── Radio ────────────────────────────────────────── */
.stRadio > div > div > label {
    font-size: 0.82rem !important;
    color: var(--text) !important;
}

/* ── Checkboxes ───────────────────────────────────── */
.stCheckbox > label {
    font-size: 0.82rem !important;
    color: var(--text) !important;
}

/* ── Slider ───────────────────────────────────────── */
.stSlider .st-bd { background: var(--accent) !important; }
.stSlider p { font-size: 0.8rem !important; color: var(--text) !important; }

/* ── Headings override ────────────────────────────── */
h1 { font-family:'DM Serif Display',serif !important; font-size:1.7rem !important; color:var(--text-hi) !important; font-weight:400 !important; }
h2 { font-family:'DM Sans',sans-serif !important; font-size:1.1rem !important; color:var(--text-hi) !important; font-weight:700 !important; }
h3 { font-family:'DM Sans',sans-serif !important; font-size:0.95rem !important; color:var(--text) !important; font-weight:600 !important; }
p, li { font-size:0.85rem !important; color:var(--text) !important; }
strong { color:var(--text-hi) !important; }
code { font-family:'IBM Plex Mono',monospace !important; font-size:0.78rem !important;
       background:var(--surface3) !important; padding:1px 5px !important; border-radius:4px !important; }

/* ── Progress bar ─────────────────────────────────── */
.stProgress > div > div > div { background: linear-gradient(90deg, var(--accent), var(--accent2)) !important; }

/* ── Caption ──────────────────────────────────────── */
.stCaption, .stCaption p { font-size:0.72rem !important; color:var(--text-lo) !important; font-style:italic !important; }

hr { border-color: var(--border2) !important; margin: 1.2rem 0 !important; }
</style>
""", unsafe_allow_html=True)

# ============================================================
# AMINO ACID / ALLERGENICITY CONSTANTS
# ============================================================
HYDROPHOBICITY = {
    'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,
    'I':4.5,'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,
    'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3
}
POLAR_RESIDUES = set("DERKHQNSTY")
AA_COLORS = {
    'A':'#4ade80','C':'#facc15','D':'#f87171','E':'#f87171','F':'#818cf8',
    'G':'#94a3b8','H':'#22d3ee','I':'#4ade80','K':'#60a5fa','L':'#4ade80',
    'M':'#a78bfa','N':'#fb923c','P':'#fbbf24','Q':'#fb923c','R':'#60a5fa',
    'S':'#34d399','T':'#34d399','V':'#4ade80','W':'#818cf8','Y':'#c084fc'
}

# ============================================================
# PLOTLY THEME
# ============================================================
PLOT_LAYOUT = dict(
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(15,22,35,0.9)',
    font=dict(family='DM Sans, sans-serif', color='#94a3b8', size=11),
    title_font=dict(family='DM Sans, sans-serif', color='#cbd5e1', size=13, ),
    xaxis=dict(gridcolor='rgba(56,189,248,0.07)', zerolinecolor='rgba(56,189,248,0.15)',
               tickfont=dict(family='IBM Plex Mono', size=10, color='#475569')),
    yaxis=dict(gridcolor='rgba(56,189,248,0.07)', zerolinecolor='rgba(56,189,248,0.15)',
               tickfont=dict(family='IBM Plex Mono', size=10, color='#475569')),
    margin=dict(t=45, b=35, l=50, r=20),
    legend=dict(font=dict(family='DM Sans', size=11, color='#94a3b8'),
                bgcolor='rgba(0,0,0,0)', bordercolor='rgba(0,0,0,0)'),
)

def sf(fig):
    fig.update_layout(**PLOT_LAYOUT)
    return fig

# ============================================================
# SHARED HELPERS
# ============================================================
def fasta_stats(records):
    lengths = [len(r.seq) for r in records]
    return dict(count=len(records), min=min(lengths), max=max(lengths),
                avg=round(np.mean(lengths),1), total=sum(lengths))

def length_hist(df, title, color="#38bdf8"):
    fig = px.histogram(df, x="length", nbins=30, title=title,
                       color_discrete_sequence=[color])
    fig.update_traces(marker_line_color='rgba(0,0,0,0.4)', marker_line_width=1)
    return sf(fig)

def aa_comp_bar(seq, title):
    from collections import Counter
    c = Counter(seq.upper())
    aas = sorted(c.keys())
    fig = go.Figure(go.Bar(
        x=aas, y=[c[a] for a in aas],
        marker_color=[AA_COLORS.get(a,'#94a3b8') for a in aas],
        hovertemplate='%{x}: %{y}<extra></extra>'
    ))
    fig.update_layout(title=title, xaxis_title="Amino Acid", yaxis_title="Count", **PLOT_LAYOUT)
    return fig

def input_viz(records, tag):
    s = fasta_stats(records)
    st.markdown(f"""
    <div class="metric-row">
        <div class="mbox"><div class="mbox-val">{s['count']}</div><div class="mbox-lbl">Sequences</div></div>
        <div class="mbox"><div class="mbox-val" style="color:#818cf8">{s['avg']}</div><div class="mbox-lbl">Avg Length</div></div>
        <div class="mbox"><div class="mbox-val" style="color:#34d399">{s['min']}</div><div class="mbox-lbl">Min aa</div></div>
        <div class="mbox"><div class="mbox-val" style="color:#fbbf24">{s['max']}</div><div class="mbox-lbl">Max aa</div></div>
        <div class="mbox"><div class="mbox-val" style="color:#38bdf8; font-size:1.25rem">{s['total']:,}</div><div class="mbox-lbl">Total Residues</div></div>
    </div>
    """, unsafe_allow_html=True)
    df_l = pd.DataFrame({"id":[r.id for r in records],"length":[len(r.seq) for r in records]})
    c1, c2 = st.columns(2)
    with c1:
        st.plotly_chart(length_hist(df_l, "Length Distribution"), use_container_width=True, key=f"lh_{tag}_{id(records)}")
    with c2:
        all_s = "".join(str(r.seq) for r in records)
        st.plotly_chart(aa_comp_bar(all_s, "AA Composition"), use_container_width=True, key=f"ac_{tag}_{id(records)}")

def step_card(num, title, desc):
    st.markdown(f"""
    <div class="step-card">
        <span class="step-badge">Step {num}</span>
        <div class="step-title">{title}</div>
        <p class="step-desc">{desc}</p>
    </div>""", unsafe_allow_html=True)

def sec(label):
    st.markdown(f'<div class="sec-label">{label}</div>', unsafe_allow_html=True)

# ============================================================
# ALLERGENICITY
# ============================================================
def allergen_score(seq):
    h = sum(HYDROPHOBICITY.get(a,0) for a in seq)/len(seq)
    p = sum(1 for a in seq if a in POLAR_RESIDUES)/len(seq)
    return round((abs(h)*0.6)+((1-p)*0.4), 3)

def allergen_label(score, thr): return "Allergen" if score>=thr else "Non-Allergen"

# ============================================================
# PROTPARAM
# ============================================================
def get_protparam(seq, max_retries=3):
    url = 'https://web.expasy.org/cgi-bin/protparam/protparam'
    for _ in range(max_retries):
        try:
            from bs4 import BeautifulSoup
            res = requests.post(url, data={'prot_id':'','mandatory':'',None:'','sequence':seq}, timeout=20)
            if res.status_code != 200: raise Exception()
            text = BeautifulSoup(res.text,'html.parser').get_text()
            def pv(pat):
                m = re.search(pat, text)
                return float(re.findall(r'-?\d+\.?\d*', m.group(0))[0]) if m else None
            return {
                'Molecular weight': pv(r'Molecular weight: -?[0-9]+\.[0-9]+'),
                'Instability index': pv(r'computed to be -?[0-9]+\.[0-9]+'),
                'Instability': 'Unstable' if 'unstable' in text.lower() else 'Stable',
                'Aliphatic index': pv(r'Aliphatic index: -?[0-9]+\.[0-9]+'),
                'GRAVY score': pv(r'Grand average of hydropathicity \(GRAVY\):\s*-?\d+\.?\d*')
            }
        except: time.sleep(0.5)
    raise RuntimeError("ProtParam failed")

# ============================================================
# ANTIGENICITY HELPERS
# ============================================================
def extract_uniprot_id(header):
    parts = header.split("|")
    return parts[1] if len(parts)>=3 else None

def map_antigenic(csv_df, fasta_records, id_col="Header", label_col="Antigenicity_Category"):
    fasta_dict = {extract_uniprot_id(r.description): r for r in fasta_records if extract_uniprot_id(r.description)}
    lines, matched = [], set()
    for _, row in csv_df.iterrows():
        if row.get(label_col,"") != "High": continue
        uid = extract_uniprot_id(str(row[id_col]))
        if uid and uid in fasta_dict:
            r = fasta_dict[uid]
            lines.append(f">{r.description}\n{r.seq}")
            matched.add(uid)
    return "\n".join(lines), matched

# ============================================================
# SESSION STATE INIT
# ============================================================
for k in [f"step{i}_done" for i in range(1,9)]:
    if k not in st.session_state: st.session_state[k] = False

# ============================================================
# HEADER
# ============================================================
st.markdown("""
<div class="pipe-header">
    <div class="pipe-title">🧬 PVC_predictor</div>
    <div class="pipe-sub">Integrated Proteome Analysis Suite For Potential Vaccine Candidate · 8-Stage Subtractive Genomics Workflow</div>
</div>""", unsafe_allow_html=True)

# Flow indicator
STEP_NAMES = ["01 Paralog","02 Non-Homology","03 Essential","04 Virulence",
              "05 Allergenicity","06 ProtParam","07 Antigenicity","08 Localization"]
DONE_KEYS  = [f"step{i}_done" for i in range(1,9)]

html = '<div class="flow-wrap">'
for i,(name,dk) in enumerate(zip(STEP_NAMES, DONE_KEYS)):
    cls = "flow-node done" if st.session_state[dk] else "flow-node"
    html += f'<span class="{cls}">{name}</span>'
    if i < len(STEP_NAMES)-1: html += '<span class="flow-arrow">›</span>'
html += '</div>'
st.markdown(html, unsafe_allow_html=True)

# ============================================================
# SIDEBAR
# ============================================================
with st.sidebar:
    st.markdown("### ⚙ Pipeline Settings")
    st.markdown("**BLAST Parameters**")
    identity_cutoff = st.slider("Min % Identity", 30, 100, 40)
    coverage_cutoff = st.slider("Min % Coverage", 30, 100, 70)
    bitscore_cutoff = st.number_input("Min Bitscore", value=50)
    evalue_cutoff   = st.selectbox("E-value", [1e-5,1e-10,1e-20,1e-50])
    threads         = st.number_input("Threads", min_value=1, value=4)
    st.markdown("---")
    st.markdown("**Allergenicity**")
    allergen_thr = st.slider("Threshold", 0.4, 0.9, 0.6, 0.01)
    st.markdown("---")
    st.markdown("**ProtParam**")
    chunk_opt = st.checkbox("Split into 30 aa chunks", value=False)
    st.markdown("---")
    st.markdown("""<div style="font-size:0.72rem;color:#475569;line-height:1.7">
    <strong style="color:#38bdf8;font-family:'IBM Plex Mono'">Workflow</strong><br>
    1. Remove paralogous proteins<br>
    2. Filter human-homologous hits<br>
    3. Identify essential proteins<br>
    4. Screen virulence factors<br>
    5. Predict allergenicity<br>
    6. Compute physicochemical params<br>
    7. Predict antigenicity (IApred)<br>
    8. Predict subcellular localization
    </div>""", unsafe_allow_html=True)

# ============================================================
# TABS
# ============================================================
tabs = st.tabs([
    "Step 1 · Paralog",
    "Step 2 · Non-Homology",
    "Step 3 · Essential",
    "Step 4 · Virulence",
    "Step 5 · Allergenicity",
    "Step 6 · ProtParam",
    "Step 7 · Antigenicity",
    "Step 8 · Localization",
    "Summary"
])

# ============================================================
#  STEP 1 — PARALOG CLUSTERING
# ============================================================
with tabs[0]:
    step_card("01","Paralog Clustering & Representative Selection",
              "All-vs-all BLASTP to detect paralogous proteins. Uses graph-based connected components (NetworkX) to cluster paralogs, then picks the longest representative per cluster.")

    f1 = st.file_uploader("Upload Proteome FASTA", type=["fa","fasta","faa"], key="up1")
    if f1:
        raw1 = f1.read()
        recs1 = list(SeqIO.parse(StringIO(raw1.decode()), "fasta"))
        sec("INPUT SEQUENCES")
        with st.expander("View Input Analysis", expanded=True):
            input_viz(recs1, "s1in")
        if st.button("Run Paralog Clustering", key="btn1"):
            t0 = time.time()
            with st.spinner("Running all-vs-all BLASTP…"):
                with tempfile.TemporaryDirectory() as d:
                    pp=os.path.join(d,"p.fasta"); db=os.path.join(d,"db"); op=os.path.join(d,"b.tsv")
                    open(pp,"wb").write(raw1)
                    pdf=pd.DataFrame({"protein_id":[r.id for r in recs1],
                                      "sequence":[str(r.seq) for r in recs1],
                                      "protein_length":[len(r.seq) for r in recs1]})
                    subprocess.run(["makeblastdb","-in",pp,"-dbtype","prot","-out",db],check=True)
                    subprocess.run(["blastp","-query",pp,"-db",db,"-evalue",str(evalue_cutoff),
                                    "-num_threads",str(threads),"-outfmt","6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore","-out",op],check=True)
                    cols=["qseqid","sseqid","pident","align_length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
                    bdf=pd.read_csv(op,sep="\t",names=cols)
                    bdf=bdf[bdf.qseqid!=bdf.sseqid]
                    bdf=bdf.merge(pdf[["protein_id","protein_length"]],left_on="qseqid",right_on="protein_id",how="left")
                    bdf["coverage"]=(bdf["align_length"]/bdf["protein_length"])*100
                    bdf=bdf[(bdf.pident>=identity_cutoff)&(bdf.coverage>=coverage_cutoff)&(bdf.bitscore>=bitscore_cutoff)]
                    G=nx.Graph()
                    for _,row in bdf.iterrows(): G.add_edge(row.qseqid,row.sseqid)
                    clusters=list(nx.connected_components(G))
                    reps=[pdf[pdf.protein_id.isin(c)].sort_values("protein_length",ascending=False).iloc[0] for c in clusters]
                    rep_df=pd.DataFrame(reps) if reps else pd.DataFrame(columns=pdf.columns)
                    clustered=set().union(*clusters) if clusters else set()
                    single_df=pdf[~pdf.protein_id.isin(clustered)]
                    final_df=pd.concat([rep_df,single_df])
                    st.session_state.update({"s1_final":final_df,"s1_total":len(pdf),
                        "s1_reps":len(rep_df),"s1_singles":len(single_df),
                        "s1_rt":round(time.time()-t0,2),"step1_done":True})
                    st.rerun()

    if st.session_state.step1_done:
        st.success(f"Paralog clustering complete — {st.session_state.s1_rt}s")
        col1,col2,col3=st.columns(3)
        col1.metric("Input Proteins",st.session_state.s1_total)
        col2.metric("Paralog Representatives",st.session_state.s1_reps)
        col3.metric("Singleton Proteins",st.session_state.s1_singles)
        sec("OUTPUT ANALYSIS")
        with st.expander("View Output Visualizations", expanded=True):
            c1,c2=st.columns(2)
            with c1:
                fig=go.Figure(go.Pie(labels=["Representatives","Singletons"],
                    values=[st.session_state.s1_reps,st.session_state.s1_singles],hole=0.55,
                    marker=dict(colors=["#38bdf8","#818cf8"])))
                fig.update_layout(title="Protein Classification",**PLOT_LAYOUT)
                st.plotly_chart(fig,use_container_width=True)
            with c2:
                dfl=st.session_state.s1_final.rename(columns={"protein_length":"length"})
                st.plotly_chart(length_hist(dfl,"Output: Length Distribution","#818cf8"),use_container_width=True)
        sec("RESULT TABLE")
        st.dataframe(st.session_state.s1_final[["protein_id","protein_length"]].rename(
            columns={"protein_id":"Protein ID","protein_length":"Length (aa)"}),use_container_width=True)
        fasta_out="".join(f">{r.protein_id}|REPRESENTATIVE\n{r.sequence}\n" for _,r in st.session_state.s1_final.iterrows())
        st.download_button("⬇ Download Non-Redundant FASTA", fasta_out, "step1_nonredundant.fasta", key="dl1")

# ============================================================
#  STEP 2 — NON-HOMOLOGY
# ============================================================
with tabs[1]:
    step_card("02","Non-Homology Filter vs. Host Proteome",
              "BLASTP query against a reference host proteome (e.g. human). Retains only proteins lacking significant homology — candidates safe from cross-reactivity.")
    c1,c2=st.columns(2)
    with c1: qf2=st.file_uploader("Query Proteome (FASTA)",type=["fa","fasta","faa"],key="up2q")
    with c2: rf2=st.file_uploader("Reference Proteome (FASTA)",type=["fa","fasta","faa"],key="up2r")
    if qf2 and rf2:
        qb2=qf2.read(); rb2=rf2.read()
        rq2=list(SeqIO.parse(StringIO(qb2.decode()),"fasta"))
        rr2=list(SeqIO.parse(StringIO(rb2.decode()),"fasta"))
        with st.expander("View Input Analysis", expanded=True):
            c1,c2=st.columns(2)
            with c1: st.markdown("**Query**"); input_viz(rq2,"s2q")
            with c2: st.markdown("**Reference**"); input_viz(rr2,"s2r")
        if st.button("Run Non-Homology Filter", key="btn2"):
            t0=time.time()
            with st.spinner("Running BLASTP…"):
                with tempfile.TemporaryDirectory() as d:
                    qp=os.path.join(d,"q.fasta");rp=os.path.join(d,"r.fasta");db=os.path.join(d,"db");op=os.path.join(d,"o.tsv")
                    open(qp,"wb").write(qb2); open(rp,"wb").write(rb2)
                    qdf=pd.DataFrame({"qseqid":[r.id for r in rq2],"sequence":[str(r.seq) for r in rq2],"protein_length":[len(r.seq) for r in rq2]})
                    subprocess.run(["makeblastdb","-in",rp,"-dbtype","prot","-out",db],check=True)
                    subprocess.run(["blastp","-query",qp,"-db",db,"-evalue",str(evalue_cutoff),"-num_threads",str(threads),
                        "-outfmt","6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore","-out",op],check=True)
                    cols=["qseqid","sseqid","pident","align_length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
                    bdf=pd.read_csv(op,sep="\t",names=cols)
                    ql=dict(zip(qdf.qseqid,qdf.protein_length))
                    bdf["query_length"]=bdf.qseqid.map(ql)
                    bdf["query_coverage"]=(bdf.align_length/bdf.query_length)*100
                    filt=bdf[(bdf.pident>=identity_cutoff)&(bdf.query_coverage>=coverage_cutoff)&(bdf.bitscore>=bitscore_cutoff)&(bdf.evalue<=evalue_cutoff)]
                    hids=set(filt.qseqid)
                    nhdf=qdf[~qdf.qseqid.isin(hids)].copy()
                    st.session_state.update({"s2_nhdf":nhdf,"s2_total":len(qdf),"s2_hom":len(hids),
                        "s2_nonhom":len(nhdf),"s2_rt":round(time.time()-t0,2),"step2_done":True})
                    st.rerun()
    if st.session_state.step2_done:
        st.success(f"Non-homology filter complete — {st.session_state.s2_rt}s")
        c1,c2,c3=st.columns(3)
        c1.metric("Input Proteins",st.session_state.s2_total)
        c2.metric("Homologous (excluded)",st.session_state.s2_hom)
        c3.metric("Non-Homologous (kept)",st.session_state.s2_nonhom)
        with st.expander("View Output Visualizations", expanded=True):
            c1,c2=st.columns(2)
            with c1:
                fig=go.Figure(go.Pie(labels=["Non-Homologous","Homologous"],
                    values=[st.session_state.s2_nonhom,st.session_state.s2_hom],hole=0.55,
                    marker=dict(colors=["#34d399","#f87171"])))
                fig.update_layout(title="Homology Filter Result",**PLOT_LAYOUT)
                st.plotly_chart(fig,use_container_width=True)
            with c2:
                dfl=st.session_state.s2_nhdf.rename(columns={"protein_length":"length"})
                st.plotly_chart(length_hist(dfl,"Non-Homologous: Length Distribution","#34d399"),use_container_width=True)
        sec("RESULT TABLE")
        st.dataframe(st.session_state.s2_nhdf[["qseqid","protein_length"]].rename(columns={"qseqid":"Protein ID","protein_length":"Length (aa)"}),use_container_width=True)
        fo="".join(f">{r.qseqid}|NON_HOMOLOGOUS\n{r.sequence}\n" for _,r in st.session_state.s2_nhdf.iterrows())
        c1,c2=st.columns(2)
        with c1: st.download_button("⬇ Download CSV",st.session_state.s2_nhdf.to_csv(index=False),"step2_nonhom.csv",key="dl2c")
        with c2: st.download_button("⬇ Download FASTA",fo,"step2_nonhom.fasta",key="dl2f")

# ============================================================
#  STEP 3 — ESSENTIAL PROTEINS
# ============================================================
with tabs[2]:
    step_card("03","Essential Protein Detection (DEG Database)",
              "BLASTP against the Database of Essential Genes (DEG). Proteins with significant hits are flagged as essential — high-priority drug targets.")
    c1,c2=st.columns(2)
    with c1: pf3=st.file_uploader("Query Proteome (FASTA)",type=["fasta","fa","faa"],key="up3p")
    with c2:
        ef3=st.file_uploader("Essential Reference FASTA",type=["fasta","fa","faa"],key="up3e")
        st.caption("Reference: DEG database — https://drive.google.com/file/d/1woTl_J3TW4y1SVaNF2IhCTNP9uwk6WxF/view")
    if pf3 and ef3:
        pb3=pf3.read(); eb3=ef3.read()
        rp3=list(SeqIO.parse(StringIO(pb3.decode()),"fasta")); re3=list(SeqIO.parse(StringIO(eb3.decode()),"fasta"))
        with st.expander("View Input Analysis", expanded=True):
            c1,c2=st.columns(2)
            with c1: st.markdown("**Query**"); input_viz(rp3,"s3q")
            with c2: st.markdown("**Essential Reference**"); input_viz(re3,"s3r")
        if st.button("Run Essential Protein Detection", key="btn3"):
            t0=time.time()
            with st.spinner("Running BLASTP…"):
                with tempfile.TemporaryDirectory() as d:
                    pp=os.path.join(d,"q.faa");ep=os.path.join(d,"e.faa");db=os.path.join(d,"db");op=os.path.join(d,"o.tsv")
                    open(pp,"wb").write(pb3); open(ep,"wb").write(eb3)
                    recs=list(SeqIO.parse(pp,"fasta"))
                    subprocess.run(["makeblastdb","-in",ep,"-dbtype","prot","-out",db],check=True)
                    subprocess.run(["blastp","-query",pp,"-db",db,"-evalue",str(evalue_cutoff),"-num_threads",str(threads),
                        "-outfmt","6 qseqid sseqid pident length qlen slen evalue bitscore","-out",op],check=True)
                    bdf=pd.read_csv(op,sep="\t",names=["qseqid","sseqid","pident","length","qlen","slen","evalue","bitscore"])
                    bdf["coverage"]=(bdf.length/bdf.qlen)*100
                    edf=bdf[(bdf.pident>=identity_cutoff)&(bdf.coverage>=coverage_cutoff)&(bdf.bitscore>=bitscore_cutoff)&(bdf.evalue<=evalue_cutoff)]
                    eids=set(edf.qseqid)
                    ess=[r for r in recs if r.id in eids]; ness=[r for r in recs if r.id not in eids]
                    st.session_state.update({"s3_edf":edf,"s3_ess":ess,"s3_ness":ness,"s3_total":len(recs),
                        "s3_rt":round(time.time()-t0,2),"step3_done":True})
                    st.rerun()
    if st.session_state.step3_done:
        st.success(f"Essential detection complete — {st.session_state.s3_rt}s")
        c1,c2,c3=st.columns(3)
        c1.metric("Total Input",st.session_state.s3_total)
        c2.metric("Essential",len(st.session_state.s3_ess))
        c3.metric("Non-Essential",len(st.session_state.s3_ness))
        with st.expander("View Output Visualizations", expanded=True):
            c1,c2=st.columns(2)
            with c1:
                fig=go.Figure(go.Pie(labels=["Essential","Non-Essential"],
                    values=[len(st.session_state.s3_ess),len(st.session_state.s3_ness)],hole=0.55,
                    marker=dict(colors=["#fbbf24","#1c2840"])))
                fig.update_layout(title="Essential vs Non-Essential",**PLOT_LAYOUT)
                st.plotly_chart(fig,use_container_width=True)
            with c2:
                fig2=px.scatter(st.session_state.s3_edf,x="pident",y="bitscore",color="coverage",
                    title="BLAST: Identity vs Bitscore",color_continuous_scale="Blues",
                    labels={"pident":"% Identity","bitscore":"Bitscore"})
                sf(fig2); st.plotly_chart(fig2,use_container_width=True)
        sec("BLAST RESULTS TABLE")
        st.dataframe(st.session_state.s3_edf,use_container_width=True)
        ef_f="".join(f">{r.id}|ESSENTIAL\n{r.seq}\n" for r in st.session_state.s3_ess)
        nef_f="".join(f">{r.id}|NON_ESSENTIAL\n{r.seq}\n" for r in st.session_state.s3_ness)
        c1,c2,c3=st.columns(3)
        with c1: st.download_button("⬇ Essential FASTA",ef_f,"step3_essential.faa",key="dl3e")
        with c2: st.download_button("⬇ Non-Essential FASTA",nef_f,"step3_nonessential.faa",key="dl3n")
        with c3: st.download_button("⬇ BLAST Table (CSV)",st.session_state.s3_edf.to_csv(index=False),"step3_blast.csv",key="dl3c")

# ============================================================
#  STEP 4 — VIRULENCE
# ============================================================
with tabs[3]:
    step_card("04","Virulence Factor Detection (VFDB)",
              "Screens the proteome against a virulence factor database to identify pathogenicity-associated proteins — key for target prioritisation.")
    c1,c2=st.columns(2)
    with c1: pf4=st.file_uploader("Query Proteome (FASTA)",type=["fasta","fa","faa"],key="up4p")
    with c2:
        vf4=st.file_uploader("Virulence Reference FASTA",type=["fasta","fa","faa"],key="up4v")
        st.caption("Reference: VFDB — https://drive.google.com/file/d/1uh5h_fngqL8bOeFRp4vPyVhBJoymsGUd/view")
    if pf4 and vf4:
        pb4=pf4.read(); vb4=vf4.read()
        rp4=list(SeqIO.parse(StringIO(pb4.decode()),"fasta")); rv4=list(SeqIO.parse(StringIO(vb4.decode()),"fasta"))
        with st.expander("View Input Analysis", expanded=True):
            c1,c2=st.columns(2)
            with c1: st.markdown("**Query**"); input_viz(rp4,"s4q")
            with c2: st.markdown("**Virulence Reference**"); input_viz(rv4,"s4r")
        if st.button("Run Virulence Detection", key="btn4"):
            t0=time.time()
            with st.spinner("Running BLASTP…"):
                with tempfile.TemporaryDirectory() as d:
                    pp=os.path.join(d,"q.faa");vp=os.path.join(d,"v.faa");db=os.path.join(d,"db");op=os.path.join(d,"o.tsv")
                    open(pp,"wb").write(pb4); open(vp,"wb").write(vb4)
                    recs=list(SeqIO.parse(pp,"fasta"))
                    subprocess.run(["makeblastdb","-in",vp,"-dbtype","prot","-out",db],check=True)
                    subprocess.run(["blastp","-query",pp,"-db",db,"-evalue",str(evalue_cutoff),"-num_threads",str(threads),
                        "-outfmt","6 qseqid sseqid pident length qlen slen evalue bitscore","-out",op],check=True)
                    bdf=pd.read_csv(op,sep="\t",names=["qseqid","sseqid","pident","length","qlen","slen","evalue","bitscore"])
                    bdf["coverage"]=(bdf.length/bdf.qlen)*100
                    vdf=bdf[(bdf.pident>=identity_cutoff)&(bdf.coverage>=coverage_cutoff)&(bdf.bitscore>=bitscore_cutoff)&(bdf.evalue<=evalue_cutoff)]
                    vids=set(vdf.qseqid)
                    virs=[r for r in recs if r.id in vids]; nvirs=[r for r in recs if r.id not in vids]
                    st.session_state.update({"s4_vdf":vdf,"s4_virs":virs,"s4_nvirs":nvirs,"s4_total":len(recs),
                        "s4_rt":round(time.time()-t0,2),"step4_done":True})
                    st.rerun()
    if st.session_state.step4_done:
        st.success(f"Virulence detection complete — {st.session_state.s4_rt}s")
        c1,c2,c3=st.columns(3)
        c1.metric("Total Input",st.session_state.s4_total)
        c2.metric("Virulent",len(st.session_state.s4_virs))
        c3.metric("Non-Virulent",len(st.session_state.s4_nvirs))
        with st.expander("View Output Visualizations", expanded=True):
            c1,c2=st.columns(2)
            with c1:
                fig=go.Figure(go.Pie(labels=["Virulent","Non-Virulent"],
                    values=[len(st.session_state.s4_virs),len(st.session_state.s4_nvirs)],hole=0.55,
                    marker=dict(colors=["#f87171","#1c2840"])))
                fig.update_layout(title="Virulence Classification",**PLOT_LAYOUT)
                st.plotly_chart(fig,use_container_width=True)
            with c2:
                fig2=px.histogram(st.session_state.s4_vdf,x="pident",nbins=20,
                    title="Identity Distribution — Virulent Hits",color_discrete_sequence=["#f87171"])
                sf(fig2); st.plotly_chart(fig2,use_container_width=True)
        sec("BLAST RESULTS TABLE")
        st.dataframe(st.session_state.s4_vdf,use_container_width=True)
        vf_f="".join(f">{r.id}|VIRULENT\n{r.seq}\n" for r in st.session_state.s4_virs)
        nvf_f="".join(f">{r.id}|NON_VIRULENT\n{r.seq}\n" for r in st.session_state.s4_nvirs)
        c1,c2,c3=st.columns(3)
        with c1: st.download_button("⬇ Virulent FASTA",vf_f,"step4_virulent.faa",key="dl4v")
        with c2: st.download_button("⬇ Non-Virulent FASTA",nvf_f,"step4_nonvirulent.faa",key="dl4n")
        with c3: st.download_button("⬇ BLAST Table (CSV)",st.session_state.s4_vdf.to_csv(index=False),"step4_blast.csv",key="dl4c")

# ============================================================
#  STEP 5 — ALLERGENICITY
# ============================================================
with tabs[4]:
    step_card("05","Allergenicity Prediction (AllerTOP-style)",
              "Offline prediction based on hydrophobicity and polarity metrics. Compatible with CSV (Epitope column) and FASTA input formats.")
    st.info(f"Current threshold: **{allergen_thr}** — adjust in sidebar")
    uf5=st.file_uploader("Upload CSV (Epitope col) or FASTA",type=["csv","fasta","fa","faa","txt"],key="up5")
    if uf5:
        raw5=uf5.read(); seqs5=[]; ids5=[]
        if uf5.name.endswith(".csv"):
            dfc=pd.read_csv(StringIO(raw5.decode()))
            if "Epitope" not in dfc.columns: st.error("CSV must contain an 'Epitope' column.")
            else: seqs5=dfc["Epitope"].astype(str).tolist(); ids5=[f"Seq_{i+1}" for i in range(len(seqs5))]
        else:
            try:
                ra=list(SeqIO.parse(StringIO(raw5.decode()),"fasta"))
                seqs5=[str(r.seq) for r in ra]; ids5=[r.id for r in ra]
            except Exception as e: st.error(str(e))
        if seqs5:
            sec("INPUT SEQUENCES")
            with st.expander("View Input Analysis", expanded=True):
                dfl=pd.DataFrame({"id":ids5,"length":[len(s) for s in seqs5]})
                c1,c2=st.columns(2)
                with c1: st.plotly_chart(length_hist(dfl,"Input: Length Distribution"),use_container_width=True,key="s5lh")
                with c2: st.plotly_chart(aa_comp_bar("".join(seqs5),"Input: AA Composition"),use_container_width=True,key="s5ac")
            if st.button("Run Allergenicity Prediction", key="btn5"):
                t0=time.time()
                res5=[{"ID":ids5[i],"Length":len(s.strip().upper()),
                       "Allergenicity_Score":allergen_score(s.strip().upper()),
                       "Allergenicity":allergen_label(allergen_score(s.strip().upper()),allergen_thr)}
                      for i,s in enumerate(seqs5) if s.strip()]
                df5=pd.DataFrame(res5).sort_values("Allergenicity_Score",ascending=False).reset_index(drop=True)
                st.session_state.update({"s5_df":df5,"s5_seqs":seqs5,"s5_ids":ids5,
                    "s5_rt":round(time.time()-t0,2),"step5_done":True})
                st.rerun()
    if st.session_state.step5_done:
        df5=st.session_state.s5_df
        na=len(df5[df5.Allergenicity=="Allergen"]); nn=len(df5[df5.Allergenicity=="Non-Allergen"])
        st.success(f"Allergenicity analysis complete — {st.session_state.s5_rt}s")
        c1,c2,c3=st.columns(3)
        c1.metric("Total Sequences",len(df5)); c2.metric("Allergens",na); c3.metric("Non-Allergens",nn)
        with st.expander("View Output Visualizations", expanded=True):
            c1,c2,c3=st.columns(3)
            with c1:
                fig=go.Figure(go.Pie(labels=["Allergen","Non-Allergen"],values=[na,nn],hole=0.55,
                    marker=dict(colors=["#f87171","#34d399"])))
                fig.update_layout(title="Classification",**PLOT_LAYOUT)
                st.plotly_chart(fig,use_container_width=True)
            with c2:
                fig2=px.histogram(df5,x="Allergenicity_Score",nbins=20,color="Allergenicity",
                    color_discrete_map={"Allergen":"#f87171","Non-Allergen":"#34d399"},title="Score Distribution")
                sf(fig2); st.plotly_chart(fig2,use_container_width=True)
            with c3:
                fig3=px.scatter(df5,x="Length",y="Allergenicity_Score",color="Allergenicity",
                    color_discrete_map={"Allergen":"#f87171","Non-Allergen":"#34d399"},title="Length vs Score")
                fig3.add_hline(y=allergen_thr,line_dash="dash",line_color="#fbbf24",
                    annotation_text=f"Threshold ({allergen_thr})")
                sf(fig3); st.plotly_chart(fig3,use_container_width=True)
        sec("RESULTS TABLE")
        st.dataframe(df5,use_container_width=True)
        nn_df=df5[df5.Allergenicity=="Non-Allergen"]
        ids_l=st.session_state.s5_ids; seqs_l=st.session_state.s5_seqs
        nn_fa="".join(f">{r.ID}|Score={r.Allergenicity_Score}\n{seqs_l[ids_l.index(r.ID)]}\n"
                      for _,r in nn_df.iterrows() if r.ID in ids_l)
        c1,c2=st.columns(2)
        with c1: st.download_button("⬇ All Results (CSV)",df5.to_csv(index=False),"step5_allergenicity.csv",key="dl5c")
        with c2: st.download_button("⬇ Non-Allergen FASTA",nn_fa,"step5_nonallergen.faa",key="dl5f")

# ============================================================
#  STEP 6 — PROTPARAM
# ============================================================
with tabs[5]:
    step_card("06","ProtParam Physicochemical Analysis",
              "Fetches molecular weight, instability index, aliphatic index and GRAVY score from ExPASy ProtParam. Requires internet connection.")
    uf6=st.file_uploader("Upload FASTA",type=["fasta","fa","faa"],key="up6")
    if uf6:
        raw6=uf6.read()
        recs6=list(SeqIO.parse(StringIO(raw6.decode()),"fasta"))
        with st.expander("View Input Analysis", expanded=True):
            input_viz(recs6,"s6in")
        st.info(f"30 aa chunk mode: **{'ON' if chunk_opt else 'OFF'}** — toggle in sidebar")
        if st.button("Compute ProtParam", key="btn6"):
            t0=time.time(); all6=[]; prog=st.progress(0,"Fetching ProtParam data…")
            try:
                for i,rec in enumerate(recs6):
                    prog.progress((i+1)/len(recs6), text=f"Processing {rec.id} ({i+1}/{len(recs6)})…")
                    if chunk_opt:
                        for ci,ch in enumerate(textwrap.wrap(str(rec.seq).upper(),30),1):
                            r=get_protparam(ch); r.update({"Protein_ID":rec.id,"Chunk":ci,"Sequence":ch}); all6.append(r)
                    else:
                        r=get_protparam(str(rec.seq).upper()); r.update({"Protein_ID":rec.id,"Chunk":1,"Sequence":str(rec.seq).upper()}); all6.append(r)
                prog.empty()
                df6=pd.DataFrame(all6)[['Protein_ID','Chunk','Sequence','Molecular weight','Instability index','Instability','Aliphatic index','GRAVY score']]
                st.session_state.update({"s6_df":df6,"s6_rt":round(time.time()-t0,2),"step6_done":True})
                st.rerun()
            except Exception as ex: st.error(f"ProtParam error: {ex}")
    if st.session_state.step6_done:
        df6=st.session_state.s6_df
        st_df=df6[df6.Instability=="Stable"]; un_df=df6[df6.Instability=="Unstable"]
        st.success(f"ProtParam complete — {st.session_state.s6_rt}s")
        c1,c2,c3,c4=st.columns(4)
        c1.metric("Total",len(df6)); c2.metric("Stable",len(st_df)); c3.metric("Unstable",len(un_df))
        c4.metric("Avg MW (Da)",round(df6["Molecular weight"].dropna().mean(),1))
        with st.expander("View Output Visualizations", expanded=True):
            c1,c2=st.columns(2)
            with c1:
                fig=go.Figure(go.Pie(labels=["Stable","Unstable"],values=[len(st_df),len(un_df)],hole=0.55,
                    marker=dict(colors=["#34d399","#f87171"])))
                fig.update_layout(title="Stability Classification",**PLOT_LAYOUT)
                st.plotly_chart(fig,use_container_width=True)
            with c2:
                fig2=px.scatter(df6.dropna(subset=["Molecular weight","GRAVY score"]),
                    x="Molecular weight",y="GRAVY score",color="Instability",hover_data=["Protein_ID"],
                    color_discrete_map={"Stable":"#34d399","Unstable":"#f87171"},title="MW vs GRAVY Score")
                sf(fig2); st.plotly_chart(fig2,use_container_width=True)
            c3,c4=st.columns(2)
            with c3:
                fig3=px.histogram(df6.dropna(subset=["Instability index"]),x="Instability index",nbins=20,color="Instability",
                    color_discrete_map={"Stable":"#34d399","Unstable":"#f87171"},title="Instability Index Distribution")
                fig3.add_vline(x=40,line_dash="dash",line_color="#fbbf24",annotation_text="Threshold (40)")
                sf(fig3); st.plotly_chart(fig3,use_container_width=True)
            with c4:
                fig4=px.box(df6.dropna(subset=["Aliphatic index"]),x="Instability",y="Aliphatic index",color="Instability",
                    color_discrete_map={"Stable":"#34d399","Unstable":"#f87171"},title="Aliphatic Index by Stability")
                sf(fig4); st.plotly_chart(fig4,use_container_width=True)
        sec("RESULTS TABLE")
        st.dataframe(df6,use_container_width=True)
        sf_fa="\n".join(f">{r.Protein_ID}_Chunk{r.Chunk}\n{r.Sequence}" for _,r in st_df.iterrows())
        c1,c2=st.columns(2)
        with c1: st.download_button("⬇ All Results (CSV)",df6.to_csv(index=False),"step6_protparam.csv",key="dl6c")
        with c2: st.download_button("⬇ Stable Sequences (FASTA)",sf_fa,"step6_stable.fasta",key="dl6f")

# ============================================================
#  STEP 7 — ANTIGENICITY  (IApred)
# ============================================================
with tabs[6]:
    step_card("07","Antigenicity Prediction (IApred)",
              "Predicts intrinsic antigenicity using IApred. Point the path below to your local IApred.py. Extracts high-antigenicity sequences for downstream use.")

    # ── IApred path configuration ──────────────────────────
    st.markdown('<div class="sec-label">IAPRED SCRIPT PATH</div>', unsafe_allow_html=True)

    # Auto-detect common locations
    _default_candidates = [
        os.path.join("IAPred-main", "IApred.py"),
        os.path.join("IApred-main", "IApred.py"),
        os.path.join("IApred", "IApred.py"),
        "IApred.py",
    ]
    _auto_path = next((p for p in _default_candidates if os.path.isfile(p)), "IAPred-main/IApred.py")

    iapred_path = st.text_input(
        "Path to IApred.py",
        value=st.session_state.get("iapred_path_val", _auto_path),
        help="Relative or absolute path. E.g.  IAPred-main/IApred.py  or  C:/tools/IAPred-main/IApred.py",
        key="iapred_path_input"
    )
    st.session_state["iapred_path_val"] = iapred_path

    # Show path status
    if os.path.isfile(iapred_path):
        st.success(f"✅ Found: `{os.path.abspath(iapred_path)}`")
    else:
        st.warning(f"⚠ Script not found at `{iapred_path}` — check the path or use the CSV upload below.")

    st.markdown("---")

    uf7=st.file_uploader("Upload FASTA file",type=["fasta","fa","faa"],key="up7")
    if uf7:
        raw7=uf7.read()
        recs7=list(SeqIO.parse(StringIO(raw7.decode()),"fasta"))
        sec("INPUT SEQUENCES")
        with st.expander("View Input Analysis", expanded=True):
            input_viz(recs7,"s7in")

        if st.button("Run IApred Prediction", key="btn7"):
            if not os.path.isfile(iapred_path):
                st.error(f"IApred.py not found at `{iapred_path}`. Please correct the path above.")
            else:
                t0=time.time()
                with st.spinner("Running IApred…"):
                    try:
                        import sys
                        # Write input FASTA to a temp file in the SAME folder as IApred.py
                        # so relative imports inside IApred.py resolve correctly
                        iapred_dir = os.path.dirname(os.path.abspath(iapred_path))
                        iapred_script = os.path.abspath(iapred_path)

                        with tempfile.NamedTemporaryFile(
                            delete=False, suffix=".fasta",
                            dir=iapred_dir, prefix="iapred_input_"
                        ) as tmp:
                            tmp.write(raw7)
                            fasta_path7 = tmp.name

                        out_csv7 = fasta_path7.replace(".fasta", "_results.csv")

                        result = subprocess.run(
                            [sys.executable, iapred_script, fasta_path7, out_csv7],
                            capture_output=True, text=True,
                            cwd=iapred_dir   # run from IApred's own directory
                        )

                        if result.returncode != 0:
                            st.error("IApred returned an error. See details below.")
                            with st.expander("🔍 Error Details (stderr / stdout)", expanded=True):
                                if result.stderr.strip():
                                    st.markdown("**stderr:**")
                                    st.code(result.stderr, language="text")
                                if result.stdout.strip():
                                    st.markdown("**stdout:**")
                                    st.code(result.stdout, language="text")
                                st.markdown(f"**Return code:** `{result.returncode}`")
                                st.markdown(f"**Script used:** `{iapred_script}`")
                                st.markdown(f"**Input FASTA:** `{fasta_path7}`")
                                st.markdown(f"**Expected output CSV:** `{out_csv7}`")
                        elif os.path.exists(out_csv7):
                            csv7=pd.read_csv(out_csv7)
                            ant_fa,matched7=map_antigenic(csv7,recs7)
                            st.session_state.update({"s7_csv":csv7,"s7_ant_fa":ant_fa,
                                "s7_matched":matched7,"s7_rt":round(time.time()-t0,2),"step7_done":True})
                            # Clean up temp files
                            try: os.remove(fasta_path7)
                            except: pass
                            st.rerun()
                        else:
                            st.warning("IApred ran but no output CSV was produced.")
                            with st.expander("🔍 Debug Info"):
                                st.markdown(f"**Expected output at:** `{out_csv7}`")
                                st.markdown(f"**IApred directory contents:**")
                                st.code("\n".join(os.listdir(iapred_dir)), language="text")
                                if result.stdout.strip():
                                    st.code(result.stdout, language="text")

                    except Exception as e:
                        st.error(f"Unexpected error: {e}")

        # ── Manual CSV upload fallback ─────────────────────
        st.markdown("---")
        st.markdown('<div class="sec-label">MANUAL CSV UPLOAD — FALLBACK</div>', unsafe_allow_html=True)
        st.caption("Already ran IApred externally? Upload its output CSV directly here.")
        man_csv7=st.file_uploader("IApred Results CSV",type=["csv"],key="up7csv")
        if man_csv7 and uf7:
            csv7=pd.read_csv(man_csv7)
            ant_fa,matched7=map_antigenic(csv7,recs7)
            st.session_state.update({"s7_csv":csv7,"s7_ant_fa":ant_fa,
                "s7_matched":matched7,"s7_rt":0,"step7_done":True})
            st.rerun()

    if st.session_state.step7_done:
        csv7=st.session_state.s7_csv
        st.success(f"Antigenicity prediction complete — {st.session_state.s7_rt}s")

        n_high=len(csv7[csv7.get("Antigenicity_Category","") == "High"]) if "Antigenicity_Category" in csv7.columns else len(st.session_state.s7_matched)
        c1,c2,c3=st.columns(3)
        c1.metric("Total Sequences",len(csv7))
        c2.metric("High Antigenicity",n_high)
        c3.metric("Matched to FASTA",len(st.session_state.s7_matched))

        with st.expander("View Output Visualizations", expanded=True):
            if "Antigenicity_Category" in csv7.columns:
                c1,c2=st.columns(2)
                with c1:
                    counts=csv7["Antigenicity_Category"].value_counts()
                    fig=go.Figure(go.Pie(labels=counts.index.tolist(),values=counts.values.tolist(),hole=0.55,
                        marker=dict(colors=["#38bdf8","#818cf8","#34d399","#fbbf24"])))
                    fig.update_layout(title="Antigenicity Categories",**PLOT_LAYOUT)
                    st.plotly_chart(fig,use_container_width=True)
                with c2:
                    fig2=px.bar(counts.reset_index(),x="Antigenicity_Category",y="count",
                        title="Category Distribution",color="Antigenicity_Category",
                        color_discrete_sequence=["#38bdf8","#818cf8","#34d399","#fbbf24"])
                    sf(fig2); st.plotly_chart(fig2,use_container_width=True)
            else:
                st.info("Visualization requires 'Antigenicity_Category' column in the IApred CSV.")

        sec("PREDICTION RESULTS TABLE")
        st.dataframe(csv7,use_container_width=True)

        if st.session_state.s7_ant_fa.strip():
            with st.expander("Antigenic FASTA Preview"):
                st.code(st.session_state.s7_ant_fa[:2000]+("…" if len(st.session_state.s7_ant_fa)>2000 else ""),language="text")

        c1,c2,c3=st.columns(3)
        with c1: st.download_button("⬇ Full Results (CSV)",csv7.to_csv(index=False),"step7_iapred_full.csv",key="dl7c")
        with c2: st.download_button("⬇ High Antigenicity Only (CSV)",
            csv7[csv7.get("Antigenicity_Category","")=="High"].to_csv(index=False) if "Antigenicity_Category" in csv7.columns else csv7.to_csv(index=False),
            "step7_antigenic_records.csv",key="dl7h")
        with c3: st.download_button("⬇ Antigenic Sequences (FASTA)",st.session_state.s7_ant_fa,"step7_antigenic.faa",key="dl7f")

# ============================================================
#  STEP 8 — SUBCELLULAR LOCALIZATION  (DeepLocPro)
# ============================================================
with tabs[7]:
    step_card("08","Subcellular Localization Prediction (DeepLocPro)",
              "Predicts prokaryotic subcellular localization using DeepLocPro. Supports Gram-negative and Gram-positive bacteria. Outputs per-localization FASTA files.")

    uf8=st.file_uploader("Upload Protein FASTA",type=["fasta","fa","faa"],key="up8")

    group_label8=st.radio("Bacterial type",
        ["Gram-negative bacteria","Gram-positive bacteria"],key="grp8")
    group_map8={"Gram-negative bacteria":"negative","Gram-positive bacteria":"positive"}
    grp8=group_map8[group_label8]

    if uf8:
        raw8=uf8.read()
        recs8=list(SeqIO.parse(StringIO(raw8.decode()),"fasta"))
        sec("INPUT SEQUENCES")
        with st.expander("View Input Analysis", expanded=True):
            input_viz(recs8,"s8in")

        if st.button("Run DeepLocPro", key="btn8"):
            t0=time.time()
            with st.spinner("Running DeepLocPro…"):
                temp8="temp_dlp_input.fasta"
                open(temp8,"wb").write(raw8)
                out8="deeploc_outputs"
                os.makedirs(out8,exist_ok=True)
                cmd8=f"deeplocpro -f {temp8} -o {out8} -g {grp8}"
                ret8=subprocess.run(cmd8,shell=True)
                if ret8.returncode!=0:
                    st.error("DeepLocPro execution failed. Ensure deeplocpro is installed.")
                else:
                    csvs=glob.glob(os.path.join(out8,"*.csv"))
                    if not csvs:
                        st.error("No output CSV found from DeepLocPro.")
                    else:
                        latest_csv=max(csvs,key=os.path.getmtime)
                        st.session_state.update({"s8_csv_path":latest_csv,"s8_raw8":raw8,
                            "s8_total":len(recs8),"s8_grp":group_label8,"s8_grpv":grp8,
                            "s8_rt":round(time.time()-t0,2),"step8_done":True})
                        st.rerun()

        # Manual CSV upload fallback
        st.markdown("---")
        st.markdown("**— or upload an existing DeepLocPro output CSV —**")
        man_csv8=st.file_uploader("DeepLocPro Results CSV",type=["csv"],key="up8csv")
        if man_csv8:
            tmp_path8="temp_dlp_manual.csv"
            open(tmp_path8,"wb").write(man_csv8.read())
            st.session_state.update({"s8_csv_path":tmp_path8,"s8_raw8":raw8,
                "s8_total":len(recs8),"s8_grp":group_label8,"s8_grpv":grp8,
                "s8_rt":0,"step8_done":True})
            st.rerun()

    if st.session_state.step8_done and "s8_csv_path" in st.session_state:
        df8=pd.read_csv(st.session_state.s8_csv_path)
        st.success(f"Localization prediction complete — {st.session_state.s8_rt}s")

        # Build acc→loc map
        acc_loc={}
        for _,row in df8.iterrows():
            acc=str(row["ACC"]); loc=str(row["Localization"]).lower().replace(" ","_")
            parts=acc.split("|"); pid=parts[1] if len(parts)>=2 else acc
            acc_loc[pid.strip()]=loc

        # Split FASTA by localization
        loc_recs={}; all_matched=[]
        raw8_use=st.session_state.s8_raw8
        for rec in SeqIO.parse(StringIO(raw8_use.decode()),"fasta"):
            parts=rec.id.split("|"); rid=parts[1] if len(parts)>=2 else rec.id
            if rid.strip() in acc_loc:
                loc=acc_loc[rid.strip()]
                loc_recs.setdefault(loc,[]).append(rec)
                all_matched.append(rec)

        if st.session_state.s8_grpv=="positive":
            INVALID={"outer_membrane","periplasmic"}
            loc_recs={k:v for k,v in loc_recs.items() if k not in INVALID}
            all_matched=[r for r in all_matched
                if acc_loc.get(r.id.split("|")[1] if "|" in r.id else r.id,"") not in INVALID]

        c1,c2=st.columns(2)
        c1.metric("Input Proteins",st.session_state.s8_total)
        c2.metric("Matched & Classified",len(all_matched))

        with st.expander("View Output Visualizations", expanded=True):
            stats8=pd.DataFrame([{"Localization":k,"Count":len(v)} for k,v in loc_recs.items()]).sort_values("Count",ascending=False)
            c1,c2=st.columns(2)
            with c1:
                colors8=["#38bdf8","#818cf8","#34d399","#fbbf24","#f87171","#a78bfa","#fb923c","#22d3ee"]
                fig=go.Figure(go.Pie(labels=stats8.Localization.tolist(),values=stats8.Count.tolist(),hole=0.5,
                    marker=dict(colors=colors8[:len(stats8)])))
                fig.update_layout(title="Localization Distribution",**PLOT_LAYOUT)
                st.plotly_chart(fig,use_container_width=True)
            with c2:
                fig2=px.bar(stats8,x="Localization",y="Count",title="Proteins per Compartment",
                    color="Localization",color_discrete_sequence=colors8)
                fig2.update_xaxes(tickangle=-30)
                sf(fig2); st.plotly_chart(fig2,use_container_width=True)

        sec("RESULTS TABLE")
        st.dataframe(df8,use_container_width=True)

        sec("LOCALIZATION SUMMARY")
        st.dataframe(stats8,use_container_width=True)

        # Build ZIP
        zip_buf=BytesIO()
        with zipfile.ZipFile(zip_buf,"w") as zf:
            for loc,recs_l in loc_recs.items():
                fa_io=StringIO(); SeqIO.write(recs_l,fa_io,"fasta"); zf.writestr(f"{loc}.fasta",fa_io.getvalue())
            all_io=StringIO(); SeqIO.write(all_matched,all_io,"fasta"); zf.writestr("all_classified.fasta",all_io.getvalue())
            zf.writestr("localization_stats.csv",stats8.to_csv(index=False))

        c1,c2=st.columns(2)
        with c1: st.download_button("⬇ Full Results (CSV)",df8.to_csv(index=False),"step8_deeploc.csv",key="dl8c")
        with c2: st.download_button("⬇ FASTAs + Stats (ZIP)",zip_buf.getvalue(),"step8_deeploc_fastas.zip",key="dl8z")

# ============================================================
#  SUMMARY TAB
# ============================================================
with tabs[8]:
    st.markdown("""
    <div class="step-card">
        <span class="step-badge">Summary</span>
        <div class="step-title">Pipeline Completion Overview</div>
        <p class="step-desc">Track the status of all 8 analysis stages and review results across the full workflow.</p>
    </div>""", unsafe_allow_html=True)

    SUMMARY = [
        ("01","Paralog Clustering","step1_done",
         lambda: f"{st.session_state.get('s1_total',0)} input → {len(st.session_state.get('s1_final',pd.DataFrame()))} representative proteins",
         "#38bdf8"),
        ("02","Non-Homology Filter","step2_done",
         lambda: f"{st.session_state.get('s2_total',0)} input → {st.session_state.get('s2_nonhom',0)} non-homologous",
         "#818cf8"),
        ("03","Essential Proteins","step3_done",
         lambda: f"{st.session_state.get('s3_total',0)} input → {len(st.session_state.get('s3_ess',[]))} essential",
         "#fbbf24"),
        ("04","Virulence Detection","step4_done",
         lambda: f"{st.session_state.get('s4_total',0)} input → {len(st.session_state.get('s4_virs',[]))} virulent",
         "#f87171"),
        ("05","Allergenicity","step5_done",
         lambda: f"{len(st.session_state.get('s5_df',pd.DataFrame()))} sequences screened",
         "#34d399"),
        ("06","ProtParam","step6_done",
         lambda: f"{len(st.session_state.get('s6_df',pd.DataFrame()))} parameters computed",
         "#a78bfa"),
        ("07","Antigenicity (IApred)","step7_done",
         lambda: f"{len(st.session_state.get('s7_csv',pd.DataFrame()))} predictions · {len(st.session_state.get('s7_matched',set()))} matched",
         "#38bdf8"),
        ("08","Localization (DeepLocPro)","step8_done",
         lambda: f"{st.session_state.get('s8_total',0)} input proteins classified",
         "#fbbf24"),
    ]

    for num,name,dk,rfn,color in SUMMARY:
        done=st.session_state.get(dk,False)
        icon="✅" if done else "⏳"
        txt=rfn() if done else "Not yet run"
        r=int(color.lstrip('#')[0:2],16); g=int(color.lstrip('#')[2:4],16); b=int(color.lstrip('#')[4:6],16)
        bg=f"rgba({r},{g},{b},0.07)" if done else "rgba(22,32,48,0.6)"
        bc=f"rgba({r},{g},{b},0.25)" if done else "rgba(71,85,105,0.2)"
        nc=color if done else "#475569"
        st.markdown(f"""
        <div style="display:flex;align-items:center;gap:1rem;padding:0.9rem 1.2rem;
                    background:{bg};border:1px solid {bc};border-radius:10px;margin-bottom:0.5rem;">
            <div style="font-size:1.1rem;line-height:1">{icon}</div>
            <div style="flex:1">
                <div style="font-family:'DM Sans',sans-serif;font-weight:700;font-size:0.88rem;color:{nc}">
                    Step {num} — {name}
                </div>
                <div style="font-family:'IBM Plex Mono',monospace;font-size:0.72rem;color:#475569;margin-top:2px">{txt}</div>
            </div>
        </div>""", unsafe_allow_html=True)

    completed=sum(1 for _,_,dk,_,_ in SUMMARY if st.session_state.get(dk,False))
    pct=int(completed/8*100)
    st.markdown(f"""
    <div style="margin:1.5rem 0 1rem;padding:1.4rem;background:rgba(56,189,248,0.04);
                border:1px solid rgba(56,189,248,0.15);border-radius:12px;text-align:center;">
        <div style="font-family:'DM Serif Display',serif;font-size:2.5rem;color:#38bdf8;line-height:1">
            {completed} <span style="font-size:1.2rem;color:#475569">/ 8</span>
        </div>
        <div style="font-family:'IBM Plex Mono',monospace;font-size:0.68rem;color:#475569;
                    text-transform:uppercase;letter-spacing:2px;margin-top:6px">
            Stages Completed · {pct}% of pipeline
        </div>
        <div style="margin-top:12px;height:5px;background:#162030;border-radius:99px;overflow:hidden">
            <div style="height:100%;width:{pct}%;background:linear-gradient(90deg,#38bdf8,#818cf8);border-radius:99px;transition:width 0.5s"></div>
        </div>
    </div>""", unsafe_allow_html=True)

    if completed >= 2:
        sec("CROSS-STEP PROTEIN RETENTION FUNNEL")
        funnel=[]
        if st.session_state.get("step1_done"): funnel.append(("01 · Paralog Output",len(st.session_state.get("s1_final",pd.DataFrame()))))
        if st.session_state.get("step2_done"): funnel.append(("02 · Non-Homologous",st.session_state.get("s2_nonhom",0)))
        if st.session_state.get("step3_done"): funnel.append(("03 · Essential",len(st.session_state.get("s3_ess",[]))))
        if st.session_state.get("step4_done"): funnel.append(("04 · Virulent",len(st.session_state.get("s4_virs",[]))))
        if funnel:
            fig_f=go.Figure(go.Funnel(
                y=[d[0] for d in funnel],x=[d[1] for d in funnel],
                marker=dict(color=["#38bdf8","#818cf8","#fbbf24","#f87171"][:len(funnel)]),
                textinfo="value+percent previous",
                textfont=dict(family="IBM Plex Mono",size=11)
            ))
            fig_f.update_layout(title="Protein Retention Through Pipeline",**PLOT_LAYOUT)
            st.plotly_chart(fig_f,use_container_width=True)
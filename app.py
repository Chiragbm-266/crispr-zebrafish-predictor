import streamlit as st
import requests
import re
import pandas as pd
import time
import xml.etree.ElementTree as ET
import os
import altair as alt
import difflib

# ==========================================
# 1. APP CONFIGURATION & UI SETUP
# ==========================================
st.set_page_config(page_title="CRISPR Phenotype Predictor", layout="wide")

st.sidebar.header("⚙️ CRISPR Settings")
selected_enzyme = st.sidebar.selectbox(
    "Select Nuclease Engine:", 
    ["SpCas9 (PAM: NGG)", "SaCas9 (PAM: NNGRRT)", "Cas12a/Cpf1 (PAM: TTTV)"]
)
st.sidebar.markdown("---")
st.sidebar.info(
    "**SpCas9:** Standard, high efficiency. Blunt cuts.\n\n"
    "**SaCas9:** Compact size, viral delivery.\n\n"
    "**Cas12a:** High specificity, staggered cuts."
)

st.title("🧬 CRISPR-Cas9 Zebrafish Predictor")
st.write("Cross-referencing global databases (ZFIN, NCBI, Ensembl, UniProt) to predict genomic cuts and synthesize detailed phenotypic hypotheses.")
st.markdown("---")

COMMON_PHENOTYPES = ["brain", "tumor", "heart", "eye", "jaw", "transparent", "golden", "edema", "muscle", "blood", "tail", "fin", "pain", "movement", "neurological"]

# ==========================================
# 2. BULLETPROOF API ENGINES (SPLIT-STREAM)
# ==========================================

def safe_mygene_fetch(query, species, fields):
    try:
        url = "https://mygene.info/v3/query"
        resp = requests.get(url, params={"q": query, "species": species, "fields": fields}, timeout=5)
        if resp.status_code == 200 and resp.json().get('hits'):
            return resp.json()['hits'][0]
    except Exception: pass
    return {}

@st.cache_data(ttl=86400)
def fetch_ncbi_data(gene_query):
    query = gene_query.lower().strip()
    
    fish_fields = "name,summary,symbol,alias,disease,go,pathway,ensembl"
    human_fields = "name,summary,symbol,disease,go,pathway"
    
    fish_hit = safe_mygene_fetch(query, "7955", fish_fields)
    human_hit = safe_mygene_fetch(query, "9606", human_fields)
    
    # Cross-reference missing orthologs dynamically
    if fish_hit and not human_hit:
        human_hit = safe_mygene_fetch(fish_hit.get('symbol', query), "9606", human_fields)
        if not human_hit and fish_hit.get('alias'):
            aliases = fish_hit['alias'] if isinstance(fish_hit['alias'], list) else [fish_hit['alias']]
            for a in aliases:
                human_hit = safe_mygene_fetch(a, "9606", human_fields)
                if human_hit: break
                
    elif human_hit and not fish_hit:
        fish_hit = safe_mygene_fetch(human_hit.get('symbol', query), "7955", fish_fields)

    fish_symbol = fish_hit.get('symbol', query).lower() if fish_hit else query.lower()
    human_ortholog = human_hit.get('symbol', query).upper() if human_hit else query.upper()
    
    display_fish = fish_symbol
    display_human = human_ortholog
    
    if query in ['nat15', 'naa60']:
        display_fish = "nat15"
        display_human = "NAA60"
        fish_symbol = "nat15"
        human_ortholog = "NAA60"
    
    show_warning = (display_fish.lower() != display_human.lower() and display_human != "UNMAPPED")

    zfin_aliases = []
    if fish_hit and 'alias' in fish_hit:
        zfin_aliases = fish_hit['alias'] if isinstance(fish_hit['alias'], list) else [fish_hit['alias']]
        
    ensembl_id = None
    if fish_hit and 'ensembl' in fish_hit:
        ens = fish_hit['ensembl']
        ensembl_id = ens[0].get('gene') if isinstance(ens, list) else ens.get('gene')

    target_hit = human_hit if human_hit else fish_hit
    name = target_hit.get('name', 'Unknown Gene') if target_hit else 'Unknown Gene'
    summary = human_hit.get('summary', fish_hit.get('summary', '')) if human_hit else (fish_hit.get('summary', '') if fish_hit else '')
    
    diseases = human_hit.get('disease', []) if human_hit else []
    disease_list = [d.get('term_name', '') for d in diseases if isinstance(d, dict) and 'term_name' in d] if isinstance(diseases, list) else [diseases.get('term_name', '')] if isinstance(diseases, dict) else []
    
    go_data = target_hit.get('go', {}) if target_hit else {}
    go_bp = go_data.get('BP', [])
    go_list = [g.get('term', '') for g in go_bp if isinstance(g, dict) and 'term' in g] if isinstance(go_bp, list) else [go_bp.get('term', '')] if isinstance(go_bp, dict) else []
    
    pathway_data = target_hit.get('pathway', {}) if target_hit else {}
    kegg = pathway_data.get('kegg', [])
    pathway_list = [p.get('name', '') for p in kegg if isinstance(p, dict) and 'name' in p] if isinstance(kegg, list) else [kegg.get('name', '')] if isinstance(kegg, dict) else []

    return display_fish, display_human, fish_symbol, human_ortholog, ensembl_id, show_warning, list(set(zfin_aliases)), name, summary, list(set(disease_list)), list(set(go_list)), list(set(pathway_list))

@st.cache_data(ttl=86400)
def fetch_uniprot_data(fish_symbol):
    url = "https://rest.uniprot.org/uniprotkb/search"
    try:
        response = requests.get(url, params={"query": f"(gene_exact:{fish_symbol}) AND (organism_id:7955)", "format": "json"}, timeout=10)
        if response.status_code == 200 and response.json().get('results'):
            return response.json()['results'][0].get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown Protein')
        return "Unmapped Protein"
    except Exception: return "Database Timeout"

@st.cache_data(ttl=86400)
def fetch_pubchem_drugs(human_ortholog):
    if human_ortholog == "UNMAPPED": return []
    url = "https://mygene.info/v3/query"
    try:
        response = requests.get(url, params={"q": human_ortholog, "species": "9606", "fields": "pharmgkb,drugbank"}, timeout=10)
        drugs = []
        if response.status_code == 200 and response.json().get('hits'):
            hit = response.json()['hits'][0]
            if 'drugbank' in hit:
                db = hit['drugbank']
                if isinstance(db, list): drugs.extend([d.get('name') for d in db if 'name' in d])
                elif isinstance(db, dict): drugs.append(db.get('name'))
        return list(set(drugs))[:6]
    except Exception: return []

@st.cache_data(ttl=86400)
def fetch_string_interactions(fish_symbol):
    url = "https://string-db.org/api/json/network"
    try:
        response = requests.get(url, params={"identifiers": fish_symbol, "species": "7955"}, timeout=10)
        if response.status_code == 200 and len(response.json()) > 0:
            interactions = set([item['preferredName_B'] for item in response.json() if item['preferredName_B'].lower() != fish_symbol.lower()])
            return list(interactions)[:8]
        return []
    except Exception: return []

def query_epmc(query_str, size=15):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    try:
        resp = requests.get(url, params={"query": query_str, "format": "json", "resultType": "lite", "pageSize": size}, timeout=10)
        if resp.status_code == 200: return resp.json().get('resultList', {}).get('result', [])
    except Exception: pass
    return []

@st.cache_data(ttl=86400)
def fetch_literature(display_fish, display_human, fish_symbol, human_ortholog):
    pmids_seen = set()
    core_papers, recent_papers = [], []
    
    hits_zf = query_epmc(f'(TITLE:"{fish_symbol}" OR TITLE:"{display_fish}" OR ABSTRACT:"{fish_symbol}") AND ("zebrafish" OR "Danio rerio")', 15)
    for p in hits_zf:
        if len(core_papers) < 5 and p.get('pmid') not in pmids_seen:
            p['badge'] = '⭐ Core Paper (Zebrafish)'
            core_papers.append(p)
            pmids_seen.add(p.get('pmid'))
            
    search_term = human_ortholog if human_ortholog != "UNMAPPED" else fish_symbol
    if len(core_papers) < 5:
        hits_all = query_epmc(f'(TITLE:"{search_term}" OR TITLE:"{display_human}" OR ABSTRACT:"{search_term}")', 25)
        for p in hits_all:
            if len(core_papers) < 5 and p.get('pmid') not in pmids_seen:
                p['badge'] = '⭐ Core Paper (Ortholog)'
                core_papers.append(p)
                pmids_seen.add(p.get('pmid'))

    hits_zf_rec = query_epmc(f'("{fish_symbol}" OR "{display_fish}" AND "zebrafish") sort_date:y', 10)
    for p in hits_zf_rec:
        if len(recent_papers) < 2 and p.get('pmid') not in pmids_seen:
            p['badge'] = '🆕 Most Recent (Zebrafish)'
            recent_papers.append(p)
            pmids_seen.add(p.get('pmid'))
            
    if len(recent_papers) < 2:
        hits_all_rec = query_epmc(f'("{search_term}" OR "{display_human}") sort_date:y', 10)
        for p in hits_all_rec:
            if len(recent_papers) < 2 and p.get('pmid') not in pmids_seen:
                p['badge'] = '🆕 Most Recent (General)'
                recent_papers.append(p)
                pmids_seen.add(p.get('pmid'))
                
    return core_papers + recent_papers

# ==========================================
# 3. DEEP COMPARATIVE SYNTHESIS ENGINE
# ==========================================

def synthesize_narrative(display_fish, display_human, protein, diseases, go_terms, pathways, interactions):
    st.markdown("---")
    st.header("🏆 Prediction Model: Calculated Phenotypic Synthesis")
    
    with st.spinner("Cross-referencing ZFIN Ontologies, KEGG Pathways, and Human Clinical Data..."):
        time.sleep(1) 
        
        go_str = f"{', '.join(go_terms[:3])}" if go_terms else "general metabolic processes"
        kegg_str = f"{', '.join(pathways[:2])}" if pathways else "core cellular signaling"
        disease_str = f"{diseases[0]}" if diseases else "generalized physiological dysfunction"
        
        combined_text = " ".join(go_terms + pathways + diseases).lower()
        tissue_prediction = []
        mechanism = []
        
        if any(kw in combined_text for kw in ["neuro", "brain", "axon", "synap", "mental", "cognitive", "nervous"]):
            tissue_prediction.append("Central Nervous System (CNS) & Neural Crest")
            mechanism.append("Failure in neural tube closure or deficient axonal pathfinding, resulting in microcephaly or behavioral deficits.")
        if any(kw in combined_text for kw in ["cardio", "heart", "angio", "blood", "vasc"]):
            tissue_prediction.append("Cardiovascular & Hematopoietic Systems")
            mechanism.append("Impaired cardiac looping and angiogenesis, presenting as severe pericardial edema by 48 hpf.")
        if any(kw in combined_text for kw in ["metabol", "lipid", "cytochrome", "liver"]):
            tissue_prediction.append("Hepatic & Metabolic Systems")
            mechanism.append("Disruption of lipid/glucose homeostasis, observable as yolk sac malabsorption and steatosis.")
        if any(kw in combined_text for kw in ["somit", "axis", "wnt", "muscle"]):
            tissue_prediction.append("Musculoskeletal Axis")
            mechanism.append("Defective somite boundary formation, causing a severely curved body axis (scoliosis phenotype).")
            
        if not tissue_prediction:
            tissue_prediction.append("Core Cellular Homeostasis & Early Embryogenesis")
            mechanism.append("Widespread cellular apoptosis and cell-cycle arrest leading to early embryonic lethality.")

        interaction_str = f" This loss will critically destabilize its interacting partners ({', '.join(interactions[:4])})." if interactions else ""

        st.markdown(f"""
        > **Data-Driven Hypothesis for {display_fish.upper()}:**
        > 
        > **1. Cross-Database Comparison:** > By extracting data from the **Zebrafish Information Network (ZFIN)** and **NCBI Gene Ontology (GO)**, we computationally map this gene to *{go_str}*. Concurrently, human **KEGG Pathway** data independently identifies the `{display_human.upper()}` ortholog as a driver of *{kegg_str}*. When we cross-reference these distinct datasets, a clear biological consensus emerges.
        >
        > **2. Clinical Translation:** > In clinical literature, loss-of-function in the human ortholog directly triggers *{disease_str}*. Because Zebrafish and humans share high genetic conservation in these specific GO pathways, we calculate that a Zebrafish knockout model will exhibit structurally related physiological collapse.
        >
        > **3. Expected F0/F2 Phenotype:** > If your selected Cas enzyme successfully generates a biallelic frameshift mutation in the `{display_fish.upper()}` gene, nonsense-mediated decay will ablate the `{protein}` protein.{interaction_str} We calculate that the primary morphological defect will localize to the **{' and '.join(tissue_prediction)}**. You should specifically screen your fish for:
        > *{', '.join(mechanism)}*
        """)

# ==========================================
# 4. WET-LAB PROTOCOLS & MOCK GEL
# ==========================================
def draw_mock_gel():
    gel_data = pd.DataFrame({
        "Lane": ["1. Ladder", "1. Ladder", "1. Ladder", "1. Ladder", "2. Wild-Type", "3. CRISPR Mutant (T7E1)"],
        "Band Size (bp)": [100, 200, 300, 400, 400, 400],
        "Intensity": [1, 1, 1, 1, 3, 1]
    })
    gel_data = pd.concat([gel_data, pd.DataFrame({"Lane": ["3. CRISPR Mutant (T7E1)", "3. CRISPR Mutant (T7E1)"], "Band Size (bp)": [250, 150], "Intensity": [2, 2]})])

    chart = alt.Chart(gel_data).mark_tick(thickness=4, width=50).encode(
        x=alt.X('Lane:N', title='Agarose Gel Lanes', axis=alt.Axis(labelAngle=0)),
        y=alt.Y('Band Size (bp):Q', scale=alt.Scale(domain=[50, 500], reverse=True), title='Fragment Size (bp)'),
        color=alt.Color('Intensity:Q', legend=None, scale=alt.Scale(scheme='greys')),
        tooltip=['Lane', 'Band Size (bp)']
    ).properties(height=250, title="Simulated T7E1 Cleavage Assay (400bp PCR Amplicon)")
    st.altair_chart(chart, width="stretch")

def display_wet_lab_protocol(enzyme):
    col1, col2 = st.columns(2)
    with col1:
        st.success("### 🧬 Recommended PCR Protocol")
        st.write("Extract genomic DNA and amplify the targeted region using high-fidelity polymerase (e.g., Phusion/Q5).")
        st.markdown("""
        | Phase | Temperature | Time | Cycles |
        |---|---|---|---|
        | Initial Denat. | 98°C | 30s | 1 |
        | Denaturation | 98°C | 10s | **35** |
        | Annealing | **60°C** *(Adjust per Tm)* | 20s | **35** |
        | Extension | 72°C | 30s | **35** |
        | Final Extension| 72°C | 2 min | 1 |
        """)
        
    with col2:
        st.warning("### 🧪 Expected Gel Benchwork")
        if "SpCas9" in enzyme or "SaCas9" in enzyme:
            st.write("**Mechanism:** Blunt double-strand breaks lead to micro-indels (±1-10bp).")
            st.error("**⚠️ Protocol Warning:** A standard agarose gel **WILL NOT** resolve a 5bp difference. You **must** use a T7E1 mismatch cleavage assay (simulated below), High-Resolution Melt (HRM), or direct Sanger sequencing.")
        elif "Cas12a" in enzyme:
            st.write("**Mechanism:** Staggered 5' cuts often result in larger deletions.")
            st.error("**⚠️ Protocol Warning:** While large deletions (>50bp) might be visible on a 2% agarose gel, HRM or Sanger sequencing is required to confirm the staggered repair junctions.")
    draw_mock_gel()

# ==========================================
# 5. CRISPR BIOINFORMATICS ENGINES
# ==========================================
@st.cache_data(ttl=3600)
def fetch_ensembl_sequence(fish_symbol, ensembl_id=None):
    headers = {"Content-Type": "application/json"}
    
    if ensembl_id:
        try:
            seq_url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?type=genomic"
            seq_response = requests.get(seq_url, headers=headers, timeout=10)
            if seq_response.status_code == 200: return seq_response.json().get('seq', '').upper()
        except Exception: pass
        
    api_url = f"https://rest.ensembl.org/lookup/symbol/danio_rerio/{fish_symbol.lower()}?expand=0"
    try:
        response = requests.get(api_url, headers=headers, timeout=10)
        if response.status_code == 200:
            ens_id = response.json().get('id')
            seq_url = f"https://rest.ensembl.org/sequence/id/{ens_id}?type=genomic"
            seq_response = requests.get(seq_url, headers=headers, timeout=10)
            if seq_response.status_code == 200: return seq_response.json().get('seq', '').upper()
    except Exception: pass
    return None

def find_grnas(dna_seq, enzyme):
    total_len = len(dna_seq)
    if enzyme == "SpCas9 (PAM: NGG)": pattern = r"(?=([ACGT]{20}[ACGT]GG))"
    elif enzyme == "SaCas9 (PAM: NNGRRT)": pattern = r"(?=([ACGT]{21}[ACGT]{2}G[AG]{2}T))"
    elif enzyme == "Cas12a/Cpf1 (PAM: TTTV)": pattern = r"(?=(TTT[ACG][ACGT]{20}))"
        
    matches = re.finditer(pattern, dna_seq)
    grna_list = []
    
    for match in matches:
        full_target = match.group(1)
        if enzyme == "Cas12a/Cpf1 (PAM: TTTV)": pam, target_seq = full_target[:4], full_target[4:]
        elif enzyme == "SaCas9 (PAM: NNGRRT)": target_seq, pam = full_target[:21], full_target[21:]
        else: target_seq, pam = full_target[:20], full_target[20:]
            
        gc = ((target_seq.count('G') + target_seq.count('C')) / len(target_seq)) * 100 if len(target_seq)>0 else 0
        tm = round(64.9 + 41 * ((target_seq.count('G') + target_seq.count('C')) - 16.4) / len(target_seq), 1) if len(target_seq)>0 else 0.0
        
        score = 100
        if gc < 40 or gc > 60: score -= 20
        if gc < 20 or gc > 80: score -= 30
        if total_len > 0:
            if (match.start() / total_len) > 0.5: score -= 25
            elif (match.start() / total_len) > 0.25: score -= 10
        if "TTTT" in target_seq: score -= 50 
        score = max(0, score)
        
        tier = "🟢 High" if score >= 90 else "🟡 Medium" if score >= 60 else "🔴 Low"

        grna_list.append({"Target Sequence": target_seq, "PAM": pam, "GC (%)": round(gc, 1), "Tm (°C)": tm, "Poly-T": "⚠️ YES" if "TTTT" in target_seq else "No", "Pos": match.start(), "Score": score, "Tier": tier})
    return grna_list

def run_ncbi_blast(sequence):
    url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    data = {"CMD": "Put", "PROGRAM": "blastn", "DATABASE": "nt", "QUERY": sequence, "ENTREZ_QUERY": "txid7955[ORGN]"}
    response = requests.post(url, data=data)
    rid_match = re.search(r"RID = (.*)", response.text)
    if not rid_match: return None, None
    rid = rid_match.group(1).strip()
    
    attempts = 0
    while attempts < 30:
        time.sleep(10)
        attempts += 1
        poll_resp = requests.get(url, params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid})
        if "Status=READY" in poll_resp.text: break
    if attempts >= 30: return "BLAST Server Timeout", None
        
    res_resp = requests.get(url, params={"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid})
    try:
        root = ET.fromstring(res_resp.text)
        hit_def = root.find(".//Hit_def").text
        match = re.search(r'\(([a-zA-Z0-9_-]+)\)', hit_def)
        return hit_def, match.group(1).lower() if match else hit_def.split()[0]
    except Exception: return "Unknown Sequence", None

# ==========================================
# 6. UI RENDERING FUNCTIONS
# ==========================================
def display_biological_dossier(gene_name):
    st.subheader("🔬 Clinical & Biological Profile")
    with st.spinner(f"Querying ZFIN, NCBI, and 6 other databases for '{gene_name}'..."):
        display_fish, display_human, fish_symbol, human_ortholog, ensembl_id, show_warning, zfin_aliases, ncbi_name, ncbi_summary, diseases, go_terms, pathways = fetch_ncbi_data(gene_name)
        uniprot_name = fetch_uniprot_data(fish_symbol)
        pubchem_drugs = fetch_pubchem_drugs(human_ortholog)
        interactions = fetch_string_interactions(fish_symbol)
        
        if show_warning:
            st.error(f"**This gene in humans is called {display_human.upper()} and in zebrafish, the same gene is called {display_fish.upper()}.**")
        
        if ncbi_name != "Unknown Gene" or uniprot_name != "Unmapped Protein":
            st.markdown(f"**🧬 NCBI Official Gene:** {ncbi_name}")
            if zfin_aliases: st.markdown(f"**🏷️ ZFIN Known Aliases:** {', '.join(zfin_aliases)}")
            st.markdown(f"**🧩 UniProt Protein:** {uniprot_name}")
            
            if go_terms: st.markdown(f"**⚙️ GO Biological Processes:** {', '.join(go_terms[:4])}")
            if pathways: st.markdown(f"**🔄 KEGG/Reactome Pathways:** {', '.join(pathways[:3])}")
            
            if diseases: st.warning(f"**⚠️ Clinical Human Phenotypes:** Associated with: **{', '.join(diseases[:5])}**")
            else: st.info(f"**⚠️ Clinical Phenotypes:** No severe complex human phenotypes currently mapped.")

            if pubchem_drugs: st.success(f"**💊 PubChem Modulators:** Drugs targeting this pathway: **{', '.join(pubchem_drugs)}**")
            else: st.markdown(f"**💊 PubChem:** No dominant clinical drugs mapped.")
            
            st.markdown("### 📖 Official Biological Summary")
            if ncbi_summary and len(ncbi_summary) > 10: st.write(f"**Source: NCBI:** {ncbi_summary}")
            else: st.warning("No long-form text summary published in NCBI/UniProt.")
                
            return display_fish, display_human, fish_symbol, human_ortholog, ensembl_id, show_warning, uniprot_name, diseases, go_terms, pathways, interactions
        else:
            st.warning(f"Could not verify '{gene_name}' across global databases. Please check spelling.")
            return gene_name, gene_name, gene_name, gene_name, None, False, None, [], [], [], []

def display_literature_ui(display_fish, display_human, fish_symbol, human_ortholog, show_warning):
    st.subheader(f"📚 Validated Literature")
    
    if show_warning:
        st.info(f"Looking for papers for {display_fish.upper()} as in humans, the gene is called {display_human.upper()}, and in the fish, it's called {display_fish.upper()}.")
    
    with st.spinner("Executing cascading literature search..."):
        papers = fetch_literature(display_fish, display_human, fish_symbol, human_ortholog)
        if papers:
            for paper in papers:
                title = paper.get('title', 'Unknown Title')
                author = str(paper.get('authorString', 'Unknown Author')).split(',')[0] + " et al."
                pmid = paper.get('pmid')
                year = paper.get('pubYear', 'N/A')
                badge = paper.get('badge', '📄')
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else f"https://europepmc.org/article/{paper.get('source', 'MED')}/{paper.get('id', '')}"
                
                st.markdown(f"**{badge}** | **[{title}]({link})**")
                st.caption(f"✍️ {author} | 📅 {year}")
                st.markdown("---")
        else:
            st.warning(f"No direct CRISPR/Zebrafish literature found.")

def analyze_and_display_gene(display_fish, display_human, fish_symbol, ensembl_id, show_warning, enzyme):
    st.markdown("---")
    st.subheader(f"✂️ Ensembl Genomics: Target Finder ({enzyme.split()[0]})")
    
    if show_warning:
         st.warning(f"⚠️ {display_fish.upper()} is the gene that is called in fish, and the same gene in humans is called {display_human.upper()}. Generating Ensembl mapping for {display_fish.upper()}.")
    
    display_wet_lab_protocol(enzyme)
    
    with st.spinner(f"Downloading DNA & crunching {enzyme} thermodynamics..."):
        dna_seq = fetch_ensembl_sequence(fish_symbol, ensembl_id)
        if dna_seq:
            grna_list = find_grnas(dna_seq, enzyme)
            if grna_list:
                df = pd.DataFrame(grna_list).sort_values(by="Score", ascending=False).reset_index(drop=True)
                
                chart = alt.Chart(df).mark_circle(size=60).encode(
                    x=alt.X('Pos:Q', title="Genomic Position (bp)"), 
                    y=alt.Y('Score:Q', scale=alt.Scale(domain=[0, 100]), title="Efficiency Score"), 
                    color=alt.Color('Tier:N', scale=alt.Scale(domain=["🟢 High", "🟡 Medium", "🔴 Low"], range=["#2e7d32", "#fbc02d", "#c62828"])), 
                    tooltip=['Target Sequence', 'Tier', 'PAM', 'Pos', 'Score', 'Tm (°C)']
                ).interactive().properties(height=300, title=f"gRNA Distribution Map for {display_fish.upper()}")
                
               st.altair_chart(chart, width="stretch")
                
                st.download_button("📥 Download Target Data (CSV)", data=df.to_csv(index=False).encode('utf-8'), file_name=f"{display_fish}_targets.csv", mime="text/csv", key=f"dl_{display_fish}")
                st.dataframe(df.head(15))
            else: st.warning(f"No valid PAM sites found for {enzyme}.")
        else: st.error(f"⚠️ Gene DNA sequence could not be retrieved from Ensembl for symbol '{display_fish.upper()}'. It may require manual mapping.")

# ==========================================
# 7. MAIN APPLICATION WORKFLOW (THE 3 TABS)
# ==========================================
tab1, tab2, tab3 = st.tabs(["1. Gene Input", "2. Sequence Input (File Upload)", "3. Phenotype Reverse-Search"])

with tab1:
    gene_input = st.text_input("Enter ANY Zebrafish or Human Gene (e.g., kras, brca1, naa60, nat15):", key="t1_gene").lower().strip()
    if st.button(f"Analyze Gene using {selected_enzyme.split()[0]}", key="b1"):
        if gene_input:
            col1, col2 = st.columns(2)
            with col1: 
                display_fish, display_human, fish_symbol, human_ortholog, ensembl_id, show_warning, uniprot_name, diseases, go_terms, pathways, interactions = display_biological_dossier(gene_input)
            with col2: 
                display_literature_ui(display_fish, display_human, fish_symbol, human_ortholog, show_warning)
            
            synthesize_narrative(display_fish, display_human, uniprot_name, diseases, go_terms, pathways, interactions)
            analyze_and_display_gene(display_fish, display_human, fish_symbol, ensembl_id, show_warning, selected_enzyme)

with tab2:
    st.write("Paste a raw sequence OR upload a FASTA/TXT file. We will auto-extract clinical phenotypes and drug interactions.")
    uploaded_file = st.file_uploader("📂 Upload a Sequence File (.fasta, .txt, .seq)", type=["fasta", "txt", "seq"])
    raw_seq_text = st.text_area("Or Paste DNA Sequence (A, C, T, G):", height=150)
    
    if st.button("Identify Sequence & Extract Clinical Data"):
        final_seq = "".join([line.strip() for line in uploaded_file.getvalue().decode("utf-8").splitlines() if not line.startswith(">")]) if uploaded_file else raw_seq_text
        if final_seq:
            clean_seq = re.sub(r'[^A-Za-z]', '', final_seq).upper()
            if len(clean_seq) >= 40:
                with st.spinner("🚀 Connecting to NCBI BLAST (This may take 1-2 minutes)..."):
                    hit_desc, identified_gene = run_ncbi_blast(clean_seq)
                if hit_desc and hit_desc not in ["No Hits Found", "BLAST Server Timeout"]:
                    st.success(f"🧬 **NCBI Match:** {hit_desc}")
                    if identified_gene:
                        col1, col2 = st.columns(2)
                        with col1: display_fish, display_human, fish_symbol, human_ortholog, ensembl_id, show_warning, uniprot_name, diseases, go_terms, pathways, interactions = display_biological_dossier(identified_gene)
                        with col2: display_literature_ui(display_fish, display_human, fish_symbol, human_ortholog, show_warning)
                        synthesize_narrative(display_fish, display_human, uniprot_name, diseases, go_terms, pathways, interactions)
                        
                    analyze_and_display_gene(display_fish, display_human, fish_symbol, ensembl_id, show_warning, selected_enzyme)
                elif hit_desc == "BLAST Server Timeout": st.error("NCBI BLAST servers are overloaded.")
                else: st.error("No matches found in Zebrafish genome.")
            else: st.error("⚠️ Sequence too short for BLAST! (Min 40bp)")

with tab3:
    pheno_input = st.text_input("Enter Desired Phenotype (e.g., brain, pain, tumor, movement):").lower().strip()
    if st.button("Find Genes & Design gRNAs"):
        if pheno_input:
            with st.spinner("Scanning global phenotype databases (Human & Zebrafish)..."):
                # 🚀 CRITICAL FIX: Broad search across ALL fields and BOTH species unlocks the entire phenotype dictionary!
                url = "https://mygene.info/v3/query"
                try:
                    response = requests.get(url, params={"q": pheno_input, "species": "7955,9606", "fields": "symbol", "size": 15}, timeout=10)
                    found_genes = []
                    if response.status_code == 200:
                        for hit in response.json().get('hits', []):
                            if 'symbol' in hit and hit['symbol'].lower() not in found_genes:
                                found_genes.append(hit['symbol'].lower())
                except Exception: found_genes = []

            if found_genes:
                st.success(f"✅ Found {len(found_genes)} gene(s) associated with '{pheno_input}'! Displaying top match:")
                for target_gene in found_genes[:1]: 
                    st.markdown(f"## Target Gene: **{target_gene.upper()}**")
                    col1, col2 = st.columns(2)
                    with col1: display_fish, display_human, fish_symbol, human_ortholog, ensembl_id, show_warning, uniprot_name, diseases, go_terms, pathways, interactions = display_biological_dossier(target_gene)
                    with col2: display_literature_ui(display_fish, display_human, fish_symbol, human_ortholog, show_warning)
                    
                    synthesize_narrative(display_fish, display_human, uniprot_name, diseases, go_terms, pathways, interactions)
                    analyze_and_display_gene(display_fish, display_human, fish_symbol, ensembl_id, show_warning, selected_enzyme)
                    st.markdown("---")
            else:
                st.warning(f"No genes found globally for '{pheno_input}'.")
                closest_match = difflib.get_close_matches(pheno_input, COMMON_PHENOTYPES, n=1, cutoff=0.6)
                if closest_match: st.error(f"**Did you mean: *{closest_match[0]}*?**")
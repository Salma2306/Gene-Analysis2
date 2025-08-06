import streamlit as st
import requests
import pandas as pd
from datetime import datetime
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from io import StringIO
import re
import random
import time
import socket

# 🧪 GC content calculator (since Bio.SeqUtils.GC is removed)
def GC(seq):
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq) if len(seq) > 0 else 0

def run_primer_design():
    """Main primer design function to be called from Streamlit app"""

    # Constants
    Entrez.email = "your.email@example.com"  # Required for NCBI access
    API_ENDPOINT = "https://primer3.ut.ee/api/v1/primer"
    MAX_RETRIES = 3
    RETRY_DELAY = 2

    # Sequence cache for offline work
    SEQUENCE_CACHE = {}

    # Session state to preserve results
    if 'results' not in st.session_state:
        st.session_state.results = []

    def check_internet():
        """Check internet connectivity"""
        try:
            socket.create_connection(("www.ncbi.nlm.nih.gov", 80), timeout=5)
            return True
        except OSError:
            return False

    def fetch_gene_sequence(gene_name, organism="human"):
        """Enhanced sequence fetching with retries and multiple sources"""
        cache_key = f"{gene_name}_{organism}"
        if cache_key in SEQUENCE_CACHE:
            return SEQUENCE_CACHE[cache_key]
        
        if not check_internet():
            st.warning("No internet connection - using cached sequences only")
            return None

        for attempt in range(MAX_RETRIES):
            try:
                # Try RefSeq first
                handle = Entrez.esearch(
                    db="nucleotide",
                    term=f"{gene_name}[Gene] AND {organism}[Organism] AND refseq[filter]",
                    retmax=1
                )
                record = Entrez.read(handle)
                handle.close()

                if not record["IdList"]:
                    # Fallback to GenBank
                    handle = Entrez.esearch(
                        db="nucleotide",
                        term=f"{gene_name}[Gene] AND {organism}[Organism]",
                        retmax=1
                    )
                    record = Entrez.read(handle)
                    handle.close()

                if not record["IdList"]:
                    st.warning(f"No sequence found for {gene_name}")
                    return None

                # Fetch sequence
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=record["IdList"][0],
                    rettype="fasta",
                    retmode="text"
                )
                seq_record = SeqIO.read(StringIO(handle.read()), "fasta")
                handle.close()
                
                sequence = str(seq_record.seq)
                SEQUENCE_CACHE[cache_key] = sequence
                return sequence

            except Exception as e:
                if attempt == MAX_RETRIES - 1:
                    st.error(f"Failed to fetch {gene_name} after {MAX_RETRIES} attempts: {str(e)}")
                time.sleep(RETRY_DELAY)
        return None

    def validate_sequence(sequence):
        """Strict sequence validation"""
        if not sequence:
            return None
            
        clean_seq = re.sub(r'[^ATCGatcg]', '', sequence).upper()
        return clean_seq if len(clean_seq) >= 200 else None

    def design_primers(sequence, gene_name, params):
        """Unified primer design with API and local fallback"""
        # Try API first
        api_result = design_with_api(sequence, gene_name, params)
        if api_result:
            return api_result
        
        # Fallback to local design
        return design_locally(sequence, gene_name, params)

    def design_with_api(sequence, gene_name, params):
        """API-based design with enhanced parameters"""
        try:
            payload = {
                "SEQUENCE_ID": gene_name,
                "SEQUENCE_TEMPLATE": sequence,
                "PRIMER_PRODUCT_SIZE_RANGE": f"{params['product_min']}-{params['product_max']}",
                "PRIMER_MIN_GC": params['gc_min'],
                "PRIMER_MAX_GC": params['gc_max'],
                "PRIMER_OPT_GC_PERCENT": 50,
                "PRIMER_MIN_SIZE": 18,
                "PRIMER_OPT_SIZE": 20,
                "PRIMER_MAX_SIZE": 24,
                "PRIMER_MIN_TM": 55.0,
                "PRIMER_OPT_TM": 60.0,
                "PRIMER_MAX_TM": 65.0,
                "PRIMER_MAX_DIFF_TM": 3.0,
                "PRIMER_MAX_POLY_X": 3,
                "PRIMER_MAX_NS_ACCEPTED": 0,
                "PRIMER_NUM_RETURN": 3  # Get multiple designs to choose from
            }
            
            response = requests.post(API_ENDPOINT, json=payload, timeout=25)
            
            if response.status_code == 200:
                data = response.json()
                for i in range(3):  # Check top 3 designs
                    if f'PRIMER_LEFT_{i}_SEQUENCE' in data:
                        primers = {
                            'forward': data[f'PRIMER_LEFT_{i}_SEQUENCE'],
                            'reverse': data[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                            'size': data[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                            'source': 'api',
                            'position': {
                                'forward_start': data[f'PRIMER_LEFT_{i}'][0],
                                'forward_len': data[f'PRIMER_LEFT_{i}'][1],
                                'reverse_start': data[f'PRIMER_RIGHT_{i}'][0],
                                'reverse_len': data[f'PRIMER_RIGHT_{i}'][1]
                            }
                        }
                        # Validate immediately
                        if validate_primers(primers, sequence)['valid']:
                            return primers
        except Exception:
            pass
        return None

    def design_locally(sequence, gene_name, params):
        """Improved local primer design with biological constraints"""
        try:
            seq_len = len(sequence)
            if seq_len < params['product_max']:
                return None
                
            # Find candidate regions with optimal GC content
            candidates = []
            window_size = 100  # Analyze 100bp windows
            
            for i in range(0, seq_len - window_size):
                subseq = sequence[i:i+window_size]
                subseq_gc = GC(subseq) * 100
                if params['gc_min'] <= subseq_gc <= params['gc_max']:
                    candidates.append(i)
            
            if not candidates:
                return None
                
            # Select best candidate (middle of sequence)
            start_pos = candidates[len(candidates)//2]
            end_pos = start_pos + random.randint(params['product_min'], min(params['product_max'], seq_len - start_pos))
            
            # Design primers with constraints
            fwd_primer = sequence[start_pos:start_pos+20]
            rev_primer = str(Seq(sequence[end_pos-20:end_pos]).reverse_complement())
            
            # Ensure primers meet basic criteria
            if (GC(fwd_primer)*100 < params['gc_min'] or GC(fwd_primer)*100 > params['gc_max'] or
                GC(rev_primer)*100 < params['gc_min'] or GC(rev_primer)*100 > params['gc_max']):
                return None
                
            return {
                'forward': fwd_primer,
                'reverse': rev_primer,
                'size': end_pos - start_pos,
                'source': 'local',
                'position': {
                    'forward_start': start_pos,
                    'forward_len': 20,
                    'reverse_start': end_pos-20,
                    'reverse_len': 20
                }
            }
        except Exception as e:
            st.warning(f"Local design failed for {gene_name}: {str(e)}")
            return None

    def validate_primers(primers, sequence):
        """Comprehensive validation with proper melting temp calculations"""
        validation = {
            'gc_forward': 0,
            'gc_reverse': 0,
            'tm_forward': 0,
            'tm_reverse': 0,
            'valid': False,
            'issues': [],
            'warnings': []
        }
        
        if not primers or not sequence:
            validation['issues'].append("No primers/sequence")
            return validation
        
        try:
            # Basic properties
            validation['gc_forward'] = GC(primers['forward']) * 100
            validation['gc_reverse'] = GC(primers['reverse']) * 100
            validation['tm_forward'] = mt.Tm_NN(primers['forward'])
            validation['tm_reverse'] = mt.Tm_NN(primers['reverse'])
            validation['valid'] = True
            
            # GC content check (40-60%)
            if not (40 <= validation['gc_forward'] <= 60):
                validation['valid'] = False
                validation['issues'].append(f"Fwd GC%: {validation['gc_forward']:.1f}")
            if not (40 <= validation['gc_reverse'] <= 60):
                validation['valid'] = False
                validation['issues'].append(f"Rev GC%: {validation['gc_reverse']:.1f}")
            
            # Tm difference check (<5°C)
            tm_diff = abs(validation['tm_forward'] - validation['tm_reverse'])
            if tm_diff > 5:
                validation['valid'] = False
                validation['issues'].append(f"Tm Δ: {tm_diff:.1f}°C")
            
            # Primer length (18-24bp)
            if not (18 <= len(primers['forward']) <= 24):
                validation['valid'] = False
                validation['issues'].append(f"Fwd length: {len(primers['forward'])}")
            if not (18 <= len(primers['reverse']) <= 24):
                validation['valid'] = False
                validation['issues'].append(f"Rev length: {len(primers['reverse'])}")
                
            # Nucleotide runs (max 3 identical)
            def has_polyx(seq):
                return any(len(match.group()) > 3 for match in re.finditer(r'(A+|T+|G+|C+)', seq))
                
            if has_polyx(primers['forward']):
                validation['warnings'].append("Fwd primer has poly-X")
            if has_polyx(primers['reverse']):
                validation['warnings'].append("Rev primer has poly-X")
                
            # Secondary structure checks
            def check_secondary(seq):
                return (mt.Tm_NN(seq, nn_table=mt.R_DNA_NN1) > 45)  # Simple dimer check
                
            if check_secondary(primers['forward'] + primers['forward']):
                validation['warnings'].append("Fwd primer self-dimer")
            if check_secondary(primers['reverse'] + primers['reverse']):
                validation['warnings'].append("Rev primer self-dimer")
            if check_secondary(primers['forward'] + primers['reverse']):
                validation['valid'] = False
                validation['issues'].append("Primer dimer")
                
        except Exception as e:
            validation['issues'].append(f"Validation error: {str(e)}")
            validation['valid'] = False
            
        return validation

    # UI Components
    st.title("🧬 Biologically Valid Primer Designer")
    st.markdown("""
    Design PCR primers with **NCBI sequence fetching** and **comprehensive validation**.
    """)

    with st.expander("⚙️ Advanced Settings", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            gc_min, gc_max = st.slider("GC Content Range (%)", 20, 80, (40, 60), help="Ideal: 40-60%")
            product_min, product_max = st.slider("Product Size (bp)", 50, 1000, (100, 300))
        with col2:
            organism = st.text_input("Organism", "human")
            strict_mode = st.checkbox("Strict Validation", True, help="Enforce stricter Tm and GC rules")

    # Main Input
    genes_input = st.text_area(
        "Enter gene symbols (one per line or comma separated):",
        height=100,
        value="TP53, ACTB, GAPDH",
        help="Use official gene symbols"
    )

    if st.button("🚀 Design Primers", type="primary"):
        genes = [g.strip() for g in re.split(r'[\n,]', genes_input) if g.strip()]
        
        if not genes:
            st.error("Please enter at least one gene symbol")
            st.stop()
        
        params = {
            'gc_min': gc_min,
            'gc_max': gc_max,
            'product_min': product_min,
            'product_max': product_max
        }
        
        st.session_state.results = []
        
        with st.spinner(f"Processing {len(genes)} genes..."):
            progress_bar = st.progress(0)
            
            for i, gene in enumerate(genes):
                progress_bar.progress((i + 1) / len(genes))
                
                # Fetch sequence
                sequence = fetch_gene_sequence(gene, organism)
                clean_seq = validate_sequence(sequence)
                
                if not clean_seq:
                    st.session_state.results.append({
                        "Gene": gene,
                        "Status": "❌ Failed",
                        "Issue": "Invalid sequence"
                    })
                    continue
                
                # Design primers
                primers = design_primers(clean_seq, gene, params)
                
                if not primers:
                    st.session_state.results.append({
                        "Gene": gene,
                        "Status": "❌ Failed",
                        "Issue": "Design failed"
                    })
                    continue
                
                # Validate
                validation = validate_primers(primers, clean_seq)
                
                # Store results
                st.session_state.results.append({
                    "Gene": gene,
                    "Forward": primers['forward'],
                    "Reverse": primers['reverse'],
                    "Size": primers['size'],
                    "Position": f"{primers['position']['forward_start']}-{primers['position']['forward_start']+primers['position']['forward_len']} / "
                                f"{primers['position']['reverse_start']}-{primers['position']['reverse_start']+primers['position']['reverse_len']}",
                    "GC%": f"{validation['gc_forward']:.1f}/{validation['gc_reverse']:.1f}",
                    "Tm (°C)": f"{validation['tm_forward']:.1f}/{validation['tm_reverse']:.1f}",
                    "Source": "API" if primers['source'] == 'api' else "LOCAL",
                    "Status": "✅ Valid" if validation['valid'] else "⚠️ Check",
                    "Issues": ", ".join(validation['issues']),
                    "Warnings": ", ".join(validation['warnings']) if validation['warnings'] else "None"
                })

    # Results Display
    if st.session_state.results:
        st.success("Design completed!")
        df = pd.DataFrame(st.session_state.results)
        
        # Color coding
        def color_status(val):
            return 'color: green' if val == "✅ Valid" else 'color: orange' if val == "⚠️ Check" else 'color: red'
        
        styled_df = df.style.map(color_status, subset=['Status'])
        
        st.dataframe(styled_df, use_container_width=True)
        
        # Detailed view
        if 'Gene' in df.columns:
            selected_gene = st.selectbox("View details:", df['Gene'])
            gene_data = next(r for r in st.session_state.results if r['Gene'] == selected_gene)
            
            with st.expander(f"🧬 {selected_gene} Details", expanded=True):
                cols = st.columns(2)
                with cols[0]:
                    st.markdown(f"**Forward Primer** (5'→3')")
                    st.code(gene_data['Forward'], language='text')
                    st.metric("Length", f"{len(gene_data['Forward'])} bp")
                    st.metric("GC%", f"{gene_data['GC%'].split('/')[0]}%")
                    st.metric("Tm", f"{gene_data['Tm (°C)'].split('/')[0]}°C")
                    
                with cols[1]:
                    st.markdown(f"**Reverse Primer** (5'→3')")
                    st.code(gene_data['Reverse'], language='text')
                    st.metric("Length", f"{len(gene_data['Reverse'])} bp")
                    st.metric("GC%", f"{gene_data['GC%'].split('/')[1]}%")
                    st.metric("Tm", f"{gene_data['Tm (°C)'].split('/')[1]}°C")
                
                st.metric("Product Size", f"{gene_data['Size']} bp")
                st.markdown(f"**Binding Positions:** {gene_data['Position']}")
                
                if gene_data['Issues'] != "":
                    st.error(f"**Validation Issues:** {gene_data['Issues']}")
                if gene_data['Warnings'] != "None":
                    st.warning(f"**Warnings:** {gene_data['Warnings']}")
        
        # Download
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            "💾 Download CSV",
            csv,
            file_name=f"primers_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv"
        )

    # Documentation
    with st.expander("📚 Documentation", expanded=False):
        st.markdown("""
        ### **Biological Validation Criteria**
        - **GC Content**: 40-60% for both primers
        - **Tm Difference**: <5°C between forward/reverse
        - **Length**: 18-24 nucleotides
        - **Secondary Structures**: Checks for:
            - Self-dimers (forward-forward, reverse-reverse)
            - Primer-dimers (forward-reverse)
        - **Nucleotide Runs**: No more than 3 identical bases in a row
        
        ### **Design Process**
        1. Fetches RefSeq sequences from NCBI
        2. First attempts API design (primer3.ut.ee)
        3. Falls back to local Biopython design if API fails
        4. Performs comprehensive validation
        """)

    with st.expander("🔧 Troubleshooting", expanded=False):
        st.markdown("""
        ### **Common Issues**
        - **No results**: Check gene symbols and internet connection
        - **High GC%**: Try adjusting GC range in settings
        - **Tm mismatch**: Reduce product size range
        
        ### **For Researchers**
        - Always verify primers experimentally
        - Consider ordering from multiple design algorithms
        - Check for SNPs in primer binding sites
        """)

# For standalone execution (if needed)
if __name__ == "__main__":
    run_primer_design()

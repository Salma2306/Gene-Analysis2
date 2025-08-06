import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
from datetime import datetime
import os
import sys
import zipfile
import io
from pathlib import Path
import base64
import re
from bs4 import BeautifulSoup
import shutil
import time
from gene_annotations import GeneAnnotationTool
from enrich_analysis_enhanced import run_enrichment_analysis
from Design_Primer import run_primer_design

# Configure the app
st.set_page_config(
    page_title="Gene Analysis Toolkit",
    layout="wide",
    page_icon="🧬"
)

# Custom CSS
st.markdown("""
<style>
    /* [Previous CSS styles remain unchanged] */
</style>
""", unsafe_allow_html=True)

def cleanup_temp_files(file_path):
    """Remove temporary files"""
    try:
        if file_path and os.path.exists(file_path):
            os.remove(file_path)
    except Exception as e:
        st.warning(f"Cleanup error: {str(e)}")

def embed_images_in_report(output_dir):
    """Create a version of the report with embedded images"""
    report_path = os.path.join(output_dir, "publication_report.html")
    if not os.path.exists(report_path):
        return None
        
    with open(report_path, 'r', encoding='utf-8') as f:
        html_content = f.read()
    
    soup = BeautifulSoup(html_content, 'html.parser')
    
    for img in soup.find_all('img'):
        if img['src'].endswith('.png'):
            img_path = os.path.join(output_dir, img['src'])
            if os.path.exists(img_path):
                with open(img_path, 'rb') as img_file:
                    img_base64 = base64.b64encode(img_file.read()).decode('utf-8')
                img['src'] = f"data:image/png;base64,{img_base64}"
    
    modified_html = str(soup)
    modified_path = os.path.join(output_dir, "report_with_embedded_images.html")
    with open(modified_path, 'w', encoding='utf-8') as f:
        f.write(modified_html)
    
    return modified_path

def create_zip_with_embedded_report(output_dir):
    """Create ZIP with embedded images report"""
    modified_report = embed_images_in_report(output_dir)
    
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(output_dir):
            for file in files:
                if file != "report_with_embedded_images.html":
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, output_dir)
                    zipf.write(file_path, arcname)
        
        if modified_report:
            zipf.write(modified_report, "report_with_embedded_images.html")
    
    zip_buffer.seek(0)
    return zip_buffer.getvalue()

def display_enrichment_results(output_dir):
    """Display all enrichment results with plots"""
    db_mapping = {
        'KEGG': 'KEGG_2021_Human',
        'Reactome': 'Reactome_2022',
        'GO_Biological_Process': 'GO_Biological_Process_2023',
        'GO_Molecular_Function': 'GO_Molecular_Function_2023',
        'GO_Cellular_Component': 'GO_Cellular_Component_2023'
    }
    
    st.subheader("Enrichment Plots")
    cols = st.columns(2)
    col_idx = 0
    
    for db_name in db_mapping.keys():
        plot_path = os.path.join(output_dir, f"{db_name}_plot.png")
        if os.path.exists(plot_path):
            with cols[col_idx % 2]:
                st.image(plot_path, 
                        use_column_width=True,
                        caption=f"{db_name.replace('_', ' ')} Enrichment")
            col_idx += 1
    
    report_path = os.path.join(output_dir, "publication_report.html")
    if os.path.exists(report_path):
        st.subheader("Analysis Report")
        with open(report_path, 'r', encoding='utf-8') as f:
            st.components.v1.html(f.read(), height=800, scrolling=True)

def safe_primer_display_and_download():
    """Safe display and download of primer results"""
    if 'results' not in st.session_state:
        return
        
    try:
        results = st.session_state.results
        if not results or not isinstance(results, list):
            return
            
        # Create DataFrame with error handling
        df = pd.DataFrame(results)
        
        # Display results if valid
        if not df.empty and 'Gene' in df.columns:
            st.dataframe(df)
            
            # Safe download button
            if {'Gene', 'Forward', 'Reverse'}.issubset(df.columns):
                csv = df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    "💾 Download Primer Results",
                    csv,
                    file_name=f"primers_{datetime.now().strftime('%Y%m%d')}.csv",
                    mime="text/csv"
                )
            else:
                st.warning("Primer data incomplete for download")
                
        # Handle individual primer display safely
        if 'Gene' in df.columns:
            selected_gene = st.selectbox("View details:", df['Gene'].unique())
            gene_data = next((r for r in results if r['Gene'] == selected_gene), None)
            
            if gene_data:
                with st.expander(f"🧬 {selected_gene} Details", expanded=True):
                    cols = st.columns(2)
                    with cols[0]:
                        st.markdown("**Forward Primer** (5'→3')")
                        if 'Forward' in gene_data:
                            st.code(gene_data['Forward'], language='text')
                        else:
                            st.warning("Forward primer not available")
                    
                    with cols[1]:
                        st.markdown("**Reverse Primer** (5'→3')")
                        if 'Reverse' in gene_data:
                            st.code(gene_data['Reverse'], language='text')
                        else:
                            st.warning("Reverse primer not available")
    except Exception as e:
        st.error(f"Error displaying primer results: {str(e)}")

# Sidebar navigation
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/2772/2772128.png", 
             use_column_width=True)
    st.title("Gene Analysis Toolkit")
    selected = option_menu(
        menu_title=None,
        options=["Home", "Gene Annotation", "Enrichment Analysis", "Primer Design"],
        icons=["house", "book", "bar-chart", "scissors"],
        default_index=0
    )

# Home Page
if selected == "Home":
    st.header("Welcome to the Gene Analysis Toolkit")
    st.markdown("""
    This web server provides three powerful tools for genetic research:
    
    1. **Gene Annotation** - Comprehensive gene information
    2. **Enrichment Analysis** - Pathway and functional analysis
    3. **Primer Design** - Biologically validated PCR primers
    """)
    st.image("https://www.genome.gov/sites/default/files/tg/en/illustration/dna_sequence.jpg", 
             use_column_width=True)

# Gene Annotation Tool
elif selected == "Gene Annotation":
    st.header("🧬 Gene Annotation Tool")
    
    with st.expander("⚙️ Settings", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            max_variants = st.number_input("Max variants per gene", 1, 50, 10)
            max_trials = st.number_input("Max clinical trials", 1, 20, 5)
        with col2:
            max_literature = st.number_input("Max literature", 1, 20, 5)
            threads = st.number_input("Threads", 1, 10, 3)
    
    genes_input = st.text_area(
        "Enter gene symbols (comma-separated):",
        value="BRCA1, TP53, MBP",
        height=100
    )
    
    if st.button("Run Annotation", type="primary"):
        if not genes_input.strip():
            st.error("Please enter gene symbols")
        else:
            genes = [g.strip().upper() for g in genes_input.split(",") if g.strip()]
            
            with st.spinner(f"Processing {len(genes)} genes..."):
                try:
                    tool = GeneAnnotationTool()
                    config = {
                        "max_variants": max_variants,
                        "max_trials": max_trials,
                        "max_literature": max_literature,
                        "MAX_WORKERS": threads
                    }
                    results = tool.process_genes(genes, config)
                    
                    if not results:
                        st.error("No results obtained")
                        st.stop()
                    
                    st.success("Analysis complete!")
                    
                    tab_names = ["Gene Info", "Variants", "Clinical Trials", "Literature"]
                    tabs = st.tabs(tab_names)
                    
                    excel_buffer = io.BytesIO()
                    with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
                        for i, tab_name in enumerate(tab_names):
                            with tabs[i]:
                                if tab_name in results and results[tab_name]:
                                    df = pd.DataFrame(results[tab_name])
                                    st.dataframe(df)
                                    df.to_excel(writer, sheet_name=tab_name[:30], index=False)
                                else:
                                    st.warning(f"No {tab_name} data available")
                                    pd.DataFrame().to_excel(writer, sheet_name=tab_name[:30])
                    
                    st.download_button(
                        "💾 Download Excel (Multi-Sheet)",
                        data=excel_buffer.getvalue(),
                        file_name=f"gene_annotations_{datetime.now().strftime('%Y%m%d')}.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                
                except Exception as e:
                    st.error(f"Error: {str(e)}")

# Enrichment Analysis Tool
elif selected == "Enrichment Analysis":
    st.header("📊 Gene Enrichment Analysis")
    
    upload_option = st.radio(
        "Input method:",
        ["Upload file", "Paste gene list"],
        horizontal=True
    )
    
    file_path = None
    if upload_option == "Upload file":
        uploaded_file = st.file_uploader("Choose a file (TXT/CSV)", type=["txt", "csv"])
        if uploaded_file:
            file_path = f"temp_{uploaded_file.name}"
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
    else:
        gene_list = st.text_area("Enter genes (one per line or comma-separated):", height=150)
        if gene_list:
            file_path = "temp_gene_list.txt"
            with open(file_path, "w") as f:
                for gene in re.split(r'[\n,]', gene_list):
                    if gene.strip():
                        f.write(gene.strip() + "\n")
    
    if st.button("Run Analysis", type="primary") and file_path:
        with st.spinner("Analyzing..."):
            try:
                result = run_enrichment_analysis(file_path)
                if not result:
                    st.error("Analysis failed")
                    st.stop()
                
                output_dir = result['output_dir']
                st.success("Analysis complete!")
                
                display_enrichment_results(output_dir)
                
                st.download_button(
                    "📦 Download All Results (with embedded plots)",
                    data=create_zip_with_embedded_report(output_dir),
                    file_name=f"enrichment_results_{datetime.now().strftime('%Y%m%d')}.zip",
                    mime="application/zip"
                )
                
            except Exception as e:
                st.error(f"Error: {str(e)}")
            finally:
                if file_path:
                    cleanup_temp_files(file_path)

# Primer Design Tool
elif selected == "Primer Design":
    st.header("✂️ PCR Primer Design")
    
    # Run primer design
    run_primer_design()
    
    # Add safe display and download
    st.markdown("---")
    safe_primer_display_and_download()

# Footer
st.markdown("---")
st.markdown("""
**Gene Analysis Toolkit v2.1**  
*Initial issues fixed including primer KeyError*  
[Report issues](mailto:Salmaloukman37@gmail.com)

""")




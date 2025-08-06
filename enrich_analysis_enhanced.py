import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
import requests
import os
from datetime import datetime
import warnings
import zipfile
import io
import json
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import time

# Configure visual settings
warnings.filterwarnings('ignore')
sns.set_style("white")
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = 10
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9

class EnrichmentAnalyzer:
    def __init__(self):
        self.ENRICHR_URL = "https://maayanlab.cloud/Enrichr"
        self.MAX_WORKERS = 3
        self.REQUEST_TIMEOUT = 30
        self.RETRY_DELAY = 2
        self.MAX_RETRIES = 3
        
        # Configure databases to analyze
        self.DATABASES = {
            'KEGG': 'KEGG_2021_Human',
            'Reactome': 'Reactome_2022',
            'GO_Biological_Process': 'GO_Biological_Process_2023',
            'GO_Molecular_Function': 'GO_Molecular_Function_2023',
            'GO_Cellular_Component': 'GO_Cellular_Component_2023'
        }
    
    def setup_output(self):
        """Create organized output directory with timestamp"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(f"enrichment_results_{timestamp}")
        output_dir.mkdir(exist_ok=True)
        return output_dir
    
    def load_genes(self, input_source):
        """
        Load gene list from either file path or direct input
        Supports both text and CSV formats
        """
        try:
            # Handle file path input
            if isinstance(input_source, (str, Path)) and os.path.exists(input_source):
                if str(input_source).endswith('.txt'):
                    with open(input_source) as f:
                        genes = [line.strip().upper() for line in f if line.strip()]
                elif str(input_source).endswith('.csv'):
                    df = pd.read_csv(input_source, header=None)
                    genes = df[0].str.upper().dropna().tolist()
                else:
                    raise ValueError("Unsupported file format")
            # Handle direct gene list input
            elif isinstance(input_source, str):
                genes = [g.strip().upper() for g in input_source.split(',') if g.strip()]
            elif isinstance(input_source, list):
                genes = [g.strip().upper() for g in input_source if g.strip()]
            else:
                raise ValueError("Invalid input source type")
            
            # Validate genes through HGNC
            valid_genes = self.validate_genes(genes)
            if not valid_genes:
                print("ERROR: No valid HGNC-approved genes found", file=sys.stderr)
                return None
                
            print(f"\nLoaded {len(valid_genes)} valid genes", file=sys.stderr)
            return valid_genes
            
        except Exception as e:
            print(f"Error loading genes: {e}", file=sys.stderr)
            return None
    
    def validate_genes(self, genes):
        """Check genes against HGNC database with retry logic"""
        try:
            url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&status=Approved&format=text"
            
            for attempt in range(self.MAX_RETRIES):
                try:
                    response = requests.get(url, timeout=self.REQUEST_TIMEOUT)
                    if response.status_code == 200:
                        approved_genes = set(line.strip().upper() for line in response.text.split('\n') if line.strip())
                        return [g for g in genes if g in approved_genes]
                    time.sleep(self.RETRY_DELAY)
                except requests.exceptions.RequestException:
                    if attempt == self.MAX_RETRIES - 1:
                        print(f"Warning: Couldn't validate genes with HGNC after {self.MAX_RETRIES} attempts", file=sys.stderr)
                        return list(set(genes))  # Deduplicate but can't validate
                    time.sleep(self.RETRY_DELAY)
        
        except Exception as e:
            print(f"Gene validation error: {str(e)}", file=sys.stderr)
            return list(set(genes))  # Fallback to deduplication only
    
    def submit_to_enrichr(self, genes, description="Gene list"):
        """Submit gene list to Enrichr with error handling"""
        submit_url = f"{self.ENRICHR_URL}/addList"
        payload = {
            'list': (None, '\n'.join(genes)),
            'description': (None, description)
        }
        
        for attempt in range(self.MAX_RETRIES):
            try:
                response = requests.post(submit_url, files=payload, timeout=self.REQUEST_TIMEOUT)
                if not response.ok:
                    print(f"Enrichr submission failed (HTTP {response.status_code})", file=sys.stderr)
                    return None
                
                data = response.json()
                if 'userListId' not in data:
                    print("Unexpected response from Enrichr", file=sys.stderr)
                    return None
                
                return data['userListId']
            
            except Exception as e:
                if attempt == self.MAX_RETRIES - 1:
                    print(f"Failed to submit to Enrichr after {self.MAX_RETRIES} attempts: {str(e)}", file=sys.stderr)
                time.sleep(self.RETRY_DELAY)
        
        return None
    
    def query_database(self, list_id, db_name, db_key):
        """Query a single Enrichr database with retry logic"""
        enrich_url = f"{self.ENRICHR_URL}/enrich"
        params = {
            'userListId': list_id,
            'backgroundType': db_key
        }
        
        for attempt in range(self.MAX_RETRIES):
            try:
                response = requests.get(enrich_url, params=params, timeout=self.REQUEST_TIMEOUT)
                
                if response.ok:
                    data = response.json()
                    if db_key in data and data[db_key]:
                        df = pd.DataFrame(data[db_key])
                        df.columns = [
                            'Rank', 'Term', 'P-value', 'Z-score', 'Combined Score',
                            'Genes', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value'
                        ]
                        df['Database'] = db_name
                        return df
                    return pd.DataFrame()  # Return empty DataFrame if no results
                
                print(f"Query failed for {db_name} (HTTP {response.status_code})", file=sys.stderr)
                time.sleep(self.RETRY_DELAY)
            
            except Exception as e:
                if attempt == self.MAX_RETRIES - 1:
                    print(f"Failed to query {db_name} after {self.MAX_RETRIES} attempts: {str(e)}", file=sys.stderr)
                time.sleep(self.RETRY_DELAY)
        
        return pd.DataFrame()  # Return empty DataFrame on failure
    
    def run_enrichment_analysis(self, genes):
        """Run parallel enrichment analysis across all databases"""
        list_id = self.submit_to_enrichr(genes)
        if not list_id:
            return {}
        
        results = {}
        
        with ThreadPoolExecutor(max_workers=self.MAX_WORKERS) as executor:
            futures = {
                executor.submit(self.query_database, list_id, name, key): name
                for name, key in self.DATABASES.items()
            }
            
            for future in tqdm(futures, total=len(self.DATABASES), desc="Querying databases"):
                db_name = futures[future]
                try:
                    df = future.result()
                    if not df.empty:
                        results[db_name] = df
                except Exception as e:
                    print(f"Error processing {db_name}: {str(e)}", file=sys.stderr)
        
        return results
    
    def create_plot(self, df, title, output_dir, db_name):
        """Generate publication-quality enrichment plot - modified to always show top terms"""
        try:
            # Process data - show top 15 terms regardless of significance
            plot_df = df.sort_values('Adjusted P-value').head(15)
            if plot_df.empty:
                return None
            
            plot_df['-log10(padj)'] = -np.log10(plot_df['Adjusted P-value'])
            plot_df['Term'] = (
                plot_df['Term']
                .str.replace('_', ' ')
                .str.replace('GO:\d+\)', '', regex=True)
                .str[:50]
            )
            
            # Create dot plot
            plt.figure(figsize=(8, 5))
            scatter = sns.scatterplot(
                data=plot_df, 
                x='-log10(padj)',
                y='Term',
                size='Combined Score',
                hue='-log10(padj)',
                palette='viridis_r',
                sizes=(50, 200),
                linewidth=0.5,
                edgecolor='gray'
            )
            
            # Add significance threshold line
            plt.axvline(-np.log10(0.1), color='red', linestyle='--', linewidth=1)
            
            plt.title(f"{db_name}: {title}", pad=20)
            plt.xlabel("-log10(Adjusted P-value)", labelpad=10)
            plt.ylabel("")
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Metrics')
            plt.tight_layout()
            
            plot_path = output_dir / f"{db_name}_plot.png"
            plt.savefig(plot_path, dpi=600, bbox_inches='tight')
            plt.close()
            
            return plot_path
        
        except Exception as e:
            print(f"Error creating plot for {db_name}: {str(e)}", file=sys.stderr)
            return None
    
    def save_results(self, results, output_dir):
        """Save results to CSV files in the output directory"""
        for db_name, df in results.items():
            csv_path = output_dir / f"Enrichr_{db_name}.csv"
            df.to_csv(csv_path, index=False)
    
    def generate_report(self, results, output_dir, gene_count):
        """Generate comprehensive HTML report with plots - modified to always show all databases"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Create plots first - now creates plots for all databases regardless of significance
        plot_paths = {}
        for db_name, df in results.items():
            plot_path = self.create_plot(
                df,
                f"Top {min(15, len(df))} Terms",
                output_dir,
                db_name
            )
            if plot_path:
                plot_paths[db_name] = plot_path.name
        
        # Generate HTML
        html = f"""
        <html>
        <head>
            <title>Gene Enrichment Analysis Report</title>
            <style>
                body {{ font-family: 'Arial', sans-serif; margin: 40px; line-height: 1.6; }}
                h1 {{ color: #2c3e50; font-size: 24px; border-bottom: 1px solid #eee; padding-bottom: 10px; }}
                h2 {{ color: #3498db; font-size: 20px; margin-top: 30px; }}
                .section {{ margin-bottom: 40px; }}
                table {{ border-collapse: collapse; width: 100%; margin: 20px 0; font-size: 14px; }}
                th {{ background-color: #f8f9fa; text-align: left; padding: 10px; border-bottom: 2px solid #ddd; }}
                td {{ padding: 10px; border-bottom: 1px solid #eee; }}
                .sig {{ color: #d62728; font-weight: bold; }}
                .plot-container {{ text-align: center; margin: 30px 0; }}
                .plot {{ max-width: 800px; width: 100%; height: auto; }}
                .summary {{ background-color: #f8f9fa; padding: 20px; border-radius: 5px; margin-bottom: 30px; }}
                .methods {{ margin-top: 40px; font-size: 13px; color: #555; }}
                .gene-list {{ font-family: monospace; background-color: #f5f5f5; padding: 5px; }}
            </style>
        </head>
        <body>
            <h1>Gene Enrichment Analysis Report</h1>
            <p>Generated on {timestamp}</p>
            
            <div class="summary">
                <h2>Analysis Summary</h2>
                <p>This report presents comprehensive enrichment analysis results for {gene_count} genes, 
                including pathway (KEGG, Reactome) and Gene Ontology (Biological Process, Molecular Function, 
                Cellular Component) analyses. All results are adjusted for multiple testing using the Benjamini-Hochberg method.</p>
                <p><strong>Significant findings:</strong> {sum(len(df[df['Adjusted P-value'] < 0.1]) for df in results.values())} significant terms found across all databases (FDR &lt; 0.1)</p>
            </div>
        """
        
        # Add results sections for ALL databases, not just significant ones
        for db_name, df in results.items():
            sig_df = df[df['Adjusted P-value'] < 0.1]
            html += f"""
            <div class="section">
                <h2>{db_name.replace('_', ' ')}</h2>
                <p>Found {len(sig_df)} significantly enriched terms (adjusted p-value &lt; 0.1) out of {len(df)} total terms</p>
                
                <div class="plot-container">
                    <img src="{plot_paths[db_name]}" class="plot" alt="{db_name} enrichment plot">
                </div>
                
                <table>
                    <tr>
                        <th>Term</th>
                        <th>P-value</th>
                        <th>Adj. P</th>
                        <th>Genes</th>
                    </tr>
            """
            
            for _, row in df.head(10).iterrows():  # Show top 10 terms regardless of significance
                genes = ', '.join(eval(row['Genes'])) if isinstance(row['Genes'], str) else str(row['Genes'])
                sig_class = "sig" if row['Adjusted P-value'] < 0.1 else ""
                html += f"""
                    <tr>
                        <td>{row['Term']}</td>
                        <td>{row['P-value']:.2e}</td>
                        <td class="{sig_class}">{row['Adjusted P-value']:.2e}</td>
                        <td class="gene-list">{genes}</td>
                    </tr>
                """
            
            html += """
                </table>
            </div>
            """
        
        # Add methods section
        html += """
            <div class="methods">
                <h2>Methods</h2>
                <p><strong>Gene Set Enrichment Analysis:</strong> Enrichment was performed using the Enrichr API 
                (Kuleshov et al., 2016) with the following databases: KEGG 2021, Reactome 2022, and Gene Ontology 
                (2023). Statistical significance was determined using Fisher's exact test followed by Benjamini-Hochberg 
                correction for multiple testing. Only terms with adjusted p-values &lt; 0.1 were considered significant.</p>
                <p><strong>Visualization:</strong> Dot plots show the -log10 transformed adjusted p-values (x-axis) 
                versus enriched terms (y-axis), with point size representing the combined enrichment score and color 
                indicating statistical significance. The red dashed line indicates the significance threshold (p=0.1).</p>
                <p><strong>Reference:</strong> Kuleshov MV, et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. 
                Nucleic Acids Res. 2016;44(W1):W90-7. doi:10.1093/nar/gkw377</p>
            </div>
        </body>
        </html>
        """
        
        report_path = output_dir / "publication_report.html"
        with open(report_path, "w") as f:
            f.write(html)
        
        return report_path
    
    def create_zip_archive(self, output_dir):
        """Create a ZIP archive of all results"""
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            for file in output_dir.glob('*'):
                if file.is_file():
                    zip_file.write(file, file.name)
        zip_buffer.seek(0)
        return zip_buffer
    
    def run_analysis(self, input_source, output_dir=None):
        """Main analysis workflow"""
        if not output_dir:
            output_dir = self.setup_output()
        
        # Load and validate genes
        genes = self.load_genes(input_source)
        if not genes:
            return None
        
        # Run enrichment analysis
        results = self.run_enrichment_analysis(genes)
        if not results:
            print("No enrichment results obtained", file=sys.stderr)
            return None
        
        # Save results and generate report
        self.save_results(results, output_dir)
        report_path = self.generate_report(results, output_dir, len(genes))
        
        # Create ZIP archive
        zip_buffer = self.create_zip_archive(output_dir)
        
        return {
            'output_dir': str(output_dir),
            'report_path': str(report_path),
            'results': {k: v.to_dict('records') for k, v in results.items()},
            'zip_data': zip_buffer.getvalue()
        }

def run_enrichment_analysis(input_source):
    """Main function to run enrichment analysis"""
    analyzer = EnrichmentAnalyzer()
    return analyzer.run_analysis(input_source)

def main():
    """Command-line interface"""
    analyzer = EnrichmentAnalyzer()
    
    # Get input file
    input_file = input("Enter path to your gene list (txt/csv): ").strip()
    while not os.path.exists(input_file):
        print("File not found. Try again.")
        input_file = input("Enter path to your gene list (txt/csv): ").strip()
    
    # Run analysis
    print("\nRunning enrichment analysis...")
    result = analyzer.run_analysis(input_file)
    
    if result:
        print(f"\nAnalysis complete! Results saved to: {result['output_dir']}")
        print(f"1. Open '{result['report_path']}' in your browser")
        print("2. Includes high-resolution plots suitable for publication")
        print("3. Contains complete methods description for your paper")
    else:
        print("\nAnalysis failed. Please check your input and try again.")

if __name__ == "__main__":
    # Check and install required packages
    try:
        import seaborn
        import matplotlib
    except ImportError:
        import subprocess
        subprocess.check_call(["pip", "install", "seaborn", "matplotlib"])
    
    main()
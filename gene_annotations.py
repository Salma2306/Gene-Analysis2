import requests
import pandas as pd
from io import BytesIO
import time
from pathlib import Path
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
import warnings
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from collections import defaultdict
import openpyxl
from datetime import datetime
import json
import sys
from Bio import Entrez
Entrez.email = "salmaloukman37@gmail.com"  # Replace with your email


# Suppress warnings
warnings.filterwarnings("ignore")

class GeneAnnotationTool:
    def __init__(self):
        # Configuration
        self.CACHE_DIR = Path("gene_data_cache")
        self.CACHE_DIR.mkdir(exist_ok=True)
        self.REQUEST_TIMEOUT = 45
        self.MAX_RETRIES = 3
        self.RETRY_DELAY = 2
        self.MAX_WORKERS = 3
        
        # Configure requests session
        self.session = self._configure_session()
        
    def _configure_session(self):
        """Configure requests session with retry logic"""
        session = requests.Session()
        retries = Retry(
            total=self.MAX_RETRIES,
            backoff_factor=1,
            status_forcelist=[500, 502, 503, 504],
            allowed_methods=["GET"]
        )
        session.mount("https://", HTTPAdapter(max_retries=retries))
        session.mount("http://", HTTPAdapter(max_retries=retries))
        return session
    
    def get_gene_info(self, gene_symbol):
        """Get comprehensive gene information from multiple sources."""
        # Try Ensembl first
        ensembl_data = self._get_ensembl_gene_info(gene_symbol)
        if ensembl_data:
            return ensembl_data
        
        # Fall back to NCBI
        ncbi_data = self._get_ncbi_gene_info(gene_symbol)
        if ncbi_data:
            return ncbi_data
        
        # Final fallback
        return {
            "Gene": gene_symbol,
            "Description": "Not available",
            "Chromosome": "N/A",
            "Genomic Location": "N/A",
            "Gene ID": "N/A",
            "Protein Name": "N/A",
            "Function": "N/A",
            "Source": "No data found",
            "Last Updated": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

    def _get_ensembl_gene_info(self, gene_symbol):
        """Get gene info from Ensembl."""
        try:
            url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
            headers = {"Content-Type": "application/json"}
            response = self.session.get(url, headers=headers, timeout=self.REQUEST_TIMEOUT)
            response.raise_for_status()
            data = response.json()
            
            return {
                "Gene": gene_symbol,
                "Description": data.get("description", "N/A").split("[")[0].strip(),
                "Chromosome": data.get("seq_region_name", "N/A"),
                "Genomic Location": f"{data.get('start', 'N/A'):,}-{data.get('end', 'N/A'):,}",
                "Gene ID": data.get("id", "N/A"),
                "Protein Name": data.get("display_name", "N/A"),
                "Function": "See UniProt for detailed function",
                "Source": "Ensembl",
                "Last Updated": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        except Exception as e:
            print(f"Ensembl failed for {gene_symbol}: {str(e)}", file=sys.stderr)
            return None

    def _get_ncbi_gene_info(self, gene_symbol):
        """Get gene info from NCBI."""
        try:
            # Search for gene ID
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            search_params = {
                "db": "gene",
                "term": f"{gene_symbol}[Gene Name] AND human[Organism]",
                "retmode": "json"
            }
            response = self.session.get(search_url, params=search_params, timeout=self.REQUEST_TIMEOUT)
            response.raise_for_status()
            search_data = response.json()
            
            if not search_data.get("esearchresult", {}).get("idlist"):
                return None
            
            gene_id = search_data["esearchresult"]["idlist"][0]
            
            # Get gene summary
            summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            summary_params = {
                "db": "gene",
                "id": gene_id,
                "retmode": "json"
            }
            response = self.session.get(summary_url, params=summary_params, timeout=self.REQUEST_TIMEOUT)
            response.raise_for_status()
            gene_data = response.json()
            
            return {
                "Gene": gene_symbol,
                "Description": gene_data.get("result", {}).get(gene_id, {}).get("summary", "N/A"),
                "Chromosome": gene_data.get("result", {}).get(gene_id, {}).get("chromosome", "N/A"),
                "Genomic Location": "N/A",
                "Gene ID": gene_id,
                "Protein Name": "N/A",
                "Function": gene_data.get("result", {}).get(gene_id, {}).get("summary", "N/A")[:200] + "..." if gene_data.get("result", {}).get(gene_id, {}).get("summary") else "N/A",
                "Source": "NCBI",
                "Last Updated": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        except Exception as e:
            print(f"NCBI failed for {gene_symbol}: {str(e)}", file=sys.stderr)
            return None

    def get_literature(self, gene_symbol, max_results=5):
        """Get recent literature from PubMed with abstracts and keywords."""
        try:
            # Search for articles
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            search_params = {
                "db": "pubmed",
                "term": f"{gene_symbol}[Title/Abstract] OR {gene_symbol}[MeSH Terms] OR {gene_symbol}[Keyword]",
                "retmax": max_results,
                "sort": "relevance",
                "retmode": "json",
                "email": "salmaloukman37@gmail.com",
                "api_key": "5cd0f5c19a3b1096410f79325b0169db7f08"  

            }
            response = self.session.get(search_url, params=search_params, timeout=self.REQUEST_TIMEOUT)
            response.raise_for_status()
            search_data = response.json()
            
            literature = []
            if search_data.get("esearchresult", {}).get("idlist"):
                pmids = search_data["esearchresult"]["idlist"]
                
                # Get detailed article information
                fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                fetch_params = {
                    "db": "pubmed",
                    "id": ",".join(pmids),
                    "retmode": "xml",
                    "rettype": "abstract"
                }
                response = self.session.get(fetch_url, params=fetch_params, timeout=self.REQUEST_TIMEOUT)
                response.raise_for_status()
                
                # Parse XML response
                from bs4 import BeautifulSoup
                soup = BeautifulSoup(response.text, 'lxml')
                
                for article in soup.find_all('PubmedArticle'):
                    pmid = article.find('PMID').text if article.find('PMID') else "N/A"
                    title = article.find('ArticleTitle').text if article.find('ArticleTitle') else "N/A"
                    
                    # Get abstract
                    abstract = ""
                    abstract_section = article.find('Abstract')
                    if abstract_section:
                        abstract = ' '.join([text.text for text in abstract_section.find_all('AbstractText')])
                    
                    # Get keywords
                    keywords = []
                    keyword_list = article.find('KeywordList')
                    if keyword_list:
                        keywords = [kw.text for kw in keyword_list.find_all('Keyword')]
                    
                    # Get authors
                    authors = []
                    author_list = article.find('AuthorList')
                    if author_list:
                        authors = [f"{auth.find('LastName').text if auth.find('LastName') else ''}, "
                                  f"{auth.find('ForeName').text if auth.find('ForeName') else ''}" 
                                  for auth in author_list.find_all('Author')]
                    first_author = authors[0] if authors else "N/A"
                    
                    # Get journal info
                    journal = article.find('Journal')
                    journal_title = journal.find('Title').text if journal and journal.find('Title') else "N/A"
                    pub_date = journal.find('PubDate') 
                    year = pub_date.find('Year').text if pub_date and pub_date.find('Year') else "N/A"
                    month = pub_date.find('Month').text if pub_date and pub_date.find('Month') else ""
                    day = pub_date.find('Day').text if pub_date and pub_date.find('Day') else ""
                    pub_date_str = f"{year}-{month}-{day}" if (year or month or day) else "N/A"
                    
                    # Get DOI
                    article_id = article.find('ArticleId', IdType="doi")
                    doi = article_id.text if article_id else "N/A"
                    
                    literature.append({
                        "Gene": gene_symbol,
                        "PMID": pmid,
                        "Title": title,
                        "First Author": first_author,
                        "Authors": "; ".join(authors) if authors else "N/A",
                        "Journal": journal_title,
                        "Publication Date": pub_date_str,
                        "Abstract": abstract[:1000] + "..." if len(abstract) > 1000 else abstract or "N/A",
                        "Keywords": "; ".join(keywords) if keywords else "N/A",
                        "DOI": doi,
                        "Source": "PubMed"
                    })
            
            return literature[:max_results]
        except Exception as e:
            print(f"PubMed failed for {gene_symbol}: {str(e)}", file=sys.stderr)
            return []

    def get_variants(self, gene_symbol, max_results=10):
        """Get genetic variants from ClinVar."""
        try:
            # Search for variants
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            search_params = {
                "db": "clinvar",
                "term": f"{gene_symbol}[gene]",
                "retmax": max_results,
                "retmode": "json"
            }
            response = self.session.get(search_url, params=search_params, timeout=self.REQUEST_TIMEOUT)
            response.raise_for_status()
            search_data = response.json()
            
            variants = []
            if search_data.get("esearchresult", {}).get("idlist"):
                variant_ids = search_data["esearchresult"]["idlist"]
                
                # Get variant details
                summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                summary_params = {
                    "db": "clinvar",
                    "id": ",".join(variant_ids),
                    "retmode": "json"
                }
                response = self.session.get(summary_url, params=summary_params, timeout=self.REQUEST_TIMEOUT)
                response.raise_for_status()
                variant_data = response.json()
                
                for variant_id, variant in variant_data.get("result", {}).items():
                    if variant_id == "uids":
                        continue
                    if isinstance(variant, dict):
                        variant_info = {
                            "Gene": gene_symbol,
                            "Variant ID": variant.get("accession", "N/A"),
                            "Source": "ClinVar"
                        }
                        
                        # Add clinical significance if available
                        clin_sig = variant.get("clinical_significance", {}).get("description")
                        if clin_sig:
                            variant_info["Clinical Significance"] = clin_sig
                        
                        # Add condition if available
                        condition = variant.get("trait_set", [{}])[0].get("trait_name")
                        if condition:
                            variant_info["Condition"] = condition
                        
                        # Add review status if available
                        review_status = variant.get("review_status")
                        if review_status:
                            variant_info["Review Status"] = review_status
                        
                        variants.append(variant_info)
            
            return variants[:max_results]
        except Exception as e:
            print(f"ClinVar failed for {gene_symbol}: {str(e)}", file=sys.stderr)
            return []

    def get_clinical_trials(self, gene_symbol, max_results=5):
        """Get clinical trials from ClinicalTrials.gov."""
        try:
            url = "https://clinicaltrials.gov/api/v2/studies"
            params = {
                "query.term": gene_symbol,
                "pageSize": max_results,
                "fields": "NCTId,BriefTitle,Condition,OverallStatus,Phase,StartDate,CompletionDate"
            }
            response = self.session.get(url, params=params, timeout=self.REQUEST_TIMEOUT)
            response.raise_for_status()
            data = response.json()
            
            trials = []
            for study in data.get("studies", [])[:max_results]:
                protocol = study.get("protocolSection", {})
                identification = protocol.get("identificationModule", {})
                status = protocol.get("statusModule", {})
                design = protocol.get("designModule", {})
                
                trials.append({
                    "Gene": gene_symbol,
                    "NCT ID": identification.get("nctId", "N/A"),
                    "Title": identification.get("briefTitle", "N/A"),
                    "Conditions": ", ".join(protocol.get("conditionsModule", {}).get("conditions", ["N/A"])),
                    "Status": status.get("overallStatus", "N/A"),
                    "Phase": design.get("phases", ["N/A"])[0],
                    "Start Date": status.get("startDate", "N/A"),
                    "Completion Date": status.get("completionDate", "N/A"),
                    "Source": "ClinicalTrials.gov"
                })
            
            return trials[:max_results]
        except Exception as e:
            print(f"ClinicalTrials.gov failed for {gene_symbol}: {str(e)}", file=sys.stderr)
            return []

    def process_gene(self, gene_symbol, config):
        """Process a single gene with configurable limits."""
        gene_data = {
            "Gene Info": self.get_gene_info(gene_symbol),
            "Variants": self.get_variants(gene_symbol, config["max_variants"]),
            "Clinical Trials": self.get_clinical_trials(gene_symbol, config["max_trials"]),
            "Literature": self.get_literature(gene_symbol, config["max_literature"])
        }
        return gene_data

    def process_genes(self, genes, config):
        """Process multiple genes in parallel."""
        all_results = defaultdict(list)
        start_time = time.time()
        
        with ThreadPoolExecutor(max_workers=config["MAX_WORKERS"]) as executor:
            futures = {executor.submit(self.process_gene, gene, config): gene for gene in genes}
            
            for future in tqdm(futures, total=len(genes), desc="Processing genes"):
                gene = futures[future]
                try:
                    result = future.result()
                    for key, value in result.items():
                        if value:
                            if isinstance(value, list):
                                all_results[key].extend(value)
                            else:
                                all_results[key].append(value)
                    time.sleep(self.RETRY_DELAY)  # Be gentle with the APIs
                except Exception as e:
                    print(f"Error processing gene {gene}: {str(e)}", file=sys.stderr)
        
        return dict(all_results)

    def save_to_excel(self, data, filename):
        """Save data to Excel with proper formatting."""
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            for sheet_name, sheet_data in data.items():
                if sheet_data:
                    df = pd.DataFrame(sheet_data)
                    
                    # Apply formatting based on sheet type
                    if sheet_name == "Gene Info":
                        df = df[["Gene", "Description", "Chromosome", "Genomic Location", 
                                "Gene ID", "Protein Name", "Function", "Source", "Last Updated"]]
                    elif sheet_name == "Variants":
                        available_cols = ["Gene", "Variant ID", "Source"]
                        if any("Clinical Significance" in d for d in sheet_data):
                            available_cols.append("Clinical Significance")
                        if any("Condition" in d for d in sheet_data):
                            available_cols.append("Condition")
                        if any("Review Status" in d for d in sheet_data):
                            available_cols.append("Review Status")
                        df = df[available_cols]
                    elif sheet_name == "Clinical Trials":
                        df = df[["Gene", "NCT ID", "Title", "Conditions", "Status", 
                                "Phase", "Start Date", "Completion Date", "Source"]]
                    elif sheet_name == "Literature":
                        df = df[["Gene", "PMID", "Title", "First Author", "Authors", 
                                "Journal", "Publication Date", "Abstract", "Keywords", 
                                "DOI", "Source"]]
                    
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
                    
                    # Auto-adjust column widths
                    worksheet = writer.sheets[sheet_name]
                    for column in worksheet.columns:
                        max_length = 0
                        column = [cell for cell in column]
                        for cell in column:
                            try:
                                if len(str(cell.value)) > max_length:
                                    max_length = len(str(cell.value))
                            except:
                                pass
                        adjusted_width = (max_length + 2)
                        worksheet.column_dimensions[column[0].column_letter].width = adjusted_width

    def run_from_web(self, genes, config):
        """Run annotation from web interface and return JSON-serializable results."""
        results = self.process_genes(genes, config)
        
        # Convert DataFrames to dict for JSON serialization
        web_results = {}
        for key, value in results.items():
            if value:
                if isinstance(value, list):
                    web_results[key] = value
                else:
                    web_results[key] = [value]
        
        return web_results

def run_gene_annotations(genes, config):
    """Main function to run gene annotations with config"""
    tool = GeneAnnotationTool()
    return tool.process_genes(genes, config)

def main():
    """Command-line interface."""
    print("Enhanced Gene Annotation Tool")
    print("----------------------------")
    
    # Get gene list from user
    genes_input = input("\nEnter gene symbols (comma-separated, e.g., 'BRCA1,TP53,MBP'): ")
    genes = [g.strip().upper() for g in genes_input.split(",") if g.strip()]
    
    if not genes:
        print("No valid genes provided. Exiting.")
        return
    
    # Default configuration
    config = {
        "max_variants": 10,
        "max_trials": 5,
        "max_literature": 5,
        "MAX_WORKERS": 3
    }
    
    # Process genes
    results = run_gene_annotations(genes, config)
    
    # Save to Excel
    if any(results.values()):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"gene_annotations_{timestamp}.xlsx"
        tool = GeneAnnotationTool()
        tool.save_to_excel(results, filename)
        
        print(f"\nResults saved to '{filename}'")
        print("Contains the following sheets:")
        for sheet in results.keys():
            if results[sheet]:
                print(f"- {sheet}: {len(results[sheet])} records")
    else:
        print("\nNo data could be retrieved for the provided genes.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {str(e)}")

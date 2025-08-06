# Gene Analysis Toolkit

![App Screenshot](https://cdn-icons-png.flaticon.com/512/2772/2772128.png)

A comprehensive web application for genetic research with three powerful tools: Gene Annotation, Enrichment Analysis, and Primer Design.

## 🌟 Key Features

### 🧬 Gene Annotation Tool
- Multi-source aggregation from Ensembl, NCBI, ClinVar, PubMed, and ClinicalTrials.gov
- Configurable limits for variants, trials, and literature
- Parallel processing for faster results
- Excel export with multiple sheets

### 📊 Enrichment Analysis
- API integration with Enrichr
- Interactive visualizations with Matplotlib/Seaborn
- Publication-ready HTML reports
- Embedded images in downloadable ZIP

### ✂️ Primer Design
- Dual-mode operation (API + local fallback)
- Biological validation with 9 quality checks
- NCBI sequence fetching with caching
- Detailed primer metrics (GC%, Tm, etc.)

## 🚀 Deployment Guide

### Recommended: Streamlit Community Cloud
```bash
1. Fork this repository
2. Go to [Streamlit Community Cloud](https://share.streamlit.io/)
3. Click "New app" → "From existing repo"
4. Select your forked repository
5. Set main file path to `app.py`
6. Click "Deploy"

# PVC_predictor: A modular pipeline for the discovery of potential vaccine candidates (PVCs) from the bacterial proteome for a vaccine construct 
![License: MIT License](https://img.shields.io/badge/License-MIT-blue.svg)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![Streamlit](https://img.shields.io/badge/Streamlit-UI-green)
![Contributors](https://img.shields.io/badge/contributors-1-blueviolet)


## Table of Contents
1. [overview](#overview)
2. [Features](#features)
3. [Prerequisites](#prerequisites)
4. [PVC_predictor Installation](#biovix-installation)
5. [Running BioVix](#running-biovix)
6. [Outputs](#outputs)
7. [Tips for Success](#tips-for-success)
8. [Methodolgy](#methodology)
9. [Reference](#reference)
10. [License and Issues](#license-and-issues)
11. [Authors and Contacts](#authors-and-contacts)

## Overview

PVC_predictor is a modular pipeline to discover the potential vaccine candidates (PVCs) from the whole bacterial proteome by subtractive proteomics for the construction of effective vaccine.
       
## Features

- **Host Homologous protein removal**
- **Removal of Paralogous Proteins**
- **Non-Virulent Proteins**
- **Non-Essential Proteins**
- **Allergenicity of Proteins**
- **Antigenicity of Proteins**,
- **Proteins Stability**
- **Proteins Subcellular Localization**
- **Data Visualization and analytics**

## Prerequisites

- [Python 3.8.5](https://www.python.org/downloads/)
- Prerequisite Tools or Software:
  - [ncbi-blast-2.17.0+-win64.exe ](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  - [deeplocpro 1.0](https://github.com/Jaimomar99/deeplocpro)
  - [IApred](https://github.com/sebamiles/IAPred)


## PVC_predictor Installation

### 1. Clone the repository

```bash
git clone https://github.com/MuhammadZain-Butt/PVC_predictor.git
cd PVC_predictor
```

### 2. Create a Virtual Environment (Recommended)
It is highly recommended to use a virtual environment to avoid dependency conflicts.

#### Windows:
```powershell
python -m venv env
env\Scripts\activate
```
Note: Replace `env_name` with your preferred name for the virtual environment.

### 4. Install prerequisite tools 
``` Install the tools by following the instructions from their official GitHub pages.```
Note: deeplocpro and IApred should be installed inside the PVC_predictor folder, or adjust the paths
       in integrated_pipeline.py, and blastp will be installed in its default path.

### 3. Install Dependencies 

```powershell
pip install -r requirements.txt
```

## Running PVC_predictor

After installing dependencies and setting up your environment, you can start BioVix using Streamlit.

```powershell
streamlit run pipeline.py
```
Once the command runs, the app will automatically open in your default browser at: `http://localhost:8501`

## Outputs

This pipeline consists of 8 steps, with the hard-coded data visualization and statistics at every step. At the end of each step, the user is provided with the downloadable FASTA file for the next step. The user can perform any module independently of other steps. After all analyses are complete, a comparative analysis funnel graph is visualized to show the reduction of the proteome or protein sequences into a few potential vaccine candidates. 

## Methodology
- HTTP post/get requests
- Tools Local installation
- Python Traditional Functions
- Files on cloud space 


## Tips for Success

- Ensure that input files are correctly formatted (e.g., fasta format) 
- Set the parameters for blastp, allergenicity, deeplocpro and so on, according to your study.
- Read the instructions for each step carefully and if you don't have required data in some steps links are provided.
- Consider the statistics and visualization after every step, that would be automatically generated.
  

## References
In Process.

## License and Issues

This PVC_predictor is licensed under the MIT License - see the [LICENSE](License) file for details.
Submit issues or contributions via [GitHub Issues](https://github.com/MuhammadZain-Butt/BioVix/issues).

## Authors and Contacts

**Mr. Muhammad Zain Butt**  
*Integrative Omics and Molecular Modeling Laboratory, Department of Bioinformatics and Biotechnology, Government College University Faisalabad (GCUF), Faisalabad, 38000, Pakistan*  
Email: [zain.202302328@gcuf.edu.pk](mailto:zain.202302328@gcuf.edu.pk)








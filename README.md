
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
7. [Deployment](#deployment)
8. [Tips for Success](#tips-for-success)
9. [Reference](#reference)
10. [License and Issues](#license-and-issues)
11. [Authors and Contacts](#authors-and-contacts)

## Overview

PVC_predictor is a modular pipeline to discover the potential vaccine candidates (PVCs) from the whole bacterial proteome by subtractive proteomics for the construction of effective vaccine.


![User Interface](https://github.com/MuhammadZain-Butt/BioVix/blob/main/overview.jpg)

       
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
- Tools or Softwares:
  - [ncbi-blast-2.17.0+-win64.exe ](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  - [deeplocpro 1.0](https://github.com/Jaimomar99/deeplocpro)
  - [IApred](https://github.com/sebamiles/IAPred)
  

## BioVix Installation

### 1. Clone the repository

```bash
git clone https://github.com/MuhammadZain-Butt/BioVix.git
cd BioVix
```

### 2. Create a Virtual Environment (Recommended)
It is highly recommended to use a virtual environment to avoid dependency conflicts.

#### Windows:
```powershell
python -m venv env
env\Scripts\activate
```
Note: Replace `env_name` with your preferred name for the virtual environment.
#### Linux / macOS:
```powershell
python3 -m venv env_name
source env_name/bin/activate
```
Note: Replace `env_name` with your preferred name for the virtual environment.

### 3. Install Dependencies

```powershell
pip install -r requirements.txt
```

### 4. Configure Environment Variables

Create a `.env` file in the project root directory (as in BioVix):

```
DEEPSEEK_API_KEY="your_deepseek_key_here"
GPT_API_KEY="your_gpt_key_here"
QWEN_API_KEY="your_qwen_key_here"
SEMANTIC_SCHOLAR_API_KEY="your_semantic_scholar_key_here"
```
Tip: If you do not have these API keys, you can create accounts here to generate them:
- [OpenRouter](https://openrouter.ai/models) 
- [Semantic Scholar](https://www.semanticscholar.org/product/api)

## Running BioVix

After installing dependencies and setting up your environment, you can start BioVix using Streamlit.

```powershell
streamlit run app.py
```
Once the command runs, the app will automatically open in your default browser at: `http://localhost:8501`

## Outputs

The following panels illustrate the outputs of BioVix across varying datasets.  **(A)** displays the raw input data, while **(B,C,D)** presents the corresponding interactive visualization rendered with Plotly. **(E)** provides the AI-generated interpretation of the graph, along with the derived search query. Finally, **(F)** lists the relevant research papers retrieved from Semantic Scholar using the formulated query: 

- **Figures**:
  1. **Gene-level Protein Expression Dataset**
      
     ![Gene-level Protein Expression](https://github.com/MuhammadZain-Butt/BioVix/blob/main/results/Figure%2004.jpg)
     
  2. **Peak Annotation dataset**
     
     ![Peak Annotation dataset](https://github.com/MuhammadZain-Butt/BioVix/blob/main/results/Figure%2005.jpg)
     
  3. **Clinical Diabetic Dataset**
     
     ![Clinical Diabetic Dataset](https://github.com/MuhammadZain-Butt/BioVix/blob/main/results/Figure06.jpg)



## Deployment

BioVix is deployed on Hugging Face and can be tested or used directly, [click here](https://huggingface.co/spaces/MuhammadZain10/BioVix)

## Tips for Success

- Ensure that input files are correctly formatted (e.g., CSV, XLSX, or TSV) and contain all information required for visualization.
- Write queries in a clear and detailed manner, and avoid using informal language.
- Use consistent naming conventions for columns and variables to improve clarity and interpretation.
  

## References
In Process.

## License and Issues

This PVC_predictor is licensed under the MIT License - see the [LICENSE](License) file for details.
Submit issues or contributions via [GitHub Issues](https://github.com/MuhammadZain-Butt/BioVix/issues).

## Authors and Contacts

**Mr. Muhammad Zain Butt**  
*Integrative Omics and Molecular Modeling Laboratory, Department of Bioinformatics and Biotechnology, Government College University Faisalabad (GCUF), Faisalabad, 38000, Pakistan*  
Email: [zain.202302328@gcuf.edu.pk](mailto:zain.202302328@gcuf.edu.pk)








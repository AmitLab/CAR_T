# Single-cell RNA sequencing (scRNA-seq) of pretreatment blood and CAR T-cell infusion products associated with treatment response in DLBCL

Analysis associated with the *Cancer Research* paper and GEO accession **GSE267097**.

This repository contains analysis code associated with:

Gurevich-Shapiro A, Zwicky P, Winter E, et al. **Distinct Cellular and Molecular Patterns in Pretreatment Peripheral Blood Are Associated with CAR T-cell Outcomes in Diffuse Large B-cell Lymphoma.** *Cancer Research* (2025). [DOI: 10.1158/0008-5472.CAN-25-3596](https://doi.org/10.1158/0008-5472.CAN-25-3596).

The study used single-cell RNA sequencing (scRNA-seq) and single-cell transcriptomics to characterize pretreatment peripheral blood and anti-CD19 CAR T-cell products from 57 patients with diffuse large B-cell lymphoma (DLBCL). The cohort included patients treated with axicabtagene ciloleucel (axi-cel) and tisagenlecleucel (tisa-cel).

The analyses focused on anti-CD19 CAR T-cell therapy, treatment response, treatment resistance, immune biomarkers, gene signatures, transcriptional programs, inflammatory signaling, and cellular features associated with clinical outcome following immunotherapy.

## Data

Raw and processed single-cell RNA-seq data are available through GEO:

- **GSE267097:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267097

The GEO dataset contains scRNA-seq data relevant to single-cell transcriptomics, anti-CD19 CAR-T response and resistance, DLBCL, axi-cel, tisa-cel, immune biomarkers, gene signatures, and inflammatory signaling.

Supplementary material is available through [AACR Figshare](https://aacr.figshare.com/collections/Data_from_Distinct_Cellular_and_Molecular_Patterns_in_Pretreatment_Peripheral_Blood_Are_Associated_with_CAR_T-cell_Outcomes_in_Diffuse_Large_B-cell_Lymphoma/8545291).

## Repository contents

The repository contains R and Python code used for:

- Single-cell transcriptomic analysis of pretreatment peripheral blood
- Analysis of anti-CD19 CAR T-cell products, including axicabtagene ciloleucel (axi-cel) and tisagenlecleucel (tisa-cel)
- Immune-cell composition and transcriptional analyses
- Gene programs and gene signatures
- Inflammatory signaling and pathway analyses
- Identification of immune biomarkers associated with treatment response and treatment resistance
- Machine-learning-based clinical prediction analyses presented in Figure 4

## Research topics

- Single-cell RNA-seq / single-cell transcriptomics
- Anti-CD19 CAR T-cell therapy
- Diffuse large B-cell lymphoma (DLBCL)
- Treatment response and resistance
- Immune biomarkers and gene signatures
- Inflammatory signaling
- Anti-CD19 CAR T-cell products, including axi-cel and tisa-cel

## Clinical prediction

### 1. Clone the repository

```bash
git clone https://github.com/AmitLab/CAR_T.git
cd CAR_T
```

### 2. Create the Python environment

The analyses were performed using **Python 3.9**.

Install the required packages from PyPI using:

```bash
pip install -r requirements
```

### 3. Run the prediction experiments

Run:

`fig4/prediction_notebooks/CAR_T_experiment_rev.ipynb`

For the pre-revision analyses, use:

`fig4/prediction_notebooks/CAR_T_experiment.ipynb`

The experiment notebook trains prediction models using different configurations and generates `.pkl` result files.

To reproduce the analyses presented in the paper:

1. Edit the general configuration to specify the repository directory.
2. Edit the experiment configuration according to the instructions in the notebook.
3. Run the notebook three times using the three configurations specified after the imports.
4. The resulting `.pkl` files will be saved in the `data` directory. **Do not delete these files**, as they are required for the visualization step.

### 4. Generate the visualizations

Run:

`fig4/prediction_notebooks/CAR_T_vis.ipynb`

Edit the repository directory in the general configuration before running the notebook.

The visualization notebook uses the saved experiment results and creates a figure directory labeled with the execution date containing the generated figures.

## Citation

If you use this code or the associated single-cell RNA-seq dataset, please cite:

> Gurevich-Shapiro A, Zwicky P, Winter E, et al. Distinct Cellular and Molecular Patterns in Pretreatment Peripheral Blood Are Associated with CAR T-cell Outcomes in Diffuse Large B-cell Lymphoma. *Cancer Research*. 2025;85(24). [https://doi.org/10.1158/0008-5472.CAN-25-3596](https://doi.org/10.1158/0008-5472.CAN-25-3596)

**Keywords:** single-cell RNA-seq, scRNA-seq, single-cell transcriptomics, anti-CD19 CAR-T, CAR T-cell therapy, axicabtagene ciloleucel, axi-cel, tisagenlecleucel, tisa-cel, diffuse large B-cell lymphoma, DLBCL, cancer immunotherapy, treatment response, treatment resistance, determinants of response, predictive biomarkers, immune biomarkers, patient stratification, precision medicine, personalized medicine, gene signatures, transcriptional programs, inflammatory signaling, pretreatment peripheral blood, CAR T-cell products, B cells, T cells, myeloid cells.

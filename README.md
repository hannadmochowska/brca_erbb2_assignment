## TCGA BRCA – ERBB2 (HER2) Differential Expression

This project contains a single Quarto file, `brca_erbb2.qmd`, which runs a full RNA-seq analysis on TCGA breast cancer (BRCA) data, focusing on ERBB2 (HER2) amplification.

The analysis:

- Loads TCGA BRCA Pan-Can Atlas 2018 data (RNA-seq, clinical, copy-number)
- Defines HER2-amplified vs non-amplified tumours from ERBB2 copy-number
- Performs differential expression (DE) with **DESeq2**
- Runs pathway enrichment with **clusterProfiler** and **ReactomePA**
- Builds an overall survival model using **glmnet** (LASSO Cox)

## 1. Files and folders

- **`brca_erbb2.qmd`** - Main analysis file. Contains all code, figures, and tables.
- **`brca_erbb2.nb.html`** - Rendered HTML report with all results.
- **`data/`**
    - `brca_tcga_pan_can_atlas_2018.tar` – downloaded TCGA Pan-Can Atlas BRCA archive
    - Extracted files inside `data/brca_tcga_pan_can_atlas_2018/`, including:
        - `data_mrna_seq_v2_rsem.txt` - RNA-seq expression matrix
        - `data_clinical_patient.txt` - clinical information
        - `data_cna.txt` - copy-number data

## 2. Required R packages

The code uses the following R packages:

- **DESeq2** - differential expression and variance-stabilising transformation
- **clusterProfiler** - functional enrichment analysis
- **ReactomePA** - Reactome pathway enrichment
- **org.Hs.eg.db** - human gene annotation (for mapping IDs)
- **ggplot2** - general plotting
- **pheatmap** - heatmaps for expression data
- **survival** - survival objects and Cox models
- **glmnet** - LASSO Cox model for survival

All are loaded in the **Setup** section at the top of `brca_erbb2.qmd`.

## 3. What each section of the code does

### 3.1 Setup

- Loads all required packages (listed above).
- Sets basic options, e.g. `stringsAsFactors = FALSE`.
- This is where you would add any global options (e.g. themes).

### 3.2 Locate and load TCGA files

- Assumes the Pan-Can Atlas .tar file is saved as
`data/brca_tcga_pan_can_atlas_2018.tar`.
- Uses `untar()` to unpack it into `data/`.
- Builds file paths for:
    - `rna_file` - RNA-seq data
    - `patient_file` - clinical data
    - `cna_file` - copy-number data
- Reads each of these into R with `read.delim()`:
    - `rna` - gene expression (rows = genes, columns = samples)
    - `patients` - one row per patient with clinical variables
    - `cna` - copy-number calls per gene per sample

### 3.3 Clean clinical data

- Cleans and reshapes the clinical data from `patients`:
    - Standardises patient IDs.
    - Keeps the relevant clinical columns (e.g. survival time/status, basic characteristics).
- Produces a clean `clin` data frame that can be merged with expression and CNA data later.

### 3.4 Define HER2 groups from ERBB2 copy-number

- Extracts ERBB2 copy-number values from the `cna` matrix.
- Aligns copy-number samples to the RNA-seq samples.
- Creates a HER2 group based on ERBB2 CNA:
    - **Amp** - ERBB2 copy-number > 0
    - **NotAmp** - ERBB2 copy-number ≤ 0
- Builds a `coldata` data frame with:
    - `patient_id`
    - `HER2` (factor with levels `NotAmp`, `Amp`)
    - row names = sample IDs
- Drops samples with missing HER2 status.

This `coldata` will be used as the sample metadata for DESeq2.

### 3.5 Build DESeq2 dataset and normalise counts

- Rebuilds a counts matrix `expr_mat` from the RNA-seq data (counts per gene per sample).
- Constructs a **DESeqDataSet**:
    
    ```r
    dds <- DESeqDataSetFromMatrix(
      countData = expr_mat,
      colData   = coldata,
      design    = ~ HER2
    )
    ```
    
- Filters out very lowly expressed genes (rows with very low counts across samples).
- Runs the DESeq2 pipeline:
    
    ```r
    dds <- DESeq(dds)
    ```
    
- Computes a **variance stabilising transformation** (VST) for downstream plots:
    
    ```r
    vst_mat <- assay(vst(dds))
    ```
    

### 3.6 QC and exploratory plots

Typical things in this section:

- **Sample counts / group sizes** (e.g. `table(coldata$HER2)` to show how many `Amp` vs `NotAmp`).
- **PCA plot** using `vst_mat`:
    - Reduces expression data to principal components.
    - Colours points by HER2 group (`Amp` vs `NotAmp`).
    - Checks if HER2 groups separate in expression space.
- **Heatmap** of the most variable genes:
    - Computes per-gene variance in `vst_mat`.
    - Selects the top N most variable genes (e.g. 500).
    - Plots a heatmap with `pheatmap()`, with columns annotated by HER2 status.

This is mainly quality control and a visual check that HER2-amplified tumours look different.

### 3.7 Differential expression results

- Extracts DESeq2 results comparing HER2 **Amp** vs **NotAmp**:
    
    ```r
    res <- results(dds, contrast = c("HER2", "Amp", "NotAmp"))
    ```
    
- Converts this to a clean data frame:
    
    ```r
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    ```
    
- Drops rows with `NA` p-values and orders by adjusted p-value (`padj`).
- Summarises:
    - Number of significantly DE genes (`padj < 0.05`)
    - Top 10 DE genes ranked by absolute log2 fold-change or by `padj`
- Shows these top genes as a table (e.g. with `knitr::kable()`).

This section answers: *“Which genes are differentially expressed between HER2-amplified and non-amplified tumours?”*

### 3.8 Pathway enrichment analysis

- Selects a list of significant genes (e.g. all genes with `padj < 0.05`).
- Converts gene identifiers as needed (e.g. to Entrez IDs) using `org.Hs.eg.db`.
- Runs enrichment using **clusterProfiler** and **ReactomePA**:
    - GO enrichment (e.g. `enrichGO()`)
    - KEGG / Reactome or other pathways (e.g. `enrichPathway()` from ReactomePA)
- Visualises results with enrichment plots (dotplots, barplots, etc.).

This section answers: *“What biological processes and pathways are enriched in HER2-amplified tumours?”*

### 3.9 Overall survival model (LASSO Cox with glmnet)

- Takes the top DE genes (e.g. top 100 by `padj`) and extracts their VST-normalised expression (`vst_mat`).
- Pulls overall survival time and status from the cleaned clinical data `clin`:
    
    ```r
    time_col   <- grep("Overall.Survival..Months|OS_MONTH", colnames(clin), value = TRUE)[1]
    status_col <- grep("Overall.Survival.Status|OS_STATUS", colnames(clin), value = TRUE)[1]
    
    ```
    
- Builds a survival object:
    
    ```r
    y <- Surv(time_raw[keep], status_num[keep])
    
    ```
    
- Creates a predictor matrix `x` where each column is one gene (VST expression), and each row is a patient.
- Fits a **LASSO Cox model** with `glmnet`:
    
    ```r
    cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1)
    coef_lasso <- coef(cvfit, s = "lambda.min")
    
    ```
    
- Extracts the non-zero coefficients (genes selected by LASSO) and saves them as a table (gene name + coefficient).

This section answers: *“Which of the DE genes are most predictive of overall survival, in a multivariable LASSO Cox model?”*

## 4. How to run the analysis

1. Place the TCGA BRCA tar file in `data/`:
    
    ```
    data/brca_tcga_pan_can_atlas_2018.tar
    ```
    
2. Open `brca_erbb2.qmd` in RStudio or use Quarto.
3. Render the report:
    - In RStudio: click **Render** (for Quarto) or **Knit** (for Rmd/Notebook).
    - Or from the terminal:
        
        ```bash
        quarto render brca_erbb2.qmd
        
        ```
        
4. Open the generated `brca_erbb2.html` (or `brca_erbb2.nb.html`) in a browser.

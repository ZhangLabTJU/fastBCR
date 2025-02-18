# fastBCR-p: Integrating Light-Chain Refinements into Heavy-Chain Clustering

</div>

## License

Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0) (see [License file](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode.zh-hans)).
This software is not to be used for commercial purposes.

## Presentation

fastBCR-p is a novel computational pipeline that incorporates public antibody filtering and strategic light-chain segregation, significantly improving both the accuracy and throughput of <strong><ins>large-scale paired BCR</ins></strong> clustering.

Existing heavy chain clonal family inference methods all suffer from the issue of confounded light chain features. To address this, based on the heavy chain clustering results from fastBCR, fastBCR-p split sequences within heavy chain clusters that have different light chain VJ genes.

* fastBCR-p provides an <strong><ins>optimized clonal family tool</ins></strong> inference for paired chain datasets.

* fastBCR-p provides an evaluation of the <strong><ins>consistency score</ins></strong> of light and heavy chains within clusters:\
   &emsp;&emsp; 1. `CDR3H`: consensus score of heavy chain junction amino acid sequence,\
   &emsp;&emsp; 2. `CDR3L`: consensus score of light chain junction amino acid sequence,\
   &emsp;&emsp; 3. `LC V gene`: the fraction of dominate V gene usage,\
   &emsp;&emsp; 4. `LC J gene`: the fraction of dominate J gene usage.

* fastBCR-p can additionally be used to predict whether a heavy chain or light chain antibody is a <strong>public antibody</strong>.\
  &emsp;&emsp; 1. <strong>public VH</strong>: binary classification of heavy chain antibodies\
  &emsp;&emsp; 2. <strong>public VL</strong>: regression prediction of light chain antibodies.

## Setup fastBCR to use fastBCR-p

fastBCR-p is an improved version of fastBCR designed for paired data. To use fastBCR-p, you need to install or update fastBCR.

### Installation of fastBCR

> :warning:  **R >= 4.1 and fastBCR >= 1.2 is required**  :warning:


Before installing fastBCR, you need to download the dependency packages
'proj4', 'msa', 'ggtree' and 'ggmsa' using Bioconductor. To install these packages, 
start R (version >= 4.1.0) and enter:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("proj4","msa","ggtree","ggmsa"))
```

Now you can install fastBCR like so:

``` r
if(!require(devtools)){
  install.packages("devtools")
}
devtools::install_github("ZhangLabTJU/fastBCR", ref = "v1.2.1")
```

## Usage

### Example pipeline

The following outlines an R-based computational pipeline meticulously crafted for the efficient analysis of bulk BCR repertoire paired chain sequencing data. This pipeline leverages the capabilities of fastBCR-p, an automated algorithm specialized in rapid clonal family inference.
You can also run the pipeline refer to [fastBCR-p_pipeline.Rmd](https://github.com/ZhangLabTJU/fastBCR/blob/v1.2.1/example/fastBCR-p_pipeline.Rmd).

### 1. Paired BCR Clonal Family Inference:

Utilizing fastBCR-p, this step efficiently identifies and categorizes paired B cell clonal families within the paired sequencing dataset.

***Example datasets***

Two real paired example AIRR Rearrangement datasets are included in the fastBCR package in the 'example_paired' folder. The datasets consist of BCR sequencing data from unsorted B Cells of a CMV patient (Run: 1287174, [Jaffe et al., 2020](https://doi.org/10.5281/zenodo.6471398)) and memory B Cells of a SARS-CoV-2 patient (Run: 1279050, [Jaffe et al., 2020](https://doi.org/10.5281/zenodo.6471398)), based on open-source antibody repertoires from the Observed Antibody Space (OAS).

***Required columns***

To infer clonal families successfully, the dataset should include essential columns:

  &emsp;&emsp; 1. `v_call_heavy`: heavy chain V gene with or without allele,\
  &emsp;&emsp; 2. `j_call_heavy`: heavy chain J gene with or without allele,\
  &emsp;&emsp; 3. `junction_aa_heavy`: amino acid translation of the heavy chain junction,\
  &emsp;&emsp; 4. `v_call_light`: light chain V gene with or without allele,\
  &emsp;&emsp; 5. `j_call_light`: light chain J gene with or without allele,\
  &emsp;&emsp; 6. `junction_aa_light`: amino acid translation of the light chain junction

***1. Data preprocessing:***

Load fastBCR package:
``` r
library(fastBCR)
```
Data preprocessing:
``` r
# load paired dataset from PAIRED_folder_path
PAIRED_folder_path <- "example/example_paired"
# the sample datasets are csv files
raw_sample_list <- data.load(folder_path = PAIRED_folder_path,
                             storage_format = "csv")
# preprocess datasets and check required columns
paired_sample_list <- paired.preprocess(raw_sample_list)

```

***2. Clonal family inferring:***
``` r
# Set 'paired = TRUE' to use fastBCR-p to infer clonal families from paired BCR sample list
cluster_lists <- data.BCR.clusters(paired_sample_list, cluster_thre = 3,
                                     overlap_thre = 0.1, consensus_thre = 0.8, 
                                     paired = TRUE)
# Set 'paired = TRUE' to use fastBCR-p to infer clonal families from paired BCR data
cluster_list_1 <- BCR.cluster(paired_sample_list[[1]], cluster_thre = 3,
                                 overlap_thre = 0.1, consensus_thre = 0.8, 
                                 paired = TRUE)
```

### 2. Public Classification Analysis:

### Setup

#### Overview
The public antibody prediction module in **fastBCR-p** leverages the pre-trained **BCR-V-BERT** model to classify heavy chain antibodies (binary classification) and predict light chain antibodies (regression). This feature allows for the identification of public antibodies using advanced transformer-based sequence analysis.

#### Python Dependency
The public antibody prediction module depends on the **BCR-V-BERT** and **PubBCRPredictor** Python package. Before running predictions, you need to install the required Python environment and dependencies as described in the [BCR-V-BERT](https://github.com/ZhangLabTJU/BCR-V-BERT) and [PubBCRPredictor](https://github.com/ZhangLabTJU/PubBCRPredictor).

***Installation***

  1. Install `reticulate` for R-Python integration:
  ```r
  install.packages("reticulate")
  library(reticulate)
  ```

  2. Create and activate the Python environment with the necessary dependencies in Anaconda Prompt:
  ```bash
  conda create --name r-py-env python=3.9
  conda activate r-py-env
  ```
  
  3. CLone repositories or download ZIP files from  [BCR-V-BERT](https://github.com/ZhangLabTJU/BCR-V-BERT) and [PubBCRPredictor](https://github.com/ZhangLabTJU/PubBCRPredictor).
  4. Install bcr_v_bert and pubbcrpredictor manually

      Setup BCR_V_BERT package manualy
      ```bash
        unzip BCR_V_BERT.zip
        cd BCR_V_BERT
        pip install -r requirements.txt
        python setup.py install
      ```

      Setup PubBCRPredictor package manualy
      ```bash
        cd ..
        unzip PubBCRPredictor.zip
        cd PubBCRPredictor
        pip install -r requirements.txt
        python setup.py install
      ```

#### Input Requirements
The public antibody prediction models require specific input formats depending on the model. Below is a summary of input requirements for each model:

| Model       | Input Requirements |
| ----------- | ----------- | 
| Heavy chain classification (all CDRs)| CDR1-3 amino acid sequences and V gene of the heavy chain antibody |
| Heavy chain classification (CDR3 only)| CDR3 amino acid sequence and V gene of the heavy chain antibody  |
| Light chain regression (all CDRs)| CDR1-3 amino acid sequences and V gene of the light chain antibody|
| Light chain regression (CDR3 only)| CDR3 amino acid sequence and V gene of the light chain antibody  |

***Example Input Format***

| Column Name | Description |
| ----------- | ----------- |
| cdr1        | CDR1 amino acid sequence       |
| cdr2        | CDR1 amino acid sequence       |
| cdr3        | CDR1 amino acid sequence       |
| vgene       | V gene identifier (e.g., “IGHV1-69”)       |

***Functions and Model Mapping***

fastBCR-p provides a unified function interface for all four models. The specific model to use depends on the model argument passed to the function:
```r
# Predict function interface
predict_public_antibody <- function(data, model = "cdrh", python_env = "r-py-env") {
  # data: Data frame with required input columns
  # model: Model type ("cdrh", "cdrh3", "cdrl", "cdrl3")
  # python_env: Python environment configured with BCR-V-BERT
  ...
}
```

| Model | model Argument | Prediction Type |
| ----------- | ----------- | ----------- |
| Heavy chain classification (V gene + all CDRs) | "cdrh"  | Binary classification |
| Heavy chain classification (V gene + CDR3 only)| "cdrh3" | Binary classification |
| Light chain regression (V gene + all CDRs)     | "cdrl"  | Regression |
| Light chain regression (V gene + CDR3 only)    | "cdrl3" | Regression |

### Usage Example
1.	Load the fastBCR library and sample data:
```r
library(fastBCR)

# Load sample data (provided in the 'example_paired' folder)
data_heavy <- read.csv("example/example_public/public_heavy_antibody.csv")
data_light <- read.csv("example/example_public/public_light_antibody.csv")
```
2. Load reticulate for R-Python integration and activate the Python environment with BCR_V_BERT and PubBCRPredictor packages
```r
library(reticulate)
# Activate the Python environment in R
use_condaenv("r-py-env", required = TRUE)
```

3.	Run public antibody prediction:
```r
# Predict public heavy chain antibodies (all CDRs)
heavy_chain_results <- predict_public_antibody(data_heavy, model = "cdrh", python_env = "r-py-env")

# Predict public light chain antibodies (all CDRs)
light_chain_results <- predict_public_antibody(data_light, model = "cdrl", python_env = "r-py-env")

# Predict public heavy chain antibodies (CDR3)
heavy_chain_results_cdr3 <- predict_public_antibody(data_heavy, model = "cdrh3", python_env = "r-py-env")

# Predict public light chain antibodies (CDR3)
light_chain_results_cdr3 <- predict_public_antibody(data_light, model = "cdrl3", python_env = "r-py-env")
```
Output:
```plaintext
> head(heavy_chain_results)
       cdr1     cdr2              cdr3    vgene label public_score
1  GYTFTGYW ILPGSGST      ARDDYDGAWFAY  IGHV1-9     1  0.999894738
2  GYTFTDYN INPNNGGT       ARGGYYYAMDY IGHV1-18     1  0.999999881
3 GYSITSGYY  ISYDGSN      ARYDYDGAWFAY  IGHV3-6     1  0.999982476
4  GFTFSSYW ISSSSSYI         ARSKQWKDY IGHV3-21     0  0.119376965
5  GGSISIYY  IYTSGST ARVIAAAGIGGGLRFDP  IGHV4-4     0  0.005968043
6  GFTFSSYA ISGHGGST  ARSNYDVWSGYPDPED IGHV3-23     0  0.443642795

> head(light_chain_results)
       cdr1 cdr2         cdr3    vgene label public_score
1    SLRSYY  GKN    QAWDSSTVV IGLV3-19     6     2.122555
2 SSNIGAGYD  GNS QSYDSSLSASVV IGLV1-40    13    22.445288
3 SSNIGAGYD  GNS QSYDSSLSGEVV IGLV1-40     8    24.287237
4  SSNIGNNY  ENN  GPWDSSLSVYV IGLV1-51     5     3.795967
5    SLRSFY  GKN  NSRDSSGNHLV IGLV3-19    14     5.469535
6 SSNIGAGSD  GNS QSYDSSLSGSWV IGLV1-40     5     5.265253
```

4.	Filter and separate high and low public antibodies:
```r

# Apply the filter to the prediction results
filtered_results <- filter_public_antibodies(heavy_chain_results, threshold = 0.5)

print(filtered_results$high_public)  # High public antibodies
print(filtered_results$low_public)  # Low public antibodies
```
Output:
```plaintext
> head(filtered_results$high_public)
        cdr1     cdr2            cdr3    vgene label public_score
1   GYTFTGYW ILPGSGST    ARDDYDGAWFAY  IGHV1-9     1    0.9998947
2   GYTFTDYN INPNNGGT     ARGGYYYAMDY IGHV1-18     1    0.9999999
3  GYSITSGYY  ISYDGSN    ARYDYDGAWFAY  IGHV3-6     1    0.9999825
7 GGSISSGSYY  IYTSGST       ATTTIRLGY IGHV4-61     1    0.5784608
8   GFTFSRYN ISYSAAGT   AKDRMGGVTFFDY IGHV3-23     0    0.5041691
9   GFTFSSYA ISGSGGST ATRPPSETYYGVLDY IGHV3-23     1    0.9698342

> head(filtered_results$low_public)
       cdr1     cdr2              cdr3    vgene label public_score
4  GFTFSSYW ISSSSSYI         ARSKQWKDY IGHV3-21     0  0.119376965
5  GGSISIYY  IYTSGST ARVIAAAGIGGGLRFDP  IGHV4-4     0  0.005968043
6  GFTFSSYA ISGHGGST  ARSNYDVWSGYPDPED IGHV3-23     0  0.443642795
10 GFTFSSYG IWYDGSNK      ARDFRVAAALDY IGHV3-33     0  0.226780578
```
***Sample Data***
Sample datasets for testing the public antibody prediction module are available in the example_paired folder. These datasets include required columns such as cdr1, cdr2, cdr3, and vgene for both heavy and light chain antibodies.

For detailed information about the BCR-V-BERT Python package and its capabilities, please refer to the [BCR-V-BERT](https://github.com/ZhangLabTJU/BCR-V-BERT).


## Issues

If you experience any issues please add an issue to the [fastBCR Issues](https://github.com/ZhangLabTJU/fastBCR/issues).

1. If you encounter the following issues during the installation of BCR_V_BERT and PubBCRPredictor packages:
    ```bash
    ERROR: pip’s dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
    xxx_b requires xxx_a>=x.xx, which is not installed.
    ```
    
   You need to install “xxx_a” and “xxx_b” package separately, and then repeat the previous steps:
    ```bash
    conda install xxx_a
    conda install xxx_b
    
    pip install -r requirements.txt
    python setup.py install
    ```

1. If you have load necessary packages but still cannot run the predict_public_antibody function, you need check if the pretrained model are loaded successfully.
    ```r
    # e.g. The cdrl3 pretrained model is not loaded into "BCR_V_BERT/model_pretrained" directory:
    > light_chain_results <- predict_public_antibody(data, model = "cdrl3", python_env = "r-py-env")
    Error in py_run_string_impl(code, local, convert) : 
      FileNotFoundError: [Errno 2] No such file or directory: '/Users/XXX/anaconda3/envs/r-py-env/lib/python3.9/site-packages/BCR_V_BERT-1.0.0-py3.9.egg/BCR_V_BERT/model_pretrained/cdrl3/v_vocab.npy'

    # e.g. The BCR_V_BERT and PubBCRPredictor packages are not installed successfully:
    > heavy_chain_results <- predict_public_antibody(data, model = "cdrh", python_env = "r-py-env")
    Error in py_run_string_impl(code, local,convert)
       NameError: name 'data' is not defined
    ```

## Reference

> Original heavy-chain clustering publication: [Wang  et al., 2023](https://doi.org/10.1016/j.crmeth.2023.100601).

> Original heavy-chain clustering protocol: [Wang  et al., 2024](https://doi.org/10.1016/j.xpro.2024.102969).

## Acknowledgements

The sample data is downloaded from Observed Antibody Space:
> Olsen, T.H., Boyles, F., and Deane, C.M. (2022). Observed Antibody Space: A diverse database of cleaned, annotated, and translated unpaired and paired antibody sequences. Protein Sci 31, 141–146. https://doi.org/10.1002/pro.4205.

<div>

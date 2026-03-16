# mitoSomatic
mitoSomatic is a highly efficient and accurate machine-learning-based tool designed to distinguish somatic mutations from germline mutations in mitochondrial DNA (mtDNA), without the need for paired control samples.

By extracting multi-dimensional biological and database features (including ANNOVAR annotations, mutation contexts, dbSNP, mitomap VAF, Phylotree, mtDB, etc.) and utilizing a pre-trained Random Forest model, mitoSomatic provides a one-step pipeline to automatically annotate and classify your mitochondrial mutation lists.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Input Format](#input-format)
- [Usage](#usage)
- [Output Files](#output-files)
- [Citation](#citation)

## Prerequisites

mitoSomatic is a hybrid pipeline built with Python and Perl. It requires minimal setup.

* **Perl**: Required for running the underlying ANNOVAR scripts (`annotate_variation.pl`, `coding_change.pl`). Pre-installed on most Linux/Mac systems.
* **Python 3.x**: Required for feature extraction and machine learning predictions.
* **Python Packages**:
  You need to install `pandas` and `scikit-learn`. You can easily install them via `pip`:
  
  ```bash
  pip install pandas scikit-learn
  ```
## Installation

**No installation is required!**
Simply download or clone the repository to your local machine, and it is ready to use out of the box. Ensure that all the `.pl` scripts, `.py` scripts, `humandb/` folder, and `.txt` database files remain in the **same directory**.

```bash
git clone https://github.com/FMMU-Xing-Lab/mitoSomatic.git
cd mitoSomatic
```

## Input Format

The input should be a plain text file (tab-separated or space-separated) containing **5 columns without a header**.

The columns must be strictly in the following order:

1. `Sample` (Sample ID)
2. `Position` (Mitochondrial position, e.g., 3469)
3. `Ref` (Reference allele)
4. `Alt` (Alternative allele)
5. `VAF` (Variant Allele Frequency, e.g., 0.005)

Example ( `input_muts.txt`):

```bash
Sample1 3469 C T 0.051
Sample1 3470 T A 0.122
Sample2 16 A T 0.880
Sample3 9999 A C 0.001
```

## Usage
You can run the entire pipeline with a single command using `main_pipeline.py`.

**Important Notes:**
1. It is highly recommended to use the **absolute path** when calling main_pipeline.py to avoid any relative path issues.
2. The output directory **will be created automatically** if it does not exist.

**Command:**
```bash
python /absolute/path/to/mitoSomatic/main_pipeline.py <Input_File> <Output_Directory>
```

**Example:**
```bash
python /home/user/software/mitoSomatic/main_pipeline.py ./data/input_muts.txt ./results/
```

## Output Files

After a successful run, two files will be generated in your specified output directory:

1. `*.features.txt` (Intermediate File)
  Contains the original 5 columns plus **9 extracted features** (`ANNOVAR`, `sub_type`, `refbase`, `Mutational_distribution`, `Mutation_assessor`, `dbSNP`, `VAF_mitomap`, `Phylotree`, `mtDB`).
2. `*.predicted.txt` (Final Result File)
   This is the final output file. It contains all the features along with a newly appended `label` column at the end. The `label` column explicitly classifies each mutation as either `somatic` or `germline`.


## Citation
If you use this software in your research, please cite our paper:
> mitoSomatic: a tool for accurate identification of mitochondrial DNA somatic mutations without paired controls

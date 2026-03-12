# Processing CellProfiler SQLite outputs with CytoTable and pycytominer

In this module, we perform preprocessing of the CellProfiler SQLite outputs to generate formatted single cell morphology data.
We take the following steps in order with each software tool.

**CytoTable**
1. Convert the SQLite output into parquet files by extracting the whole image (including metadata) table and nuclei table which both include morphology features. After conversion, metadata from the image file are added to the nuclei file to create a complete merged single cell dataset.

**pycytominer**
2. Perform annotation on the converted single cell dataset where perturbation and cell line metadata are added to each row based on well.
3. Normalize the annotated data using the standard scalar method.
4. Perform feature selection on the normalized data which performs multiple operations, including dropping features with NaNs, removing features that are highly correlated, etc.

## Perform preprocessing of morphology features

To perform the above preprocessing steps, run the below command:

```bash
# Make sure you are in the 4.preprocess_features directory in terminal
source preprocess_features.sh
```

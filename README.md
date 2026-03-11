# Nuclear speckle image-based profiling

In this repository, image analysis and image-based profiling are performed on a nuclear speckle dataset to extract CellProfiler morphology features from each channel using nucleus segmentation.

We have 3 channels in this assay:

- **DAPI** - nucleus stain
- **A647** - SON protein marker for nuclear speckles
- **Cy5/GOLD** - SRRM2 protein marker for nuclear speckles

Both of these proteins are essential for nuclear speckle formation.

![ex_image_montage](./examples/ex_image_montage.png)
> This montage shows an example image set from one site after maximum projection of 9 z-slices and illumination correction.

There are two waves of data in this project:

## Wave 1

There are 2 cell lines in this project, `786O` which is treated with siRNA and `293T` which is not treated.
Four slides/plates of data have been collected, with 8 wells each.
Of the plates, there are two layouts with two replicates each.

This dataset will be used for generating two different machine learning models to predict each nuclear speckle protein morphology:

1. Traditional regression machine learning model
2. Deep learning model

These pipelines can be found in the analysis repo called [`nuclear_speckles_analysis`](https://github.com/WayScience/nuclear_speckles_analysis).

## Wave 2

There is one cell line in this dataset: `U2OS` (bone cancer).
Three slides/plates of data was collected, with 4 wells each.

## Environments

There are three conda environments we use in this repository:

1. [`cellprofiler_env.yml`](./cellprofiler_env.yml) is used for running CellProfiler in parallel across plates.
2. [`python_env.yml`](./python_env.yml) is used for notebooks focused on either preprocessing steps or other analysis.
3. [`r_env.yml`](./r_env.yml) is used for visualization of data and generating plots.

These environments can be created using the command:

```bash
# can use either conda or mamba
conda env create -f ...
```

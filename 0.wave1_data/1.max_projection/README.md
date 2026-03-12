# Perform max projection per FOV

In this module, we use a CellProfiler pipeline per slide to generate max-projected images per image set (FOV per well).
Since the image data is in `nd2` file format, the relevant metadata for each image set is embedded in the file.
Metadata can be extracted from the file in CellProfiler, but only through the GUI at this time.
To process each plate and extract max-projected TIFF files, we must run each CellProfiler pipeline per slide in the GUI sequentially.

## Steps to process the data

1. Activate the conda environment for CellProfiler called `nuclear_speckle_cp_env`.
2. Initiate the CellProfiler GUI by using the `cellprofiler` command.
3. Once in the GUI, drag the pipeline and respective slide `nd2` file into the GUI.
4. Before running, go to the `Metadata` module and press the `Extract metadata`. This is required for CellProfiler to extract metadata from this file type.
5. Press `Run analysis` and repeat the steps for each slide/plate.

## TIFF output

In this module, we extract the maximum projected images as TIFF files per site. 
The naming convention of the files are as shown below:

**Example:**

`slide2_B4_M55_CH2_Z09.tiff`

- ***slide2*** = slide/plate
- ***B4*** = well
- ***M55*** = site
- ***CH2*** = channel/stain
- ***Z09*** = number of z-slice projected

# Download nuclear speckle dataset

In this module, we will include instructions to download the dataset once it has become publicly available. 
The current format of the raw files is `nd2`. 
There is one large file with all the images per "slide" or plate. 
There are four total slides in this dataset. 
Each come with a `positions.txt` which doubles as a plate map file given that there is well-level metadata, but it can not be used as a platemap file for CellProfiler when extracting images in [the next module](../1.max_projection/).

## Generate CSV platemap files

We convert the `txt` files into `CSV` file format for use downstream.
During this process, we also fix columns to improve format or to fix a bug that impacts downstream analysis.

To generate the platemap files, run the command below:

```bash
source generate_platemap_CSVs.sh
```

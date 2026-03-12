suppressPackageStartupMessages(library(ggplot2)) #plotting
suppressPackageStartupMessages(library(dplyr)) #work with data frames
suppressPackageStartupMessages(library(ggExtra)) #marginal density plot to add to UMAP

# Set directory and file structure
umap_dir <- file.path("results")
umap_files <- list.files(umap_dir, pattern = "\\.tsv$", full.names = TRUE)
print(umap_files)

output_fig_dir <- file.path("figures")
umap_prefix <- "UMAP_"
plate_suffix <- "_sc_feature_selected.tsv"

# Define output figure paths as a dictionary where each plate has a figure output path
output_umap_files <- list()
for (umap_file in umap_files) {
    # Use the file name to extract plate
    plate <- stringr::str_remove(
        stringr::str_remove(
            unlist(
                strsplit(umap_file, "/")
            )[2],
            umap_prefix
        ),
        plate_suffix
    )
    output_umap_files[plate] <- file.path(
        output_fig_dir,
        paste0(umap_prefix, plate)
    )
}
        
print(output_umap_files)


# Load data
umap_cp_df <- list()
for (plate in names(output_umap_files)) {
    # Find the umap file associated with the plate
    umap_file <- umap_files[stringr::str_detect(umap_files, plate)]
    
    # Load in the umap data
    df <- readr::read_tsv(
        umap_file,
        col_types = readr::cols(
            .default = "d",
            "Metadata_Plate" = "c",
            "Metadata_Well" = "c",
            "Metadata_Site" = "c",
            "Metadata_CellLine" = "c",
            "Metadata_Condition" = "c",
            "Metadata_Nuclei_Site_Count" = "d"
        )
    )

    # Append the data frame to the list
    umap_cp_df[[plate]] <- df 
}

for (plate in names(umap_cp_df)) {
    # Genotype UMAP file path
    condition_output_file <- paste0(output_umap_files[[plate]], "_condition.png")

    # UMAP labeled with condition
    condition_gg <- (
        ggplot(umap_cp_df[[plate]], aes(x = UMAP0, y = UMAP1))
        + geom_point(
            aes(color = Metadata_Condition), size = 0.8, alpha = 0.6
        )
        + theme_bw()
        + scale_color_brewer(palette = "Dark2")
        + labs(color = "siRNA Treatment")  # Change legend title
        + theme(
            legend.position = "left"
        )
    )

    # Add marginal density plots
    condition_gg_marginal <- ggMarginal(
        condition_gg,
        type = "density",
        margins = "both",
        size = 5,
        groupColour = TRUE,
        groupFill = TRUE,
    )
    
    ggsave(condition_output_file, condition_gg_marginal, dpi = 500, height = 6, width = 6)

    # UMAP labeled with cell count
    cell_count_output_file <- paste0(output_umap_files[[plate]], "_cell_count.png")
    
    umap_cell_count_gg <- (
        ggplot(umap_cp_df[[plate]], aes(x = UMAP0, y = UMAP1))
        + geom_point(
            aes(color = Metadata_Nuclei_Site_Count), size = 1.2, alpha = 0.6
        )
        + theme_bw()
        + theme(
            strip.background = element_rect(colour = "black", fill = "#fdfff4")
        )
        + scale_color_viridis_c(name = "Number of\nsingle cells\nper FOV", direction = -1, option = "inferno")
    )

    ggsave(cell_count_output_file, umap_cell_count_gg, dpi = 500, height = 6, width = 6)
}




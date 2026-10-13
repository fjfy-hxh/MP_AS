#https://github.com/stephenwilliams22/seurat/blob/stephen/spaceranger-4.0/R/preprocessing.R
library(sf)
Load10X_Spatial_seg <- function (
    data.dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1",
    bin.size = NULL,
    load.segmentations = FALSE,
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = NULL,
    ...
) {
  # if more than one directory is passed in
  if (length(x = data.dir) > 1) {
    # party on with the first value
    data.dir <- data.dir[1]
    # but also raise a warning
    warning(
      paste0(
        "`data.dir` expects a single value but received multiple - ",
        "continuing using the first: '",
        data.dir,
        "'."
      ),
      immediate. = TRUE,
    )
  }
  # if the specified directory does not exist
  if (!file.exists(data.dir)) {
    # raise an error
    stop(paste0("No such file or directory: ", "'", data.dir, "'"))
  }
  
  # if `bin.size` is not set but `data.dir` points to a folder with binned data
  if (is.null(bin.size) & file.exists(paste0(data.dir, "/binned_outputs"))) {
    # point `bin.size` to the "standard" set - i.e. everything in the default
    # output except the 8 um binning because it's a memory hog
    bin.size <- c(16, 8)
  }
  # if `bin.size` is specified
  if(!is.null(bin.size)) {
    # convert `bin.size` to a character vector and pad values to three digits
    bin.size.pretty <- paste0(sprintf("%03d", bin.size), "um")
    # point `data.dirs` to the specified binnings
    data.dirs <- paste0(
      data.dir,
      "/binned_outputs/",
      "square_",
      bin.size.pretty
    )
    # suffix assay/slice names with each bin size
    assay.names <- paste0(assay, ".", bin.size.pretty)
    slice.names <- paste0(slice, ".", bin.size.pretty)
  } else {
    # otherwise just hold onto the top-level directory
    data.dirs <- data.dir
    # and keep the assay/slice names unchanged
    assay.names <- assay
    slice.names <- slice
  }
  
  # read in counts matrices from specified h5 files
  counts.paths <- lapply(data.dirs, file.path, filename)
  counts.list <- lapply(counts.paths, Read10X_h5, ...)
  # maybe convert Cell identifiers to uppercase
  if (to.upper) {
    rownames(counts) <- lapply(rownames(counts), toupper)
  }
  
  if (is.null(image)) {
    # read in the corresponding images and coordinate mappings
    image.list <- mapply(
      Read10X_Image,
      file.path(data.dirs, "spatial"),
      assay = assay.names,
      slice = slice.names,
      MoreArgs = list(filter.matrix = filter.matrix)
    )
  } else {
    # make sure any passed images are in a vector
    image.list <- c(image)
  }
  
  # check that for each counts matrix there is a corresponding image
  if (length(image.list) != length(counts.list)) {
    stop(
      paste0(
        "The number of images does not match the number of counts matrices. ",
        "Ensure each spatial dataset has a corresponding image."
      )
    )
  }
  
  # for each counts matrix, build a Seurat object
  object.list <- mapply(CreateSeuratObject, counts.list, assay = assay.names)
  # associate each counts matrix with its corresponding image
  object.list <- mapply(
    function(
    .object,
    .image,
    .assay,
    .slice
    ) {
      # align the image's identifiers with the object's
      .image <- .image[Cells(.object)]
      # add the image to the corresponding Seurat instance
      .object[[.slice]] <- .image
      return (.object)
    },
    object.list,
    image.list,
    assay.names,
    slice.names
  )
  # merge the Seurat instances - each assay should have unique Cell identifiers
  object <- merge(
    object.list[[1]],
    y = object.list[-1]
  )
  
  # Possibly make this into a new function
  if(load.segmentations) {
    segmentation.assay.name <- "Segmentations"
    seg.data.dir <- paste0(data.dir, "/segmented_outputs/")
    seg.counts.path <- paste0(seg.data.dir, "filtered_feature_cell_matrix.h5")
    if (!file.exists(seg.counts.path)) {
      stop("No cell segmentations found. Is your Visium sample HD?")
    }
    segmentation.counts <- Read10X_h5(seg.counts.path)
    segmentation.count.cellIDs <- colnames(segmentation.counts)
    #segmentation.image <- Read10X_Image(image.dir = paste0(seg.data.dir, "spatial"), assay = "Segmentations", slice = "slice1.segmentations")
    segmentation.object <- CreateSeuratObject(segmentation.counts, assay = segmentation.assay.name)
    #segmentation.object[["slice1.segmentation.image"]] <- segmentation.image
    segmentation.fov <- Read10X_Segmentations(image.dir = paste0(seg.data.dir, "spatial"), data.dir = data.dir, cell.names = segmentation.count.cellIDs)
    segmentation.object[["slice1.polygons"]] <- segmentation.fov
    object <- merge(x = object, y = segmentation.object)
    DefaultAssay(object = object) <- segmentation.assay.name
  }
  
  return(object)
}

Read10X_Segmentations <- function (image.dir,
                                   data.dir,
                                   image.name = "tissue_lowres_image.png",
                                   assay = "Segmentations",
                                   slice = "slice1.segmentations",
                                   segmentation.type = "cell",
                                   cell.names
)
{
  
  image <- png::readPNG(source = file.path(image.dir, image.name))
  scale.factors <- Read10X_ScaleFactors(filename = file.path(image.dir,
                                                             "scalefactors_json.json"))
  coordinates <- Read10X_HD_GeoJson(data.dir = data.dir, image.dir = image.dir)
  key <- Key(slice, quiet = TRUE)
  
  
  st_coords_df <- as.data.frame(sf::st_coordinates(coordinates))
  st_coords_df$index <- 1:nrow(st_coords_df)
  st_coords <- merge(
    st_coords_df,
    data.frame(
      cell_id = coordinates$cell_id,
      cell = coordinates$barcodes
    ),
    by.x = 'L2',
    by.y = 'cell_id',
    all = TRUE
  )
  st_coords <- st_coords[complete.cases(st_coords), ]
  st_coords <- st_coords[order(st_coords$index), ]
  st_coords <- st_coords[st_coords$cell %in% cell.names, ]
  image.scale <- ifelse(grepl("lowres", image.name, ignore.case = TRUE), "lowres", "highres")
  fov <- CreateFOV(st_coords,
                   type = "segmentation",
                   #radius = scale.factors[[image.scale]],
                   assay = assay, key = key)
  return(fov)
}

#' Format 10X Genomics GeoJson cell IDs
#'
#' A helper function to format cell IDs from the segmentation GeoJson to the same type as in the matrix.h5
#' The GeoJson has cell IDs as integers (eg 1). They need to be in the format cellid_000000001-1
#'
#' @return Vector of formatted cell IDs
Read10X_ScaleFactors <- function(filename) {
  raw.data <- jsonlite::fromJSON(file.path(filename))
  
  scale.factors <- scalefactors(
    spot = raw.data$spot_diameter_fullres,
    fiducial = raw.data$fiducial_diameter_fullres,
    hires = raw.data$tissue_hires_scalef,
    lowres = raw.data$tissue_lowres_scalef
  )
  
  return (scale.factors)
}

Read10X_HD_GeoJson <- function(data.dir, image.dir, segmentation.type = "cell", scale.factor = NULL) {
  segmentation_polygons <- read_sf(file.path(data.dir,"segmented_outputs", paste0(segmentation.type, "_segmentations.geojson")))
  if (!is.null(scale.factor)) {
    scale.factors <- Read10X_ScaleFactors(
      filename = file.path(image.dir, "scalefactors_json.json")
    )
    segmentation_polygons$geometry <- segmentation_polygons$geometry*scale.factors[[scale.factor]]
  }
  
  segmentation_polygons$barcodes <- Format10X_GeoJson_CellID(segmentation_polygons$cell_id)
  segmentation_polygons
}


Format10X_GeoJson_CellID <- function(ids, prefix = "cellid_", suffix = "-1", digits = 9) {
  format_string <- paste0("%0", as.integer(digits), "d")
  
  formatted_ids <- sapply(ids, function(id) {
    numeric_part <- sprintf(format_string, as.integer(id))
    paste0(prefix, numeric_part, suffix)
  })
  
  return(formatted_ids)
}

#' Extract a subset based on spatial coordinates
#'
#' @param seurat_obj Seurat object (must contain spatial coordinates)
#' @param x_prop X axis proportion threshold (0-1)
#' @param y_prop Y axis proportion threshold (0-1)
#' @param x_dir X axis direction ("<" or ">")
#' @param y_dir Y axis direction ("<" or ">")
#' @param coord_slot Where the coordinates are stored (default "spatial")
#' @param x_col X coordinate column name (default "x")
#' @param y_col Y coordinate column name (default "y")
#' @param cell_id_col Cell ID column name (default "cell")
#'
#' @return Seurat subset containing the spots in the specified region









#' Spatial subset extraction function based on GetTissueCoordinates
#'
#' @return Seurat subset containing the spots in the specified region
ExtractByCoordinates <- function(seurat_obj, 
                                 x_prop = 0.5, 
                                 y_prop = 0.5,
                                 x_dir = "<",
                                 y_dir = "<",
                                 x_col = "x",
                                 y_col = "y") {
  
  # Get the coordinate data (force using GetTissueCoordinates)
  coords <- GetTissueCoordinates(seurat_obj)
  # Compute the coordinate bounds
  x_range <-range(coords$x)
  y_range <- range(coords$y)
  
  # Compute the dynamic thresholds
  x_threshold <- x_range[1] + x_prop * diff(x_range)
  y_threshold <- y_range[1] + y_prop * diff(y_range)
  
  # Build the filter conditions
  x_condition <- switch(x_dir,
                        ">" = coords$x > x_threshold,
                        "<" = coords$x < x_threshold,
                        stop("x_dir must be '<' or '>'"))
  
  y_condition <- switch(y_dir,
                        ">" = coords$y > y_threshold,
                        "<" = coords$y < y_threshold,
                        stop("y_dir must be '<' or '>'"))
  
  # Get the cell barcodes (assume the coordinate row names are barcodes)
  selected_cells <- coords$cell[x_condition & y_condition]
  
  # Create the subset
  subset_obj <- subset(seurat_obj, cells = selected_cells)
  return(subset_obj)
}

##############################################################
###########################################################
Cal_avg_coord_all_images <- function(seurat_obj, 
                                     cell_col = "cell",
                                     x_col = "x", 
                                     y_col = "y") {
  
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  # Check whether images exist
  if (length(seurat_obj@images) == 0) {
    stop("No images data found in the Seurat object")
  }
  
  # Initialize the list storing all image results
  all_avg_coords <- list()
  
  # Iterate over each image
  for (n in seq_along(seurat_obj@images)) {
    # Extract the coordinate data of the current image
    coords <- GetTissueCoordinates(seurat_obj@images[[n]])
    
    # Check whether the required columns exist
    required_cols <- c(cell_col, x_col, y_col)
    missing_cols <- setdiff(required_cols, colnames(coords))
    if (length(missing_cols) > 0) {
      warning(paste0("Image ", n, " missing columns: ", paste(missing_cols, collapse = ", "), ", skipping"))
      next
    }
    
    # Extract unique cell IDs
    unique_cells <- unique(coords[[cell_col]])
    
    # Compute the average coordinates of each cell in the current image
    avg_coords <- do.call(rbind, lapply(unique_cells, function(cell) {
      cell_data <- coords[coords[[cell_col]] == cell, ]
      data.frame(
        cell = cell,
        x_avg = mean(cell_data[[x_col]], na.rm = TRUE),
        y_avg = mean(cell_data[[y_col]], na.rm = TRUE),
        image = names(seurat_obj@images)[n],  # record the image name
        stringsAsFactors = FALSE
      )
    }))
    
    # Store the current image result
    all_avg_coords[[n]] <- avg_coords
  }
  
  # Merge the results of all images
  combined_avg_coords <- do.call(rbind, all_avg_coords)
  rownames(combined_avg_coords) <- paste(combined_avg_coords$cell, combined_avg_coords$image, sep = "_")
  
  return(combined_avg_coords)
}


#' @title Perform RCTD Deconvolution and Add Annotations to Seurat Object
#' @description This function performs RCTD deconvolution, adds the results as metadata to a Seurat object,
#'              and returns both the modified Seurat object and the full RCTD results.
#' @param seurat_obj A Seurat object containing spatial transcriptomics data
#' @param reference_RCTD A pre-built RCTD reference object
#' @param assay_name Name of the assay to use for counts data (default: "Segmentations")
#' @param slot_name Name of the slot to use (default: "counts")
#' @param cell_col Name of the column containing cell IDs in metadata (default: "cell")
#' @param x_col Name of the column containing x coordinates (default: "x")
#' @param y_col Name of the column containing y coordinates (default: "y")
#' @param UMI_min Minimum UMI count threshold for cell filtering (default: 10)
#' @param doublet_mode RCTD doublet mode ("doublet", "full", or "multi") (default: "doublet")
#' @param max_cores Maximum number of cores to use (default: 20)
#' @param annotation_col Name for the new annotation column in metadata (default: "Cell_Anno")
#' @return Returns a list containing:
#'         - $seurat_obj: The Seurat object with added RCTD annotations
#'         - $rctd_results: The complete RCTD results object
RunRCTD_Annotations <- function(seurat_obj, 
                                reference_RCTD,
                                data_type = "xenium",
                                assay_name = "Segmentations",
                                slot_name = "counts",
                                cell_col = "cell",
                                x_col = "x",
                                y_col = "y",
                                UMI_min = 10,
                                doublet_mode = "doublet",
                                max_cores = 20,
                                annotation_col = "Cell_Anno") {
  
  if (!requireNamespace("spacexr", quietly = TRUE)) {
    stop("Please install the spacexr package: devtools::install_github('dmcable/spacexr')")
  }
  
  counts <- tryCatch({
    # Try Seurat v5 syntax first
    GetAssayData(seurat_obj, assay = assay_name, layer = slot_name)
  }, error = function(e) {
    # If it fails, fall back to the Seurat v4 slot syntax
    GetAssayData(seurat_obj, assay = assay_name, slot = slot_name)
  })
  
  if (data_type == "visiumHD") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Please install Seurat package for visiumHD coordinate calculation")
    }
    if (!exists("Cal_avg_coord")) {
      stop("Cal_avg_coord function is not defined. Please define it for visiumHD data.")
    }
    rctd_coords <- Cal_avg_coord(
      seurat_obj,
      cell_col = cell_col,
      x_col = x_col,
      y_col = y_col
    )
    
  } else if (data_type == "xenium") {
    if (is.null(seurat_obj@images) || length(seurat_obj@images) == 0) {
      stop("No image data found in seurat_obj@images. Check if fov exists.")
    }
    
    image_name <- names(seurat_obj@images)[1]
    if (is.null(seurat_obj@images[[image_name]]$centroids)) {
      stop(paste("No centroids found in", image_name, 
                 "Try checking seurat_obj@images[[image_name]] structure"))
    }
    
    rctd_coords <- as.data.frame(seurat_obj@images[[image_name]]$centroids@coords)
    seurat_cells <- colnames(seurat_obj)
    
    if (nrow(rctd_coords) != length(seurat_cells)) {
      stop(paste("Number of coordinates (", nrow(rctd_coords), 
                 ") does not match number of cells in seurat object (", 
                 length(seurat_cells), ")"))
    }
    
    rownames(rctd_coords) <- seurat_cells
    cat("Set coordinates rownames from seurat object cell names: ", 
        length(seurat_cells), " cells\n")
    
  } else if (data_type == "pcf") {
    if (!"spatial" %in% names(seurat_obj@reductions)) {
      stop("No 'spatial' reduction found. Please create it first.")
    }
    rctd_coords <- as.data.frame(Embeddings(seurat_obj, reduction = "spatial"))
    
  } else {
    stop("Invalid data_type. Must be 'xenium', 'visiumHD', or 'pcf'")
  }
  
  if (!is.data.frame(rctd_coords)) {
    rctd_coords <- as.data.frame(rctd_coords)
  }
  
  if (ncol(rctd_coords) == 2) {
    colnames(rctd_coords) <- c("xcoord", "ycoord")
  }
  
  common_cells <- intersect(rownames(rctd_coords), colnames(counts))
  if (length(common_cells) == 0) {
    stop("No common cells found between coordinates and counts matrix")
  }
  
  counts <- counts[, common_cells, drop = FALSE]
  rctd_coords <- rctd_coords[common_cells, , drop = FALSE]
  nUMI <- colSums(counts)
  names(nUMI) <- common_cells
  
  cat("Creating SpatialRNA object with", nrow(rctd_coords), "cells and", 
      nrow(counts), "genes\n")
  
  spatial_data <- spacexr::SpatialRNA(
    coords = rctd_coords,
    counts = counts,
    nUMI = nUMI
  )
  
  myRCTD <- spacexr::create.RCTD(
    spatial_data,
    reference_RCTD,
    max_cores = max_cores,
    UMI_min = UMI_min
  )
  
  myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = doublet_mode)
  
  rctd_cells <- rownames(myRCTD@results$results_df)
  seurat_cells <- colnames(seurat_obj)
  common_cells_result <- intersect(rctd_cells, seurat_cells)
  
  cat("Cell matching statistics:\n",
      "Total cells in Seurat object:", length(seurat_cells), "\n",
      "Cells retained by RCTD:", length(rctd_cells), "\n",
      "Successfully matched cells:", length(common_cells_result), "\n",
      "Percentage matched:", round(length(common_cells_result)/length(seurat_cells)*100, 1), "%\n")
  
  rctd_annot <- myRCTD@results$results_df$first_type
  names(rctd_annot) <- rownames(myRCTD@results$results_df)
  
  cell_decon <- rep(NA, ncol(seurat_obj))
  names(cell_decon) <- colnames(seurat_obj)
  cell_decon[common_cells_result] <- rctd_annot[common_cells_result]
  
  cell_decon_labels <- ifelse(
    is.na(cell_decon) | cell_decon > length(levels(rctd_annot)),
    NA,
    levels(rctd_annot)[cell_decon]
  )
  
  seurat_obj <- AddMetaData(
    seurat_obj,
    metadata = cell_decon_labels,
    col.name = annotation_col
  )
  
  return(seurat_obj)
}

#' @title Generate a barcode information data frame with a custom suffix
#' @description Extract barcodes (with a custom suffix appended) and sample origin information from multiple Seurat objects
#' @param seurat_list Named list containing Seurat objects
#' @param suffix Suffix to append (default "_slice2")
#' @return Returns a data frame with the following columns:
#'         - original_barcode: original barcode
#'         - barcode_with_suffix: barcode with the suffix appended
#'         - group: sample name of origin (the original list name)
#' @examples
#' # Example usage:
#' # my_list <- list(s1 = af.5w_1, s2 = af.5w_2)
#' # df <- generate_barcode_df(my_list, suffix = "_slice2")
generate_barcode_df <- function(seurat_list, suffix = "_slice2") {
  # Input validation
  if (!is.list(seurat_list) || is.null(names(seurat_list))) {
    stop("Input must be a named list")
  }
  
  if (!is.character(suffix) || length(suffix) != 1) {
    stop("Suffix must be a single string")
  }
  
  result_df <- data.frame()
  
  for (obj_name in names(seurat_list)) {
    seurat_obj <- seurat_list[[obj_name]]
    
    # Skip non-Seurat objects
    if (!inherits(seurat_obj, "Seurat")) {
      warning("Skipping non-Seurat object: ", obj_name)
      next
    }
    
    # Get the barcodes
    barcodes <- colnames(seurat_obj)
    if (is.null(barcodes)) {
      warning("Object is missing barcode information: ", obj_name)
      next
    }
    
    # Add to the result data frame (rename the sample column to group)
    result_df <- rbind(result_df, data.frame(
      original_barcode = barcodes,
      barcode_with_suffix = paste0(barcodes, suffix),
      group = obj_name,  # key change: rename the column to group
      stringsAsFactors = FALSE
    ))
  }
  
  rownames(result_df) <- NULL
  return(result_df)
}
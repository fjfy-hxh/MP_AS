my_CreateScRank <- function(
    input, 
    meta = NULL, 
    cell_type = "celltype", 
    species, 
    drug = NULL, 
    target = NULL, 
    type = "antagonist", 
    if_cluster = FALSE, 
    var.gene = NULL
) {
  if (is(input, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Please confirm \"Seurat\" R package is installed correctly")
    }
    seuratObj <- input
    meta <- input@meta.data %>% droplevels()
    cell_types <- meta[[cell_type]]
    
    # 🌟 core fix: compatible with Seurat v5 low-level counts extraction
    assay_name <- "RNA"
    if (!assay_name %in% names(input@assays)) {
      assay_name <- Seurat::DefaultAssay(input)
    }
    
    # completely bypass the slot error, safely extract the counts matrix and convert to matrix
    raw_counts <- tryCatch({
      SeuratObject::LayerData(object = input, layer = "counts", assay = assay_name)
    }, error = function(e) {
      assay_obj <- input[[assay_name]]
      if (.hasSlot(assay_obj, "layers") && "counts" %in% names(assay_obj@layers)) {
        mat <- assay_obj@layers[["counts"]]
        rownames(mat) <- rownames(input)
        colnames(mat) <- colnames(input)
        return(mat)
      } else {
        return(slot(assay_obj, "counts"))
      }
    })
    
    data <- as.matrix(raw_counts)
    
  } else {
    if (all(is.null(meta), if_cluster)) {
      message(crayon::cyan("Note: Metadata will be created after clustering"))
    } else if (is.null(meta)) {
      stop("Please provide a meta if not input a Seurat object!")
    } else if (all(is.data.frame(meta), colnames(meta) %in% cell_type)) {
      stop("Please provide a correct meta formatted with data.frame!")
    }
    
    seuratObj <- NULL
    valid_format <- is(input, "dgCMatrix") || tester::is_numeric_matrix(input) || tester::is_numeric_dataframe(input)
    if (valid_format) {
      data <- methods::as(input, "dgCMatrix")
    } else {
      data <- input
    }
    
    if (all(is.null(meta), if_cluster)) {
      if (requireNamespace("scSHC", quietly = TRUE)) {
        message(crayon::cyan("Running clustering by scSHC..."))
        if (is.null(var.gene)) {
          clusters <- scSHC::scSHC(data)
        } else {
          commone_features <- intersect(rownames(data), var.gene)
          clusters <- scSHC::scSHC(data[commone_features, ], num_features = length(commone_features))
        }
        cells <- names(clusters[[1]])
        labels <- unname(clusters[[1]])
        meta <- data.frame(celltype = labels)
        rownames(meta) <- cells
      } else {
        stop("scSHC is required for this functionality. Please install it using devtools::install_github(\"igrabski/sc-SHC\").")
      }
    }
    
    if (!identical(ncol(data), nrow(meta))) {
      stop(paste0("The cells in meta and data are mismatched, where meta with ", 
                  nrow(meta), " cells and data with ", ncol(data), " cells."))
    }
    if (!identical(colnames(data), rownames(meta))) {
      stop("colnames(data) is not consistent with rownames(meta)")
    }
  }
  
  if (is.null(type)) {
    stop("Error! Please provide correct type of drug (antagonist or agnoist).")
  } else if (!type %in% c("antagonist", "agonist")) {
    stop("Error! Please provide correct type of drug (antagonist or agnoist).")
  }
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Please confirm \"Seurat\" R package is installed correctly", .call = FALSE)
  }
  
  # 🌟 fix: call the private internal function .hvg of scRank
  hvg_func <- get(".hvg", envir = asNamespace("scRank"))
  utile_db <- get("utile_database", envir = asNamespace("scRank"))
  
  if (is(input, "Seurat")) {
    if (nrow(input) > 2000) {
      hvg <- hvg_func(input)
    } else {
      hvg <- hvg_func(input, nfeatures = nrow(input))
    }
  } else {
    seuratObj <- Seurat::CreateSeuratObject(data, meta.data = meta)
    if (nrow(seuratObj) > 2000) {
      hvg <- hvg_func(seuratObj)
    } else {
      hvg <- hvg_func(seuratObj, nfeatures = nrow(seuratObj))
    }
  }
  
  # database field retrieval
  if (species == "human") {
    drug_gene <- utile_db$Drug_Target$human$Symbol
    tf_gene <- utile_db$Gene_TF$human$Symbol
  } else if (species == "mouse") {
    drug_gene <- utile_db$Drug_Target$mouse$mousegene
    tf_gene <- utile_db$Gene_TF$mouse$Symbol
  } else {
    stop("species must be either 'human' or 'mouse'")
  }
  
  if (is.null(drug)) {
    if (is.null(target)) {
      stop("Please input the target gene of inhibitor or inhibitor name")
    } else {
      if (all(target %in% rownames(data))) {
        target_gene <- target
      } else {
        stop("Please check if the target gene is in the gene expression profile.")
      }
    }
  } else {
    drug_name <- toupper(drug)
    if (drug_name %in% utile_db$Drug_Info$drug_name) {
      target <- utile_db$Drug_Info %>% 
        dplyr::filter(drug %in% drug_name) %>% 
        dplyr::top_n(1, score) %>% 
        dplyr::pull(gene_name)
    } else {
      stop("Drug name is not in our list, please assign the known inhibited\n target gene name in to CreateScRank(target=`known target`)")
    }
  }
  
  gene4use <- unique(c(target, hvg, tf_gene, drug_gene))
  
  # remove mitochondrial and ribosomal genes
  drop_idx <- grep(pattern = "^RP[[:digit:]]+|^RPL|^RPS|^MT-", toupper(gene4use))
  if (length(drop_idx) > 0) {
    gene4use <- gene4use[-drop_idx]
  }
  gene4use <- gene4use[gene4use %in% rownames(data)]
  
  object <- new("scRank", 
                data = list(data = data, seuratObj = seuratObj), 
                meta = list(rawmeta = meta), 
                para = list(species = species, 
                            cell_type = cell_type, 
                            drug = drug, 
                            target = target, 
                            gene4use = gene4use, 
                            type = type))
  return(object)
}

BIGprofiler <- function(gene_list = NULL, gene_df = NULL, ID = "SYMBOL",
                        species = "human",
                        collection = NULL, subcollection = NULL,
                        db = NULL,
                        minGSSize = 10,
                        maxGSSize = 500,
                        category = NULL,
                        subcategory = NULL){
  db_join <- pathway_GOID <- gs_exact_source <- BgRatio <- Description <- FDR <- GeneRatio <- ensembl_gene <- entrez_gene <- geneID <- gene_symbol <- genes <- group <- n_query_genes_in_pathway <- gs_collection <- gs_name <- gs_subcollection <- `k/K` <- p.adjust <- pathway <- pval <- pvalue <- qvalue <- n_background_genes <- n_query_genes <- n_pathway_genes <- gene_list_overlap <- db_species <- NULL
  
  if(!is.null(category)){stop("category is deprecated. Please use collection.")}
  
  ##### Database #####
  #Load gene ontology
  if(!is.null(collection)){
    #Check that collection exists in msigdb
    all_cat <- msigdbr::msigdbr_collections() %>%
      dplyr::pull(gs_collection) %>% unique()
    if(!collection %in% all_cat){
      stop("collection does not exist. Use msigdbr::msigdbr_collections() to see options.") }
    
    #Recode species
    if(species == "human"){
      species <- "Homo sapiens"
      db_species <- "HS"
    }
    if(species == "mouse"){
      species <- "Mus musculus"
      db_species <- "MM"
    }
    db.format <- msigdbr::msigdbr(species, db_species, collection=collection)
    #Subset subcollection if selected
    if(!is.null(subcollection)){
      #Check that subcollection exists in msigdb
      all_subcat <- msigdbr::msigdbr_collections() %>%
        dplyr::pull(gs_subcollection) %>% unique()
      if(!subcollection %in% all_subcat){
        stop("Subcollection does not exist. Use msigdbr::msigdbr_collections() to see options.") }
      
      db.format <- db.format %>%
        dplyr::filter(grepl(paste0("^",subcollection), gs_subcollection))
    }
  } else if(!is.null(db)){
    db.format <- db %>%
      dplyr::mutate(gs_collection = "custom", gs_subcollection=NA)
    #Name columns to match mgsibdbr
    colnames(db.format)[1] <- "gs_name"
    if(ID == "SYMBOL"){
      colnames(db.format)[2] <- "gene_symbol"
    } else if(ID == "ENSEMBL"){
      colnames(db.format)[2] <- "ensembl_gene"
    } else if(ID == "ENTREZ"){
      colnames(db.format)[2] <- "entrez_gene"
    } else{
      stop("Please like ID from SYMBOL, ENSEMBL, or ENTREZ.")
    }
  } else {
    stop("Please provide gene set information as Broad collection/subcollection or in a data frame as db.")
  }
  
  #Get gene ID
  if(ID == "SYMBOL"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, gene_symbol, gs_collection, gs_subcollection)
  } else if(ID == "ENSEMBL"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, ensembl_gene, gs_collection, gs_subcollection)
  } else if(ID == "ENTREZ"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, entrez_gene, gs_collection, gs_subcollection)
  } else{
    stop("Please like ID from SYMBOL, ENSEMBL, or ENTREZ.")
  }
  
  ##### Format data ####
  if(!is.null(gene_df)){
    gene_list_format <- list()
    col1 <- colnames(gene_df)[1]
    col2 <- colnames(gene_df)[2]
    for(g in unique(unlist(gene_df[,1]))){
      gene_list_format[[g]] <- gene_df %>%
        dplyr::filter(get(col1) == g) %>%
        dplyr::pull(get(col2)) %>% unique()
    }
  } else if(!is.null(gene_list)){
    gene_list_format <- gene_list
    
  } else{
    stop("Please provide either gene_list or gene_df.")
  }
  
  ##### Loop through gene df #####
  #Blank holders
  all.results <- list()
  
  for(g in names(gene_list_format)){
    print(g)
    
    #Check that genes exist in db
    gene_list_overlap <- base::intersect(gene_list_format[[g]],
                                         unlist(db.format2[,2]))
    if(length(gene_list_overlap)>0){
      #run enrichment on gene list
      enrich.result <- clusterProfiler::enricher(gene=gene_list_overlap,
                                                 TERM2GENE=db.format2[,1:2],
                                                 minGSSize = minGSSize,
                                                 maxGSSize = maxGSSize)
      
      #handle no enrichment results
      if(is.null(enrich.result)){
        if(!is.null(db)){
          result.clean <- data.frame(
            group=g,
            gs_collection="custom",
            pathway="No enriched terms")
        } else{
          result.clean <- data.frame(
            group=g,
            gs_collection=collection,
            gs_subcollection=subcollection,
            pathway="No enriched terms")
        }
        
      } else{
        #Format collection labels
        db.species.clean <- db.format2 %>%
          dplyr::distinct(gs_collection, gs_subcollection, gs_name) %>%
          dplyr::rename(pathway=gs_name)
        
        #Format results
        #Format gene column to vector
        result.clean <- enrich.result@result %>%
          tibble::remove_rownames() %>%
          dplyr::rename(pathway=Description, FDR=p.adjust, pval=pvalue) %>%
          dplyr::mutate(genes = strsplit(geneID, split="/")) %>%
          #Extract values from ratios
          tidyr::separate_wider_delim(BgRatio, names=c("n_pathway_genes","n_background_genes"), delim="/") %>%
          tidyr::separate_wider_delim(GeneRatio, names=c("n_query_genes_in_pathway",
                                                         "group_in_cat.subcat"),
                                      delim="/") %>%
          dplyr::mutate_at(dplyr::vars("n_pathway_genes","n_background_genes",
                                       "n_query_genes_in_pathway","group_in_cat.subcat"),
                           as.numeric) %>%
          #Calculate k/K
          dplyr::mutate("k/K"=n_query_genes_in_pathway/n_pathway_genes) %>%
          
          #Add ID columns for database names
          dplyr::left_join(db.species.clean, by = "pathway") %>%
          #Add columns for group info
          dplyr::mutate(group=g, n_query_genes = length(gene_list_format[[g]])) %>%
          #Reorder variables
          dplyr::select(group, n_query_genes,
                        gs_collection, gs_subcollection, n_background_genes,
                        # group_in_cat.subcat,
                        pathway, n_pathway_genes, n_query_genes_in_pathway, `k/K`,
                        pval, FDR, qvalue, genes) %>%
          dplyr::arrange(FDR)
        
        # add GO term reference ID to results
        if(!is.null(collection)){
          if(collection == "C5"){
            db_join <- db.format %>%
              dplyr::select(c("gs_name", "gs_exact_source")) %>%
              dplyr::distinct()
            result.clean <- result.clean %>%
              dplyr::left_join(db_join, by = c("pathway" = "gs_name")) %>%
              dplyr::rename(pathway_GOID = gs_exact_source) %>%
              dplyr::relocate(pathway_GOID, .after = pathway)
          }
        }
        
        #Run enrich and save to results list
        all.results[[g]] <- result.clean
      }
    }else {
      all.results[[g]] <- tibble::tibble(
        group=g, n_query_genes = length(gene_list_format[[g]]),
        gs_collection=collection, gs_subcollection=subcollection,
        n_background_genes=NA, #group_in_cat.subcat=0,
        pathway="No overlap of query genes and specified database.")
    }
  }
  
  ##### Save results #####
  #combine list of df results
  results.all.df <- dplyr::bind_rows(all.results)
  return(results.all.df)
}
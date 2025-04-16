#' Filter rare and low abundance genes
#'
#' Filter genes at a minimum counts per million (CPM) in a minmum number or percent of total samples.
#'
#' @param dat DGEList output by edgeR::DEGList( )
#' @param min_CPM numeric minimum counts per million (CPM)
#' @param gene_var character name for column with gene names in dat$genes that matches names in expression data dat$E. Default "ensembl_gene_id"
#' @param min_sample numeric minimum number of samples
#' @param min_pct numeric minimum percent of samples (0-100)
#' @param plot logical if should plot mean variance trends
#'
#' @param min.CPM Deprecated form of min_CPM
#' @param gene.var Deprecated form of gene_var
#' @param min.sample Deprecated form of min_sample
#' @param min.pct Deprecated form of min_pct
#'
#' @return DGEList object filtered to not rare genes
#' @export
#'
#' @examples
#' dat.filter <- filter_rare(dat = example.dat, min_CPM = 0.1, min_sample = 3,
#'                           gene_var="geneName")
#' dat.filter <- filter_rare(dat = example.dat, min_CPM = 0.1, min_pct = 10,
#'                           plot = TRUE, gene_var="geneName")

filter_rare <- function(dat, min_CPM, gene_var="ensembl_gene_id",
                        min_sample=NULL, min_pct=NULL, plot=FALSE,
                        #Deprecated
                        min.CPM=NULL, gene.var=NULL,
                        min.sample=NULL, min.pct=NULL){
  # BAck compatibility
  if(!is.null(min.CPM)){min_CPM <- min.CPM}
  if(!is.null(min.sample)){min_sample <- min.sample}
  if(!is.null(min.pct)){min_pct <- min.pct}
  if(!is.null(gene.var)){gene_var <- gene.var}
  
  x <- y <- linex <- liney <- NULL
  
  ##### Check parameters #####
  #Correct input object type?
  if(!isa(dat, "DGEList")){ stop("dat object must be a DGEList object from edgeR") }
  #Set at least one min sample param
  if(is.null(min_sample) & is.null(min_pct)){ stop("Please provide one of min_sample or min_pct") }
  #min samples or percent only?
  if(!is.null(min_sample) & !is.null(min_pct)){ stop("Please provide only one of min_sample or min_pct") }
  #percent as non-decimal?
  if(!is.null(min_pct)){
    if(min_pct<1){ warning("min_pct should be percentage between 0 and 100. Your value appears to be a proportion less than 1. Please verify.") }}
  
  ##### Define min number of samples #####
  #Calculate min samples based on percent if provided
  if(!is.null(min_pct)){
    min_sample <- round(nrow(dat$samples)*min_pct/100, digits=0)
  }
  
  ##### List not rare genes #####
  # Convert counts to counts per million
  dat.cpm <- edgeR::cpm(dat$counts)
  # Calculate number of samples that meet the cutoff per gene
  not.rare.samples <- rowSums(dat.cpm >= min_CPM)
  # List not rare genes to be RETAINED
  not.rare.genes <- names(not.rare.samples[not.rare.samples >= min_sample])
  
  ##### Filter data to remove rare genes #####
  dat.filter <- dat
  # Filter counts
  dat.filter$counts <- dat.filter$counts[rownames(dat.filter$counts) %in% not.rare.genes,]
  
  #If gene info exists, filter as well
  if(!is.null(dat.filter$genes)){
    #Gene name column in gene info table?
    if(!(gene_var %in% colnames(dat.filter$genes))){
      stop("Gene name variable not present in gene info (dat$genes)") }
    
    # Filter gene key
    dat.filter$genes <- dat.filter$genes[dat.filter$genes[,gene_var] %in% not.rare.genes,]
  }
  
  ##### plot #####
  if(plot){
    temp <- limma::voom(dat, plot=FALSE, save.plot = TRUE)
    temp2 <- limma::voom(dat.filter, plot=FALSE, save.plot = TRUE)
    
    plot <- data.frame(
      group = c(rep("All genes",length(temp$voom.xy$x)),
                rep("Filtered genes",length(temp2$voom.xy$x))),
      x = c(temp$voom.xy$x, temp2$voom.xy$x),
      y = c(temp$voom.xy$y, temp2$voom.xy$y),
      linex = c(temp$voom.line$x, temp2$voom.line$x),
      liney = c(temp$voom.line$y, temp2$voom.line$y)) %>%
      
      ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x=x, y=y), size=0.5) +
      ggplot2::geom_path(ggplot2::aes(x=linex, y=liney), color="red") +
      ggplot2::theme_classic() +
      ggplot2::facet_wrap(~group, nrow=1) +
      ggplot2::labs(x="log2( count size + 0.5 )", y="Sqrt (stdev)")
    print(plot)
  }
  
  ##### message #####
  #calculate total genes and removed
  tot.genes <- nrow(dat)
  filter.tot <- tot.genes-length(not.rare.genes)
  filter.perc <- round(filter.tot/tot.genes*100, digits=2)
  message(paste0(filter.tot, " (", filter.perc,  "%) of ", tot.genes, " genes removed."))
  
  return(dat.filter)
}
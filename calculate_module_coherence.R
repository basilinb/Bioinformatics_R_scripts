calculate_module_coherence<-
  function(mods,
           mods_title="STUDY1",
           mods_var = "module.char",
           remove_mods = c(0,"0","00","grey"),
           dat,
           gene_var="geneName",
           dat_title="STUDY2",
           return_plot=TRUE,
           r_cutoff = 0.3,
           p_cutoff = 0.01){
    
    # set undefined global variables to null
    X1 <- X2 <- allCorVals<- allPVals <- V3<- Cor<- P <- Set <- gene1<- gene2<- P.est<- h_label<- h_pos<- negLogP<- value<- h_line<- NULL
    
    # convenience define module definition dataframe and check data formatting
    moduleSets<- mods
    if(class(dat)[1]== "EList"){voomData<-dat$E}else{voomData<-dat}
    if(!is.numeric(voomData)){
      print("Warning, count/expression matrix is non-numeric.
        This usually arises when a non-expression column has been retained from a metadata object (eg due to grouping).
        Be advised to check for stray metadata in the expression matrix.")}
    
    
    #-----------------------------------------------------------------
    #  Calculate the gene set means (module expression)
    #  [this is currently not used, but may be incorporated later]
    #-----------------------------------------------------------------
    uniGS<-moduleSets%>%
      dplyr::pull(mods_var)%>%
      unique()%>%as.character()%>%gtools::mixedsort() # Pull unique modules & order correctly
    
    # establish storage
    GSsub<-c()
    GSabsent<-list()
    
    # for each geneset pull corresponding gene indices and compute means
    for (i in uniGS){
      curGS<-i
      curIDs<-as.character(moduleSets[which(moduleSets[, mods_var]==curGS), gene_var])
      
      matchIndex<-match(curIDs, rownames(voomData))
      if (any(is.na(matchIndex))){
        matchIndex<-matchIndex[-which(is.na(matchIndex))]}
      GSabsent[[i]]<-curIDs[which(!(curIDs %in% rownames(voomData)))]
      if (length(curIDs[which(curIDs %in% rownames(voomData))]) > 1){
        GSsub<-rbind(GSsub, apply(voomData[matchIndex,], 2, mean))
      }else{
        GSsub<-rbind(GSsub, voomData[matchIndex,])}
    }
    
    # print message if genes are missing from genesets.
    if(length(unlist(GSabsent))>0){
      print(paste0(length(unique(unlist(GSabsent))), " Module Genes absent in voom object. Missing genes accesible in .$GSabsent"))}
    
    #-------------------------------------------
    #  Set function for correlation p values
    #-------------------------------------------
    
    # Correlation p value function
    cor2pvalue = function(r, n) {
      t <- (r*sqrt(n-2))/sqrt(1-r^2)
      p <- 2*stats::pt(abs(t),(n-2), lower.tail=FALSE)
      se <- sqrt((1-r*r)/(n-2))
      out <- list(r, n, t, p, se)
      names(out) <- c("r", "n", "t", "p", "se")
      return(out)
    }
    
    
    # remove gene sets that are zero
    if(!is.null(remove_mods)){uniGS<-uniGS[-which(uniGS %in% remove_mods)]}
    
    # Setup storage
    SubGeneCorDF<-list()
    SubGeneCorMedian<-c()
    SubGenePMedian<-c()
    
    # for each geneset pull the correct gene IDs and indices, calculate correlation matrices, and format results.
    for (i in uniGS){
      curGS<-i
      curIDs<-as.character(moduleSets[which(moduleSets[, mods_var]==curGS), gene_var])
      matchIndex<-match(curIDs, rownames(voomData))
      if (any(is.na(matchIndex))){
        matchIndex<-matchIndex[-which(is.na(matchIndex))]}
      
      if(length(matchIndex)>1){ # run correlation calculations as long as modules contain multiple genes
        setCor<-stats::cor(t(voomData[matchIndex,])) # Calculate pairwise gene correlations from expression data
        n<-t(!is.na(t(voomData[matchIndex,]))) %*% (!is.na(t(voomData[matchIndex,]))) # create matrix with sample size to calculate correlation p (dim=curIDs*curIDs, value = sample size)
        
        allCorInfo<-cor2pvalue(setCor, n) # calculate correlation p values from matrix.
        setP<-allCorInfo$p
        
        SubGeneCorDF[[i]]<- # Assemble df of gene pairs, correlations, and p values
          data.frame(t(utils::combn(colnames(setCor), 2)),
                     allCorVals=setCor[t(utils::combn(colnames(setCor), 2))],
                     allPVals = setP[t(utils::combn(colnames(setP), 2))])%>%
          dplyr::rename(gene1 = X1, gene2 = X2)%>%
          dplyr::mutate(V3=i)} else {print(paste0("Faulty gene-set:", curGS, " only contains 1 gene present in expression data. Skipping."))}
    }
    
    
    # Assemble and format results dataframe
    SubGeneCorDF<-
      dplyr::bind_rows(SubGeneCorDF)%>%
      dplyr::rename(Cor = allCorVals, P = allPVals, Set = V3)%>%
      dplyr::mutate(Cor = as.numeric(as.character(Cor)))%>%
      dplyr::mutate(P = as.numeric(as.character(P)))%>%
      dplyr::mutate(Set = factor(Set, levels =uniGS))%>%
      dplyr::select(Set, gene1, gene2, Cor, P)
    
    #Warning if 0 exist in p-values
    if(min(SubGeneCorDF$P) ==0){
      print(paste(length(SubGeneCorDF$P[SubGeneCorDF$P==0]),
                  "p-values = 0. In plots, these will be replaced with the lowest non-zero value",
                  formatC(sort(unique(SubGeneCorDF$P))[2], digits=2, format="e")))
    }
    
    
    #Boxplot of all pairwise correlations per module
    coherence_boxplot_cor<-
      SubGeneCorDF%>%
      ggplot2::ggplot(ggplot2::aes(y=Cor, x=Set))+
      ggplot2::geom_hline(yintercept =0, color="black")+
      ggplot2::geom_boxplot(outlier.shape = NA)+
      ggplot2::geom_hline(yintercept = r_cutoff, linetype = "dashed", color = "red")+
      ggplot2::scale_y_continuous(breaks=c(seq(-1,1,0.2), r_cutoff))+
      ggplot2::ylab("Correlation Strength")+
      ggplot2::xlab(paste0(mods_title, " Modules"))+
      ggplot2::labs(title = paste0(mods_title, " Module Coherence in ", dat_title, " Samples"),
                    subtitle = "Inter-Gene PearsonCorrelation")+
      ggplot2::theme_bw()+
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust=0.5),
                     plot.subtitle = ggplot2::element_text(hjust=0.5),
                     axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1))
    
    # define plot labels
    labels_df_P<-data.frame(name=rep("negLogP", length(p_cutoff)),
                            h_pos = -log10(p_cutoff)-0.2,
                            h_label = paste("p=",p_cutoff,sep=""))
    
    #Boxplot of all pairwise correlations pvals per module
    coherence_boxplot_p<-
      SubGeneCorDF%>%
      dplyr::mutate(P.est = ifelse(P==0, sort(unique(SubGeneCorDF$P))[2], P)) %>% #Fill in true 0 with lowest P-value in dataset
      ggplot2::ggplot(ggplot2::aes(y=-log10(P.est), x=Set))+
      ggplot2::geom_hline(yintercept =0, color="black")+
      ggplot2::geom_boxplot(outlier.shape = NA)+
      ggplot2::geom_hline(yintercept = -log10(p_cutoff), color="red", linetype="dashed")+
      ggplot2::geom_text(data=labels_df_P, ggplot2::aes(label = h_label, x=3, y=h_pos), color="red")+
      ggplot2::ylab("-log10(Correlation P Value)")+
      ggplot2::xlab(paste0(mods_title, " Modules"))+
      ggplot2::labs(title = paste0(mods_title, " Module Coherence in ", dat_title, " Samples"),
                    subtitle = "Inter-Gene PearsonCorrelation")+
      ggplot2::theme_bw()+
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust=0.5),
                     plot.subtitle = ggplot2::element_text(hjust=0.5),
                     axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1))
    
    # Plot equivalent with facets
    # define plot labels
    labels_df<-data.frame(name=c(rep("Cor", length(r_cutoff)), rep("negLogP", length(p_cutoff))),
                          h_line = c(r_cutoff, -log10(p_cutoff)),
                          h_pos = c(r_cutoff-0.02, -log10(p_cutoff)-0.2),
                          h_label = c(paste("r=",r_cutoff,sep=""), paste("p=",p_cutoff,sep="")))
    
    #Boxplot of all pairwise correlations & pvals per module
    coherence_boxplot_faceted<-
      SubGeneCorDF%>%
      dplyr::mutate(P.est = ifelse(P==0, sort(unique(SubGeneCorDF$P))[2], P)) %>% #Fill in true 0 with lowest P-value in dataset
      dplyr::mutate(negLogP = -log10(P.est))%>%
      tidyr::pivot_longer(cols = c(Cor, negLogP))%>%
      ggplot2::ggplot(ggplot2::aes(y=value, x=Set))+
      ggplot2::geom_hline(yintercept =0, color="black")+
      ggplot2::geom_boxplot(outlier.shape = NA)+
      ggplot2::geom_hline(data=labels_df, ggplot2::aes(yintercept = h_line), color = "red", linetype="dashed")+
      ggplot2::geom_text(data=labels_df, ggplot2::aes(label = h_label, x=3, y=h_pos), color="red")+
      ggplot2::xlab(paste0(mods_title, " Modules"))+
      ggplot2::labs(title = paste0(mods_title, " Module Coherence in ", dat_title, " Samples"),
                    subtitle = "Inter-Gene PearsonCorrelation")+
      ggplot2::theme_bw()+
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        strip.placement = "outside",
        strip.background = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust=0.5),
        plot.subtitle = ggplot2::element_text(hjust=0.5),
        axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1))+
      ggplot2::facet_wrap(~name, scales = "free_y", strip.position = "left",
                          labeller = ggplot2::labeller(name = c("Cor" = "Correlation Strength", "negLogP" = "-log10(Correlation P Value)")))
    
    # Format results for return
    if(return_plot){print(coherence_boxplot_faceted)}
    
    if(length(unlist(GSabsent))>0){
      return(list(coherence_boxplot_combined = coherence_boxplot_faceted,
                  coherence_boxplot_cor = coherence_boxplot_cor,
                  coherence_boxplot_p = coherence_boxplot_p,
                  subgene_correlation_df = SubGeneCorDF,
                  GSabsent = GSabsent))} else
                    
                    return(list(coherence_boxplot_combined = coherence_boxplot_faceted,
                                coherence_boxplot_cor = coherence_boxplot_cor,
                                coherence_boxplot_p = coherence_boxplot_p,
                                subgene_correlation_df = SubGeneCorDF))
    
  }
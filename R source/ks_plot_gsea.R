
ks_plot_gsea <- function(data,#result file from clusterProfiler GSEA analysis, a Large gseaResult
                         term,#terms to display
                         rankColor,#line color
                         labels_size=T#if T, show NES and pvalue as bubbles; otherwise show as numbers
){
  
  require(ggplot2)
  require(patchwork)
  require(clusterProfiler)
  require(viridis)
  require(grDevices)
  
  #clusterprofile
  gseaScores <- getFromNamespace("gseaScores", "DOSE")
  
  # define function
  gsInfo <- function(object, terms) {
    geneList <- object@geneList
    
    gsea_result_table <- object@result
    site=which(gsea_result_table$Description==terms, arr.ind = TRUE)
    genesetid=gsea_result_table$ID[site]
    
    geneSetID <- object@result[genesetid, "ID"]
    
    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    
    df$Description <- object@result[geneSetID, "Description"]
    return(df)
  }
  
  
  
  
  if (length(term) == 1) {
    
    gsea_plotData <- gsInfo(data, terms = term)
    
  } else {
    
    gsea_plotData <- do.call(rbind, lapply(term, gsInfo, object = data))
    
  }
  
  #get gsea rank data
  
  gsea_plotData <- gsea_plotData %>%
    mutate("gene_name" = rep(data@gene2Symbol, length(term))) %>%
    filter(position == 1)
  
  test = gsea_plotData
  
  
  #process gsea rank data for plot
  test$xend = test$x+1
  test$xend2 = ''
  test$xend3 = ''
  
  for (i in 1:nrow(test)) {
    test$xend2[i] = test$xend[i+1] - test$xend[i]
  }
  
  for (i in 1:ncol(test)){
    test[,i][is.na(test[,i])] <- 5
  } 
  
  test$xend2 <- as.numeric(test$xend2)
  
  for (i in 1:nrow(test)) {
    test$xend3[i] = (test$x[i] + test$xend2[i])
  }
  
  test$xend3 <- as.numeric(test$xend3)
  
  
  test$y = rep(1:length(term), as.numeric(table(test$Description)))
  
  
  label_breaks <- seq(1:length(term))+0.5
  
  p=ggplot(test) + 
    # geom_rect(aes(xmin = x-1,xmax = length(data@gene2Symbol), ymin = y-0.01 , ymax = y+0.8, fill=x), lwd=4)+
    # scale_fill_gradientn(colours = colorRampPalette(c("#90191B","white","#3F90C9"))(100))+
    geom_rect(aes(xmin = x,xmax = xend, ymin = y-0.01, ymax = y+0.8, color=x), color=rankColor, size=0.5)+
    scale_y_continuous(expand = c(0,0),labels=term, breaks = label_breaks)+
    scale_x_continuous(expand = c(0,0))+
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 15,colour = 'black'),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  
  
  
  #NES and adjust Pvalue
  gesa_table <- data@result
  gsea_table = gesa_table[gesa_table$Description %in% term, ]
  
  
  if(labels_size==T){
    p1 = ggplot(gsea_table, aes(x=1, y=label_breaks))+
      geom_point(aes(fill=-log(p.adjust),size=setSize), shape=21)+
      scale_fill_gradientn(colors=viridis(256))+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.title= element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.key.height = unit(0.5, "cm"))+
      scale_size(limits = c(0,400), breaks = c(100, 200, 300))
    
    
    
    
    p2 = ggplot(gsea_table, aes(x=1, y=label_breaks))+
      geom_point(aes(fill=NES,size=setSize), shape=21)+
      scale_fill_gradientn(colours = colorRampPalette(c("#3F90C9","white","#90191B"))(100))+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.title= element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.key.height = unit(0.5, "cm"))+
      scale_size(limits = c(0,400), breaks = c(100, 200, 300))
    
  }else{
    
    
    format_number_in_scientific_notation <- function(number) {
      # use sprintf to format numbers in scientific notation
      formatted_number <- sprintf("%e", number)
      return(formatted_number)
    }
    
    
    p1 = ggplot(gsea_table, aes(x=1, y=label_breaks))+
      geom_text(aes(label=format_number_in_scientific_notation(p.adjust)))+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.title= element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank())+
      ggtitle("p.adjust")+
      theme(plot.title = element_text(hjust = 0.5))
    
    
    
    p2 = ggplot(gsea_table, aes(x=1, y=label_breaks))+
      geom_text(aes(label=round(NES,4)))+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.title= element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank())+
      ggtitle("NES")+
      theme(plot.title = element_text(hjust = 0.5))
    
    
    
  }
  
  
  # library(patchwork)
  
  if(labels_size==T){
    
    plot_fig <- p+p1+p2+plot_layout(ncol = 3,width=c(2,0.2,0.2),guides = 'collect')&
      theme(legend.position='bottom')
    
  }else{
    
    plot_fig <- p+p1+p2+plot_layout(ncol = 3,width=c(2,0.5,0.5))
    
  }
  
  
  
  return(plot_fig)
  
}


KS_GSEA <- function(gene,
                    LogFC,
                    analysis=c('GO',"KEGG"),
                    package=c('clusterProfiler','fgsea'),
                    OrgDb=c("org.Hs.eg.db", "org.Mm.eg.db")){
  if(OrgDb=="org.Hs.eg.db"){
    organism = "hsa"
    species = "Homo sapiens"
  }
  
  if(OrgDb=="org.Mm.eg.db"){
    organism = "mmu"
    species = "Mus musculus"
  }
  
  if(package=="clusterProfiler"){
    entrezID <- bitr(gene, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = OrgDb)#convert gene symbol to ID
    
    genelist <- LogFC
    names(genelist) <- gene
    genelist <- genelist[names(genelist) %in% entrezID[,1]]
    names(genelist) <- entrezID[match(names(genelist),entrezID[,1]),2]
    genelist <- sort(genelist,decreasing = T)
    
    #install.packages('R.utils')
    R.utils::setOption( "clusterProfiler.download.method",'auto')
    
    if(analysis=="KEGG"){
      KEGG_gesa <- gseKEGG(geneList = genelist,
                           organism = organism,#query for different species: https://www.genome.jp/kegg/catalog/org_list.html
                           minGSSize = 10,
                           maxGSSize = 500,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           verbose = FALSE,
                           eps = 0)
      
      KEGG_gesa_result = setReadable(KEGG_gesa, OrgDb = OrgDb, keyType = "ENTREZID")
      
      KEGG_gesa_table <- KEGG_gesa_result@result
      write.csv(KEGG_gesa_table, file = './KEGG_gesa_table.csv')
      
      
      return(KEGG_gesa_result)
      
    }
    else{
      GO_gesa <- gseGO(geneList = genelist,
                       ont = "BP",
                       OrgDb=OrgDb,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       verbose = FALSE)
      
      GO_gesa_result = setReadable(GO_gesa,OrgDb = OrgDb, keyType = "ENTREZID")
      
      GO_gesa_table <- GO_gesa_result@result 
      write.csv(GO_gesa_table, file = './GO_gesa_table.csv')
      
      return(GO_gesa_result)
      
    }
    
    
  }
  
  if(package=="fgsea"){
    if(analysis=="KEGG"){
      
      geneset_KEGG = msigdbr(species = species,#mouse:Mus musculus
                             category = "C2", 
                             subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
      geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)#remove KEGG_ prefix
      geneset_KEGG$gs_name <- tolower(geneset_KEGG$gs_name)#convert uppercase to lowercase
      geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)#convert _ to space
      library(Hmisc)
      geneset_KEGG$gs_name <- capitalize(geneset_KEGG$gs_name)#capitalize the first letter
      GSEA_geneset <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
      
    }else{
      
      geneset_GO = msigdbr(species = species,#mouse:Mus musculus
                           category = "C5", 
                           subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
      
      geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#remove KEGG_ prefix
      geneset_GO$gs_name <- tolower(geneset_GO$gs_name)#convert uppercase to lowercase
      geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#convert _ to space
      library(Hmisc)
      geneset_GO$gs_name <- capitalize(geneset_GO$gs_name)#capitalize the first letter
      GSEA_geneset <- geneset_GO %>% split(x = .$gene_symbol, f = .$gs_name)
    }
    
    df <- data.frame(logFC = LogFC, gene = gene)
    df <- df[order(df$logFC,decreasing = T),]
    ranks <- df$logFC
    names(ranks) <- df$gene
    
    ## GSEA analysis
    GSEA_df <- fgsea(pathways = GSEA_geneset, 
                     stats = ranks,
                     minSize=10,
                     maxSize=500,
                     eps=0.0)
    library(data.table)
    fwrite(GSEA_df, file="./GSEA_df.txt", sep="\t", sep2=c("", " ", ""))
    
    return(GSEA_df)
    
  }
  
}

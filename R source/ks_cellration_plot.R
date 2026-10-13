ks_scFraction <- function(group,#single-cell grouping corresponding to cells
                          celltype,#celltype annotation corresponding to cells
                          Total=T, #whether to add total statistics
                          pal,#celltype annotation colors
                          colorLab){#fill colors of the stacked bar plot
  
  require(tidyr)
  require(ggplot2)
  require(dplyr)
  require(aplot)
  
  #count the number of each celltype per sample
  cell_counts = data.frame("Group" = group, "Cells" = celltype) %>% 
    group_by(Group, Cells) %>%
    dplyr::summarize(n = n()) %>% 
    ungroup() %>% 
    pivot_wider(names_from=Cells, values_from=n, values_fill= list(n=0)) %>% 
    as.data.frame()
  
  
  rownames(cell_counts) <- cell_counts[,1]
  cell_counts <- as.matrix(cell_counts[,-1]) %>% t()
  
  #calculate proportions
  cell_ratio <- (cell_counts/rowSums(cell_counts))%>% reshape2::melt()
  colnames(cell_ratio) <- c("cells","group","ratio")
  
  if(Total){
    
    cell_counts <- rbind(cell_counts, colSums(cell_counts))
    rownames(cell_counts)[nrow(cell_counts)] <- "total"
    
    #calculate proportions
    cell_ratio_total <- (cell_counts/rowSums(cell_counts))%>% reshape2::melt()
    colnames(cell_ratio_total) <- c("cells","group","ratio")
    
    total_data <- subset(cell_ratio_total, cells=="total")
    
  }
  
  
  
  #plot
  p <- ggplot(cell_ratio, aes(x=cells, y=ratio,fill=group))+
    geom_bar(stat = "identity", width = 0.9,size = 0.5,colour = '#222222')+
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          panel.background = element_blank(),
          plot.margin = unit(c(0.25,1,0.25,1), "cm"),
          axis.text.x = element_text(angle = 90,  
                                     color = 'black',
                                     hjust = 1,
                                     vjust = 0.5,
                                     margin = margin(0.5,0,0,0,'cm')),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(colour = 'black'),
          axis.text.y = element_text(colour = "black")) +
    scale_y_continuous(expand = c(0, 0))+
    scale_fill_manual(values = colorLab)+
    xlab("")+
    ylab("Proportion of clusters")+
    coord_cartesian(clip = 'off')
  
  
  
  if(Total){
    
    p_total <- ggplot(total_data, aes(x=cells, y=ratio,fill=group))+
      geom_bar(stat = "identity", width = 11,size = 0.5,colour = '#222222')+
      scale_x_discrete(expand = c(0, 0))+
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            panel.background = element_blank(),
            plot.margin = unit(c(0.25,1,0.25,1), "cm"),
            axis.text.x = element_text(angle = 90,  
                                       color = 'black',
                                       hjust = 1,
                                       vjust = 0.5,
                                       margin = margin(0.5,0,0,0,'cm')),
            axis.ticks = element_blank(),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_blank()) +
      scale_y_continuous(expand = c(0, 0))+
      scale_fill_manual(values = colorLab)
    
    p = p+theme(legend.position = "none")
    p_none = annoPoint(object = p,
                       annoPos = 'botomn',
                       xPosition = c(1:length(unique(celltype))),
                       ptSize=1,
                       yPosition=-0.05,
                       pCol = pal)
    
    p1 <- p_none%>%insert_right(p_total, width = 1/length(unique(celltype)))
    
    
    
  }else{
    
    p1 = annoPoint(object = p,
                   annoPos = 'botomn',
                   xPosition = c(1:length(unique(celltype))),
                   ptSize=1,
                   yPosition=-0.05,
                   pCol = pal)
    
  }
  
  
  return(p1)
  
  
  
  
}
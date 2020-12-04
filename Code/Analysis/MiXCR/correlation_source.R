#----------------------------------------------------------------------
#Sarah Russell
#correlation_source.R
#load data functions to make correlation plots with MiXCR + Bisque data
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#compare all data
#----------------------------------------------------------------------

cors <- function(dat) {
  df=merge(cells,dat,by="sample",all=T)
  # turn all three matrices (r, n, and P into a data frame)
  M <- Hmisc::rcorr(as.matrix(select_if(df,is.numeric)),type=c("spearman"))
  # return the three data frames in a list return(Mdf)
  Mdf <- map(M, ~data.frame(.x)) %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA))
    return(Mdf)
}

formatted_cors <- function(data1,data2){
pdf(paste("/Users/sarahrussell/R/plots/", date, "BISQUE_correlations.pdf", sep="_"), width=8,height=6)
#make plot
g1=cors(data1) %>% filter(., grepl(paste(celltype, collapse="|"), measure1) & !grepl(paste(celltype, collapse="|"), measure2)) %>%
      ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
      geom_tile() +
      labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation",
           subtitle="Only significant Spearman's correlation coefficients shown") + scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text() +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0))


g2=cors(data2) %>% filter(., grepl(paste(celltype, collapse="|"), measure1) & !grepl(paste(celltype, collapse="|"), measure2)) %>%
      ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
      geom_tile() +
      labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation",
      subtitle="Only significant Spearman's correlation coefficients shown") + scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text() +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0))


  print(g1)
  print(g2)

dev.off()
}


#----------------------------------------------------------------------
#compare groups of data
#----------------------------------------------------------------------
cors_gr <- function(dat,
                    comparison) {
  if (comparison=="STAGE") {
    df1=merge(cells,dat,by="sample",all=T) %>%
      filter(., sample %in% labels$sample[labels$STAGE=="ADVANCED"])
      gr1="ADVANCED"
    df2=merge(cells,dat,by="sample",all=T) %>%
      filter(., sample %in% labels$sample[labels$STAGE=="LIMITED"])
      gr2="LIMITED"
  } else if (comparison=="SNFClust") {
    df1=merge(cells,dat,by="sample",all=T) %>%
      filter(., sample %in% labels$sample[labels$SNFClust==1])
      gr1="SNFClust1"
    df2=merge(cells,dat,by="sample",all=T) %>%
      filter(., sample %in% labels$sample[labels$SNFClust==2])
      gr2="SNFClust2"
  } else {
    print("error")
  }
  # turn all three matrices (r, n, and P into a data frame)
  M1 <- Hmisc::rcorr(as.matrix(select_if(df1,is.numeric)),type=c("spearman"))
  M2 <- Hmisc::rcorr(as.matrix(select_if(df2,is.numeric)),type=c("spearman"))
  # return the three data frames in a list return(Mdf)
  Mdf1 <- map(M1, ~data.frame(.x))
  Mdf2 <- map(M2, ~data.frame(.x))
  y=list(Mdf1,Mdf2)
  names(y)=c(gr1,gr2)
  return(y)
}


make_plot <- function(df,gr){
  g=df %>% filter(., grepl(paste(celltype, collapse="|"), measure1) & !grepl(paste(celltype, collapse="|"), measure2)) %>%
    ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation", title=gr,
         subtitle="Only significant Spearman's correlation coefficients shown") + scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))
}


formatted_cors_gr <- function(data1,
                              data2,
                          comparison
                          ){
pdf(paste("/Users/sarahrussell/R/plots/", date, comparison, "correlations_groups.pdf", sep="_"), width=8,height=6)
  dat1=cors_gr(data1,comparison)[1]
  df1=dat1[[1]] %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA))

  dat2=cors_gr(data1,comparison)[2]
  df2=dat2[[1]] %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA))

  dat3=cors_gr(data2,comparison)[1]
  df3=dat3[[1]] %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA))

  dat4=cors_gr(data2,comparison)[2]
  df4=dat4[[1]] %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA))

print(make_plot(df1,names(dat1)))
print(make_plot(df2,names(dat2)))
print(make_plot(df3,names(dat3)))
print(make_plot(df4,names(dat4)))
dev.off()
}

shell_dtw<- function(D18o,MgCa,show.dtw.plot=FALSE) {
  
  require(tidyverse)
  require(dtw)
  require(ggplotify)
  
  # if (missing(limit)) {
  #   limit <- 1
  # }
  # else {
  #   limit <- limit
  # }
  
  if(is.null(show.dtw.plot))
    show.dtw.plot <- FALSE
  
  suppressMessages({
    suppressWarnings({
      
      D18o_warp <- D18o %>%
        mutate(d18o=scale(d18o)) %>%
        pull(d18o)
      
      MgCa_warp <- MgCa %>%
        mutate(mg_ca=scale(mg_ca)) %>%
        pull(mg_ca)
      
      alignment <- dtw(D18o_warp*-1,MgCa_warp,
                       keep=TRUE,
                       step.pattern=symmetric2)
      
      
      
      
      
      Aligned_d18o <- D18o %>%
        mutate(x=1:length(D18o$d18o),# you need to change this number "19" to the length of your stable isotope dataset
               proxy="d18o") %>%
        rename(value=d18o) %>%
        select(x,value,proxy)
      
      Aligned_mg_ca <- MgCa %>%
        # filter(dist<5.5) %>%
        pull(mg_ca) %>%
        .[alignment$index2] %>%
        as.data.frame() %>%
        clean_names() %>%
        rename(value=x) %>%
        cbind(x=alignment$index1) %>%
        mutate(proxy="mg_ca")
      
      
      Data_aligned <- rbind(Aligned_d18o,Aligned_mg_ca) %>%
        group_by(x,proxy) %>%
        summarise(value=mean(value)) %>% #this unifies points that were given the same location
        ungroup() %>%
        mutate(value=round(value+0.0001*row_number(),5)) %>% #this adds some error to data points that are similar but far away from each other and prevents them to be grouped as one in the grouping step below
        group_by(proxy,value) %>% # this unifies extra points that were added to fill gaps in the record
        summarise(x=mean(x))
      
      
      lmp <- function (modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
      }
      
      
      
      corr_data <- Data_aligned %>%
        pivot_wider(names_from = "proxy",
                    values_from = "value") %>%
        select(d18o,mg_ca)
      
      LM <- lm(d18o ~ mg_ca,
               data=corr_data,
      )
      
      R2 <- round(summary(LM)[[8]],2)
      
      P_value <- ifelse(lmp(LM)<0.001,"p<0.001",
                        ifelse(lmp(LM)<0.01,"p<0.01",
                               ifelse(lmp(LM)<0.05,"p<0.05","ns!!!")))
      
      Stats <- paste0("R²=", R2,", ",P_value)
      
      P1 <- Data_aligned %>%
        pivot_wider(names_from = "proxy",
                    values_from = "value") %>%
        ggplot() +
        aes(x = d18o,
            y = mg_ca) +
        geom_point() +
        # xlim(0.5, 1.8) +
        # ylim(0.7, 1) +
        geom_smooth(method = "lm",
                    se = FALSE,formula=y ~ x,
                    col="black",
                    linetype="dashed",
                    size=0.5) +
        ylab("Mg/Ca Ratio") +
        xlab(expression("δ" ^ "18" * "O [VPDB]"))+
        ggtitle("Correlation of proxies",#"Correlation graph for C16 d18O and Mg/Ca",
                subtitle = Stats)
      
      
      if (show.dtw.plot==TRUE) {
        print("show.dtw.plot does not work yet")
        # as.ggplot(~plot(alignment,
        #                 type = "twoway",
        #                 main="DTW Plot",
        #                 adj = 0,
        #                 line=1,
        #                 offset = 10)
        # )/P1
      }
      else {
        P1
      }
      
      
    })
  })
  
  
  
  
  
  
}

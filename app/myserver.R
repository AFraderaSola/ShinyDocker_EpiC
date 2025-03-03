############################
######## Libraries #########
############################

library(shinythemes)
library(mailtoR)
library(DT)
library(shinyBS)
library(plotly)
library(tidyverse)
library(rMQanalysis)
set.seed(666)

############################
####### Variables ##########
############################

load("./InputFiles/00_DNA_5mCG.RData")

############################
####### Functions ##########
############################

Enriched2 <- function (x, y = NULL, c = 1, s0 = 1, pvalue = 0.05, ...) 
{
  if (class(x) == "data.frame") {
    if (length(x) == 2) {
      y <- x[[2]]
      x <- x[[1]]
    }
    else {
      stop("x hast to be a data frame with length 2!")
    }
  }
  results <- list()
  results$curve <- createEnrichmentCurve2(c, s0, pvalue, ...)
  results$positive <- (c/(x - s0) - y + -log10(pvalue) <= 0 & 
                         x >= s0)
  results$enriched <- results$positive 
  results$background <- !results$enriched
  class(results) <- append(class(results), "Enriched")
  return(results)
}

createEnrichmentCurve2 <-  function (c, s0, pvalue = 1, step = 0.01, xlim = 50, linetype = 2, 
                                     size = 0.1) 
{
  x <- seq(s0, xlim, step)
  y <- c/(x - s0) + -log10(pvalue)
  return(annotate("path", x = c(NA, x), y = c(NA, y), 
                  linetype = linetype, size = size))
}

createLink <- function(val) {
  sprintf('<a href="https://www.uniprot.org/uniprotkb?query=%s" target="_blank" class="btn btn-link" role="button">UniProt</a>',val)
}

############################
######### Script ###########
############################

server <- function(input, output, session) {
  
  #### Data functions ####
  
  # Modify the different data sets in function of
  # user's input
  
  ##### DNA Modifications ####
  
  #### 5mCG ####
  
  ## Quantified Proteins (QP) data table (DT)
  
  DNA_5mCG_QP_DT <- reactive({

    QP_DT <- DNA_5mCG_quantified_proteins_df[[grep(pattern = gsub(pattern = " ",
                                                                  replacement = "_",
                                                                  x = input$DNA_5mCG_species),
                                                   x = names(DNA_5mCG_quantified_proteins_full_df))]]
    QP_DT

  })
  
  ## Full (F) data table (DT)
  
  DNA_5mCG_F_DT <- reactive({
    
    QP_DT <- DNA_5mCG_quantified_proteins_full_df[[grep(pattern = gsub(pattern = " ",
                                                                  replacement = "_",
                                                                  x = input$DNA_5mCG_species),
                                                   x = names(DNA_5mCG_quantified_proteins_full_df))]] 
    QP_DT
    
  })
  
  #### Output functions ####
  
  #### DNA Modifications #### 
  
  #### 5mCG ####
  
  ## Quantified Proteins (QP) data table (DT)
  
  output$DNA_5mCG_QP_DT <- renderDataTable({
    
    df <- DNA_5mCG_QP_DT()
    
    df$Database <- createLink(df$`Protein ID`)

    cols_to_round <- colnames(df)[grep(pattern = "Mean|log2|p-value", colnames(df))]

    df <- df %>%
      datatable(options = list(pageLength = 15, autoWidth = TRUE),
                rownames = F,escape = FALSE) %>%
      formatRound(columns = cols_to_round, digits = 3)

    return(df)
    
  })
  
  DNA_5mCG_QP_DT_proxy <- dataTableProxy("DNA_5mCG_QP_DT")
  
  ## Overview barplot
  
  output$DNA_5mCG_Overview_Barplot <- renderPlotly({
    
    overview_bp_df <- lapply(DNA_5mCG_quantified_proteins_full_df, function(df_list){
      
      df_list %>% 
        group_by(enriched_mCG_CG) %>%
        tally()
      
    })
    
    overview_bp_df <- bind_rows(overview_bp_df, .id = "Species")
    
    overview_bp_df$Species <- gsub(pattern = "_",replacement = " ",x = overview_bp_df$Species)
    
    overview_bp_df <- overview_bp_df %>%
      group_by(Species) %>%
      mutate(N = sum(n)) %>%
      arrange(N)
    
    overview_bp_df$Species <- factor(x = overview_bp_df$Species, levels = unique(overview_bp_df$Species))
    
    overview_bp_df$enriched_mCG_CG <- gsub(pattern = "TRUE",replacement = "Enriched",x = overview_bp_df$enriched_mCG_CG)
    
    overview_bp_df$enriched_mCG_CG <- gsub(pattern = "FALSE",replacement = "Not enriched",x = overview_bp_df$enriched_mCG_CG)
    
    overview_bp_df$enriched_mCG_CG <- factor(x = overview_bp_df$enriched_mCG_CG, levels = rev(unique(overview_bp_df$enriched_mCG_CG)))
    
    DNA_5mCG_Overview_Barplot <- overview_bp_df %>%
                                  ggplot(aes(x=Species, y=n, fill = enriched_mCG_CG)) +
                                  geom_bar(stat="identity", width=.6) +
                                  coord_flip()+
                                  scale_fill_manual(labels = c("Enriched", "Not enriched"), values = c("#1B9E77", "#7570B3")) +
                                  labs(fill = "")+
                                  theme_minimal()+
                                  ylab("Quantified Proteins")+
                                  theme(legend.text = element_text(size = 12),
                                        legend.key.size = unit(1, 'cm'),
                                        legend.position = "top")+
                                  theme(axis.text=element_text(size=12),
                                        axis.text.y = element_text(face = "italic"),
                                        axis.title=element_text(size=14,face="bold"))
    
    print(DNA_5mCG_Overview_Barplot)
    
  })
  
  ## Overview barplot
  
  output$DNA_5mCG_Volcanoplot <- renderPlotly({

    #Input variables
    
    ID <- input$DNA_5mCG_QP_DT_ui_selectedIDs
    enrich_s0 <- log2(2^as.numeric(input$DNA_5mCG_LFC))
    enrich_pvalue <- as.numeric(input$DNA_5mCG_Pvalue)

    # Fixed variables
    
    EpiMod <- "DNA_<sup>5</sup>mCG"
    color_intercept <- "#56B4E9"
    exp_a <- "mCG"
    exp_b <- "CG"
    difference <- paste('difference', exp_a, exp_b, sep='_')
    pvalue <- paste('log10.pvalue', exp_a, exp_b, sep='_')
    val_count_a <- paste0('value_count_',exp_a)
    val_count_b <- paste0('value_count_',exp_b)
    enriched_col <- paste('enriched', exp_a, exp_b, sep='_')
    techname <- c("CG", "mCG")
    nicename <- c("CG", "<sup>5</sup>mCG")
    expdf <- data.frame(techname, nicename)

    min_quant_events <- 0
    replicate_count <- 4
    enrich_c <- .05
    volcano_threshold_linetype <- 2
    volcano_threshold_linesize <- 0.5
    difference_volcano_column <- 'mean2_naimp_meas_'
    volcano_column_type <- sub('([^[:digit:]]+).*', '\\1', difference_volcano_column)

    pg_quant <- DNA_5mCG_F_DT()
    
    pg_quant$highlight <- FALSE
    
    if (length(ID) != 0) {
      
      ids_to_highlight <- gsub(pattern = "\\|",replacement = "\\\\\\|", x = ID)
      
      ids_searchstring <- paste(ids_to_highlight, collapse='|')
      
      matching_rows <- grep(ids_searchstring, as.character(pg_quant$my_label))
      
      pg_quant$highlight[matching_rows] <- TRUE
      
      
    }
    
    
    ## Filter count values. Standard set to 0, so they can check the dot all along the time course.
    ## In Mario's standard script is usually set to 2. But in this case, since it shown in the
    ## box plot and in the dot plot, would be confusing not to show it here as well.
    ## Can be adjusted on the fixed variables section of this function.
    
    filter_min_quant_events <-
      pg_quant[val_count_a] >= min_quant_events |
      pg_quant[val_count_b] >= min_quant_events
    
    pg_quant <- pg_quant[filter_min_quant_events,]
    
    ## Select and highlight the enriched dots.
    
    my_enriched <- Enriched2(pg_quant[c(difference, pvalue)],
                            c=enrich_c,
                            s0=enrich_s0,
                            pvalue=enrich_pvalue,
                            linetype=volcano_threshold_linetype,
                            size=volcano_threshold_linesize)
    pg_quant$enriched <- as.factor(my_enriched$enriched)
    
    subtitle <- sprintf('enriched: %d; background: %d; total: %d',
                        sum(my_enriched$enriched), sum(my_enriched$background),
                        nrow(pg_quant))
    
    
    my_points_nohighlight <- interaction(pg_quant[!pg_quant$highlight, c('enriched', 'highlight')])
    my_points_highlight <- interaction(pg_quant[pg_quant$highlight, c('enriched', 'highlight')])
    
    niceexp_a <- expdf[expdf$techname == exp_a,]$nicename
    niceexp_b <- expdf[expdf$techname == exp_b,]$nicename
    
    ## Generate a user friendly plotly string for the hoover.
    
    plotly_string <- paste0(sprintf('<b>Protein ID</b>: %s<br /><b>log<sub>2</sub>(Fold change):</b> %.2f<br /><b>-log<sub>10</sub>(pvalue):</b> %.2f<br /><b>Value counts:</b> %s; %s',
                                    sub('(.{30}).*', '\\1...', pg_quant$my_label),
                                    pg_quant[[difference]], pg_quant[[pvalue]],
                                    paste(pg_quant[[paste0('value_count_', exp_a)]], 'in', niceexp_a), paste(pg_quant[[paste0('value_count_', exp_b)]], 'in', niceexp_b)
                                    ),
                            switch('Gene.names' %in% names(pg_quant), 
                                   sprintf('<br />Gene.names: %s', pg_quant$Gene.names), 
                                   NULL),
                            switch('Protein.names' %in% names(pg_quant), 
                                   sprintf('<br />Protein.names: %s', pg_quant$Protein.names), 
                                   NULL),
                            switch('description' %in% names(pg_quant), 
                                   sprintf('<br />description: %s', pg_quant$description), 
                                   NULL))
    
    pg_quant$my_text <- plotly_string
    
    ## Subset the data for the volcano plot.
    
    pg_quant_subset <-
      pg_quant %>%
      rowwise() %>%
      mutate(imputation=ifelse(max(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count &
                                 min(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count,
                               'measured',
                               ifelse(min(range(.data[[val_count_a]], .data[[val_count_b]])) > 0,
                                      'some imputed',
                                      'ref imputed'))) %>%
      select(Protein.IDs, my_label, matches(difference), matches(pvalue),
             my_text, highlight, enriched, imputation)
    
    ## Plot
    
    ## Select maximum x and y values of your volcano plot.
    
    max_x_value <- max(abs(pg_quant[grep('^difference', names(pg_quant))]),
                       na.rm=TRUE)
    max_y_value <- max(c(unlist(pg_quant[grep('^log10.pvalue', names(pg_quant))]),
                         10^enrich_pvalue + 0.5))
    
    volcanoplot <- 
      ggplot(pg_quant_subset, aes_string(x=difference, y=pvalue, text='my_text')) +   
      geom_hline(yintercept=0, color= color_intercept, na.rm=TRUE) +
      geom_vline(xintercept=0, color= color_intercept, na.rm=TRUE) +
      plot(my_enriched) + # adding the threshold line of our Enriched object
      geom_point(data=subset(pg_quant_subset, !highlight),
                 aes(alpha=my_points_nohighlight, 
                     color=my_points_nohighlight,
                     fill=my_points_nohighlight, 
                     # shape=imputation,
                     size=my_points_nohighlight), 
                 na.rm=TRUE) + 
      geom_point(data=subset(pg_quant_subset, highlight),
                 aes(alpha=my_points_highlight, 
                     color=my_points_highlight,
                     fill=my_points_highlight, 
                     # shape=imputation,
                     size=my_points_highlight), 
                 na.rm=TRUE) + 
      # we have combined ENRICHED.HIGHLIGHT
      scale_alpha_manual(values=c(`TRUE.FALSE`=.7, `FALSE.FALSE`=.4, `TRUE.TRUE`=.7, `FALSE.TRUE`=.7),
                         guide="none") + 
      scale_color_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#008837', `FALSE.TRUE`='#7B3294'),
                         guide="none") +  
      scale_fill_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#008837', `FALSE.TRUE`='#7B3294'),
                        guide="none")+
      # scale_shape_manual(values=c(`measured`=21, `ref imputed`=25, `some imputed`=24))+
      scale_size_manual(values=c(`TRUE.FALSE`=1, `FALSE.FALSE`=.5, `TRUE.TRUE`=1.2, `FALSE.TRUE`=1),
                        guide="none")
    
    volcanoplot <- volcanoplot +
      ggtitle(sprintf('%s experiment\nVolcano plot of %s against %s\n%s',gsub(pattern = "_",replacement = " ",x = EpiMod), niceexp_a, niceexp_b, subtitle)) +
      coord_cartesian(xlim=c(-max_x_value,max_x_value), ylim=c(-.8,max_y_value)) +
      ylab("-log<sub>10</sub>(p-value)") + 
      xlab("log<sub>2</sub>(Fold change)") +
      annotate('text', label=niceexp_a, x=max_x_value / 2, y=-.4,size = 6) + 
      annotate('text', label=niceexp_b, x=-max_x_value / 2, y=-.4,size = 6)+
      theme_minimal()+
      # theme(panel.grid.major = element_line(linetype = 'solid', colour = colour), 
      #       panel.grid.minor = element_line(linetype = 'solid',colour = colour))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
    
    DNA_5mCG_Volcanoplot <-
      plotly::ggplotly(volcanoplot,tooltip='text') %>%
      plotly::layout(showlegend = FALSE)
    
    print(DNA_5mCG_Volcanoplot)

  })
  
  #### Observe functions ####
  
  #### DNA Modifications #### 
  
  #### 5mCG ####
  
  observe({
    updateSelectizeInput(
      session,
      'DNA_5mCG_QP_DT_ui_selectedIDs',
      choices=DNA_5mCG_QP_DT()[['Protein ID']][as.numeric(input$DNA_5mCG_QP_DT_rows_selected)],
      selected=DNA_5mCG_QP_DT()[['Protein ID']][as.numeric(input$DNA_5mCG_QP_DT_rows_selected)]
    )
  })

  observeEvent(input$update_DNA_5mCG_QP_DT, {
    rows <-
      match(input$DNA_5mCG_QP_DT_ui_selectedIDs,
            DNA_5mCG_QP_DT()[['Protein ID']])
    selectRows(DNA_5mCG_QP_DT_proxy,
               as.numeric(rows))
  })
  
}


#SECTION 1: Import, Clean, and Prep Data ---------------------------------------

#Install packages
  install.packages('tidyverse')
  install.packages('vegan')
  install.packages('ggplot2')
  install.packages('dplyr')
  install.packages('ggridges') 
  install.packages('readxl')
  install.packages("hrbrthemes")
  install.packages("systemfonts")
  install.packages("gridExtra")
  install.packages("grid")
  install.packages("abind")
  install.packages("here") 
  install.packages("openxlsx")
  install.packages("tictoc")
  install.packages("lme4")
  install.packages("lmerTest") 
  install.packages("emmeans") 
  install.packages("pbkrtest")
  install.packages("viridis") 
  install.packages("reshape2")
  install.packages("betapart")

#call packages
  library(tidyverse)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(ggridges) 
  library(readxl) 
  library(systemfonts)
  #library(hrbrthemes)
  library(gridExtra)
  library(grid)
  library(abind) 
  library(here)
  library(openxlsx)
  library(tictoc) 
  library(lme4)
  library(lmerTest) 
  library(emmeans) 
  library(pbkrtest) 
  library(viridis)
  library(reshape2)
  library(betapart)

#Import, prep data
  #set working directory
  setwd(here())
  
  ##1.1: Seeds Data ------------------------------------------------------------
  #Import data
    WideCounts <- read_excel(here("Data/All_Data_Seeds-7_26.xlsx"), sheet = "COUNTS")
    WideWeights <- read_excel(here("Data/All_Data_Seeds-7_26.xlsx"), sheet = "WEIGHTS")
    
  #Add column combining Site and Quad
    WideCounts$SiteQuad <- paste(WideCounts$Site, WideCounts$Quad, sep = "_")
    WideCounts <- WideCounts |> 
      relocate(SiteQuad, .before = Site)
  #remove unwanted data 
    #remove NA site (probably my averages line I forgot to delete from the excel)
    WideCounts <- WideCounts |> 
      filter_out(SiteQuad == "NA_NA") 
    #remove unknowns and other unneeded columns 
    WideCounts <- WideCounts[, !grepl("Unknown|Uknown|Notes|134", names(WideCounts))]
    
    
  #Flip data
    LongCounts <- gather(data = WideCounts, 
                         key = "Species", 
                         value = "Count", "Acer_platanoides_count":"Helianthus_annus_count")
    LongWeights <- gather(data = WideWeights, 
                         key = "Species", 
                         value = "Weight", "Acer_platanoides_weight":"Helianthus_annus_weight")
    
  #expand long data format so every seed observation is its own line
    #convert count to integer 
    LongCounts$Count <- as.integer(LongCounts$Count)
    #convert NAs in count  
    LongCounts$Count[is.na(LongCounts$Count)] <- 0
    #expand counts 
    ExpCounts <- uncount(LongCounts, weights = Count)
    
    LongCounts <- LongCounts |> 
      mutate(SiteQuad = case_when(
        SiteQuad == "RBP_101" ~ "RBP_0101",
        SiteQuad == "RBP_103" ~ "RBP_0103",
        SiteQuad == "RBP_105" ~ "RBP_0105",
        SiteQuad == "RBP_107" ~ "RBP_0107",
        SiteQuad == "RBP_300" ~ "RBP_0300",
        SiteQuad == "RBP_302" ~ "RBP_0302",
        SiteQuad == "TFO_5" ~ "TFO_0005",
        SiteQuad == "TFO_102" ~ "TFO_0102",
        SiteQuad == "TFO_107" ~ "TFO_0107",
        SiteQuad == "TFO_112" ~ "TFO_0112",
        SiteQuad == "TFO_117" ~ "TFO_0117",
        SiteQuad == "TFO_15" ~ "TFO_0015",
        SiteQuad == "TFO_17" ~ "TFO_0017",
        SiteQuad == "TFO_204" ~ "TFO_0204",
        SiteQuad == "TFO_209" ~ "TFO_0209",
        SiteQuad == "TFO_214" ~ "TFO_0214",
        SiteQuad == "TFO_306" ~ "TFO_0306",
        SiteQuad == "TFO_311" ~ "TFO_0311",
        SiteQuad == "TFO_408" ~ "TFO_0408",
        SiteQuad == "TFO_413" ~ "TFO_0413",
        SiteQuad == "TFO_418" ~ "TFO_0418",
        SiteQuad == "TFO_505" ~ "TFO_0505",
        SiteQuad == "TFO_510" ~ "TFO_0510",
        SiteQuad == "TFO_607" ~ "TFO_0607",
        SiteQuad == "TFO_709" ~ "TFO_0709",
        SiteQuad == "TFO_908" ~ "TFO_0908",
        TRUE ~ SiteQuad
      )) |> 
      filter(!SiteQuad %in% c("TFO_0418", "TFO_0017"))
    
  ##1.2: Canopy Loss Data 
    #import data
    CanopyLoss <- read_excel(here("Data/Canopy_Loss.xlsx"), sheet = "merged")

#SECTION 2: Goofin' around -----------------------------------------------------
    
Bar_Data <- ExpCounts |> 
      count(Species)
    
All_Bar <- ggplot(Bar_Data, aes(x = Species, y = n)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10()

#diversity stuff
WCna <- WideCounts |> 
  mutate(across(12:107, as.numeric)) |> 
  #use mutate across everything to turn all instances of "NA" into a 0
  mutate(across(c(Acer_platanoides_count:Helianthus_annus_count), ~replace_na(., 0))) 
  #convert all data columns into numeric 

div <- diversity(WCna[,12:107], index = "shannon")
even <- div/log(specnumber(WCna[,12:107]))

WCna$div <- div
WCna$even <- even

WCna.div <- WCna |> 
  group_by(Set, SiteQuad) |> 
  summarise(
    diversity = mean(div),
    .groups = "drop"
  )

WC.div.avg <- WCna.div |> 
  group_by(Set) |> 
  summarize(
    avg_div = mean(diversity),
    .groups = "drop"
  )

ggplot(WC.div.avg, aes(x = Set, y = avg_div)) +
  geom_point() +
  geom_line()

#try some beta-div stuff? 
#   #reset WCna
#       WCna <- WideCounts |> 
#         mutate(across(12:107, as.numeric)) |> 
#         #use mutate across everything to turn all instances of "NA" into a 0
#         mutate(across(c(Acer_platanoides_count:Helianthus_annus_count), ~replace_na(., 0))) 
# #convert all data columns into numeric 
#   #overall beta-div of dataset: 
#   #NOTE: Using most basic definition of beta-div here: beta = gamma/alpha
#     #find gamma-div
#       spTotals <- colSums(WCna[, 12:ncol(WCna)], na.rm = TRUE)
#       gammaDiv <- diversity(spTotals, index = "shannon")
#     #find alpha-div (div of each quad) 
#       div <- diversity(WCna[,12:107], index = "shannon")
#       even <- div/log(specnumber(WCna[,12:107]))
#       
#       WCna$div <- div
#       WCna$even <- even
#       
#       WCna.betadiv <- WCna |> 
#         group_by(SiteQuad) |> 
#         summarise(
#           alpha_div = mean(div),
#           .groups = "drop"
#         )
#     #add beta-div column 
#       WCna.betadiv <- WCna.betadiv |> 
#         mutate(beta_div = gammaDiv/alpha_div)
#       
#     #bar plot 
#       betaBar <- ggplot(WCna.betadiv, aes(x = SiteQuad, y = beta_div)) + 
#         geom_col() +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
#   #another beta div calculation 
#     #make community matrix with site/date as row identifiers
#       comm_matrix <- dcast(LongCounts, SiteQuad + Set ~ Species, 
#                            value.var = "Count", 
#                            fun.aggregate = sum,
#                            fill = 0)
#       
#     #separate metadata and species data
#       meta <- comm_matrix[, c("SiteQuad", "Set")] 
#       species_mat <- comm_matrix[, -(1:2)] 
#       
#     
#     #turnover/nestedness partitioning 
#       #convert abundance to presence/absence bc this turn/nest partitioning 
#       #needs incidence  
#       
#       #we also gotta remove rows with only 0s because that doesn't work for BC diss 
#         #there's only 52 of the 1000+ rows with 0s so it's not a huge deal
#         keep <- rowSums(species_mat) > 0
#         comm_mat_clean <- rowSums(comm_matrix) > 0
#         spmat_clean <- species_mat[keep, ]
#         meta_clean <- meta[keep, ]
#         
#         pa_mat <- ifelse(spmat_clean > 0, 1, 0)
#       #calculate turnover, nestedness, and Jaccard dissimilarity (combo of both) 
#         beta_parts <- beta.pair(pa_mat, index.family = "jaccard")
#         
#       #Bray-Curtis dissimilarity on species matrix 
#         bc_dist <- vegdist(spmat_clean, method = "bray")
#     
#     #nMDS time
#       # nmds <-  metaMDS(spmat_clean, distance = "bray")
#       # plot(nmds)   
#       #(that didn't work at all lmao) 
#       
#   #permanova 
#     adonis2(bc_dist ~ SiteQuad * Set, data = meta_clean, by = NULL, permutations = 99)
#     
#     write.xlsx(comm_matrix, here("comm_matrix.xlsx")) 
    

#SECTION 3: Beta-Analysis and Canopy Cover -------------------------------------
  ##3.1: get pairwise beta between all subplots --------------------------------
    
    LCsum <- LongCounts                                                |> 
      select("Date_Collected", "SiteQuad", "Plot", "Species", "Count") |> 
      rename(subplot = Plot)                                           |> 
      group_by(Species, SiteQuad, subplot)                             |> 
      summarise(
        Count   = sum(Count), 
        .groups = "drop")                                              |> 
      rename(plot    = SiteQuad, 
             species = Species, 
             count   = Count)
    
  #rename following plots: 
  #RBP_101, RBP_103, RBP_105, RBP_107, RBP_300, RBP_302, TFO_102
  
 

  LCsum$count <- as.integer(LCsum$count)
  
    # make a function to compute pairwise dissimilarities for each plot
   beta_subplot <- function(df) {
     
     #pivot plot data to subplot x species matrix
     comm <- df                        |> 
       select(subplot, species, count) |> 
       pivot_wider(names_from = species, 
                   values_from = count, 
                   values_fill = 0)    |> 
       column_to_rownames("subplot")
       
      #pairwise Bray-Curtis among 5 subplots
       bc <- vegdist(comm, method = "bray")
       
       #convert dist object to tidy long format df
       bc_df <- as.matrix(bc)                    |> 
         as.data.frame()                         |> 
         rownames_to_column("subplot_1")         |> 
         pivot_longer(-subplot_1, 
                      names_to = "subplot_2", 
                      values_to = "bray_curtis") |> 
      filter(subplot_1 < subplot_2) #keeps only unique pairs, dropping diagonal/dupes
       
      bc_df
       
   }
   
  pairwise_beta_results <- LCsum                                   |> 
     group_by(plot)                                                |> 
    #create new data frame for each plot
     group_split()                                                 |> 
    #run the function on every subplot individually
     map_dfr(~ beta_subplot(.x) |> mutate(plot = unique(.x$plot))) |> 
     relocate(plot)
    
  ##3.2: get delta in canopy cover for every pair of subplots ------------------
    
    CanopyLoss <- CanopyLoss                  |> 
    rename(Percent_Lost = `canopy loss`)      |> 
    select(Percent_Lost, Site, Plot, Subplot) |> 
    mutate(SitePlot = paste(Site, Plot, sep = "_"))
  
    pairwise_canopy <- CanopyLoss                |> 
      rename(subplot1     = Subplot, 
             canopy_loss1 = Percent_Lost)        |> 
      inner_join(
        CanopyLoss                               |> 
          rename(subplot2     = Subplot, 
                 canopy_loss2 = Percent_Lost),
        by = "SitePlot"
      )                                          |> 
      filter(subplot1 < subplot2)                |> 
      mutate(diff = canopy_loss1 - canopy_loss2) |> 
      select(SitePlot, subplot1, subplot2, 
             canopy_loss1, canopy_loss2, diff)   |> 
      mutate(diff = abs(diff))
  
    join_canopy <- pairwise_canopy                            |> 
      mutate(plot_pair = paste(subplot1, subplot2, sep = "")) |> 
      select(SitePlot, plot_pair, diff)
  
  ##3.3: add canopy diff data to beta data -------------------------------------
      All_Pairs <- pairwise_beta_results |> 
        mutate(canopy_diff = pairwise_canopy$diff)
    
  ByPairDiff.plot <-  ggplot(All_Pairs, aes(x = canopy_diff, y = bray_curtis)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = "Difference in Canopy Loss (%)", 
      y = "Pairwise Beta-Diversity (Bray-Curtis Dissimilarity)"
    )
    
    View(LCsum)
 
  LCsum <- LCsum |> 
    mutate(plotsp = paste(plot, subplot, sep = "_")) |> 
    select(plotsp, species, count)
       
  WCsum <- pivot_wider(
    data        = LCsum, 
    names_from  = species, 
    values_from = count, 
    values_fill = 0
  ) 
  
  CanopyLoss <- CanopyLoss |> 
    mutate(plotsp = paste(SitePlot, Subplot, sep = '_'))
  
  WCsum <- WCsum |> 
    separate(plotsp, into = c("site", "plot", "subplot"), sep = "_", remove = FALSE) |> 
    mutate(Percent_Loss = CanopyLoss$Percent_Lost) |> 
    relocate(Percent_Loss, .after = subplot)
  
  write.xlsx(WCsum, here("primere_output_1.xlsx")) 
  
  ## 3.4: repeat 3.2 but for average canopy loss between the two --------------- 
  
  canopy_pavg <- CanopyLoss |> 
    rename(subplot1     = Subplot, 
           canopy_loss1 = Percent_Lost)        |> 
    inner_join(
      CanopyLoss                               |> 
        rename(subplot2     = Subplot, 
               canopy_loss2 = Percent_Lost),
      by = "SitePlot"
    )                                          |> 
    filter(subplot1 < subplot2)                |> 
    mutate(avg_loss = (canopy_loss1 + canopy_loss2)/2) |> 
    select(SitePlot, subplot1, subplot2, 
           canopy_loss1, canopy_loss2, avg_loss)
  
  join_canopy <- pairwise_canopy                            |> 
    mutate(plot_pair = paste(subplot1, subplot2, sep = "")) |> 
    select(SitePlot, plot_pair, diff)
    
  All_Pairs <- pairwise_beta_results |> 
    mutate(canopy_avg = canopy_pavg$avg_loss)
    
  ByPairAverage.plot <-  ggplot(All_Pairs, aes(x = canopy_avg, y = bray_curtis)) + 
    geom_point() + 
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      x = "Average Canopy Loss (%)", 
      y = "Pairwise Beta-Diversity (Bray-Curtis Dissimilarity)"
    )
  
  NMDS_ord <- read_excel(here("Data/NMDS_ordination.xlsx"))
  
  View(CanopyLoss)
  
  NMDSxCLoss.df <- data.frame(
    nmds1 = NMDS_ord$nmds1, 
    nmds2 = NMDS_ord$nmds2,
    PercentLost = CanopyLoss$Percent_Lost
  )  
  
  NMDS1vsCLoss.plot <-  ggplot(NMDSxCLoss.df, aes(x = PercentLost, y = nmds1)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = "Canopy Loss (%)", 
      y = "NMDS1"
    )
  
  NMDS2vsCLoss.plot <-  ggplot(NMDSxCLoss.df, aes(x = PercentLost, y = nmds2)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = "Canopy Loss (%)", 
      y = "NMDS2"
    )
    
    
    
    
    

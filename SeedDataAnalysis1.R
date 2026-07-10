
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
  #Import data
    WideCounts <- read_excel(here("All_Data_Seeds-7_26.xlsx"), sheet = "COUNTS")
    WideWeights <- read_excel(here("All_Data_Seeds-7_26.xlsx"), sheet = "WEIGHTS")
    
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
  #reset WCna
      WCna <- WideCounts |> 
        mutate(across(12:107, as.numeric)) |> 
        #use mutate across everything to turn all instances of "NA" into a 0
        mutate(across(c(Acer_platanoides_count:Helianthus_annus_count), ~replace_na(., 0))) 
#convert all data columns into numeric 
  #overall beta-div of dataset: 
  #NOTE: Using most basic definition of beta-div here: beta = gamma/alpha
    #find gamma-div
      spTotals <- colSums(WCna[, 12:ncol(WCna)], na.rm = TRUE)
      gammaDiv <- diversity(spTotals, index = "shannon")
    #find alpha-div (div of each quad) 
      div <- diversity(WCna[,12:107], index = "shannon")
      even <- div/log(specnumber(WCna[,12:107]))
      
      WCna$div <- div
      WCna$even <- even
      
      WCna.betadiv <- WCna |> 
        group_by(SiteQuad) |> 
        summarise(
          alpha_div = mean(div),
          .groups = "drop"
        )
    #add beta-div column 
      WCna.betadiv <- WCna.betadiv |> 
        mutate(beta_div = gammaDiv/alpha_div)
      
    #bar plot 
      betaBar <- ggplot(WCna.betadiv, aes(x = SiteQuad, y = beta_div)) + 
        geom_col() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

  #another beta div calculation 
    #make community matrix with site/date as row identifiers
      comm_matrix <- dcast(LongCounts, SiteQuad + Set ~ Species, 
                           value.var = "Count", 
                           fun.aggregate = sum,
                           fill = 0)
      
    #separate metadata and species data
      meta <- comm_matrix[, c("SiteQuad", "Set")] 
      species_mat <- comm_matrix[, -(1:2)] 
      
    
    #turnover/nestedness partitioning 
      #convert abundance to presence/absence bc this turn/nest partitioning 
      #needs incidence  
      
      #we also gotta remove rows with only 0s because that doesn't work for BC diss 
        #there's only 52 of the 1000+ rows with 0s so it's not a huge deal
        keep <- rowSums(species_mat) > 0
        comm_mat_clean <- rowSums(comm_matrix) > 0
        spmat_clean <- species_mat[keep, ]
        meta_clean <- meta[keep, ]
        
        pa_mat <- ifelse(spmat_clean > 0, 1, 0)
      #calculate turnover, nestedness, and Jaccard dissimilarity (combo of both) 
        beta_parts <- beta.pair(pa_mat, index.family = "jaccard")
        
      #Bray-Curtis dissimilarity on species matrix 
        bc_dist <- vegdist(spmat_clean, method = "bray")
    
    #nMDS time
      # nmds <-  metaMDS(spmat_clean, distance = "bray")
      # plot(nmds)   
      #(that didn't work at all lmao) 
      
  #permanova 
    adonis2(bc_dist ~ SiteQuad * Set, data = meta_clean, by = NULL, permutations = 99)
    
    write.xlsx(comm_matrix, here("comm_matrix.xlsx")) 
    
    

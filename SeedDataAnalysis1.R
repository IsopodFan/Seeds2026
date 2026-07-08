
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
  #Flip data
    LongCounts <- gather(data = WideCounts, 
                         key = "Species", 
                         value = "Count", "Acer_platanoides_count":"Unknown_40_count")
    LongWeights <- gather(data = WideWeights, 
                         key = "Species", 
                         value = "Weight", "Acer_platanoides_weight":"Unknown_40_weight")
    
  #expand long data format so every seed observation is its own line
    #convert count to integer 
    LongCounts$Count <- as.integer(LongCounts$Count)
    #convert NAs in count  
    LongCounts$Count[is.na(LongCounts$Count)] <- 0
    #expand counts 
    LongCounts <- uncount(LongCounts, weights = Count)
    
  #don't want: unks 6, 10, 12, 13, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 28
  #remove unknowns
    WideCounts <- WideCounts[, !grepl("Unknown|Uknown|Notes|134", names(WideCounts))]
    LongCounts <- LongCounts |> 
        filter(!str_detect(LongCounts$Species, "Unknown")) |> 
        #remove uknowns
        filter(!str_detect(LongCounts$Species, "Uknown"))
    LongWeights <- LongWeights |> 
      filter(!str_detect(LongWeights$Species, "Unknown"))
    

#SECTION 2: Goofin' around -----------------------------------------------------
    
Bar_Data <- LongCounts |> 
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

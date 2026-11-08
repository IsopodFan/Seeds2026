
#SECTION 1: Import, Clean, and Prep Data ---------------------------------------

#call packages
  library(tidyverse)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(ggridges) 
  library(readxl) 
  library(systemfonts)
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
  library(adespatial)
  library(mgcv)

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
    
    #expand counts 
    ExpCounts <- uncount(LongCounts, weights = Count)
    
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

## 2.2: Pie Chart --------------------------------------------------------------

View(LongCounts)
Pie.df <- ExpCounts |> 
  mutate(species_lumped = fct_lump_prop(Species, prop = .005, other_level = "Other")) |> 
  count(species_lumped, name = "Abundance") %>%
  bind_rows(
    filter(., species_lumped %in% c("Rosaceae-1_count", "Rosaceae_sp._count")) |> 
      summarise(species_lumped = "Rosaceae_sp", across(where(is.numeric), sum))
  ) %>% 
  bind_rows(
    filter(., species_lumped %in% c("Rubus_bifrons_count", "Rubus_sp._count", "Rubus_phoenicolasius_count")) |> 
      summarise(species_lumped = "Rubus_sp", across(where(is.numeric), sum))
  ) |> 
  filter(!species_lumped %in% c("Rubus_sp._count", 
                                "Rubus_phoenicolasius_count", 
                                "Rubus_bifrons_count", 
                                "Rosaceae-1_count", 
                                "Rosaceae_sp._count")) |>
  mutate(pct = round((Abundance / sum(Abundance))*100, 3)) |> 
  mutate(species_pct = paste(species_lumped, pct, sep = "_")) |> 
  arrange(desc(pct)) %>%
  { bind_rows(filter(., species_pct != "Other_7.295"), filter(., species_pct == "Other_7.295")) } |> 
  mutate(species_pct = factor(species_pct, levels = species_pct))

  

  pie.plot <- ggplot(Pie.df, aes(x = "", y = Abundance, fill = species_pct)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c(
      "#902A02",
      "#F6873F",
      "#C7A897", 
      "#893345", 
      "#5c3f20", 
      "#3a130A",
      "#b5543c",
      "#da8437", 
      "#ae2b3c",
      "#643744"
    )) +
    theme_void() +
    labs(title = "Proportional Species Abdundance", fill = "Species") #+
    #geom_text(aes(label = species_lumped), position = position_stack(vjust = 0.5), size = 3) 
  

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
 
  
  ##3.4: using NMDS as beta-diversity ------------------------------------------
  
  NMDS_ord <- read_excel(here("Data/NMDS_ordination.xlsx"))
  
  View(CanopyLoss)
  
  NMDSxCLoss.df <- data.frame(
    nmds1 = NMDS_ord$nmds1, 
    nmds2 = NMDS_ord$nmds2,
    PercentLost = CanopyLoss$Percent_Lost
  )  
  
  NMDS1vsCLoss.plot <-  ggplot(NMDSxCLoss.df, aes(x = PercentLost, y = nmds1)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
    labs(
      x = "Canopy Loss (%)", 
      y = "Beta-Diversity (NMDS1)", 
      title = "Canopy Loss (%) vs Beta-Diversity"
    ) + 
    theme_minimal()
  
  summary(lm(NMDSxCLoss.df$nmds1 ~ poly(NMDSxCLoss.df$PercentLost, 2)))
  #r2 = 0.6671
  #p < 2.2 x 10^-16
  
  
  
  NMDS2vsCLoss.plot <-  ggplot(NMDSxCLoss.df, aes(x = PercentLost, y = nmds2)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
    labs(
      x = "Canopy Loss (%)", 
      y = "Beta-Diversity (NMDS2)"
    )
   
  #AIC?
  
  x <-  NMDSxCLoss.df$PercentLost
  y <- NMDSxCLoss.df$nmds2
  
  linear    <- lm(y ~ x)
  quadratic <- lm(y ~ x + I(x^2))
  cubic     <- lm(y ~ x + I(x^2) + I(x^3))
  poly_fit  <- lm(y ~ poly(x, 4))
  AIC(linear, quadratic, cubic, poly_fit) 
  
  
  ##3.6: alpha-diversity 
  ExpComm <- ExpCounts |> 
    mutate(SQP = paste(SiteQuad, Plot, sep = "_")) |> 
    select(SQP, Species) |> 
    count(SQP, Species) |> 
    pivot_wider(names_from = Species, values_from = n, values_fill = 0) |>
    column_to_rownames("SQP")
  
  ShannonExp <- data.frame(
    plot = rownames(ExpComm), 
    shannon = diversity(ExpComm, index = "shannon")
  )
  
  ShannonExp <- ShannonExp |> 
    mutate(Percent_Lost = CanopyLoss$Percent_Lost)
  
  AlpDiv.plot <- ggplot(ShannonExp, aes(x = Percent_Lost, y = shannon)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = "Canopy Loss (%)", 
      y = "Individual Subplot Alpha-Diversity (shannon)",
      title = "Canopy Loss (%) vs Subplot Alpha-Diversity"
    ) + 
    theme_minimal() 
  
  summary(lm(ShannonExp$shannon ~ ShannonExp$Percent_Lost))
  #r2 = .0026
  #p = .25
  
  
  
  #AIC
  
  x <-  ShannonExp$Percent_Lost
  y <- ShannonExp$shannon
  
  linear    <- lm(y ~ x)
  quadratic <- lm(y ~ x + I(x^2))
  cubic     <- lm(y ~ x + I(x^2) + I(x^3))
  poly_fit  <- lm(y ~ poly(x, 4))
  AIC(linear, quadratic, cubic, poly_fit) 
  
#SECTION 4: Alpha-div Over Time ------------------------------------------------ 

  Alp.df <- ExpCounts |> 
    mutate(SQP = paste(SiteQuad, Plot, sep = "_")) |> 
    mutate(Jday = as.numeric(difftime(as.Date(ExpCounts$Date_Collected), 
                                      as.Date("2023-01-01"), units = "days"))) |> 
    filter(Jday < 1000) |> 
    select(SQP, Jday, Species) |>
    count(SQP, Jday, Species) |>
    pivot_wider(names_from = Species, values_from = n, values_fill = 0)
    
  Alp.df <- data.frame(
    plot = Alp.df$SQP, 
    Jday = Alp.df$Jday,
    shannon = diversity(Alp.df[3:25], index = "shannon")
  )
  
  AlpOverTime.plot <- ggplot(Alp.df, aes(x = Jday, y = shannon)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = "Days Since January 1st, 2023", 
      y = "Individual Subplot Alpha-Diversity (shannon)",
      title = "Subplot Alpha-Diversity Over Time"
    ) + 
    theme_minimal() 
  
  summary(lm(Alp.df$shannon ~ Alp.df$Jday))
  #r2 = .02531 
  #p < 2.2 x 10^16
  
#SECTION 5: fuckin around and vibe coding while procrastinating ------------------------------------------
  
## 5.1: Long Counts to subplot x species matrix ---------------------------------
 
# check to be sure all the subplots are labelled correctly
stopifnot(all(grepl("^(RBP|TFO)_[0-9]{4}$", LongCounts$SiteQuad)))
  
# create key column 
  LongCounts <- LongCounts |> 
    mutate( 
      key = paste(SiteQuad, Plot, sep = "_")
      )
  
# take only the species columns 
  sp_cols <- names(WideCounts)[match("Acer_platanoides_count", names(WideCounts)):
                                 match("Helianthus_annus_count", names(WideCounts))]
  
# separate herbaceous from woody for later analysis 
  herb_start <- match("Herbaceous_seed_non-focal", sp_cols)
  woody_cols <- sp_cols[seq_len(herb_start - 1)]
  herb_cols  <- sp_cols[herb_start:length(sp_cols)]
  cat("woody:", length(woody_cols), "| herbaceous:", length(herb_cols), "\n")
  
# summarize across all collections to give one row per subplot, one column per species 
  sp_mat <- LongCounts |> 
    group_by(key, Species) |> 
    summarise(Count = sum(Count), .groups = "drop") |> 
    pivot_wider(names_from  = Species, 
                values_from = Count, 
                values_fill = 0)
  
## 5.2: Canopy Loss and Merge --------------------------------------------------
  
  # make canopy formatting match seed data
  canopy <- CanopyLoss |> 
    mutate(
      Plot        = str_pad(as.character(Plot), 4, pad = "0"), 
      key         = paste(Site, Plot, Subplot, sep = "_"), 
      canopy_loss = as.numeric(`canopy loss`), 
      plot_id     = paste(Site, Plot, sep = "_")
    )
  
  # join seed data and canopy data
  dat <- inner_join(canopy, sp_mat, by = "key") 
  cat("merged rows:", nrow(dat), "(expect 120)\n")
  print(anti_join(canopy, sp_mat, by = "key")$key)   #needs to be empty 
  
  # add canopy loss buckets 
  
  dat <- dat |>
    mutate(class   = cut(canopy_loss, c(-Inf, 10, 20, 30, Inf),
                         labels = c("0-10", "10-20", "20-30", "30+")),
           plot_id = factor(plot_id),
           Site    = factor(Site))
  
  ## 5.3: Analysis Function ---------------------------------------------------- 
  
  run_analysis <- function(d, cols, label, trans = sqrt) {
    
    ## output: 
      # Mean distance to centroid for each canopy loss class 
      # permutation test to see if spreads differ 
      # tukeyHSD test for which sets of classes differ 
      # spearman tests - negative rho means the subplots at higher loss are closer 
        # to center and contribute less to unique comp 
      # mean pairwise BC betadiv within class, scaled from 0 to 1
      # NMDS block
    
    
    
    cat("\n=====", label, "=====\n")
    Y <- as.matrix(d[, cols]); rownames(Y) <- d$key
    keep <- rowSums(Y) > 0 
    Y <- Y[keep, , drop = FALSE]; d <- d[keep, ]
    cls <- droplevels(d$class)
    pid <- droplevels(d$plot_id)
    
    bc <- vegdist(trans(Y), method = "bray")
    
    # (A) dispersion by canopy class 
    bd <- betadisper(bc, cls, type = "centroid") 
    cat("Mean distance to class centroid: \n"); 
    print(round(tapply(bd$distances, cls, mean), 3))
    
    ## permute whole plots because canopy loss is not independent between them 
    
    ctrl <- how(
      nperm  = 9999, 
      plots  = Plots(strata = pid, type = "free"), 
      within = Within(type = "none")
    )
    print(permutest(bd, permutations = ctrl))
    print(TukeyHSD(bd))
    
    # (B) continuous measures 
    
    d$dist_global <- betadisper(bc, factor(rep("all", nrow(Y))))$distances
    d$LCBD <- beta.div(Y, method = "hellinger", nperm = 999)$LCBD
    print(cor.test(d$canopy_loss, d$dist_global, method = "spearman", exact = FALSE))
    print(cor.test(d$canopy_loss, d$LCBD,        method = "spearman", exact = FALSE))
    
    ## (C) mean pairwise Bray-Curtis within each class (direct beta diversity)
    bcm <- as.matrix(bc)
    pw  <- sapply(levels(cls), function(k) {
      i <- which(cls == k)
      if (length(i) < 3) return(NA_real_)
      m <- bcm[i, i]; mean(m[lower.tri(m)])
    })
    cat("Mean pairwise BC within class:\n"); print(round(pw, 3))
    
    ## (D) NMDS, to check what Figure 4's NMDS1 actually captures
    nm <- metaMDS(bc, k = 2, trymax = 50, trace = 0)
    d$NMDS1 <- nm$points[, 1]; d$NMDS2 <- nm$points[, 2]
    cat("NMDS stress:", round(nm$stress, 3), "\n")
    cat("Spearman NMDS1~canopy:", round(cor(d$canopy_loss, d$NMDS1, method = "spearman"), 2),
        "| NMDS2~canopy:", round(cor(d$canopy_loss, d$NMDS2, method = "spearman"), 2), "\n")
    cat("SD of NMDS1 by class:\n"); print(round(tapply(d$NMDS1, cls, sd), 3))
    
    invisible(list(data = d, bc = bc, bd = bd))
  }
  
## 5.4: Run the function ------------------------------------------------------- 
  res_woody <- run_analysis(dat, woody_cols, "WOODY, all sites")
  res_herb  <- run_analysis(dat, herb_cols,  "HERBACEOUS (unknowns excluded), all sites")
  res_all <- run_analysis(dat, c(woody_cols, herb_cols), "All species")
  

## 5.5: Confound checks -------------------------------------------------------- 
  # I'm tired of typing this out

  Yw <- as.matrix(dat[, woody_cols])
  dat$total_woody    <- rowSums(Yw)
  dat$richness_woody <- rowSums(Yw > 0)
  print(dat %>% group_by(class) %>%
          summarise(n = n(), seeds = mean(total_woody), richness = mean(richness_woody)))
  
  ## Presence/absence: if dispersion is flat here but not with abundances,
  ## the pattern is driven by abundance, not species turnover
  bc_pa <- vegdist(Yw > 0, method = "jaccard", binary = TRUE)
  bd_pa <- betadisper(bc_pa, dat$class, type = "centroid")
  cat("\nJaccard (presence/absence) distance to centroid:\n")
  print(round(tapply(bd_pa$distances, dat$class, mean), 3))
  
  
## 5.6: Modelling continuous gradient ------------------------------------------ 
  
  tab <- res_woody$data
  tab$plot_id <- factor(tab$plot_id)
  m1 <- mgcv::gam(dist_global ~ s(canopy_loss, k = 5) + Site + s(plot_id, bs = "re"),
            data = tab, method = "REML")
  summary(m1)       # plot random effect handles subplots nested in plots
  m2 <- mgcv::gam(LCBD ~ s(canopy_loss, k = 5) + Site + s(plot_id, bs = "re"),
            data = tab, method = "REML")
  summary(m2)
  
  ## PERMANOVA for centroid shifts (report alongside betadisper)
  bc_w <- res_woody$bc
  ctrl <- how(nperm = 9999, plots = Plots(strata = tab$plot_id, type = "free"),
              within = Within(type = "none"))
  adonis2(bc_w ~ Site + canopy_loss, data = tab, permutations = ctrl, by = "terms")
  
## 5.7: Figures ---------------------------------------------------------------- 
  
  p1 <- ggplot(tab, aes(canopy_loss, dist_global, colour = Site)) +
    geom_point(alpha = .7) +
    geom_smooth(aes(group = 1), method = "gam", formula = y ~ s(x, k = 5), colour = "black") +
    labs(x = "Canopy loss (%)", y = "Distance to centroid (Bray-Curtis, sqrt counts)") +
    theme_classic()
  #ggsave("fig_dispersion_vs_canopy.png", p1, width = 6, height = 4, dpi = 300)
  
  p2 <- ggplot(tab, aes(class, dist_global)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = Site), width = .15, alpha = .6) +
    labs(x = "Canopy loss class (%)", y = "Distance to centroid") + theme_classic()
  #ggsave("fig_centroid_boxplot.png", p2, width = 6, height = 4, dpi = 300)
  
  p3 <- ggplot(tab, aes(NMDS1, NMDS2, colour = class)) +
    geom_point() + stat_ellipse(type = "t", level = .8) + theme_classic()
  #ggsave("fig_nmds_by_class.png", p3, width = 6, height = 5, dpi = 300)
  
## 5.8: test for continuous dispersion without the bins ------------------------
  
  ## Idea: fit the "expected" composition at each canopy-loss level (plus Site), then
  ## measure how far each subplot sits from its own expected position. Those residual
  ## distances are a continuous version of betadisper that no longer mixes up
  ## "composition shifts along the gradient" with "sites get more/less variable".
  ## (dist_global in section 6 does mix them: a subplot can be far from the global
  ## centroid just because the mean moved.)
  resid_dispersion <- function(d, cols, trans = sqrt, deg = 2) {
    Y <- as.matrix(d[, cols]); rownames(Y) <- d$key
    keep <- rowSums(Y) > 0
    Y <- Y[keep, , drop = FALSE]; d <- d[keep, ]
    n  <- nrow(Y)
    bc <- vegdist(trans(Y), method = "bray")
    ## PCoA coordinates with a Cailliez correction so every axis is real. Bray-Curtis
    ## produces many negative eigenvalues; just dropping those axes gave a noticeably
    ## more extreme result for woody seeds, so don't do that here.
    P  <- suppressWarnings(cmdscale(bc, k = n - 2, add = TRUE)$points)
    rhs <- if (nlevels(droplevels(d$Site)) > 1) "Site + poly(canopy_loss, deg)" else "poly(canopy_loss, deg)"
    fit <- lm(as.formula(paste("P ~", rhs)), data = d)
    R   <- residuals(fit)
    p   <- nrow(coef(fit))                                   # number of fitted terms
    d$resid_dist <- sqrt(rowSums(R^2)) * sqrt(n / (n - p))   # small-sample correction
    list(data = d, bc = bc)
  }
  
  rd_w <- resid_dispersion(dat, woody_cols)$data
  rd_h <- resid_dispersion(dat, herb_cols)$data
  rd_t <- resid_dispersion(filter(dat, Site == "TFO"), woody_cols)$data
  
  for (nm in c("rd_w", "rd_h", "rd_t")) {
    x <- get(nm)
    cat("\n", nm, "\n")
    print(cor.test(x$canopy_loss, x$resid_dist, method = "spearman", exact = FALSE))
  }
  
  ## Formal test: linear slope + smooth, with plot as a random effect
  ## (a negative canopy_loss coefficient = dispersion falls with disturbance)
  m_lin <- mgcv::gam(resid_dist ~ canopy_loss + Site + s(plot_id, bs = "re"),
               family = Gamma(link = "log"), data = rd_w, method = "REML")
  summary(m_lin)
  m_sm  <- mgcv::gam(resid_dist ~ s(canopy_loss, k = 5) + Site + s(plot_id, bs = "re"),
               family = Gamma(link = "log"), data = rd_w, method = "REML")
  summary(m_sm); plot(m_sm, select = 1, shade = TRUE)
  ## repeat m_lin with data = rd_h for herbaceous
  
  ## Moving-window beta diversity (descriptive only): mean pairwise Bray-Curtis among
  ## the k subplots closest in canopy loss to each focal subplot. Windows overlap,
  ## so don't put a p-value on this -- use it to see the SHAPE of the relationship.
  window_bd <- function(d, bc, k = 20) {
    bcm <- as.matrix(bc)
    sapply(seq_len(nrow(d)), function(i) {
      idx <- order(abs(d$canopy_loss - d$canopy_loss[i]))[seq_len(k)]
      m <- bcm[idx, idx]; mean(m[lower.tri(m)])
    })
  }
  tab$window_bc <- window_bd(res_woody$data, res_woody$bc, k = 20)
  p4 <- ggplot(tab, aes(canopy_loss, window_bc)) +
    geom_point(aes(colour = Site)) + geom_smooth(method = "loess", se = TRUE, colour = "black") +
    labs(x = "Canopy loss (%)", y = "Mean pairwise Bray-Curtis in 20-subplot window") +
    theme_classic()
  #ggsave("fig_moving_window_bd.png", p4, width = 6, height = 4, dpi = 300)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    

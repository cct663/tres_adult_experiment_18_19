# Analysis for two tree swallow experiments conducted in 2018 & 2019 in Ithaca NY.
# 
# Experiment 1: 2018, 3x2 factorial with first a challenge then a color manipulation.
# Experiment 2: 2019, 2x2 factorial with first a color manipulation then a challenge.
# 
# Written by Conor Taff
# Last updated 8/31/2022
# Run under R Studio 1.1.463 on Mac OSX 10.11.6


## Load packages ----
  pacman::p_load(tidyverse, here, lme4, scales, lmerTest, sjPlot, gridExtra, MASS, rethinking)
    ## not loading plyr because it conflicts with here and tidyverse
  
## Load data ----
  # the data files include similar data from each year but experiments will mainly be analyzed separately
  # data by female: one row per female with measures throughout the breeding season
    d_fem <- read.delim(here::here("1_raw_data", "data_by_female.txt"))
    
    # Remove records that were noted as exclude. These are mostly cases where either a nest failed so early 
      # that no treatment was delivered or something weird happened (e.g., treatments mixed up).
    d_fem <- subset(d_fem, d_fem$exclude != "Yes")
    
  # Make separate data frames for each year
    d_fem18 <- subset(d_fem, d_fem$year == "2018")
    d_fem19 <- subset(d_fem, d_fem$year == "2019")
    
      # Standardize initial breast color in each year
        d_fem18$bbright_s <- scale(d_fem18$fbbright)
        d_fem19$bbright_s <- scale(d_fem19$fbbright)
        d_fem18$fmass1_s <- scale(d_fem18$fmass1)
        d_fem19$fmass1_s <- scale(d_fem19$fmass1)
    
  # data by nestling: one row per nestling with some more individual level nestling info
    d_nestling <- read.delim(here::here("1_raw_data", "data_by_nestling.txt"))
    
  # provisioning data: one row per nest per day with provisiong info from rfid sensors
    d_provision <- read.delim(here::here("1_raw_data", "daily_provision_data.txt"))
    
  # visitation data: one row per nest per day with social visits and trips for each nest
    d_social <- read.delim(here::here("1_raw_data", "nest_visit_data.txt"))
    
  # rfid reads: total rfid reads at each box on each day (used to filter out malfunctioning rfid boards)
    d_tot_rfid <- read.delim(here::here("1_raw_data", "total_rfid_reads.txt"))
  
  # # ACTH validation experiment data  # Reported in Taff et al. 2022 Journal of Experimental Biology
  #   d_acth_fem <- read.delim(here::here("1_raw_data", "adult_acth_validation.txt"))
  #   d_acth_nest <- read.delim(here::here("1_raw_data", "nestling_acth_validation.txt"))
    
## Set colors for different treatment groups ----
    col_dull <- "slateblue"
    col_sham <- "orange"
    ch_pred <- "chartreuse3"
    ch_con <- "orange"
    ch_tape <- "coral3"
    
    # Full treatments
      col_CC <- "gray50"      # 2018 & 2019, full control
      col_CS <- "coral3"      # 2018, sham color then taping challenge
      col_CP <- "orange"      # 2018, sham color then predator
      col_DC <- "slateblue"   # 2018, dulling folowed by sham coloring
      col_DS <- "violet"      # 2018, dulling folowed by taping challenge
      col_DP <- "chartreuse3" # 2018, dulling followed by predator
        
      col_PC <- "coral3"      # 2019, predator first then sham color
      col_PD <- "chartreuse3" # 2019, predator first then dulling
      col_CD <- "slateblue"   # 2019, control first, then dulling
    
## Female Mass Models ----
      
  # Model for mass at time 2 and 3 in relation to treatment
   
    # Transform data frame to what I want separately for 2018 & 2019     
      # Make 2018 data into longer format with mass at 2nd and 3rd capture in separate rows
        d18_mass <- d_fem18 %>%
          pivot_longer(cols = c("fmass2", "fmass3"), names_to = "cap_num",
                       names_transform = list(
                         cap_num = ~ readr::parse_factor(.x, levels = c("fmass2", "fmass3"),
                                                         ordered = TRUE)),
                       values_to = "mass", values_drop_na = TRUE)
        d18_mass <- as.data.frame(d18_mass)
        
      # Add a 'stage' variable to distinguish 2nd and 3rd captures based on which parts of treatment received so far
        stage <- data.frame(cap_num = c("fmass2", "fmass3"),
                            stage = c("stage1", "stage2"))
        d18_mass <- plyr::join(d18_mass, stage, "cap_num", "left", "first")
      
        
      # Make 2019 data into longer format with mass at 2nd and 3rd capture in separate rows 
        d19_mass <- d_fem19 %>%
          pivot_longer(cols = c("fmass2", "fmass3"), names_to = "cap_num",
                       names_transform = list(
                         cap_num = ~ readr::parse_factor(.x, levels = c("fmass2", "fmass3"),
                                                         ordered = TRUE)),
                       values_to = "mass", values_drop_na = TRUE)
        d19_mass <- as.data.frame(d19_mass)
        
      # Add a 'stage' variable to distinguish 2nd and 3rd captures based on which parts of treatment received so far
        d19_mass <- plyr::join(d19_mass, stage, "cap_num", "left", "first")
      
      # Fit models for 2018. Initial brightness doesn't matter so it is dropped.  
        m_mass_18 <- lmer(mass ~ color*challenge*stage + color*bbright_s + fmass1_s + (1|fRFID), data = d18_mass)
        m_mass_18b <- lmer(mass ~ color*challenge*stage + fmass1_s + (1|fRFID), data = d18_mass)
      
      # Make a table and save it to a folder to use in markdown appendix  
        m_mass_18b_t <- tab_model(m_mass_18b,
                  pred.labels = c("Intercept", "Signal Dulled", "Predator", "Flight Reduction",
                                  "Stage", "Pre-Treatment Mass", "Dulled * Predator", "Dulled * Flight", 
                                  "Dulled * Stage", "Predator * Stage", "Flight * Stage", 
                                  "Dulled * Predator * Stage", "Dulled * Flight * Stage"),
                  dv.labels = "Post-treatment Mass")
        saveRDS(m_mass_18b, here::here("5_other_outputs/m_mass_18b.RDS"))
        saveRDS(m_mass_18b_t, here::here("5_other_outputs/m_mass_18b_t.RDS"))
      
      # Fit models for 2018. Initial brightness doesn't matter so it is dropped.  
        m_mass_19 <- lmer(mass ~ color*challenge*stage*bbright_s + fmass1_s + (1|fRFID), data = d19_mass)
        m_mass_19b <- lmer(mass ~ color*challenge*stage + fmass1_s + (1|fRFID), data = d19_mass)
      
      # Make a table and save it to a folder to use in markdown appendix
        m_mass_19b_t <- tab_model(m_mass_19b,
                  pred.labels = c("Intercept", "Predator", "Signal Dulled", "Stage", 
                                  "Pre-Treatment Mass", "Dulled * Predator",
                                  "Predator * Stage", "Dulled * Stage", "Dulled * Predator * Stage"),
                  dv.labels = "Post-treatment Mass")
        saveRDS(m_mass_19b_t, here::here("5_other_outputs/m_mass_19b_t.RDS"))

## Female Mass Plots ----      
  # Figures showing adult mass at each timepoint by treatment split by year (experiment)
      # Pivot mass into long format with one row for each mass measurement
        long_mass <- d_fem %>%
            pivot_longer(cols = c("fmass1", "fmass2", "fmass3"), names_to = "cap_num", 
                         names_transform = list(
                           cap_num = ~ readr::parse_factor(.x, levels = c("fmass1", "fmass2", "fmass3"),
                                                           ordered = TRUE)),
                         values_to = "mass", values_drop_na = TRUE)  
        long_mass <- as.data.frame(long_mass)
        long_mass$xpos <- as.numeric(long_mass$cap_num)  # capture number was treated as a factor
          
      # Pivot and summarize slightly differently to calculate mean and SE for each treatment/time point    
        sum_mass <- d_fem %>%
          pivot_longer(cols = c("fmass1", "fmass2", "fmass3"), names_to = "cap_num", 
                       names_transform = list(
                         cap_num = ~ readr::parse_factor(.x, levels = c("fmass1", "fmass2", "fmass3"),
                                                         ordered = TRUE)),
                       values_to = "mass", values_drop_na = TRUE) %>% 
                group_by(as.factor(year), full_treatment, cap_num) %>%
                summarise(n = n(), mu = mean(mass, na.rm = TRUE), se = sd(mass, na.rm = TRUE) / sqrt(n()))
        sum_mass$x_pos <- as.numeric(sum_mass$cap_num)
      
      # Make some values to use for plotting. Specifies point shapes, colors, horizontal jittering.  
        for_plot <- data.frame(full_treatment = unique(sum_mass$full_treatment),
                              p_shape = c(21, 24, 22, 21, 24, 22, 21, 24, 24),
                              p_col = c(rep(col_sham, 3), rep(col_dull, 3), col_dull, col_sham, col_dull),
                              p_dodge = c(-.25, -.15, -.05, .05, .15, .25, .0833, -.0833, .25))
      
      # Make a few adjustments for plotting and join to the dataframe with the plotting parameters   
        sum_mass2 <- as.data.frame(sum_mass)   
        colnames(sum_mass2)[1] <- "year"
        sum_mass2 <- plyr::join(sum_mass2, for_plot, "full_treatment", "left", "first")
       
      # Join the initial long form object above to the plotting parameters data frame 
        long_mass2 <- plyr::join(long_mass, for_plot, "full_treatment", "left", "first")
        
      # Split both the points and the mu + se data frames into separate 2018 & 2019 subsets
        lm2_18 <- subset(long_mass2, long_mass2$year == "2018")
        lm2_19 <- subset(long_mass2, long_mass2$year == "2019")
        
        sm2_18 <- subset(sum_mass2, sum_mass2$year == "2018")
        sm2_19 <- subset(sum_mass2, sum_mass2$year == "2019")
      
      # Set up plotting paramters  
        par(mfrow = c(1, 2))
     
      # Initiate plot for 2018   
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Dulling then Challenge", xlim = c(0.3, 3.7), ylim = c(15.5, 25),
                 xlab = "Capture Number", ylab = "Female Mass (g)")
            axis(1, seq(0, 4, 1))
            axis(2, seq(5, 30, 1), las = 2)
          
          # Add faded raw data points  
            points(lm2_18$xpos + lm2_18$p_dodge, lm2_18$mass, pch = 16, col = alpha("black", 0.3))
          
          # Add mu + se points for each treatment group  
            for(i in 1:nrow(sm2_18)){
              lines(rep(sm2_18$x_pos[i] + sm2_18$p_dodge[i], 2), 
                    c(sm2_18$mu[i] - sm2_18$se[i], sm2_18$mu[i] + sm2_18$se[i]), lwd = 2)
            }
            points(sm2_18$x_pos + sm2_18$p_dodge, sm2_18$mu, pch = sm2_18$p_shape, cex = 1.4,
                   bg = as.character(sm2_18$p_col))
            #No legend in this one since legend is shard across the two panels
            #legend("topright", c("No Color", "Dull Color", "Control", "Predator", "Handicap"),
            #       pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex = 1.2, bty = "n")
           
          # Add labeling rectangles and text for when treatments happened 
            rect(1.4, 15.7, 3.5, 16.2, col = alpha("chartreuse3", 0.3))
            text(2.15, 15.95, "Dulling", pos = 4)
            
            rect(2.2, 16.3, 3.5, 16.8, col = alpha("coral3", 0.2))
            text(2.4, 16.55, "Challenge", pos = 4)
        
        #Initiate Plot for 2018
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Challenge then Dulling", xlim = c(0.3, 3.7), ylim = c(15.5, 25),
                 xlab = "Capture Number", ylab = "Female Mass (g)")
            axis(1, seq(0, 4, 1))
            axis(2, seq(5, 30, 1), las = 2)
           
          # Add faded raw data points  
            points(lm2_19$xpos + lm2_19$p_dodge, lm2_19$mass, pch = 16, col = alpha("black", 0.3))
          
          # Add mu + se points for each treatment group   
            for(i in 1:nrow(sm2_19)){
              lines(rep(sm2_19$x_pos[i] + sm2_19$p_dodge[i], 2), 
                    c(sm2_19$mu[i] - sm2_19$se[i], sm2_19$mu[i] + sm2_19$se[i]), lwd = 2)
            }
            points(sm2_19$x_pos + sm2_19$p_dodge, sm2_19$mu, pch = sm2_19$p_shape, cex = 1.4,
                   bg = as.character(sm2_19$p_col))
            legend("topright", c("Sham", "Dulled", "Control", "Predator", "Handicap"),
                   pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex = 1, bty = "n")
          
          # Add labeling rectangles and text for when treatments happened   
            rect(2.3, 15.7, 3.6, 16.2, col = alpha("chartreuse3", 0.3))
            text(2.65, 15.95, "Dulling", pos = 4)
            
            rect(1.4, 16.3, 2.5, 16.8, col = alpha("coral3", 0.2))
            text(1.5, 16.55, "Challenge", pos = 4)
    
    
## Female Corticosterone Models ----
    
    # Scale continuous predictors to mean 0 and sd 1
        d_fem18$fbase1_s <- scale(d_fem18$fbase1)
        d_fem18$fstr1_s <- scale(d_fem18$fstr1)
        d_fem18$fdex1_s <- scale(d_fem18$fdex1)
        
    # Fit simple models for each cort measure in 2018 experiment
        m18_base2 <- lm(fbase2 ~ color + fbase1_s, data = d_fem18)
        m18_base3 <- lm(fbase3 ~ color*challenge + fbase1_s, data = d_fem18)
        m18_str3 <- lm(fstr3 ~ color*challenge + fstr1_s, data = d_fem18)
        m18_dex3 <- lm(fdex3 ~ color*challenge + fdex1_s, data = d_fem18)
   
    # Make a table with cleaned up names from all cort models 2018     
        m_18_cort_t <- tab_model(m18_base2, m18_base3, m18_str3, m18_dex3,
                                 pred.labels = c("Intercept", "Signal Dulled", "Initial Base Cort", 
                                                 "Predator", "Flight Reduction", "Dulled * Predator",
                                                 "Dulled * Stress", "Initial Stress-Induced Cort",
                                                 "Initial Post-Dex Cort"),
                                 dv.labels = c("2nd Capture Baseline", "3rd Capture Baseline",
                                               "3rd Capture Stress-Induced", "3rd Capture Post-Dexamethasone"))
        
    # Save the table to file so I can add it to the markdown appendix
        saveRDS(m_18_cort_t, here::here("5_other_outputs/m_18_cort_t.RDS"))
    
    # Scale continous predictors for 2019    
        d_fem19$fbase1_s <- scale(d_fem19$fbase1)
        d_fem19$fstr1_s <- scale(d_fem19$fstr1)
        d_fem19$fdex1_s <- scale(d_fem19$fdex1)
        
    # Fit simple models for 2019 cort. Note that 'fex3' here is actually post-cortrosyn injection
        m19_base2 <- lm(fbase2 ~ color + fbase1_s, data = d_fem19)
        m19_base3 <- lm(fbase3 ~ color*challenge + fbase1_s, data = d_fem19)
        m19_str3 <- lm(fstr3 ~ color*challenge + fstr1_s, data = d_fem19)
        m19_dex3 <- lm(fdex3 ~ color*challenge, data = d_fem19)
     
    # Make a cleaned up table for 2019 experiment models   
        m_19_cort_t <- tab_model(m19_base2, m19_base3, m19_str3, m19_dex3,
                  pred.labels = c("Intercept", "Predator", "Initial Base Cort",
                                  "Signal Dulled", "Dulled * Predator", "Initial Stress-Induced Cort"),
                  dv.labels = c("2nd Capture Baseline", "3rd Capture Baseline",
                                "3rd Capture Stress-Induced", "3rd Capture Post-Cortrosyn"))
        
    # Save the 2019 table to add to markdown appendix
        saveRDS(m_19_cort_t, here::here("5_other_outputs/m_19_cort_t.RDS"))
 
## Female Corticosterone Plots ----
        
    ## manipulate cort into long format for plotting
          long_cort <- d_fem %>%
            pivot_longer(cols = c("fbase1", "fstr1", "fdex1", "fbase2", "fbase3", "fstr3", "fdex3"), names_to = "cort_type", 
                         names_transform = list(
                           cap_num = ~ readr::parse_factor(.x, levels = c("fbase1", "fstr1", "fdex1", "fbase2", "fbase3", "fstr3", "fdex3"),
                                                           ordered = TRUE)),
                         values_to = "cort", values_drop_na = TRUE)  
          long_cort <- as.data.frame(long_cort)
          
    ## Add some plotting positions and merge with long form cort file
          cort_pos <- data.frame(cort_type = c("fbase1", "fstr1", "fdex1", "fbase2", "fbase3", "fstr3", "fdex3"),
                                 xpos = c(1, 2, 3, 5, 7, 8, 9))
          long_cort <- plyr::join(long_cort, cort_pos, "cort_type", "left", "first")
          
    ## Create a seperate data frame that has mean and se for each treatment group      
          sum_cort <- d_fem %>%
            pivot_longer(cols = c("fbase1", "fstr1", "fdex1", "fbase2", "fbase3", "fstr3", "fdex3"), names_to = "cort_type", 
                         names_transform = list(
                           cap_num = ~ readr::parse_factor(.x, levels = c("fbase1", "fstr1", "fdex1", "fbase2", "fbase3", "fstr3", "fdex3"),
                                                           ordered = TRUE)),
                         values_to = "cort", values_drop_na = TRUE) %>% 
            group_by(as.factor(year), full_treatment, cort_type) %>%
            summarise(n = n(), mu = mean(cort, na.rm = TRUE), se = sd(cort, na.rm = TRUE) / sqrt(n()))
          sum_cort <- as.data.frame(sum_cort)
          sum_cort <- plyr::join(sum_cort, cort_pos, "cort_type", "left", "first")
     
    ## Make plotting parameters for point shapes, colors, and horizontal jitter     
          for_plot <- data.frame(full_treatment = unique(sum_cort$full_treatment),
                                 p_shape = c(21, 24, 22, 21, 24, 22, 21, 24, 24),
                                 p_col = c(rep(col_sham, 3), rep(col_dull, 3), col_dull, col_sham, col_dull),
                                 p_dodge = c(-.25, -.15, -.05, .05, .15, .25, .0833, -.0833, .25))
    
    ## Join the long objects to the plotting parameters            
          sum_cort <- plyr::join(sum_cort, for_plot, "full_treatment")
          long_cort <- plyr::join(long_cort, for_plot, "full_treatment")
          colnames(sum_cort)[2] <- "year"
   
    ## Create separate data frames for each year       
          lc2_18 <- subset(long_cort, long_cort$year == "2018")
          lc2_19 <- subset(long_cort, long_cort$year == "2019")
          
          sc2_18 <- subset(sum_cort, sum_cort$year == "2018")
          sc2_19 <- subset(sum_cort, sum_cort$year == "2019")
          
    # Make 2018 plot
        #Initiate the plot  
          plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
               main = "Dulling then Challenge", xlim = c(0.3, 9.7), ylim = c(-20, 80),
               xlab = "", ylab = "Corticosterone (ng/ul)")
          axis(1, c(-1, 1, 2, 3, 5, 7, 8, 9, 15), c("", "B 1", "S 1", "D 1", "B 2", "B 3", "S 3", "D 3", ""))
          axis(2, c(-50, seq(0, 200, 10)), las = 2)
       
        # Add raw points   
          points(lc2_18$xpos + lc2_18$p_dodge, lc2_18$cort, pch = 16, col = alpha("black", 0.3))
          abline(h = 0)
        
        # Add group means and standard error lines  
          for(i in 1:nrow(sc2_18)){
            lines(rep(sc2_18$xpos[i] + sc2_18$p_dodge[i], 2), 
                  c(sc2_18$mu[i] - sc2_18$se[i], sc2_18$mu[i] + sc2_18$se[i]), lwd = 2)
          }
          points(sc2_18$xpos + sc2_18$p_dodge, sc2_18$mu, pch = sc2_18$p_shape, cex = 1.4,
                 bg = as.character(sc2_18$p_col))
        
        # Add labels for when manipulations occurred  
          rect(4.4, -15, 9.6, -10, col = alpha("chartreuse3", 0.3))
          text(6.5, -12.5, "Dulling", pos = 4)
          
          rect(6.4, -9, 9.6, -4, col = alpha("coral3", 0.2))
          text(7.3, -6.5, "Challenge", pos = 4)
        
        # Add legend  
          legend("bottomleft", c("No Color", "Dull Color", "Control", "Predator", "Handicap"),
                 pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex =0.9, bty = "n")
          
    #Make 2019 plot
         
        # Initiate plot 
          plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
               main = "Challenge then Dulling", xlim = c(0.3, 9.7), ylim = c(-20, 80),
               xlab = "", ylab = "Corticosterone (ng/ul)")
          axis(1, c(-1, 1, 2, 3, 5, 7, 8, 9, 15), c("", "B 1", "S 1", "D 1", "B 2", "B 3", "S 3", "A 3", ""))
          axis(2, c(-50, seq(0, 200, 10)), las = 2)
        
        # Add raw points  
          points(lc2_19$xpos + lc2_19$p_dodge, lc2_19$cort, pch = 16, col = alpha("black", 0.3))
          abline(h = 0)
        
        # Add group means and standard errors  
          for(i in 1:nrow(sc2_19)){
            lines(rep(sc2_19$xpos[i] + sc2_19$p_dodge[i], 2), 
                  c(sc2_19$mu[i] - sc2_19$se[i], sc2_19$mu[i] + sc2_19$se[i]), lwd = 2)
          }
          points(sc2_19$xpos + sc2_19$p_dodge, sc2_19$mu, pch = sc2_19$p_shape, cex = 1.4,
                 bg = as.character(sc2_19$p_col))
        
        # Add labeling rectangles and text  
          rect(6.4, -15, 9.6, -10, col = alpha("chartreuse3", 0.3))
          text(7, -12.5, "Dulling", pos = 4)
          
          rect(4.4, -9, 6.6, -4, col = alpha("coral3", 0.2))
          text(5, -6.5, "Challenge", pos = 4)
        
        # Add legend  
          legend("bottomleft", c("No Color", "Dull Color", "Control", "Predator"),
                 pch = c(21, 21, 21, 24), pt.bg = c(col_sham, col_dull, rep("white", 2)), cex =0.9, bty = "n")

## Reproductive Success Models ----
      #Models for 2018
          m_clutch_18 <- lm(clutch ~ color*challenge, data = d_fem18)
          m_brood_18 <- lm(maxbrood ~ color*challenge, data = d_fem18)
          m_numd6_18 <- lm(numd6 ~ color*challenge, data = d_fem18)
          m_numband_18 <- lm(numband ~ color*challenge, data = d_fem18)
          m_num_d15_18 <- lm(num_d15 ~ color*challenge, data = d_fem18)
          m_numfled_18 <- lm(numfled ~ color*challenge, data = d_fem18)
      
      # Overall effects reported from ANOVA. Group estimates from model summary tables.
          
        # Make table for 2018
          m_RS_18_t <- tab_model(m_clutch_18, m_brood_18, m_numd6_18, m_numband_18, m_num_d15_18, m_numfled_18,
                    pred.labels = c("Intercept", "Signal Dulled", "Predator", "Flight Reduction",
                                    "Dulled * Predator", "Dulled * Flight"),
                    dv.labels = c("Clutch Size", "Num. Hatched", "Num. Day 6", "Num. Day 12", "Num. Day 15", "Num. Fledged"))
          saveRDS(m_RS_18_t, here::here("5_other_outputs/m_RS_18_t.RDS"))
          
      #Models for 2019
          m_clutch_19 <- lm(clutch ~ color*challenge, data = d_fem19)
          m_brood_19 <- lm(maxbrood ~ color*challenge, data = d_fem19)
          m_numd6_19 <- lm(numd6 ~ color*challenge, data = d_fem19)
          m_numband_19 <- lm(numband ~ color*challenge, data = d_fem19)
          m_num_d15_19 <- lm(num_d15 ~ color*challenge, data = d_fem19)
          m_numfled_19 <- lm(numfled ~ color*challenge, data = d_fem19)
          
       # Overall p-values from anova. Estimates from model summary table  
          
          # Make table for 2019
          m_RS_19_t <- tab_model(m_clutch_19, m_brood_19, m_numd6_19, m_numband_19, m_num_d15_19, m_numfled_19,
                    pred.labels = c("Intercept", "Predator", "Signal Dulled", "Predator * Dulled"),
                    dv.labels = c("Clutch Size", "Num. Hatched", "Num. Day 6", "Num. Day 12", "Num. Day 15", "Num. Fledged"))
          saveRDS(m_RS_19_t, here::here("5_other_outputs/m_RS_19_t.RDS"))
          
## Reproductive Success Plots ----
          
      # Make all reproductive success measures into a long format data frame
        long_rs <- d_fem %>%
          pivot_longer(cols = c("clutch", "maxbrood", "numd6", "numband", "num_d15", "numfled"), names_to = "time_point", 
                       names_transform = list(
                         time_point = ~ readr::parse_factor(.x, levels = c("clutch", "maxbrood", "numd6", "numband", "num_d15", "numfled"),
                                                         ordered = TRUE)),
                       values_to = "count", values_drop_na = TRUE)  
        long_rs <- as.data.frame(long_rs)
        rs_pos <- data.frame(time_point = c("clutch", "maxbrood", "numd6", "numband", "num_d15", "numfled"),
                               xpos = seq(1, 6, 1))
        long_rs <- plyr::join(long_rs, rs_pos, "time_point", "left", "first")
        
     # Make a sumarized data frame that has means and errors for each group   
        sum_rs <- d_fem %>%
          pivot_longer(cols = c("clutch", "maxbrood", "numd6", "numband", "num_d15", "numfled"), names_to = "time_point", 
                       names_transform = list(
                         time_point = ~ readr::parse_factor(.x, levels = c("clutch", "maxbrood", "numd6", "numband", "num_d15", "numfled"),
                                                         ordered = TRUE)),
                       values_to = "count", values_drop_na = TRUE) %>% 
          group_by(as.factor(year), full_treatment, time_point) %>%
          summarise(n = n(), mu = mean(count, na.rm = TRUE), se = sd(count, na.rm = TRUE) / sqrt(n()))
        sum_rs <- as.data.frame(sum_rs)
        sum_rs <- plyr::join(sum_rs, rs_pos, "time_point", "left", "first")
        colnames(sum_rs)[2] <- "year"
     
      # Add plotting parameters and join to the objects   
        for_plot <- data.frame(full_treatment = unique(sum_rs$full_treatment),
                               p_shape = c(21, 24, 22, 21, 24, 22, 21, 24, 24),
                               p_col = c(rep(col_sham, 3), rep(col_dull, 3), col_dull, col_sham, col_dull),
                               p_dodge = c(-.25, -.15, -.05, .05, .15, .25, .0833, -.0833, .25))
        
        sum_rs <- plyr::join(sum_rs, for_plot, "full_treatment")
        long_rs <- plyr::join(long_rs, for_plot, "full_treatment")
      
      # Split into separate data frames for each experiment  
        rs2_18 <- subset(long_rs, long_rs$year == "2018")
        rs2_19 <- subset(long_rs, long_rs$year == "2019")
        
        srs2_18 <- subset(sum_rs, sum_rs$year == "2018")
        srs2_19 <- subset(sum_rs, sum_rs$year == "2019")
        
      # Make the 2018 plot
        # Initiate Plot
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Dulling then Challenge", xlim = c(0.5, 6.5), ylim = c(-0.1, 7.1),
                 xlab = "", ylab = "Count")
            axis(1, c(-1, 1, 2, 3, 4, 5, 6, 10), c("", "Clutch", "Hatch", "Day 6", "Day 12", "Day 15", "Fledged", ""))
            axis(2, seq(-1, 10, 1), las = 2)
        
        # Add raw points    
            points(rs2_18$xpos + rs2_18$p_dodge, rs2_18$count, pch = 16, col = alpha("black", 0.3))
            #abline(h = 0)
        
        # Add group means and errors    
            for(i in 1:nrow(srs2_18)){
              lines(rep(srs2_18$xpos[i] + srs2_18$p_dodge[i], 2), 
                    c(srs2_18$mu[i] - srs2_18$se[i], srs2_18$mu[i] + srs2_18$se[i]), lwd = 2)
            }
            points(srs2_18$xpos + srs2_18$p_dodge, srs2_18$mu, pch = srs2_18$p_shape, cex = 1.4,
                   bg = as.character(srs2_18$p_col))
          
        # Add lines to separate each time point  
            abline(v = 1.5, lty = 2)
            abline(v = 2.5, lty = 2)
            abline(v = 3.5, lty = 2)
            abline(v = 4.5, lty = 2)
            abline(v = 5.5, lty = 2)
         
        # Add labels and text
            rect(0, 6.6, 3.5, 7, col = alpha("chartreuse3", 0.3))
            text(1.7, 6.8, "Dulling", pos = 4)
            
            rect(2, 6, 3, 6.4, col = alpha("coral3", 0.2))
            text(2.1, 6.2, "Challenge", pos = 4)
          
        #   Add legend
            legend("bottomleft", c("No Color", "Dull Color", "Control", "Predator", "Handicap"),
                   pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex =0.9, bty = "n")
        
      #Make 2019 plot
        # Initiate the plot
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Challenge then Dulling", xlim = c(0.5, 6.5), ylim = c(-0.1, 7.1),
                 xlab = "", ylab = "Count")
            axis(1, c(-1, 1, 2, 3, 4, 5, 6, 10), c("", "Clutch", "Hatch", "Day 6", "Day 12", "Day 15", "Fledged", ""))
            axis(2, seq(-1, 10, 1), las = 2)
       
        # Add raw points     
            points(rs2_19$xpos + rs2_19$p_dodge, rs2_19$count, pch = 16, col = alpha("black", 0.3))
            #abline(h = 0)
        
        # Add group means and errors    
            for(i in 1:nrow(srs2_19)){
              lines(rep(srs2_19$xpos[i] + srs2_19$p_dodge[i], 2), 
                    c(srs2_19$mu[i] - srs2_19$se[i], srs2_19$mu[i] + srs2_19$se[i]), lwd = 2)
            }
            points(srs2_19$xpos + srs2_19$p_dodge, srs2_19$mu, pch = srs2_19$p_shape, cex = 1.4,
                   bg = as.character(srs2_19$p_col))
       
        # Add lines to separate time points     
            abline(v = 1.5, lty = 2)
            abline(v = 2.5, lty = 2)
            abline(v = 3.5, lty = 2)
            abline(v = 4.5, lty = 2)
            abline(v = 5.5, lty = 2)
        
        # Add labels and text    
            rect(0, 6.6, 3.5, 7, col = alpha("chartreuse3", 0.3))
            text(1.7, 6.8, "Dulling", pos = 4)
            
            text(.9, 6.2, "Challenge", pos = 4)
            arrows(0.5, 6.2, 0.9, 6.2, code = 1)
         
        # Add legend   
            legend("bottomleft", c("No Color", "Dull Color", "Control", "Predator"),
                   pch = c(21, 21, 21, 24), pt.bg = c(col_sham, col_dull, rep("white", 2)), cex =0.9, bty = "n")
    
## Female Glucose Models ----
            
        # Scale continuous predictors to mean 0 and sd 1
            d_fem18$fbgluc1_s <- scale(d_fem18$fbgluc1)
            d_fem18$fsgluc1_s <- scale(d_fem18$fsgluc1)
            
        # Fit simple models for each cort measure in 2018 experiment
            m18_bgluc2 <- lm(fbgluc2 ~ color*challenge + fbgluc1_s, data = d_fem18)
            m18_bgluc3 <- lm(fbgluc3 ~ color*challenge + fbgluc1_s, data = d_fem18)
            m18_sgluc3 <- lm(fsgluc3 ~ color*challenge + fsgluc1_s, data = d_fem18)
            
        # Make a table with cleaned up names from all cort models 2018     
            m_18_gluc_t <- tab_model(m18_bgluc2, m18_bgluc3, m18_sgluc3,
                                     pred.labels = c("Intercept", "Signal Dulled", "Predator",
                                                     "Flight Reduction", "Initial Baseline Glucose",
                                                     "Dulled * Predator",
                                                     "Dulled * Flight", "Initial Stress-Induced Glucose"),
                                     dv.labels = c("2nd Capture Base Glucose", "3rd Capture Base Glucose",
                                                   "3rd Capture Stress-Induced Glucose"))
  
            
        # Save the table to file so I can add it to the markdown appendix
            saveRDS(m_18_gluc_t, here::here("5_other_outputs/m_18_gluc_t.RDS"))
            
        # Scale continuous predictors to mean 0 and sd 1
            d_fem19$fbgluc1_s <- scale(d_fem19$fbgluc1)
            d_fem19$fsgluc1_s <- scale(d_fem19$fsgluc1)
            
        # Fit simple models for each cort measure in 2018 experiment
            m19_bgluc2 <- lm(fbgluc2 ~ color*challenge + fbgluc1_s, data = d_fem19)
            m19_bgluc3 <- lm(fbgluc3 ~ color*challenge + fbgluc1_s, data = d_fem19)
            m19_sgluc3 <- lm(fsgluc3 ~ color*challenge + fsgluc1_s, data = d_fem19)
            
        # Make a table with cleaned up names from all cort models 2018     
            m_19_gluc_t <- tab_model(m19_bgluc2, m19_bgluc3, m19_sgluc3,
                                     pred.labels = c("Intercept", "Predator", "Signal Dulled", 
                                                     "Initial Baseline Glucose",
                                                     "Dulled * Predator",
                                                     "Initial Stress-Induced Glucose"),
                                     dv.labels = c("2nd Capture Base Glucose", "3rd Capture Base Glucose",
                                                   "3rd Capture Stress-Induced Glucose"))
            
            
        # Save the table to file so I can add it to the markdown appendix
            saveRDS(m_19_gluc_t, here::here("5_other_outputs/m_19_gluc_t.RDS")) 
    
## Female Glucose Plots ----

  # no plots made of glucose data at present  
    
## Nestling Morphology Models ----
  # Join nestlings to adult table to get treatment info
      d_fem2 <- d_fem[, c("soc_uby", "color", "challenge")]
      d_nestling2 <- plyr::join(d_nestling, d_fem2, "soc_uby", "left", "first")
      d_nestling2 <- subset(d_nestling2, d_nestling2$raised_nest == "home" | 
                              d_nestling2$raised_nest == "cross")
      d_nestling2$raised_nest = relevel(d_nestling2$raised_nest, ref = "home")
            
  # Make separate dataframes for each year
      d_nestling18 <- subset(d_nestling2, d_nestling2$year == "2018")
      d_nestling19 <- subset(d_nestling2, d_nestling2$year == "2019")
      
  # Fit models for mass at day 12 & 15 plus headbill & wing at day 12 for 2018
      m_d6_mass_18 <- lm(d6_av_mass ~ color*challenge, data = d_fem18)
      m_d12_mass_18 <- lmer(d12_mass ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling18)
      m_d15_mass_18 <- lmer(d15_mass ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling18)
      m_d12_hb_18 <- lmer(d12_head ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling18)
      m_d12_wing_18 <- lmer(d12_wing ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling18)
  
  # Make a table of the models
      m18_nest_morph_t <- tab_model(m_d6_mass_18, m_d12_hb_18, m_d12_wing_18, m_d12_mass_18, m_d15_mass_18,
                pred.labels = c("Intercept", "Signal Dulled", "Predator", "Flight Reduction",
                                "Dulled * Predator", "Dulled * Flight", "Cross-Fostered"),
                dv.labels = c("Average Mass Day 6", "Head + Bill Day 12", 
                              "Wing Length Day 12", "Mass Day 12", "Mass Day 15"))
  
  # Save table
      saveRDS(m18_nest_morph_t, here::here("5_other_outputs/m18_nest_morph_t.RDS"))
      
  # Fit models for mass at day 12 & 15 plus headbill & wing at day 12 for 2018
      m_d6_mass_19 <- lm(d6_av_mass ~ color*challenge, data = d_fem19)
      m_d12_mass_19 <- lmer(d12_mass ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling19)
      m_d15_mass_19 <- lmer(d15_mass ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling19)
      m_d12_hb_19 <- lmer(d12_head ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling19)
      m_d12_wing_19 <- lmer(d12_wing ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling19)
      
  # Make a table of the models
      m19_nest_morph_t <- tab_model(m_d6_mass_19, m_d12_hb_19, m_d12_wing_19, m_d12_mass_19, m_d15_mass_19,
                                    pred.labels = c("Intercept", "Predator", "Signal Dulled", 
                                                    "Predator * Dulled", "Cross-Fostered"),
                                    dv.labels = c("Average Mass Day 6", "Head + Bill Day 12", 
                                                  "Wing Length Day 12", "Mass Day 12", "Mass Day 15"))
      
  # Save table
      saveRDS(m19_nest_morph_t, here::here("5_other_outputs/m19_nest_morph_t.RDS"))
    
## Nestling Morphology Plots ----
    
      # Make a long format object for the different time point nestling mass measures
          long_nm <- d_fem %>%
            pivot_longer(cols = c("d6_av_mass", "d12_av_mass", "d15_av_mass"), names_to = "time_point", 
                         names_transform = list(
                           time_point = ~ readr::parse_factor(.x, levels = c("d6_av_mass", "d12_av_mass", "d15_av_mass"),
                                                              ordered = TRUE)),
                         values_to = "mass", values_drop_na = TRUE)  
          long_nm <- as.data.frame(long_nm)
          nm_pos <- data.frame(time_point = c("d6_av_mass", "d12_av_mass", "d15_av_mass"),
                               xpos = seq(1, 3, 1))
          long_nm <- plyr::join(long_nm, nm_pos, "time_point", "left", "first")
          
       # Make a long format version that has mean and satndard errors for each group   
          sum_nm <- d_fem %>%
            pivot_longer(cols = c("d6_av_mass", "d12_av_mass", "d15_av_mass"), names_to = "time_point", 
                         names_transform = list(
                           time_point = ~ readr::parse_factor(.x, levels = c("d6_av_mass", "d12_av_mass", "d15_av_mass"),
                                                              ordered = TRUE)),
                         values_to = "mass", values_drop_na = TRUE) %>% 
            group_by(as.factor(year), full_treatment, time_point) %>%
            summarise(n = n(), mu = mean(mass, na.rm = TRUE), se = sd(mass, na.rm = TRUE) / sqrt(n()))
          sum_nm <- as.data.frame(sum_nm)
          sum_nm <- plyr::join(sum_nm, nm_pos, "time_point", "left", "first")
          colnames(sum_nm)[2] <- "year"
       
      # Create some plotting parameters and join to the objects above   
          for_plot <- data.frame(full_treatment = unique(sum_nm$full_treatment),
                                 p_shape = c(21, 24, 22, 21, 24, 22, 21, 24, 24),
                                 p_col = c(rep(col_sham, 3), rep(col_dull, 3), col_dull, col_sham, col_dull),
                                 p_dodge = c(-.25, -.15, -.05, .05, .15, .25, .0833, -.0833, .25))
          
          sum_nm <- plyr::join(sum_nm, for_plot, "full_treatment")
          long_nm <- plyr::join(long_nm, for_plot, "full_treatment")
      
      # Split into separate objects for each year    
          nm2_18 <- subset(long_nm, long_nm$year == "2018")
          nm2_19 <- subset(long_nm, long_nm$year == "2019")
          
          snm2_18 <- subset(sum_nm, sum_nm$year == "2018")
          snm2_19 <- subset(sum_nm, sum_nm$year == "2019")
     
      # Make a long format table like above but for headbill and wing     
          long_hbw <- d_fem %>%
            pivot_longer(cols = c("d12_av_hb", "d12_av_wing"), names_to = "measure", 
                         names_transform = list(
                           measure = ~ readr::parse_factor(.x, levels = c("d12_av_hb", "d12_av_wing"),
                                                              ordered = TRUE)),
                         values_to = "length", values_drop_na = TRUE)  
          long_hbw <- as.data.frame(long_hbw)
          hbw_pos <- data.frame(measure = c("d12_av_hb", "d12_av_wing"),
                               xpos = seq(1, 2, 1))
          long_hbw <- plyr::join(long_hbw, hbw_pos, "measure", "left", "first")
          
       # Make hte long format headbill and wing mean and standard error by each group   
          sum_hbw <- d_fem %>%
            pivot_longer(cols = c("d12_av_hb", "d12_av_wing"), names_to = "measure", 
                         names_transform = list(
                           measure = ~ readr::parse_factor(.x, levels = c("d12_av_hb", "d12_av_wing"),
                                                              ordered = TRUE)),
                         values_to = "length", values_drop_na = TRUE) %>% 
            group_by(as.factor(year), full_treatment, measure) %>%
            summarise(n = n(), mu = mean(length, na.rm = TRUE), se = sd(length, na.rm = TRUE) / sqrt(n()))
          sum_hbw <- as.data.frame(sum_hbw)
          sum_hbw <- plyr::join(sum_hbw, hbw_pos, "measure", "left", "first")
          colnames(sum_hbw)[2] <- "year"
      
      # Plotting parameters for headbill and wing    
          for_plot <- data.frame(full_treatment = unique(sum_hbw$full_treatment),
                                 p_shape = c(21, 24, 22, 21, 24, 22, 21, 24, 24),
                                 p_col = c(rep(col_sham, 3), rep(col_dull, 3), col_dull, col_sham, col_dull),
                                 p_dodge = c(-.25, -.15, -.05, .05, .15, .25, .0833, -.0833, .25))
          
          sum_hbw <- plyr::join(sum_hbw, for_plot, "full_treatment")
          long_hbw <- plyr::join(long_hbw, for_plot, "full_treatment")
     
      # Split headbill and wing into separate dataframes for each year     
          hbw2_18 <- subset(long_hbw, long_hbw$year == "2018")
          hbw2_19 <- subset(long_hbw, long_hbw$year == "2019")
          
          shbw2_18 <- subset(sum_hbw, sum_hbw$year == "2018")
          shbw2_19 <- subset(sum_hbw, sum_hbw$year == "2019")
          
          
      # Initiate the plot for 2018 nestling mass
          par(mfrow = c(1, 1))
          plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
               main = "Dulling then Challenge", xlim = c(0.5, 3.5), ylim = c(-0.1, 25),
               xlab = "", ylab = "Average Nestling Mass (g)")
          axis(1, c(-1, 1, 2, 3, 10), c("", "Day 6", "Day 12", "Day 15", ""))
          axis(2, seq(-30, 30, 5), las = 2)
          
          points(nm2_18$xpos + nm2_18$p_dodge, nm2_18$mass, pch = 16, col = alpha("black", 0.3))
          #abline(h = 0)
          
          for(i in 1:nrow(snm2_18)){
            lines(rep(snm2_18$xpos[i] + snm2_18$p_dodge[i], 2), 
                  c(snm2_18$mu[i] - snm2_18$se[i], snm2_18$mu[i] + snm2_18$se[i]), lwd = 2)
          }
          points(snm2_18$xpos + snm2_18$p_dodge, snm2_18$mu, pch = snm2_18$p_shape, cex = 1.4,
                 bg = as.character(snm2_18$p_col))
          
          abline(v = 1.5, lty = 2)
          abline(v = 2.5, lty = 2)
          
          legend("bottomright", c("No Color", "Dull Color", "Control", "Predator", "Handicap"),
                 pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex =0.9, bty = "n")
          
      # Plot for 2018 nestling headbill 
          par(mfrow = c(1, 2))
          plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
               main = "Dulling then Challenge", xlim = c(1.5, 2.5), ylim = c(15, 60),
               xlab = "", ylab = "Length on Day 12 (mm)")
          axis(1, c(-1, 1, 2, 10), c("", "Head + Bill", "Wing", ""))
          axis(2, seq(-30, 200, 5), las = 2)
          
          points(hbw2_18$xpos + hbw2_18$p_dodge, hbw2_18$length, pch = 16, col = alpha("black", 0.3))
          #abline(h = 0)
          
          for(i in 1:nrow(shbw2_18)){
            lines(rep(shbw2_18$xpos[i] + shbw2_18$p_dodge[i], 2), 
                  c(shbw2_18$mu[i] - shbw2_18$se[i], shbw2_18$mu[i] + shbw2_18$se[i]), lwd = 2)
          }
          points(shbw2_18$xpos + shbw2_18$p_dodge, shbw2_18$mu, pch = shbw2_18$p_shape, cex = 1.4,
                 bg = as.character(shbw2_18$p_col))
          
      # Plot for 2018 wing
          plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
               main = "Dulling then Challenge", xlim = c(0.5, 1.5), ylim = c(19, 29),
               xlab = "", ylab = "Length on Day 12 (mm)")
          axis(1, c(-1, 1, 2, 10), c("", "Head + Bill", "Wing", ""))
          axis(2, seq(-30, 200, 2), las = 2)
          
          points(hbw2_18$xpos + hbw2_18$p_dodge, hbw2_18$length, pch = 16, col = alpha("black", 0.3))
          #abline(h = 0)
          
          for(i in 1:nrow(shbw2_18)){
            lines(rep(shbw2_18$xpos[i] + shbw2_18$p_dodge[i], 2), 
                  c(shbw2_18$mu[i] - shbw2_18$se[i], shbw2_18$mu[i] + shbw2_18$se[i]), lwd = 2)
          }
          points(shbw2_18$xpos + shbw2_18$p_dodge, shbw2_18$mu, pch = shbw2_18$p_shape, cex = 1.4,
                 bg = as.character(shbw2_18$p_col))
          
          legend("bottomright", c("No Color", "Dull Color", "Control", "Predator", "Handicap"),
                 pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex =0.9, bty = "n")
      
      #2019 nestling mass plot
          par(mfrow = c(1, 1))
          
          plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
               main = "Challenge then Dulling", xlim = c(0.5, 3.5), ylim = c(-0.1, 25),
               xlab = "", ylab = "Average Nestling Mass (g)")
          axis(1, c(-1, 1, 2, 3, 10), c("", "Day 6", "Day 12", "Day 15", ""))
          axis(2, seq(-30, 30, 5), las = 2)
          
          points(nm2_19$xpos + nm2_19$p_dodge, nm2_19$mass, pch = 16, col = alpha("black", 0.3))
          #abline(h = 0)
          
          for(i in 1:nrow(snm2_19)){
            lines(rep(snm2_19$xpos[i] + snm2_19$p_dodge[i], 2), 
                  c(snm2_19$mu[i] - snm2_19$se[i], snm2_19$mu[i] + snm2_19$se[i]), lwd = 2)
          }
          points(snm2_19$xpos + snm2_19$p_dodge, snm2_19$mu, pch = snm2_19$p_shape, cex = 1.4,
                 bg = as.character(snm2_19$p_col))
          
          abline(v = 1.5, lty = 2)
          abline(v = 2.5, lty = 2)
          
          
          legend("bottomright", c("No Color", "Dull Color", "Control", "Predator", "Handicap"),
                 pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex =0.9, bty = "n")
          
      ## 2019 nestling headbill plot
          par(mfrow = c(1, 2))
          plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
               main = "Challenge then Dulling", xlim = c(1.5, 2.5), ylim = c(15, 60),
               xlab = "", ylab = "Length on Day 12 (mm)")
          axis(1, c(-1, 1, 2, 10), c("", "Head + Bill", "Wing", ""))
          axis(2, seq(-30, 200, 5), las = 2)
          
          points(hbw2_19$xpos + hbw2_19$p_dodge, hbw2_19$length, pch = 16, col = alpha("black", 0.3))
          #abline(h = 0)
          
          for(i in 1:nrow(shbw2_19)){
            lines(rep(shbw2_19$xpos[i] + shbw2_19$p_dodge[i], 2), 
                  c(shbw2_19$mu[i] - shbw2_19$se[i], shbw2_19$mu[i] + shbw2_19$se[i]), lwd = 2)
          }
          points(shbw2_19$xpos + shbw2_19$p_dodge, shbw2_19$mu, pch = shbw2_19$p_shape, cex = 1.4,
                 bg = as.character(shbw2_19$p_col))
          
      # 2019 nestling wing plot    
          plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
               main = "Challenge then Dulling", xlim = c(0.5, 1.5), ylim = c(19, 29),
               xlab = "", ylab = "Length on Day 12 (mm)")
          axis(1, c(-1, 1, 2, 10), c("", "Head + Bill", "Wing", ""))
          axis(2, seq(-30, 200, 2), las = 2)
          
          points(hbw2_19$xpos + hbw2_19$p_dodge, hbw2_19$length, pch = 16, col = alpha("black", 0.3))
          #abline(h = 0)
          
          for(i in 1:nrow(shbw2_19)){
            lines(rep(shbw2_19$xpos[i] + shbw2_19$p_dodge[i], 2), 
                  c(shbw2_19$mu[i] - shbw2_19$se[i], shbw2_19$mu[i] + shbw2_19$se[i]), lwd = 2)
          }
          points(shbw2_19$xpos + shbw2_19$p_dodge, shbw2_19$mu, pch = shbw2_19$p_shape, cex = 1.4,
                 bg = as.character(shbw2_19$p_col))
          
          legend("bottomright", c("No Color", "Dull Color", "Control", "Predator"),
                 pch = c(21, 21, 21, 24), pt.bg = c(col_sham, col_dull, rep("white", 2)), cex =0.9, bty = "n")
    
## Nestling Corticosterone Models ----
    # Join nestlings to adult table to get treatment info
      d_fem2 <- d_fem[, c("soc_uby", "color", "challenge")]
      d_nestling2 <- plyr::join(d_nestling, d_fem2, "soc_uby", "left", "first")
      d_nestling2 <- subset(d_nestling2, d_nestling2$raised_nest == "home" | 
                              d_nestling2$raised_nest == "cross")
      d_nestling2$raised_nest = relevel(d_nestling2$raised_nest, ref = "home")
    
    # Make separate dataframes for each year
      d_nestling18 <- subset(d_nestling2, d_nestling2$year == "2018")
      d_nestling19 <- subset(d_nestling2, d_nestling2$year == "2019")  
      
    # Fit models for mass at day 12 & 15 plus headbill & wing at day 12 for 2018
      m_bcort_18 <- lmer(d12_base ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling18)
      m_scort_18 <- lmer(d12_stress ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling18)
      m_dex_18 <- lmer(d12_dex ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling18)
      m_acth_18 <- lmer(d15_acth ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling18)
      
    # Make a table of the models
      m18_nest_cort_t <- tab_model(m_bcort_18, m_scort_18, m_dex_18, m_acth_18,
                                    pred.labels = c("Intercept", "Signal Dulled", "Predator", "Flight Reduction",
                                                    "Cross-Fostered",
                                                    "Dulled * Predator", "Dulled * Challenge"),
                                    dv.labels = c("Base Cort Day 12", "Stress-Induced Cort Day 12", "Post-Dex Cort Day 12",
                                                  "Post-Cortrosyn Cort Day 15"))
      
    # Save table
      saveRDS(m18_nest_cort_t, here::here("5_other_outputs/m18_nest_cort_t.RDS"))
      
    # Fit models for mass at day 12 & 15 plus headbill & wing at day 12 for 2019
      m_bcort_19 <- lmer(d12_base ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling19)
      m_scort_19 <- lmer(d12_stress ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling19)
      m_dex_19 <- lmer(d12_dex ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling19)
      m_acth_19 <- lmer(d15_acth ~ color*challenge + raised_nest + (1|soc_uby), data = d_nestling19)
      
    # Make a table of the models
      m19_nest_cort_t <- tab_model(m_bcort_19, m_scort_19, m_dex_19, m_acth_19,
                                   pred.labels = c("Intercept", "Predator", "Signal Dulled",
                                                   "Cross-Fostered", "Predator * Dulled"),
                                   dv.labels = c("Base Cort Day 12", "Stress-Induced Cort Day 12", "Post-Dex Cort Day 12",
                                                 "Post-Cortrosyn Cort Day 15"))
      
    # Save table
      saveRDS(m19_nest_cort_t, here::here("5_other_outputs/m19_nest_cort_t.RDS"))
    
## Nestling Corticosterone Plots ----
    # no nestling corticosterone plots made at present
## Provisioning Models ----   
      # Joining female data to provisioning 
        d_fem$uby <- paste(d_fem$unitbox, d_fem$year, sep = "_")
        d_prov2 <- plyr::join(d_provision, d_fem, "uby", "left", "first")
        d_prov3 <- subset(d_prov2, d_prov2$offset < 16 & d_prov2$doy < d_prov2$end_soc_date)
        
      # Save just the columns needed for the model
        d_prov3 <- d_prov3[, c("offset", "f_feed", "m_feed", "f_reads", "m_reads", "year", 
                               "color", "challenge", "full_treatment", "experiment",
                               "maxbrood", "doy", "f_rfid", "uby")]
        
      # Subset to days with > 10 reads to remove days of abandonment or problems with RFID board
        d_prov3 <- subset(d_prov3, d_prov3$f_reads > 10)
      
      # Subset for each year
        d_prov3_18 <- subset(d_prov3, d_prov3$year == 2018)
        d_prov3_19 <- subset(d_prov3, d_prov3$year == 2019)
        
      # Add initial brightness to each of the objects
        df18 <- d_fem18[, c("soc_uby", "bbright_s")]
        df19 <- d_fem19[, c("soc_uby", "bbright_s")]
        colnames(df18)[1] <- "uby"
        colnames(df19)[1] <- "uby"
        d_prov3_18 <- plyr::join(d_prov3_18, df18, "uby", "left", "first")
        d_prov3_19 <- plyr::join(d_prov3_19, df19, "uby", "left", "first")
        
      # Fit models for 2018
        m_feed18 <- lmer(f_feed ~ maxbrood + I(maxbrood^2) + offset + I(offset^2) + color*challenge*bbright_s + (1|doy) + (1|uby), data = d_prov3_18)
        m_feed18b <- lmer(f_feed ~ maxbrood + I(maxbrood^2) + offset + I(offset^2) + color*challenge + (1|doy) + (1|uby), data = d_prov3_18)
        m_feed18c <- lmer(f_feed ~ maxbrood + I(maxbrood^2) + offset + I(offset^2) + color*challenge + scale(m_feed) + (1|doy) + (1|uby), data = d_prov3_18)
        
        
      # Make 2018 table
        m_feed18_t <- tab_model(m_feed18, m_feed18b, m_feed18c,
                  pred.labels = c("Intercept", "Brood Size", "Brood Size^2", "Nestling Age", "Nestling Age^2", "Signal Dulled", "Predator", "Flight Reduction",
                                  "Initial Brightness", "Dulled * Predator", "Dulled * Flight", "Dulled * Brightness", "Predator * Brightness", "Flight * Brightness",
                                  "Dulled * Predator * Brightness", "Dulled * Flight * Brightness", "Male Daily Provisioning"),
                  dv.labels = rep("Daily Female Provisioning Trips", 3))
        saveRDS(m_feed18_t, here::here("5_other_outputs/m_feed18_t.RDS"))
        
      # Fit models for 2019
        m_feed19 <- lmer(f_feed ~ maxbrood + I(maxbrood^2) + offset + I(offset^2) + color*challenge*bbright_s + (1|doy) + (1|uby), data = d_prov3_19)
        m_feed19b <- lmer(f_feed ~ maxbrood + I(maxbrood^2) + offset + I(offset^2) + color*challenge + (1|doy) + (1|uby), data = d_prov3_19)
        m_feed19c <- lmer(f_feed ~ maxbrood + I(maxbrood^2) + offset + I(offset^2) + color*challenge + scale(m_feed) + (1|doy) + (1|uby), data = d_prov3_19)
        
      # Make 2019 table
        m_feed19_t <- tab_model(m_feed19, m_feed19b, m_feed19c,
                  pred.labels = c("Intercept", "Brood Size", "Brood Size^2", "Nestling Age", "Nestling Age^2", "Predator", "Signal Dulled", "Initial Brightness",
                                  "Predator * Dulled", "Predator * Brightness", "Dulled * Brightness", "Predator * Dulled * Brightness", "Male Daily Provisioning"),
                  dv.labels = rep("Daily Female Provisioning Trips", 3))
        saveRDS(m_feed19_t, here::here("5_other_outputs/m_feed19_t.RDS"))
    
## Provisioning Plots ----
  # Figure showing nestling provisioning
      # These plots are not in the style of the other plots above but this was faster to make
        # and I'm not sure that any of them will actually end up staying in the paper anyway.
    
      # Plot female provisioning from 2018 experiment
        ggplot(d_prov3_18, aes(x = offset, y = f_feed, col = full_treatment)) +
          geom_point() + geom_smooth() + xlim(0, 15) + xlab("Days after hatching") +
          ylab("Daily female provisioning trips") + ggtitle("Dulling then Challenge") +
          theme_bw()
        
      # Plot male provisioning from 2018 experiment  
        ggplot(d_prov3_18, aes(x = offset, y = m_feed, col = full_treatment)) +
          geom_point() + geom_smooth() + xlim(0, 15) + xlab("Days after hatching") +
          ylab("Daily male provisioning trips") + ggtitle("Dulling then Challenge") +
          theme_bw()
        
      # Plot female provisioning from 2019 experiment  
        ggplot(d_prov3_19, aes(x = offset, y = f_feed, col = full_treatment)) +
          geom_point() + geom_smooth() + xlim(0.6, 14.5) + xlab("Days after hatching") +
          ylab("Daily female provisioning trips") + ggtitle("Challenge then Dulling") +
          theme_bw()
      
      # Plot male provisioning from 2019 experiment  
        ggplot(d_prov3_19, aes(x = offset, y = m_feed, col = full_treatment)) +
          geom_point() + geom_smooth() + xlim(0.6, 14.5) + xlab("Days after hatching") +
          ylab("Daily male provisioning trips") + ggtitle("Challenge then Dulling") +
          theme_bw()
  
    
## Old Approach to Nest Visits Models [Not Included] ----
    
    ## Make a new data frame that has one row for each day of potential observation between the first capture and 'end_soc_date' for each bird
          # make an object with one row for each unique unit box year combination
            soc_vis <- data.frame(uby = unique(d_fem$uby))
    
          # Loop through each of these unique nest attempts
            for(i in 1:nrow(soc_vis)){
              
              # pick out the female record that matches the unit, box, year
                sub <- subset(d_fem, d_fem$uby == soc_vis$uby[i])
                
                # Extend dates basd on capture date and set 'end_soc_date' so there is one row for each day rather than one for each nest
                  # Conditionals deal with cases where nests failed early.
                    if(sub$end_soc_date[1] - sub$cap1[1] + 1 > 0){
                      dats <- seq(from = sub$cap1[1] + 1, to = sub$end_soc_date[1], by = 1)
                      temp <- data.frame(uby = rep(sub$uby, length(dats)),
                                         doy = dats)
                    
                  # For first loop, make a new object. For all other loops, bind the temperorary object onto the end of the first one
                      if(i == 1){
                        soc_visits <- temp
                      }
                      if(i > 1){
                        soc_visits <- rbind(soc_visits, temp)
                      }
                    }
                    print(i)
                  }
    
    ## Join that object back to the main female data  
        soc_visits <- plyr::join(soc_visits, d_fem, "uby", "left", "first")
    
    ## Loop through and determine unique visits and trips for each day of observation
        for(i in 1:nrow(soc_visits)){
            # Subset just non box owner trips to the focal box
              sub1 <- subset(d_social, as.character(d_social$ubox) == as.character(soc_visits$unitbox[i]) 
                             & d_social$doy == soc_visits$doy[i] & d_social$year == soc_visits$year[i])
              
            # separately subset from that male and female visitors
              subf <- subset(sub1, sub1$sex == "Female")
              subm <- subset(sub1, sub1$sex == "Male")
           
            # use those subsets to fill in the number of total and unique visitors each day 
              soc_visits$tot_f_vis[i] <- nrow(subf)
              soc_visits$tot_m_vis[i] <- nrow(subm)
              soc_visits$uni_f_vis[i] <- length(unique(subf$rfid))
              soc_visits$uni_m_vis[i] <- length(unique(subm$rfid))
              
            # Make a second subset that is all trips to other boxes by the focal female on this day
              sub2 <- subset(d_social, as.character(d_social$rfid) == as.character(soc_visits$fRFID[i])
                             & d_social$doy == soc_visits$doy[i] & d_social$year == soc_visits$year[i])
             
            # Use that subset to fill in trips made 
              soc_visits$tot_trips[i] <- nrow(sub2)
              soc_visits$uni_trips[i] <- length(unique(sub2$ubox))
            
            # print status. this takes a couple minutes to run
              print(paste(i, " of ", nrow(soc_visits), sep = ""))
        
        }
        
    ## Make a 'stage' column that indicates stage 1 (1-2 capture) vs. stage 2 (>2 capture). 
        # Note stage 2 extends past 3rd capture since rfid keeps recording.
        for(i in 1:nrow(soc_visits)){
          if(is.na(soc_visits$cap2[i]) == TRUE){soc_visits$stage[i] <- "stage1"}
          if(is.na(soc_visits$cap2[i]) == FALSE){
            if(soc_visits$doy[i] < soc_visits$cap2[i]){soc_visits$stage[i] <- "stage1"}
            if(soc_visits$doy[i] > soc_visits$cap2[i] - 1){soc_visits$stage[i] <- "stage2"}
          }
        }
        soc_visits$stage <- as.factor(soc_visits$stage)
        
    ## Cut down to the columns that are actually used for these model
        soc_models <- soc_visits[, c("uby", "doy", "year", "stage", "color", "challenge", "hatch_date", "full_treatment",
                                     "tot_f_vis", "tot_m_vis", "uni_f_vis", "uni_m_vis", "tot_trips", "uni_trips")]
        
    # Add in a number of offset days relative to hatching
        soc_models$offset <- soc_models$doy - soc_models$hatch_date
        soc_models <- subset(soc_models, soc_models$offset > -6 & soc_models$offset < 15)
        
    ggplot(data = soc_models, mapping = aes(x = offset, y = uni_m_vis, col = full_treatment)) + 
      geom_jitter(width = 0.1, col = "slateblue", size = 0.7, alpha = 0.6, height = 0) +
      geom_smooth(method = "loess") + facet_wrap(~ year)
        
    ## Split into separate frames for each year
        soc_mod18 <- subset(soc_models, soc_models$year == "2018")
        soc_mod19 <- subset(soc_models, soc_models$year == "2019")
        

        
    ## Pull out the standardized breast brightness from 2018 & 19 objects and add here. Can't standardize directly
        # because there are multiple observations per individual in this object.
        df18 <- d_fem18[, c("soc_uby", "bbright_s")]
        colnames(df18)[1] <- "uby"
        df19 <- d_fem19[, c("soc_uby", "bbright_s")]
        colnames(df19)[1] <- "uby"
        soc_mod18 <- plyr::join(soc_mod18, df18, "uby", "left", "first")
        soc_mod19 <- plyr::join(soc_mod19, df19, "uby", "left", "first")
        

        
    ## Make separate frames for stage 1 and stage 2 to allow simpler models. The reasoning for this is that the second treatment
        # isn't applied until after stage 1, so it doesn't make any sense to include it as a predictor of stage 1 behavior.
        soc_mod18s1 <- subset(soc_mod18, soc_mod18$stage == "stage1")
        soc_mod18s2 <- subset(soc_mod18, soc_mod18$stage == "stage2")
        soc_mod19s1 <- subset(soc_mod19, soc_mod19$stage == "stage1")
        soc_mod19s2 <- subset(soc_mod19, soc_mod19$stage == "stage2")
        
    ## Fit models for each response variable 2018   
        #STAGE 1
          m_mvis_18_s1 <- glmer(uni_m_vis ~ color*bbright_s + (1|uby), family = "poisson", data = soc_mod18s1)
          m_mvis_18_s1b <- glmer(uni_m_vis ~ color + (1|uby), family = "poisson", data = soc_mod18s1)
          m_fvis_18_s1 <- glmer(uni_f_vis ~ color*bbright_s + (1|uby), family = "poisson", data = soc_mod18s1)
          m_fvis_18_s1b <- glmer(uni_f_vis ~ color + (1|uby), family = "poisson", data = soc_mod18s1,
                                 control = glmerControl(optimizer = "bobyqa"))
          m_trip_18_s1 <- glmer(uni_trips ~ color*bbright_s + (1|uby), family = "poisson", data = soc_mod18s1)
          m_trip_18_s1b <- glmer(uni_trips ~ color + (1|uby), family = "poisson", data = soc_mod18s1)
        #STAGE 2
          m_mvis_18_s2 <- glmer(uni_m_vis ~ color*challenge*bbright_s + (1|uby), family = "poisson", data = soc_mod18s2,
                                control = glmerControl(optimizer = "bobyqa"))
          m_mvis_18_s2b <- glmer(uni_m_vis ~ color*challenge + (1|uby), family = "poisson", data = soc_mod18s2)
          m_fvis_18_s2 <- glmer(uni_f_vis ~ color*challenge*bbright_s + (1|uby), family = "poisson", data = soc_mod18s2,
                                control = glmerControl(optimizer = "bobyqa"))
          m_fvis_18_s2b <- glmer(uni_f_vis ~ color*challenge + (1|uby), family = "poisson", data = soc_mod18s2,
                                 control = glmerControl(optimizer = "bobyqa"))
          m_trip_18_s2 <- glmer(uni_trips ~ color*challenge*bbright_s + (1|uby), family = "poisson", data = soc_mod18s2,
                                control = glmerControl(optimizer = "bobyqa"))
          m_trip_18_s2b <- glmer(uni_trips ~ color*challenge + (1|uby), family = "poisson", data = soc_mod18s2)
          
        #STAGE 1 Table
          m_soc18_st1_t <- tab_model(m_mvis_18_s1b, m_fvis_18_s1b, m_trip_18_s1b, 
                    pred.labels = c("Intercept", "Signal Dulled"),
                    dv.labels = c("Daily Unique Male Visitors", "Daily Unique Female Visitors", "Daily Unique Trips to Other Boxes"))
          saveRDS(m_soc18_st1_t, here::here("5_other_outputs/m_soc_18_st1_t.RDS"))
          
        #STAGE 2 Table
          m_soc18_st2_t <- tab_model(m_mvis_18_s2b, m_fvis_18_s2, m_trip_18_s2b,
                    pred.labels = c("Intercept", "Signal Dulled", "Predator", "Flight Reduction",
                                    "Dulled * Predator", "Dulled * Flight", "Initial Brightness", "Dulled * Initial Bright", 
                                    "Predator * Initial Bright", "Flight * Initial Bright", "Dulled * Predator * Bright",
                                    "Dulled * Flight * Bright"),
                    dv.labels = c("Daily Unique Male Visitors", "Daily Unique Female Visitors", "Daily Unique Trips to Other Boxes"))
          saveRDS(m_soc18_st2_t, here::here("5_other_outputs/m_soc_18_st2_t.RDS"))
        
    ## Fit models for each response variable 2018   
        #STAGE 1
          m_mvis_19_s1 <- glmer(uni_m_vis ~ color*bbright_s + (1|uby), family = "poisson", data = soc_mod19s1)
          m_mvis_19_s1b <- glmer(uni_m_vis ~ color + (1|uby), family = "poisson", data = soc_mod19s1)
          m_fvis_19_s1 <- glmer(uni_f_vis ~ color*bbright_s + (1|uby), family = "poisson", data = soc_mod19s1)
          m_fvis_19_s1b <- glmer(uni_f_vis ~ color + (1|uby), family = "poisson", data = soc_mod19s1)
          m_trip_19_s1 <- glmer(uni_trips ~ color*bbright_s + (1|uby), family = "poisson", data = soc_mod19s1)
          m_trip_19_s1b <- glmer(uni_trips ~ color + (1|uby), family = "poisson", data = soc_mod19s1)
        #STAGE 2
          #m_mvis_19_s2 <- glmer(uni_m_vis ~ color*challenge*bbright_s + (1|uby), family = "poisson", data = soc_mod19s2)
          m_mvis_19_s2b <- glmer(uni_m_vis ~ color*challenge + (1|uby), family = "poisson", data = soc_mod19s2,
                                 control = glmerControl(optimizer = "bobyqa"))
          m_fvis_19_s2 <- glmer(uni_f_vis ~ color*challenge*bbright_s + (1|uby), family = "poisson", data = soc_mod19s2,
                                control = glmerControl(optimizer = "bobyqa"))
          m_fvis_19_s2b <- glmer(uni_f_vis ~ color*challenge + (1|uby), family = "poisson", data = soc_mod19s2,
                                 control = glmerControl(optimizer = "bobyqa"))
          m_trip_19_s2 <- glmer(uni_trips ~ color*challenge*bbright_s + (1|uby), family = "poisson", data = soc_mod19s2)
          m_trip_19_s2b <- glmer(uni_trips ~ color*challenge + (1|uby), family = "poisson", data = soc_mod19s2)
          
        #STAGE 1 Table
          m_soc19_st1_t <- tab_model(m_mvis_19_s1, m_fvis_19_s1b, m_trip_19_s1b, 
                    pred.labels = c("Intercept", "Predator", "Initial Brightness", "Predator * Brightness"),
                    dv.labels = c("Daily Unique Male Visitors", "Daily Unique Female Visitors", "Daily Unique Trips to Other Boxes"))
          saveRDS(m_soc19_st1_t, here::here("5_other_outputs/m_soc_19_st1_t.RDS"))
          
        #STAGE 2 Table
          m_soc19_st2_t <- tab_model(m_mvis_19_s2b, m_fvis_19_s2b, m_trip_19_s2b,
                    pred.labels = c("Intercept", "Predator", "Signal Dulled", "Predator * Dulled"),
                    dv.labels = c("Daily Unique Male Visitors", "Daily Unique Female Visitors", "Daily Unique Trips to Other Boxes"))
          saveRDS(m_soc19_st2_t, here::here("5_other_outputs/m_soc_19_st2_t.RDS"))
        
        
    
## Nest Visit Models ----
        ## Note: may depend on some objects created in old nest visits model section above
          
          ## Determine how many days of rfid records there are for each female between 1-2 capture and from 2 until end date
            for(i in 1:nrow(d_fem)){
              d_fem$days_in_1to2[i] <- d_fem$cap2[i] - d_fem$cap1[i]
              d_fem$days_in_2plus[i] <- d_fem$end_soc_date[i] + 1 - d_fem$cap2[i]
              if(is.na(d_fem$days_in_1to2[i]) == TRUE){d_fem$days_in_1to2[i] <- 0}
              if(is.na(d_fem$days_in_2plus[i]) == TRUE){d_fem$days_in_2plus[i] <- 0}
              if(d_fem$days_in_2plus[i] < 0){d_fem$days_in_2plus[i] <- 1}
              if(d_fem$days_in_2plus[i] > 17){d_fem$days_in_2plus[i] <- 17}
            }
          
          ## Calculate number of visitors or trips of different types in each stage for each female
          for(i in 1:nrow(d_fem)){
            sub <- subset(d_social, as.character(d_social$ubox) == as.character(d_fem$unitbox[i]))
            subf <- subset(sub, sub$sex == "Female")
            subm <- subset(sub, sub$sex == "Male")
            
            subf1 <- subset(subf, subf$doy >= d_fem$cap1[i] & subf$doy < d_fem$cap2[i])
            subf2 <- subset(subf, subf$doy >= d_fem$cap2[i] & subf$doy <= d_fem$cap2[i] + 17)
            
            subm1 <- subset(subm, subm$doy >= d_fem$cap1[i] & subm$doy < d_fem$cap2[i])
            subm2 <- subset(subm, subm$doy >= d_fem$cap2[i] & subm$doy <= d_fem$cap2[i] + 17)
            
            d_fem$tot_f_vis[i] <- nrow(subf)
            d_fem$tot_m_vis[i] <- nrow(subm)
            d_fem$uni_f_vis[i] <- length(unique(subf$rfid))
            d_fem$uni_m_vis[i] <- length(unique(subm$rfid))
            
            d_fem$tot_f_vis1[i] <- nrow(subf1)
            d_fem$tot_m_vis1[i] <- nrow(subm1)
            d_fem$uni_f_vis1[i] <- length(unique(subf1$rfid))
            d_fem$uni_m_vis1[i] <- length(unique(subm1$rfid))
            
            d_fem$tot_f_vis2[i] <- nrow(subf2)
            d_fem$tot_m_vis2[i] <- nrow(subm2)
            d_fem$uni_f_vis2[i] <- length(unique(subf2$rfid))
            d_fem$uni_m_vis2[i] <- length(unique(subm2$rfid))
            
            d_fem$uni_f_vis1s[i] <- length(unique(subf1$rfid)) / d_fem$days_in_1to2[i]
            d_fem$uni_m_vis1s[i] <- length(unique(subm1$rfid)) / d_fem$days_in_1to2[i]
            d_fem$uni_f_vis2s[i] <- length(unique(subf2$rfid)) / d_fem$days_in_2plus[i]
            d_fem$uni_m_vis2s[i] <- length(unique(subm2$rfid)) / d_fem$days_in_2plus[i]
            
            subt <- subset(d_social, as.character(d_social$rfid) == as.character(d_fem$fRFID[i]))
            subt1 <- subset(subt, subt$doy >= d_fem$cap1[i] & subt$doy < d_fem$cap2[i])
            subt2 <- subset(subt, subt$doy >= d_fem$cap2[i] & subt$doy <= d_fem$cap2[i] + 17)
            
            d_fem$uni_trip1[i] <- length(unique(subt1$ubox))
            d_fem$tot_trip1[i] <- nrow(subt1)
            d_fem$uni_trip2[i] <- length(unique(subt2$ubox))
            d_fem$tot_trip2[i] <- nrow(subt2)
            
          }
          
      ## Make subsets for each year and standardize plumage brightnesss
          fem18 <- subset(d_fem, d_fem$year == 2018)
          fem18$bbright_s <- scale(fem18$fbbright)
          fem19 <- subset(d_fem, d_fem$year == 2019)
          fem19$bbright_s <- scale(fem19$fbbright)
          
      ## Models for 2018
          # Fit models for each response variable 2018   
          # STAGE 1
            m_mvisu_18_s1 <- glm(uni_m_vis1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem18)
            m_mvisu_18_s1b <- glm(uni_m_vis1 ~ color + days_in_1to2, family = "quasipoisson", data = fem18)
            m_mvist_18_s1 <- glm(tot_m_vis1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem18)
            m_mvist_18_s1b <- glm(tot_m_vis1 ~ color + days_in_1to2, family = "quasipoisson", data = fem18)
            m_fvisu_18_s1 <- glm(uni_f_vis1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem18)
            m_fvisu_18_s1b <- glm(uni_f_vis1 ~ color + days_in_1to2, family = "quasipoisson", data = fem18)
            m_fvist_18_s1 <- glm(tot_f_vis1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem18)
            m_fvist_18_s1b <- glm(tot_f_vis1 ~ color + days_in_1to2, family = "quasipoisson", data = fem18)
            m_tripu_18_s1 <- glm(uni_trip1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem18)
            m_tripu_18_s1b <- glm(uni_trip1 ~ color + days_in_1to2, family = "quasipoisson", data = fem18)
            m_tript_18_s1 <- glm(tot_trip1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem18)
            m_tript_18_s1b <- glm(tot_trip1 ~ color + days_in_1to2, family = "quasipoisson", data = fem18)
          
          # STAGE 2
            m_mvisu_18_s2 <- glm(uni_m_vis2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem18)
            m_mvisu_18_s2b <- glm(uni_m_vis2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem18)
            m_mvist_18_s2 <- glm(tot_m_vis2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem18)
            m_mvist_18_s2b <- glm(tot_m_vis2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem18)
            m_fvisu_18_s2 <- glm(uni_f_vis2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem18)
            m_fvisu_18_s2b <- glm(uni_f_vis2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem18)
            m_fvist_18_s2 <- glm(tot_f_vis2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem18)
            m_fvist_18_s2b <- glm(tot_f_vis2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem18)
            m_tripu_18_s2 <- glm(uni_trip2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem18)
            m_tripu_18_s2b <- glm(uni_trip2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem18)
            m_tript_18_s2 <- glm(tot_trip2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem18)
            m_tript_18_s2b <- glm(tot_trip2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem18)
          
          # STAGE 1 Table
            m_soc18_st1_t <- tab_model(m_mvisu_18_s1b, m_fvisu_18_s1, m_tripu_18_s1b,
                         pred.labels = c("Intercept", "Signal Dulled", "Observation Effort", "Initial Brightness", "Dulled * Initial"),
                         dv.labels = c("Unique Male Visitors",
                                       "Unique Female Visitors", "Unique Trips"),
                         title = "Experiment One: Stage 1")
            saveRDS(m_soc18_st1_t, here::here("5_other_outputs/m_soc_18_st1_t.RDS"))
            
          # STAGE 2 Table
            m_soc18_st2_t <- tab_model(m_mvisu_18_s2b, m_fvisu_18_s2, m_tripu_18_s2b,
                         pred.labels = c("Intercept", "Signal Dulled", "Predator", "Flight Reduction", "Observation Effort",
                                         "Dulled * Predator", "Dulled * Flight", "Initial Brightness",
                                         "Dulled * Initial", "Predator * Initial", "Flight * Initial",
                                         "Dulled * Flight * Initial", "Dulled * Flight * Initial"),
                         dv.labels = c("Unique Male Visitors", "Unique Female Visitors", "Unique Trips"),
                         title = "Experiment One: Stage 2")
            saveRDS(m_soc18_st2_t, here::here("5_other_outputs/m_soc_18_st2_t.RDS"))

          
      ## Models for 2019
          # Fit models for each response variable 2019 
          # STAGE 1
            m_mvisu_19_s1 <- glm(uni_m_vis1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem19)
            m_mvisu_19_s1b <- glm(uni_m_vis1 ~ color + days_in_1to2, family = "quasipoisson", data = fem19)
            m_mvist_19_s1 <- glm(tot_m_vis1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem19)
            m_mvist_19_s1b <- glm(tot_m_vis1 ~ color + days_in_1to2, family = "quasipoisson", data = fem19)
            m_fvisu_19_s1 <- glm(uni_f_vis1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem19)
            m_fvisu_19_s1b <- glm(uni_f_vis1 ~ color + days_in_1to2, family = "quasipoisson", data = fem19)
            m_fvist_19_s1 <- glm(tot_f_vis1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem19)
            m_fvist_19_s1b <- glm(tot_f_vis1 ~ color + days_in_1to2, family = "quasipoisson", data = fem19)
            m_tripu_19_s1 <- glm(uni_trip1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem19)
            m_tripu_19_s1b <- glm(uni_trip1 ~ color + days_in_1to2, family = "quasipoisson", data = fem19)
            m_tript_19_s1 <- glm(tot_trip1 ~ color*bbright_s + days_in_1to2, family = "quasipoisson", data = fem19)
            m_tript_19_s1b <- glm(tot_trip1 ~ color + days_in_1to2, family = "quasipoisson", data = fem19)
            
          # STAGE 2
            m_mvisu_19_s2 <- glm(uni_m_vis2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem19)
            m_mvisu_19_s2b <- glm(uni_m_vis2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem19)
            m_mvist_19_s2 <- glm(tot_m_vis2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem19)
            m_mvist_19_s2b <- glm(tot_m_vis2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem19)
            m_fvisu_19_s2 <- glm(uni_f_vis2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem19)
            m_fvisu_19_s2b <- glm(uni_f_vis2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem19)
            m_fvist_19_s2 <- glm(tot_f_vis2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem19)
            m_fvist_19_s2b <- glm(tot_f_vis2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem19)
            m_tripu_19_s2 <- glm(uni_trip2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem19)
            m_tripu_19_s2b <- glm(uni_trip2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem19)
            m_tript_19_s2 <- glm(tot_trip2 ~ color*challenge*bbright_s + days_in_2plus, family = "quasipoisson", data = fem19)
            m_tript_19_s2b <- glm(tot_trip2 ~ color*challenge + days_in_2plus, family = "quasipoisson", data = fem19)
            
          # STAGE 1 Table
            m_soc19_st1_t <- tab_model(m_mvisu_19_s1, m_fvisu_19_s1, m_tripu_19_s1,
                                       pred.labels = c("Intercept", "Predator", "Initial Brightness", "Observation Effort", "Predator * Initial"),
                                       dv.labels = c("Unique Male Visitors",
                                                     "Unique Female Visitors", "Unique Trips"),
                                       title = "Experiment Two: Stage 1")
            saveRDS(m_soc19_st1_t, here::here("5_other_outputs/m_soc_19_st1_t.RDS"))
            
          # STAGE 2 Table
            m_soc19_st2_t <- tab_model(m_mvisu_19_s2b, m_fvisu_19_s2b, m_tripu_19_s2b,
                                       pred.labels = c("Intercept", "Predator", "Signal Dulled", "Observation Effort", "Dulled * Predator"),
                                       dv.labels = c("Unique Male Visitors", "Unique Female Visitors", "Unique Trips"),
                                       title = "Experiment Two: Stage 2")
            saveRDS(m_soc19_st2_t, here::here("5_other_outputs/m_soc_19_st2_t.RDS"))  
            
          
          
## Nest Visits Plots ----
    ##visitors to the box
        # Pivot the table longer so that unique male, female visitors by stage are each in their own row
          long_vis <- d_fem %>%
            pivot_longer(cols = c("uni_f_vis1", "uni_m_vis1", "uni_f_vis2", "uni_m_vis2"), names_to = "type", 
                         names_transform = list(
                           type = ~ readr::parse_factor(.x, levels = c("uni_f_vis1", "uni_m_vis1", "uni_f_vis2", "uni_m_vis2"),
                                                           ordered = TRUE)),
                         values_to = "count", values_drop_na = TRUE)  
          long_vis <- as.data.frame(long_vis)
          vis_pos <- data.frame(type = c("uni_f_vis1", "uni_m_vis1", "uni_f_vis2", "uni_m_vis2"),
                                xpos = seq(1, 4, 1))
          long_vis <- plyr::join(long_vis, vis_pos, "type", "left", "first")
          
        # Summarize the long table by group to get means and standard errors  
          sum_vis <- d_fem %>%
            pivot_longer(cols = c("uni_f_vis1", "uni_m_vis1", "uni_f_vis2", "uni_m_vis2"), names_to = "type", 
                         names_transform = list(
                           type = ~ readr::parse_factor(.x, levels = c("uni_f_vis1", "uni_m_vis1", "uni_f_vis2", "uni_m_vis2"),
                                                           ordered = TRUE)),
                         values_to = "count", values_drop_na = TRUE) %>% 
            group_by(as.factor(year), full_treatment, type) %>%
            summarise(n = n(), mu = mean(count, na.rm = TRUE), se = sd(count, na.rm = TRUE) / sqrt(n()))
          sum_vis <- as.data.frame(sum_vis)
          sum_vis <- plyr::join(sum_vis, vis_pos, "type", "left", "first")
          colnames(sum_vis)[2] <- "year"
        
        # Make plotting parameters for each group  
          for_plot <- data.frame(full_treatment = unique(sum_vis$full_treatment),
                                 p_shape = c(21, 24, 22, 21, 24, 22, 21, 24, 24),
                                 p_col = c(rep(col_sham, 3), rep(col_dull, 3), col_dull, col_sham, col_dull),
                                 p_dodge = c(-.25, -.15, -.05, .05, .15, .25, .0833, -.0833, .25))
         
        # Join plotting parameters back to long dataframe and make subsets for each year 
          sum_vis <- plyr::join(sum_vis, for_plot, "full_treatment")
          long_vis <- plyr::join(long_vis, for_plot, "full_treatment")
          
          vis2_18 <- subset(long_vis, long_vis$year == "2018")
          vis2_19 <- subset(long_vis, long_vis$year == "2019")
          
          svis2_18 <- subset(sum_vis, sum_vis$year == "2018")
          svis2_19 <- subset(sum_vis, sum_vis$year == "2019")
    
    ## Visitors
    
        # Plot visitors by stage for 2018
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Dulling then Challenge", xlim = c(0.5, 4.5), ylim = c(-2, 25),
                 xlab = "", ylab = "Number of Unique Visitors")
            axis(1, c(-1, 1, 2, 3, 4, 10), c("", "Females 1", "Males 1", "Females 2", "Males 2", ""))
            axis(2, seq(-30, 200, 5), las = 2)
            abline(h = 0)
            
            points(vis2_18$xpos + vis2_18$p_dodge, vis2_18$count, pch = 16, col = alpha("black", 0.3))
            #abline(h = 0)
            
            for(i in 1:nrow(svis2_18)){
              lines(rep(svis2_18$xpos[i] + svis2_18$p_dodge[i], 2), 
                    c(svis2_18$mu[i] - svis2_18$se[i], svis2_18$mu[i] + svis2_18$se[i]), lwd = 2)
            }
            points(svis2_18$xpos + svis2_18$p_dodge, svis2_18$mu, pch = svis2_18$p_shape, cex = 1.4,
                   bg = as.character(svis2_18$p_col))
            
            legend(.5, 23, c("No Color", "Dull Color", "Control", "Predator", "Handicap"),
                   pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex =0.9, bty = "n")
            rect(0.5, 23.5, 4.5, 25, col = alpha("chartreuse3", 0.3))
            text(2.3, 24.25, "Dulling", pos = 4)
            
            rect(2.5, 21.5, 4.5, 23, col = alpha("coral3", 0.2))
            text(3.2, 22.25, "Challenge", pos = 4)
            
            lines(c(2.5, 2.5), c(-5, 23.5), lty = 2)
            
            text(1.5, -1, "First to second capture")
            text(3.5, -1, "After second capture")
        
        # Plot visitors by stage for 2019    
              plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                   main = "Challenge then Dulling", xlim = c(0.5, 4.5), ylim = c(-2, 25),
                   xlab = "", ylab = "Number of Unique Visitors")
              axis(1, c(-1, 1, 2, 3, 4, 10), c("", "Females 1", "Males 1", "Females 2", "Males 2", ""))
              axis(2, seq(-30, 200, 5), las = 2)
              abline(h = 0)
              
              points(vis2_19$xpos + vis2_19$p_dodge, vis2_19$count, pch = 16, col = alpha("black", 0.3))
              #abline(h = 0)
              
              for(i in 1:nrow(svis2_19)){
                lines(rep(svis2_19$xpos[i] + svis2_19$p_dodge[i], 2), 
                      c(svis2_19$mu[i] - svis2_19$se[i], svis2_19$mu[i] + svis2_19$se[i]), lwd = 2)
              }
              points(svis2_19$xpos + svis2_19$p_dodge, svis2_19$mu, pch = svis2_19$p_shape, cex = 1.4,
                     bg = as.character(svis2_19$p_col))
              
              legend(.5, 22, c("No Color", "Dull Color", "Control", "Predator"),
                     pch = c(21, 21, 21, 24), pt.bg = c(col_sham, col_dull, rep("white", 2)), cex =0.9, bty = "n")
              rect(2.5, 23.5, 4.5, 25, col = alpha("chartreuse3", 0.3))
              text(3.3, 24.25, "Dulling", pos = 4)
              
              rect(0.5, 21.5, 2.5, 23, col = alpha("coral3", 0.2))
              text(1.2, 22.25, "Challenge", pos = 4)
              
              lines(c(2.5, 2.5), c(-5, 23.5), lty = 2)
              
              text(1.5, -1, "First to second capture")
              text(3.5, -1, "After second capture")
              
    ## Plotting model for unique female visitors from 2018 per day models
              
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Dulling then Challenge", xlim = c(32, 72), ylim = c(0, 3),
                 xlab = "Initial Plumage Brightness", ylab = "Unique Female Visitors Per Day")
            axis(1, seq(-10, 110, 5))
            axis(2, seq(-2, 50, 1), las = 2)
            
            r <- seq(min(na.omit(soc_mod18s2$bbright_s)), 
                     max(na.omit(soc_mod18s2$bbright_s)),
                     0.1)
            
            bb_mu18 <- mean(na.omit(d_fem18$fbbright))
            bb_sd18 <- sd(na.omit(d_fem18$fbbright))
            r2 <- (r * bb_sd18) + bb_mu18
            
            mu_con_con <- sapply(r, function(z)exp(fixef(m_fvis_18_s2)[1] +
                                    fixef(m_fvis_18_s2)[5]*z))
            mu_con_pred <- sapply(r, function(z)exp(fixef(m_fvis_18_s2)[1] +
                                    fixef(m_fvis_18_s2)[3] +
                                    fixef(m_fvis_18_s2)[5]*z +
                                    fixef(m_fvis_18_s2)[9]*z))
            mu_con_str <- sapply(r, function(z)exp(fixef(m_fvis_18_s2)[1] +
                                    fixef(m_fvis_18_s2)[4] +
                                    fixef(m_fvis_18_s2)[5]*z +
                                    fixef(m_fvis_18_s2)[10]*z))
            mu_dull_con <- sapply(r, function(z)exp(fixef(m_fvis_18_s2)[1] +
                                    fixef(m_fvis_18_s2)[2] +
                                    fixef(m_fvis_18_s2)[5]*z +
                                    fixef(m_fvis_18_s2)[8]*z))
            mu_dull_pred <- sapply(r, function(z)exp(fixef(m_fvis_18_s2)[1] +
                                    fixef(m_fvis_18_s2)[2] +
                                    fixef(m_fvis_18_s2)[3] +
                                    fixef(m_fvis_18_s2)[5]*z +
                                    fixef(m_fvis_18_s2)[6] +
                                    fixef(m_fvis_18_s2)[8]*z +
                                    fixef(m_fvis_18_s2)[9]*z +
                                    fixef(m_fvis_18_s2)[11]*z))
            mu_dull_str <- sapply(r, function(z)exp(fixef(m_fvis_18_s2)[1] +
                                    fixef(m_fvis_18_s2)[2] +
                                    fixef(m_fvis_18_s2)[4] +
                                    fixef(m_fvis_18_s2)[5]*z +
                                    fixef(m_fvis_18_s2)[7] +
                                    fixef(m_fvis_18_s2)[8]*z +
                                    fixef(m_fvis_18_s2)[10]*z +
                                    fixef(m_fvis_18_s2)[12]*z))
            
            
            m_fvis_18_s2 <- glmer(uni_f_vis ~ color*challenge*bbright_s + (1|uby), family = "poisson", data = soc_mod18s2,
                                  control = glmerControl(optimizer = "bobyqa"))
            
            
            lines(r2, mu_con_con, lwd = 2, col = col_sham)
            lines(r2, mu_con_pred, lwd = 2, col = col_sham, lty = 2)
            lines(r2, mu_con_str, lwd = 2, col = col_sham, lty = 3)
            lines(r2, mu_dull_con, lwd = 2, col = col_dull)
            lines(r2, mu_dull_pred, lwd = 2, col = col_dull, lty = 2)
            lines(r2, mu_dull_str, lwd = 2, col = col_dull, lty = 3)
            
      
            
          ## Plotting model for unique female visitors from 2018 total by stage
            
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Dulling then Challenge", xlim = c(32, 72), ylim = c(0, 20),
                 xlab = "Initial Plumage Brightness", ylab = "Unique Female Visitors Stage 1")
            axis(1, seq(-10, 110, 5))
            axis(2, seq(-5, 50, 5), las = 2)
            
            r <- seq(min(na.omit(soc_mod18s1$bbright_s)), 
                     max(na.omit(soc_mod18s1$bbright_s)),
                     0.1)
            
            bb_mu18 <- mean(na.omit(d_fem18$fbbright))
            bb_sd18 <- sd(na.omit(d_fem18$fbbright))
            r2 <- (r * bb_sd18) + bb_mu18
            
            mu_con <- sapply(r, function(z)exp(coef(m_fvisu_18_s1)[1] +
                                                 coef(m_fvisu_18_s1)[3]*z + coef(m_fvisu_18_s1)[4]*8))
            mu_dull <- sapply(r, function(z)exp(coef(m_fvisu_18_s1)[1] +
                                                  coef(m_fvisu_18_s1)[2] +
                                                  coef(m_fvisu_18_s1)[3]*z +
                                                  coef(m_fvisu_18_s1)[4]*8 +
                                                  coef(m_fvisu_18_s1)[5]*z))
            
            post <- mvrnorm(n = 1e5, mu = coef(m_fvisu_18_s1), Sigma = vcov(m_fvisu_18_s1))
            
            ci_con <- sapply(r, function(z)HPDI(exp(post[, 1] +
                                                      post[, 3]*z + post[, 4]*8)))
            ci_dull <- sapply(r, function(z)HPDI(exp(post[, 1] +
                                                       post[, 2] +
                                                       post[, 3]*z +
                                                       post[, 4]*8 +
                                                       post[, 5]*z)))
            
            
            shade(ci_con, r2, col = alpha(col_sham, 0.4))
            shade(ci_dull, r2, col = alpha(col_dull, 0.4))
            
            lines(r2, mu_con, lwd = 2, col = col_sham)
            lines(r2, mu_dull, lwd = 2, col = col_dull)
            
        ## Plotting model for unique female visitors from 2018 total by stage
            
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Dulling then Challenge", xlim = c(32, 72), ylim = c(0, 40),
                 xlab = "Initial Plumage Brightness", ylab = "Unique Female Visitors Stage 2")
            axis(1, seq(-10, 110, 5))
            axis(2, seq(-5, 50, 5), las = 2)
            
            r <- seq(min(na.omit(soc_mod18s2$bbright_s)), 
                     max(na.omit(soc_mod18s2$bbright_s)),
                     0.1)
            
            bb_mu18 <- mean(na.omit(d_fem18$fbbright))
            bb_sd18 <- sd(na.omit(d_fem18$fbbright))
            r2 <- (r * bb_sd18) + bb_mu18
            
            mu_con_con <- sapply(r, function(z)exp(coef(m_fvisu_18_s2)[1] +
                                                     coef(m_fvisu_18_s2)[5]*z + coef(m_fvisu_18_s2)[6]*13))
            mu_con_con <- sapply(r, function(z)exp(coef(m_fvisu_18_s2)[1] +
                                                     coef(m_fvisu_18_s2)[3] +
                                                     coef(m_fvisu_18_s2)[5]*z + 
                                                     coef(m_fvisu_18_s2)[6]*13 +
                                                      coef(m_fvisu_18_s2)[10]*z))
            mu_con_con <- sapply(r, function(z)exp(coef(m_fvisu_18_s2)[1] +
                                                     coef(m_fvisu_18_s2)[4] +
                                                     coef(m_fvisu_18_s2)[5]*z + 
                                                     coef(m_fvisu_18_s2)[6]*13 +
                                                      coef(m_fvisu_18_s2)[11]*z))
            mu_dul_con <- sapply(r, function(z)exp(coef(m_fvisu_18_s2)[1] +
                                                     coef(m_fvisu_18_s2)[2] +
                                                     coef(m_fvisu_18_s2)[5]*z +
                                                     coef(m_fvisu_18_s2)[6]*13 +
                                                     coef(m_fvisu_18_s2)[9]*z))
            mu_dul_pred <- sapply(r, function(z)exp(coef(m_fvisu_18_s2)[1] +
                                                      coef(m_fvisu_18_s2)[2] +
                                                      coef(m_fvisu_18_s2)[3] +
                                                      coef(m_fvisu_18_s2)[5]*z +
                                                      coef(m_fvisu_18_s2)[6]*13 +
                                                      coef(m_fvisu_18_s2)[7] +
                                                      coef(m_fvisu_18_s2)[9]*z +
                                                      coef(m_fvisu_18_s2)[10]*z +
                                                      coef(m_fvisu_18_s2)[12]*z))
            mu_dul_str <- sapply(r, function(z)exp(coef(m_fvisu_18_s2)[1] +
                                                     coef(m_fvisu_18_s2)[2] +
                                                     coef(m_fvisu_18_s2)[4] +
                                                     coef(m_fvisu_18_s2)[5]*z +
                                                     coef(m_fvisu_18_s2)[6]*13 +
                                                     coef(m_fvisu_18_s2)[8] +
                                                     coef(m_fvisu_18_s2)[9]*z +
                                                     coef(m_fvisu_18_s2)[11]*z +
                                                     coef(m_fvisu_18_s2)[13]*z))
            
            post <- mvrnorm(n = 1e5, mu = coef(m_fvisu_18_s2), Sigma = vcov(m_fvisu_18_s2))
            
            ci_con_con <- sapply(r, function(z)HPDI(exp(post[, 1] +
                                                      post[, 5]*z + post[, 6]*13)))
            ci_dul_con <- sapply(r, function(z)HPDI(exp(post[, 1] +
                                                       post[, 2] +
                                                       post[, 5]*z +
                                                       post[, 6]*13 +
                                                       post[, 9]*z)))
            
            
            shade(ci_con_con, r2, col = alpha(col_sham, 0.2))
            shade(ci_dul_con, r2, col = alpha(col_dull, 0.2))
            
            lines(r2, mu_con_con, lwd = 2, col = col_sham, lty = 1)
            lines(r2, mu_con_pred, lwd = 2, col = col_sham, lty = 2)
            lines(r2, mu_con_str, lwd = 2, col = col_sham, lty = 3)
            lines(r2, mu_dul_con, lwd = 2, col = col_dull, lty = 1)
            lines(r2, mu_dul_pred, lwd = 2, col = col_dull, lty = 2)
            lines(r2, mu_dul_str, lwd = 2, col = col_dull, lty = 3)
    
    
    #### Plotting same as above but per day
    
    ##visitors to the box
    
            long_vis <- d_fem %>%
              pivot_longer(cols = c("uni_f_vis1s", "uni_m_vis1s", "uni_f_vis2s", "uni_m_vis2s"), names_to = "type", 
                           names_transform = list(
                             type = ~ readr::parse_factor(.x, levels = c("uni_f_vis1s", "uni_m_vis1s", "uni_f_vis2s", "uni_m_vis2s"),
                                                          ordered = TRUE)),
                           values_to = "count", values_drop_na = TRUE)  
            long_vis <- as.data.frame(long_vis)
            vis_pos <- data.frame(type = c("uni_f_vis1s", "uni_m_vis1s", "uni_f_vis2s", "uni_m_vis2s"),
                                  xpos = seq(1, 4, 1))
            long_vis <- plyr::join(long_vis, vis_pos, "type", "left", "first")
            
            
            sum_vis <- d_fem %>%
              pivot_longer(cols = c("uni_f_vis1s", "uni_m_vis1s", "uni_f_vis2s", "uni_m_vis2s"), names_to = "type", 
                           names_transform = list(
                             type = ~ readr::parse_factor(.x, levels = c("uni_f_vis1s", "uni_m_vis1s", "uni_f_vis2s", "uni_m_vis2s"),
                                                          ordered = TRUE)),
                           values_to = "count", values_drop_na = TRUE) %>% 
              group_by(as.factor(year), full_treatment, type) %>%
              summarise(n = n(), mu = mean(count, na.rm = TRUE), se = sd(count, na.rm = TRUE) / sqrt(n()))
            sum_vis <- as.data.frame(sum_vis)
            sum_vis <- plyr::join(sum_vis, vis_pos, "type", "left", "first")
            colnames(sum_vis)[2] <- "year"
            
            for_plot <- data.frame(full_treatment = unique(sum_vis$full_treatment),
                                   p_shape = c(21, 24, 22, 21, 24, 22, 21, 24, 24),
                                   p_col = c(rep(col_sham, 3), rep(col_dull, 3), col_dull, col_sham, col_dull),
                                   p_dodge = c(-.25, -.15, -.05, .05, .15, .25, .0833, -.0833, .25))
            
            sum_vis <- plyr::join(sum_vis, for_plot, "full_treatment")
            long_vis <- plyr::join(long_vis, for_plot, "full_treatment")
            
            vis2_18 <- subset(long_vis, long_vis$year == "2018")
            vis2_19 <- subset(long_vis, long_vis$year == "2019")
            
            svis2_18 <- subset(sum_vis, sum_vis$year == "2018")
            svis2_19 <- subset(sum_vis, sum_vis$year == "2019")
    
    ## Visitors
    
    
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Dulling then Challenge", xlim = c(0.5, 4.5), ylim = c(-0.5, 4),
                 xlab = "", ylab = "Unique Visitors / Days Recorded")
            axis(1, c(-1, 1, 2, 3, 4, 10), c("", "Females 1", "Males 1", "Females 2", "Males 2", ""))
            axis(2, seq(-30, 200, 1), las = 2)
            abline(h = 0)
            
            points(vis2_18$xpos + vis2_18$p_dodge, vis2_18$count, pch = 16, col = alpha("black", 0.3))
            #abline(h = 0)
            
            for(i in 1:nrow(svis2_18)){
              lines(rep(svis2_18$xpos[i] + svis2_18$p_dodge[i], 2), 
                    c(svis2_18$mu[i] - svis2_18$se[i], svis2_18$mu[i] + svis2_18$se[i]), lwd = 2)
            }
            points(svis2_18$xpos + svis2_18$p_dodge, svis2_18$mu, pch = svis2_18$p_shape, cex = 1.4,
                   bg = as.character(svis2_18$p_col))
            
            legend(.5, 3.7, c("No Color", "Dull Color", "Control", "Predator", "Handicap"),
                   pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex =0.9, bty = "n")
            rect(0.5, 3.7, 4.5, 4, col = alpha("chartreuse3", 0.3))
            text(2.3, 3.85, "Dulling", pos = 4)
            
            rect(2.5, 3.3, 4.5, 3.6, col = alpha("coral3", 0.2))
            text(3.2, 3.45, "Challenge", pos = 4)
            
            lines(c(2.5, 2.5), c(-5, 3.7), lty = 2)
            
            text(1.5, -.3, "First to second capture")
            text(3.5, -.3, "After second capture")
            
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Challenge then Dulling", xlim = c(0.5, 4.5), ylim = c(-0.5, 4),
                 xlab = "", ylab = "Unique Visitors / Days Recorded")
            axis(1, c(-1, 1, 2, 3, 4, 10), c("", "Females 1", "Males 1", "Females 2", "Males 2", ""))
            axis(2, seq(-30, 200, 1), las = 2)
            abline(h = 0)
            
            points(vis2_19$xpos + vis2_19$p_dodge, vis2_19$count, pch = 16, col = alpha("black", 0.3))
            #abline(h = 0)
            
            for(i in 1:nrow(svis2_19)){
              lines(rep(svis2_19$xpos[i] + svis2_19$p_dodge[i], 2), 
                    c(svis2_19$mu[i] - svis2_19$se[i], svis2_19$mu[i] + svis2_19$se[i]), lwd = 2)
            }
            points(svis2_19$xpos + svis2_19$p_dodge, svis2_19$mu, pch = svis2_19$p_shape, cex = 1.4,
                   bg = as.character(svis2_19$p_col))
            
            legend(.5, 22, c("No Color", "Dull Color", "Control", "Predator"),
                   pch = c(21, 21, 21, 24), pt.bg = c(col_sham, col_dull, rep("white", 2)), cex =0.9, bty = "n")
            rect(2.5, 23.5, 4.5, 25, col = alpha("chartreuse3", 0.3))
            text(3.3, 24.25, "Dulling", pos = 4)
            
            rect(0.5, 21.5, 2.5, 23, col = alpha("coral3", 0.2))
            text(1.2, 22.25, "Challenge", pos = 4)
            
            lines(c(2.5, 2.5), c(-5, 23.5), lty = 2)
            
            text(1.5, -1, "First to second capture")
            text(3.5, -1, "After second capture")
    
  
## Trips to other boxes Models ----
    # Fit above in the visits tab since I put them into the same output table
    
## Trips to other boxes Plots ----    
    ## trips to other boxes
    
    
    for(i in 1:nrow(d_fem)){
      sub <- subset(d_social, as.character(d_social$rfid) == as.character(d_fem$fRFID[i]))
      
      sub1 <- subset(sub, sub$doy >= d_fem$cap1[i] & sub$doy < d_fem$cap2[i])
      sub2 <- subset(sub, sub$doy >= d_fem$cap2[i])
      
      d_fem$tot_trip1[i] <- nrow(sub1)
      d_fem$tot_trip2[i] <- nrow(sub2)
      d_fem$uni_trip1[i] <- length(unique(sub1$ubox))
      d_fem$uni_trip2[i] <- length(unique(sub2$ubox))
      
    }
    
    ##trips to other boxes
    
    long_vis <- d_fem %>%
      pivot_longer(cols = c("uni_trip1", "uni_trip2"), names_to = "type", 
                   names_transform = list(
                     type = ~ readr::parse_factor(.x, levels = c("uni_trip1", "uni_trip2"),
                                                  ordered = TRUE)),
                   values_to = "count", values_drop_na = TRUE)  
    long_vis <- as.data.frame(long_vis)
    vis_pos <- data.frame(type = c("uni_trip1", "uni_trip2"),
                          xpos = seq(1, 4, 1))
    long_vis <- plyr::join(long_vis, vis_pos, "type", "left", "first")
    
    
    sum_vis <- d_fem %>%
      pivot_longer(cols = c("uni_trip1", "uni_trip2"), names_to = "type", 
                   names_transform = list(
                     type = ~ readr::parse_factor(.x, levels = c("uni_trip1", "uni_trip2"),
                                                  ordered = TRUE)),
                   values_to = "count", values_drop_na = TRUE) %>% 
      group_by(as.factor(year), full_treatment, type) %>%
      summarise(n = n(), mu = mean(count, na.rm = TRUE), se = sd(count, na.rm = TRUE) / sqrt(n()))
    sum_vis <- as.data.frame(sum_vis)
    sum_vis <- plyr::join(sum_vis, vis_pos, "type", "left", "first")
    colnames(sum_vis)[2] <- "year"
    
    for_plot <- data.frame(full_treatment = unique(sum_vis$full_treatment),
                           p_shape = c(21, 24, 22, 21, 24, 22, 21, 24, 24),
                           p_col = c(rep(col_sham, 3), rep(col_dull, 3), col_dull, col_sham, col_dull),
                           p_dodge = c(-.25, -.15, -.05, .05, .15, .25, .0833, -.0833, .25))
    
    sum_vis <- plyr::join(sum_vis, for_plot, "full_treatment")
    long_vis <- plyr::join(long_vis, for_plot, "full_treatment")
    
    vis2_18 <- subset(long_vis, long_vis$year == "2018")
    vis2_19 <- subset(long_vis, long_vis$year == "2019")
    
    svis2_18 <- subset(sum_vis, sum_vis$year == "2018")
    svis2_19 <- subset(sum_vis, sum_vis$year == "2019")
    
    ## Trips
    
    
    plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
         main = "Dulling then Challenge", xlim = c(0.5, 2.5), ylim = c(-.1, 25),
         xlab = "", ylab = "Trips to Unique Boxes")
    axis(1, c(-1, 1, 2, 10), c("", "Capture 1-2", "After Capture 2", ""))
    axis(2, seq(-30, 200, 5), las = 2)
    
    points(vis2_18$xpos + vis2_18$p_dodge, vis2_18$count, pch = 16, col = alpha("black", 0.3))
    #abline(h = 0)
    
    for(i in 1:nrow(svis2_18)){
      lines(rep(svis2_18$xpos[i] + svis2_18$p_dodge[i], 2), 
            c(svis2_18$mu[i] - svis2_18$se[i], svis2_18$mu[i] + svis2_18$se[i]), lwd = 2)
    }
    points(svis2_18$xpos + svis2_18$p_dodge, svis2_18$mu, pch = svis2_18$p_shape, cex = 1.4,
           bg = as.character(svis2_18$p_col))
    
    legend(.5, 23, c("No Color", "Dull Color", "Control", "Predator", "Handicap"),
           pch = c(21, 21, 21, 24, 22), pt.bg = c(col_sham, col_dull, rep("white", 3)), cex =0.9, bty = "n")
    rect(0.5, 23.5, 2.5, 25, col = alpha("chartreuse3", 0.3))
    text(1.3, 24.25, "Dulling", pos = 4)
    
    rect(1.5, 21.5, 2.5, 23, col = alpha("coral3", 0.2))
    text(1.9, 22.25, "Challenge", pos = 4)
    
    lines(c(1.5, 1.5), c(-5, 23.5), lty = 2)
    
  
    
    plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
         main = "Challenge then Dulling", xlim = c(0.5, 2.5), ylim = c(-.1, 25),
         xlab = "", ylab = "Trips to Unique Boxes")
    axis(1, c(-1, 1, 2, 10), c("", "Capture 1-2", "After Capture 2", ""))
    axis(2, seq(-30, 200, 5), las = 2)
    
    points(vis2_19$xpos + vis2_19$p_dodge, vis2_19$count, pch = 16, col = alpha("black", 0.3))
    #abline(h = 0)
    
    for(i in 1:nrow(svis2_19)){
      lines(rep(svis2_19$xpos[i] + svis2_19$p_dodge[i], 2), 
            c(svis2_19$mu[i] - svis2_19$se[i], svis2_19$mu[i] + svis2_19$se[i]), lwd = 2)
    }
    points(svis2_19$xpos + svis2_19$p_dodge, svis2_19$mu, pch = svis2_19$p_shape, cex = 1.4,
           bg = as.character(svis2_19$p_col))
    
    legend(.5, 21.5, c("No Color", "Dull Color", "Control", "Predator"),
           pch = c(21, 21, 21, 24), pt.bg = c(col_sham, col_dull, rep("white", 2)), cex =0.9, bty = "n")
    rect(1.5, 23.5, 2.5, 25, col = alpha("chartreuse3", 0.3))
    text(1.9, 24.25, "Dulling", pos = 4)
    
    rect(0.5, 21.5, 1.5, 23, col = alpha("coral3", 0.2))
    text(0.9, 22.25, "Challenge", pos = 4)
    
    lines(c(1.5, 1.5), c(-5, 23.5), lty = 2)

## ACTH Validation Models ----
    
    # Nestlings
      nest_l <- d_acth_nest %>%
        pivot_longer(cols = c("cort1", "cort2", "cort3"), names_to = "timepoint",
                     values_to = "cort", values_drop_na = TRUE)
      nest_l <- as.data.frame(nest_l)
      nest_l$treatment <- as.factor(nest_l$treatment)
      nest_l$treatment <- relevel(nest_l$treatment, ref = "Saline")
      m_n_acth <- lmer(cort ~ timepoint*treatment + (1|band) + (1|unit_box), data = nest_l)
      
    # Adults
      fem_l <- d_acth_fem %>%
        pivot_longer(cols = c("cort1", "cort2", "cort3"), names_to = "timepoint",
                     values_to = "cort", values_drop_na = TRUE)
      fem_l <- as.data.frame(fem_l)
      fem_l$treatment <- gsub("BCC", "Saline", fem_l$treatment)
      fem_l$treatment <- gsub("BCA", "ACTH", fem_l$treatment)
      fem_l$treatment <- as.factor(fem_l$treatment)
      fem_l$treatment <- relevel(fem_l$treatment, ref = "Saline")
      m_a_acth <- lmer(cort ~ timepoint*treatment + (1|band), data = fem_l)
      
    # Make a table
      tab_model(m_n_acth, m_a_acth, 
                pred.labels = c("Intercept (Baseline)", "Timepoint 2", "Timepoint 3",
                                "Cortrosyn at Baseline", "Cortrosyn at Timepoint 2", "Cortrosyn at Timepoint 3"),
                dv.labels = c("Nestling Corticosterone", "Adult Corticosterone"))
      
      # Note that I want to put this table in the pdf output supplementary materials, but there is no direct way to 
      # use sjplot to put html tables into a markdown pdf. I manually saved the html as a pdf to put it in. That
      # means if this code is modified the saved pdf file needs to be overwritten with a new version.
    
## ACTH Validation Plots ----
    # These are plotted and saved to file here then pasted into the rmarkdown file for the supplemental materials
    
    # Nestlings
      nest_a <- d_acth_nest %>%
        pivot_longer(cols = c("cort1", "cort2", "cort3"), names_to = "timepoint",
                     values_to = "cort", values_drop_na = TRUE) %>%
        ggplot(aes(x = timepoint, y = cort, fill = treatment)) + 
            geom_line(mapping = aes(x = timepoint, y = cort, group = band, color = treatment), alpha = 0.45) +
            geom_boxplot(width = 0.25) + theme_classic() + xlab("Minutes After Capture") + ylab(expression(paste("Corticosterone ng/", mu, "l"))) +
            scale_x_discrete(labels = c("<3", "15", "30")) + geom_vline(xintercept = 1.15, lty = 2, col = "gray40") + 
            annotate("text", x = 1.1, y = 40, label = "Cortrosyn or Saline Injection", angle = 90) + labs(fill = "Treatment") + 
            ggtitle("15 Day Old Nestlings") +
            scale_fill_discrete(name = "Treatment", labels = c("Cortrosyn", "Saline")) + guides(color = FALSE) + 
            theme(legend.position = c(0.12, 0.9))
      ggsave(here::here("3_r_scripts/cortrosyn_nestlings.pdf"), plot = nest_a, width = 6, height = 5, units = "in", device = "pdf")
    
    # Adults
      fem_a <- d_acth_fem %>%
        pivot_longer(cols = c("cort1", "cort2", "cort3"), names_to = "timepoint",
                     values_to = "cort", values_drop_na = TRUE) %>%
        ggplot(aes(x = timepoint, y = cort, fill = treatment)) + 
        geom_line(mapping = aes(x = timepoint, y = cort, group = band, color = treatment), alpha = 0.45) +
        geom_boxplot(width = 0.25, outlier.size = 0, outlier.stroke = 0) + theme_classic() + xlab("Minutes After Capture") + 
        ylab(expression(paste("Corticosterone ng/", mu, "l"))) +
        scale_x_discrete(labels = c("<3", "30", "60")) + geom_vline(xintercept = 1.15, lty = 2, col = "gray40") + 
        annotate("text", x = 1.1, y = 45, label = "Saline Injection", angle = 90) + 
        geom_vline(xintercept = 2.15, lty = 2, col = "gray40") +
        annotate("text", x = 2.1, y = 50, label = "Cortrosyn or Saline Injection", angle = 90) +
        labs(fill = "Treatment") + ggtitle("Adult Females") +
        scale_fill_discrete(name = "Treatment", labels = c("Cortrosyn", "Saline")) + guides(color = FALSE) + 
        theme(legend.position = c(0.12, 0.9))
      ggsave(here::here("3_r_scripts/cortrosyn_adults.pdf"), plot = fem_a, width = 6, height = 5, units = "in", device = "pdf")
      

    

    
    
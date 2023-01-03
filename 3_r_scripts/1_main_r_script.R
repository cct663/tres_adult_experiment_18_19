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
    d_nestling <- read.delim(here::here("1_raw_data", "data_by_nestling.txt"), header = TRUE)
    
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
        png(here::here("3_r_scripts/mass_figure.png"), width = 10, height = 5.7, units = "in", res = 300)
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
                 main = "Experiment One", xlim = c(0.3, 3.7), ylim = c(15.5, 25),
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
            rect(1.4, 15.7, 3.5, 16.2, col = alpha("coral3", 0.2))
            text(1.7, 15.95, "Signal Manipulation", pos = 4, cex = 0.9)
            
            rect(2.2, 16.3, 3.5, 16.8, col = alpha("chartreuse3", 0.3))
            text(2.4, 16.55, "Challenge", pos = 4, cex = 0.9)
        
        #Initiate Plot for 2018
            plot(1, 1, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
                 main = "Experiment Two", xlim = c(0.3, 3.7), ylim = c(15.5, 25),
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
            rect(2.3, 15.7, 3.6, 16.2, col = alpha("coral3", 0.2))
            text(2.25, 15.95, "Signal Manipulation", pos = 4, cex = 0.9)
            
            rect(1.4, 16.3, 2.5, 16.8, col = alpha("chartreuse3", 0.3))
            text(1.5, 16.55, "Challenge", pos = 4, cex = 0.9)
            
    dev.off() # plot is saved to folder        
    
    
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
        # not included
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
          rect(4.4, -15, 9.6, -10, col = alpha("coral3", 0.2))
          text(6.5, -12.5, "Signal Manipulation", pos = 4)
          
          rect(6.4, -9, 9.6, -4, col = alpha("chartreuse3", 0.3))
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
          rect(6.4, -15, 9.6, -10, col = alpha("coral3", 0.2))
          text(7, -12.5, "Signal Manipulation", pos = 4)
          
          rect(4.4, -9, 6.6, -4, col = alpha("chartreuse3", 0.3))
          text(5, -6.5, "Challenge", pos = 4)
        
        # Add legend  
          legend("bottomleft", c("No Color", "Dull Color", "Control", "Predator"),
                 pch = c(21, 21, 21, 24), pt.bg = c(col_sham, col_dull, rep("white", 2)), cex =0.9, bty = "n")

## Reproductive Success Models ----
          d_fem18$numdied <- d_fem18$clutch - d_fem18$numfled 
      #Models for 2018
          m_clutch_18 <- lm(clutch ~ color*challenge, data = d_fem18)
          m_brood_18 <- lm(maxbrood ~ color*challenge, data = d_fem18)
          m_numd6_18 <- lm(numd6 ~ color*challenge, data = d_fem18)
          m_numband_18 <- lm(numband ~ color*challenge, data = d_fem18)
          m_num_d15_18 <- lm(num_d15 ~ color*challenge, data = d_fem18)
          m_numfled_18 <- lm(numfled ~ color*challenge, data = d_fem18)
          
          # model to report
            m_fled18 <- glmer(cbind(numfled, numdied) ~ color*challenge + bbright_s + (1|soc_uby),
                              data = d_fem18, family = "binomial")
            m_fled18b <- glmer(cbind(numfled, numdied) ~ color*challenge + (1|soc_uby),
                               data = d_fem18, family = "binomial")
      
      # Overall effects reported from ANOVA. Group estimates from model summary tables.
          
        # Make table for 2018
          m_RS_18_t <- tab_model(m_clutch_18, m_brood_18, m_numd6_18, m_numband_18, m_num_d15_18, m_numfled_18,
                    pred.labels = c("Intercept", "Signal Dulled", "Predator", "Flight Reduction",
                                    "Dulled * Predator", "Dulled * Flight"),
                    dv.labels = c("Clutch Size", "Num. Hatched", "Num. Day 6", "Num. Day 12", "Num. Day 15", "Num. Fledged"))
          saveRDS(m_RS_18_t, here::here("5_other_outputs/m_RS_18_t.RDS"))
          
      #Models for 2019
          d_fem19$numdied <- d_fem19$clutch - d_fem19$numfled
          m_clutch_19 <- lm(clutch ~ color*challenge, data = d_fem19)
          m_brood_19 <- lm(maxbrood ~ color*challenge, data = d_fem19)
          m_numd6_19 <- lm(numd6 ~ color*challenge, data = d_fem19)
          m_numband_19 <- lm(numband ~ color*challenge, data = d_fem19)
          m_num_d15_19 <- lm(num_d15 ~ color*challenge, data = d_fem19)
          m_numfled_19 <- lm(numfled ~ color*challenge, data = d_fem19)
          
      # make a table
          fled19_tab <- d_fem %>%
            group_by(year, full_treatment) %>%
            summarise(n = n(), fled_mu = mean(numfled), fled_sd = sd(numfled))
          
          # model to report
            m_fled19 <- glmer(cbind(numdied, numfled) ~ color*challenge + bbright_s + (1|soc_uby),
                              data = d_fem19, family = "binomial")
          
       # Overall p-values from anova. Estimates from model summary table  
          
          # Make table for 2019
          m_RS_19_t <- tab_model(m_clutch_19, m_brood_19, m_numd6_19, m_numband_19, m_num_d15_19, m_numfled_19,
                    pred.labels = c("Intercept", "Predator", "Signal Dulled", "Predator * Dulled"),
                    dv.labels = c("Clutch Size", "Num. Hatched", "Num. Day 6", "Num. Day 12", "Num. Day 15", "Num. Fledged"))
          saveRDS(m_RS_19_t, here::here("5_other_outputs/m_RS_19_t.RDS"))
          
## Reproductive Success Plots ----
          
      # alternate plot for manuscript
          dfem2 <- d_fem[, c("soc_uby", "fband", "year", "color", "challenge", "full_treatment",
                             "numfled")]
          for(i in 1:nrow(dfem2)){
            if(dfem2$year[i] == 2018){
              dfem2$var1[i] <- dfem2$color[i]
              dfem2$var2[i] <- dfem2$challenge[i]
            }
            if(dfem2$year[i] == 2019){
              dfem2$var1[i] <- dfem2$challenge[i]
              dfem2$var2[i] <- dfem2$color[i]
            }
          }
          
          dfem2$var1 <- gsub("Dull", "Dulled", dfem2$var1)
          dfem2$var1 <- gsub("Control", "Sham", dfem2$var1)
          dfem2$var2 <- gsub("Stress", "Handicap", dfem2$var2)
          dfem2$fullt2 <- paste(dfem2$full_treatment, dfem2$year, sep = "_")
          
          dfem2$year <- gsub("2018", "Experiment One", dfem2$year)
          dfem2$year <- gsub("2019", "Experiment Two", dfem2$year)
          g <- ggplot(dfem2, mapping = aes(x = fullt2, y = numfled,
                                      color = var1, shape = var2, fill = var1)) +
            geom_boxplot(alpha = 0.25, outlier.shape = NA) +
            facet_grid(~year, scales = 'free') +
            scale_color_manual(values = c("orange", "slateblue")) +
            scale_fill_manual(values = c("orange", "slateblue")) +
            geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
            scale_shape_manual(values = c(21, 22, 24)) +
            theme_bw() +
            guides(fill = "none", color = "none", shape = "none") +
            theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                  axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
            xlab("") + ylab("Number fledged") +
            theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
            scale_x_discrete(labels = c('Control_Control_2018' = 'Sham-Control',
                                        'Control_Predator_2018' = "Sham-Predator",
                                        'Control_Stress_2018' = 'Sham-Handicap',
                                        'Dull_Control_2018' = 'Dulled-Control',
                                        'Dull_Predator_2018' = 'Dulled-Predator',
                                        'Dull_Stress_2018' = 'Dulled-Handicap',
                                        'Control_Control_2019' = 'Control-Sham',
                                        'Control_Dull_2019' = 'Control-Dulled',
                                        'Predator_Control_2019' = 'Predator-Sham',
                                        'Predator_Dull_2019' = 'Predator-Dulled'))
          
          gt <- ggplot_gtable(ggplot_build(g))
          gt$widths[5] <- 1.5*gt$widths[5]
          grid::grid.draw(gt)
          
          saveRDS(gt, here::here("5_other_outputs/rs_plot.rds"))
          
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
        
        lrsf <- subset(long_rs, long_rs$time_point == "numfled")
        

        
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
        m_feed18 <- lmer(f_feed ~ maxbrood + offset*color*challenge + color*challenge*bbright_s + (1|doy) + (1|uby), data = d_prov3_18)
        m_feed18b <- lmer(f_feed ~ maxbrood + offset*color*challenge + 
                          color*challenge + (1|doy) + (1|uby), data = d_prov3_18)
        m_feed18c <- lmer(f_feed ~ maxbrood + I(maxbrood^2) + offset + I(offset^2) + color*challenge + scale(m_feed) + (1|doy) + (1|uby), data = d_prov3_18)
        
      # get model predicted values for plotting  
        post18 <- mvrnorm(n = 1e5, mu = fixef(m_feed18b), Sigma = vcov(m_feed18b))
        mu_brd <- mean(d_prov3_18$maxbrood)
        mu_brd2 <- mu_brd*mu_brd
        r <- seq(1, 14, 1)
        mu_cd18 <- sapply(r, function(z)mean(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 4] + post18[, 7]*z))
        mu_cs18 <- sapply(r, function(z)mean(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z))
        mu_hd18 <- sapply(r, function(z)mean(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 4] + post18[, 6] +
                                               post18[, 7]*z + post18[, 9]*z + post18[, 11] + post18[, 13]*z))
        mu_hs18 <- sapply(r, function(z)mean(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 6] +
                                               post18[, 9]*z))
        mu_pd18 <- sapply(r, function(z)mean(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 4] + post18[, 5] +
                                               post18[, 7]*z + post18[, 8]*z + post18[, 10] + post18[, 12]*z))
        mu_ps18 <- sapply(r, function(z)mean(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 4] + post18[, 5] +
                                               post18[, 8]*z))
        
        ci_cd18 <- sapply(r, function(z)rethinking::HPDI(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 4] + post18[, 7]*z))
        ci_cs18 <- sapply(r, function(z)rethinking::HPDI(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z))
        ci_hd18 <- sapply(r, function(z)rethinking::HPDI(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 4] + post18[, 6] +
                                               post18[, 7]*z + post18[, 9]*z + post18[, 11] + post18[, 13]*z))
        ci_hs18 <- sapply(r, function(z)rethinking::HPDI(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 6] +
                                               post18[, 9]*z))
        ci_pd18 <- sapply(r, function(z)rethinking::HPDI(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 4] + post18[, 5] +
                                               post18[, 7]*z + post18[, 8]*z + post18[, 10] + post18[, 12]*z))
        ci_ps18 <- sapply(r, function(z)rethinking::HPDI(post18[, 1] + post18[, 2]*mu_brd + post18[, 3]*z +
                                               post18[, 4] + post18[, 5] +
                                               post18[, 8]*z))
        
        prov18 <- data.frame(feed = c(mu_cs18, mu_cd18, mu_hs18, mu_hd18, mu_ps18, mu_pd18),
                             lo = c(ci_cs18[1, ], ci_cd18[1, ], ci_hs18[1, ], ci_hd18[1, ], ci_ps18[1, ], ci_pd18[1, ]),
                             hi = c(ci_cs18[2, ], ci_cd18[2, ], ci_hs18[2, ], ci_hd18[2, ], ci_ps18[2, ], ci_pd18[2, ]),
                             age = rep(r, 6),
                             color = rep(c(rep("Sham", length(r)), rep("Dulled", length(r))), 3), 
                             challenge = c(rep("Control", length(r)*2), rep("Handicap", length(r)*2), rep("Predator", length(r)*2)))
        prov18$shape <- prov18$challenge
        

        
        
      # Make 2018 table
        m_feed18_t <- tab_model(m_feed18, m_feed18b, m_feed18c,
                  pred.labels = c("Intercept", "Brood Size", "Brood Size^2", "Nestling Age", "Nestling Age^2", "Signal Dulled", "Predator", "Flight Reduction",
                                  "Initial Brightness", "Dulled * Predator", "Dulled * Flight", "Dulled * Brightness", "Predator * Brightness", "Flight * Brightness",
                                  "Dulled * Predator * Brightness", "Dulled * Flight * Brightness", "Male Daily Provisioning"),
                  dv.labels = rep("Daily Female Provisioning Trips", 3))
        saveRDS(m_feed18_t, here::here("5_other_outputs/m_feed18_t.RDS"))
        
      # Fit models for 2019
        m_feed19 <- lmer(f_feed ~ maxbrood + offset*color*challenge + color*challenge*bbright_s + (1|doy) + (1|uby), data = d_prov3_19)
        m_feed19b <- lmer(f_feed ~ maxbrood + offset*color*challenge + color*challenge + (1|doy) + (1|uby), data = d_prov3_19)
        m_feed19b2 <- lmer(f_feed ~ maxbrood + offset*color*challenge + color*challenge +
                             (1|doy) + (1|uby), data = d_prov3_19)
        m_feed19c <- lmer(f_feed ~ maxbrood + I(maxbrood^2) + offset + I(offset^2) + color*challenge + scale(m_feed) + (1|doy) + (1|uby), data = d_prov3_19)
        
      # Make 2019 table
        m_feed19_t <- tab_model(m_feed19, m_feed19b, m_feed19c,
                  pred.labels = c("Intercept", "Brood Size", "Brood Size^2", "Nestling Age", "Nestling Age^2", "Predator", "Signal Dulled", "Initial Brightness",
                                  "Predator * Dulled", "Predator * Brightness", "Dulled * Brightness", "Predator * Dulled * Brightness", "Male Daily Provisioning"),
                  dv.labels = rep("Daily Female Provisioning Trips", 3))
        saveRDS(m_feed19_t, here::here("5_other_outputs/m_feed19_t.RDS"))
        
      # get posterior estimates from 2019 table
        post19 <- mvrnorm(n = 1e5, mu = fixef(m_feed19b), Sigma = vcov(m_feed19b))
        mu_brood19 <- mean(d_prov3_19$maxbrood)
        r2 <- seq(1, 14, 1)
        mu_sc19 <- sapply(r2, function(z) mean(post19[, 1] + post19[, 2]*mu_brood19 + post19[, 3]*z))
        mu_dc19 <- sapply(r2, function(z) mean(post19[, 1] + post19[, 2]*mu_brood19 + post19[, 3]*z +
                                                 post19[, 5] + post19[, 7]*z))
        mu_sp19 <- sapply(r2, function(z) mean(post19[, 1] + post19[, 2]*mu_brood19 + post19[, 3]*z +
                                                 post19[, 4] + post19[, 6]*z))
        mu_dp19 <- sapply(r2, function(z) mean(post19[, 1] + post19[, 2]*mu_brood19 + post19[, 3]*z +
                                                 post19[, 4] + post19[, 5] + post19[, 6]*z + post19[, 7]*z + post19[, 8] +
                                                 post19[, 9]*z))
        
        ci_sc19 <- sapply(r2, function(z) rethinking::HPDI(post19[, 1] + post19[, 2]*mu_brood19 + post19[, 3]*z))
        ci_dc19 <- sapply(r2, function(z) rethinking::HPDI(post19[, 1] + post19[, 2]*mu_brood19 + post19[, 3]*z +
                                                 post19[, 5] + post19[, 7]*z))
        ci_sp19 <- sapply(r2, function(z) rethinking::HPDI(post19[, 1] + post19[, 2]*mu_brood19 + post19[, 3]*z +
                                                 post19[, 4] + post19[, 6]*z))
        ci_dp19 <- sapply(r2, function(z) rethinking::HPDI(post19[, 1] + post19[, 2]*mu_brood19 + post19[, 3]*z +
                                                 post19[, 4] + post19[, 5] + post19[, 6]*z + post19[, 7]*z + post19[, 8] +
                                                 post19[, 9]*z))
        
        prov19 <- data.frame(feed = c(mu_sc19, mu_dc19, mu_sp19, mu_dp19),
                             lo = c(ci_sc19[1, ], ci_dc19[1, ], ci_sp19[1, ], ci_dp19[1, ]),
                             hi = c(ci_sc19[2, ], ci_dc19[2, ], ci_sp19[2, ], ci_dp19[2, ]),
                             age = rep(r2, 4),
                             challenge = rep(c(rep("Sham", length(r)), rep("Dulled", length(r))), 2),
                             color = c(rep("Control", length(r)*2), rep("Predator", length(r)*2)))
        prov19$shape <- prov19$color
    
## Provisioning Plots ----
  # Figure showing nestling provisioning
      # These plots are not in the style of the other plots above but this was faster to make
        # and I'm not sure that any of them will actually end up staying in the paper anyway.
        
        dp3 <- d_prov3 %>%
          dplyr::group_by(full_treatment, year, offset, color, challenge) %>%
          dplyr::summarise(n = n(), f_feed2 = mean(f_feed, na.rm = TRUE), sd = sd(f_feed, na.rm = TRUE))
        dp3$se <- dp3$sd / sqrt(dp3$n)
        dp3_18 <- subset(dp3, dp3$year == 2018)
        dp3_18$f_feed <- dp3_18$f_feed2
        
        dp3_18$challenge <- gsub("Stress", "Handicap", dp3_18$challenge)
        dp3_18$color <- gsub("Control", "Sham", dp3_18$color)
        dp3_18$color <- gsub("Dull", "Dulled", dp3_18$color)
        
        d_prov3_18$challenge <- gsub("Stress", "Handicap", d_prov3_18$challenge)
        d_prov3_18$color <- gsub("Control", "Sham", d_prov3_18$color)
        d_prov3_18$color <- gsub("Dull", "Dulled", d_prov3_18$color)
        d_prov3_18$shape <- d_prov3_18$challenge
    
      # Plot female provisioning from 2018 experiment
       p1 <- ggplot(d_prov3_18, aes(x = offset, y = f_feed, shape = challenge,
                                    fill = color)) +
          geom_jitter(alpha = 0.3, width = 0.1, mapping = aes(color = color)) + 
          #geom_smooth(se = FALSE) + 
          xlim(0, 15) + xlab("Days after hatching") +
          ylab("Daily female provisioning trips") + ggtitle("Experiment One") +
          theme_bw() + theme(legend.position = c(0.15, 0.8)) +
          scale_color_manual(values = c("slateblue", "orange")) +
          scale_fill_manual(values = c("slateblue", "orange")) +
          geom_line(data = dp3_18, mapping = aes(color = color)) +
          geom_point(data = dp3_18, size = 2) +
          guides(fill = "none") +
          scale_shape_manual(values = c(21, 22, 24)) +
         theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
               axis.text = element_text(size = 12), axis.title = element_text(size = 14),
               legend.title = element_blank(), legend.background = element_rect(fill = alpha("white", 0)))
       
      # plot the data from the model rather than just raw
       p18prov <- ggplot(data = prov18, mapping = aes(x = age, y = feed, color = color, linetype = challenge)) +
         geom_jitter(data = d_prov3_18, mapping = aes(x = offset, y = f_feed, shape = shape, fill = color), alpha = 0.2) +
         #geom_ribbon(mapping = aes(ymin = lo, ymax = hi, fill = color), alpha = 0.2, color = NA) +
         geom_line(size = 1.1) +
         xlim(0, 15) + xlab("Days after hatching") +
         ylab("Daily female provisioning trips") + ggtitle("Experiment One") +
         theme_bw() + #theme(legend.position = c(0.15, 0.8)) +
         annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "A", size = 6) +
         scale_color_manual(values = c("slateblue", "orange")) +
         scale_fill_manual(values = c("slateblue", "orange")) +
         #guides(fill = "none", color = "none", linetype = "none", shape = "none") +
         theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
               axis.text = element_text(size = 12), axis.title = element_text(size = 14),
               legend.title = element_blank(), legend.background = element_rect(fill = alpha("white", 0))) +
         scale_shape_manual(values = c(21, 22, 24)) +
         guides(shape = guide_legend(override.aes = list(fill = "black", alpha = 1)),
                color = guide_legend(override.aes = list(alpha = 1)),
                linetype = guide_legend(override.aes = list(size = 0.9))) +
         scale_linetype_manual(values = c("solid", "dotted", "dashed"))
        
      # Plot male provisioning from 2018 experiment  
        # ggplot(d_prov3_18, aes(x = offset, y = m_feed, col = full_treatment)) +
        #   geom_point() + geom_smooth() + xlim(0, 15) + xlab("Days after hatching") +
        #   ylab("Daily male provisioning trips") + ggtitle("Dulling then Challenge") +
        #   theme_bw()
       
       # same for 2019
       # note that 'challenge' and 'color' columns are reversed in this dataset to allow merging with 2018
       
       dp3_19 <- subset(dp3, dp3$year == 2019)
       dp3_19$f_feed <- dp3_19$f_feed2
       
       dp3_19$challenge <- gsub("Control", "Sham", dp3_19$challenge)
       dp3_19$challenge <- gsub("Dull", "Dulled", dp3_19$challenge)
       dp3_19 <- subset(dp3_19, dp3_19$offset < 15)
       
       d_prov3_19$challenge <- gsub("Control", "Sham", d_prov3_19$challenge)
       d_prov3_19$challenge <- gsub("Dull", "Dulled", d_prov3_19$challenge)
       d_prov3_19 <- subset(d_prov3_19, d_prov3_19$offset < 15) 
       
      # Plot female provisioning from 2019 experiment  
       p2 <- ggplot(d_prov3_19, aes(x = offset, y = f_feed, shape = color,
                                    fill = challenge)) +
         geom_jitter(alpha = 0.3, width = 0.1, mapping = aes(color = challenge)) + 
         #geom_smooth(se = FALSE) + 
         xlim(0, 15) + xlab("Days after hatching") +
         ylab("Daily female provisioning trips") + ggtitle("Experiment Two") +
         theme_bw() + theme(legend.position = c(0.16, 0.75)) +
         scale_color_manual(values = c("slateblue", "orange")) +
         scale_fill_manual(values = c("slateblue", "orange")) +
         geom_line(data = dp3_19, mapping = aes(color = challenge)) +
         geom_point(data = dp3_19, size = 2) +
         guides(fill = "none", shape = "none", color = "none") +
         scale_shape_manual(values = c(21, 24)) +
         theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
               axis.text = element_text(size = 12), axis.title = element_text(size = 14),
               legend.title = element_blank(), legend.background = element_rect(fill = alpha("white", 0)))
       
      # plot the data from the model rather than just raw for 2019
       p19prov <- ggplot(data = prov19, mapping = aes(x = age, y = feed, color = challenge, linetype = color)) +
         geom_jitter(data = d_prov3_19, mapping = aes(x = offset, y = f_feed, shape = color, fill = challenge), alpha = 0.2) +
         #geom_ribbon(mapping = aes(ymin = lo, ymax = hi, fill = challenge), alpha = 0.2, color = NA) +
         geom_line(size = 1.1) +
         xlim(0, 15) + xlab("Days after hatching") +
         ylab("Daily female provisioning trips") + ggtitle("Experiment Two") +
         theme_bw() + #theme(legend.position = c(0.15, 0.8)) +
         annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "B", size = 6) +
         scale_color_manual(values = c("slateblue", "orange")) +
         scale_fill_manual(values = c("slateblue", "orange")) +
         guides(color = "none", fill = "none", shape = "none", linetype = "none") +
         theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
               axis.text = element_text(size = 12), axis.title = element_text(size = 14),
               legend.title = element_blank(), legend.background = element_rect(fill = alpha("white", 0))) +
         scale_shape_manual(values = c(21, 22)) +
         scale_linetype_manual(values = c("solid", "dashed"))
      
      # Plot male provisioning from 2019 experiment  
        # ggplot(d_prov3_19, aes(x = offset, y = m_feed, col = full_treatment)) +
        #   geom_point() + geom_smooth() + xlim(0.6, 14.5) + xlab("Days after hatching") +
        #   ylab("Daily male provisioning trips") + ggtitle("Challenge then Dulling") +
        #   theme_bw()
       
       
       prov_plot <- ggpubr::ggarrange(p18prov, p19prov, widths = c(1, 0.75))
       #prov_plot <- ggpubr::ggarrange(p1, p2)
        saveRDS(prov_plot, here::here("5_other_outputs/prov_plot.rds"))
    
## Social Interaction Models ----
    
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
        
    # ggplot(data = soc_models, mapping = aes(x = offset, y = uni_f_vis, col = full_treatment)) + 
    #   geom_jitter(width = 0.1, col = "slateblue", size = 0.7, alpha = 0.6, height = 0) +
    #   geom_smooth(method = "loess") + facet_wrap(~ year)
        
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
        soc_mod18s1 <- subset(soc_mod18, soc_mod18$stage == "stage1" & soc_mod18$offset < 2)
        soc_mod18s2 <- subset(soc_mod18, soc_mod18$stage == "stage2" & soc_mod18$offset > 1)
        soc_mod19s1 <- subset(soc_mod19, soc_mod19$stage == "stage1" & soc_mod19$offset < 0)
        soc_mod19s2 <- subset(soc_mod19, soc_mod19$stage == "stage2" & soc_mod19$offset > -1)
        
    ## Fit models for each response variable 2018   
        #STAGE 1
          m_fvis_18_s1 <- glmer(tot_f_vis ~ color*bbright_s + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod18s1)
          m_fvis_18_s1b <- glmer(tot_f_vis ~ color + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod18s1)
          m_trip_18_s1 <- glmer(tot_trips ~ color*bbright_s + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod18s1)
          m_trip_18_s1b <- glmer(tot_trips ~ color + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod18s1)
        #STAGE 2
          m_fvis_18_s2 <- glmer(tot_f_vis ~ color*challenge*bbright_s + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod18s2,
                                control = glmerControl(optimizer = "bobyqa"))
          m_trip_18_s2 <- glmer(tot_trips ~ color*challenge*bbright_s + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod18s2,
                                control = glmerControl(optimizer = "bobyqa"))
          m_trip_18_s2b <- glmer(tot_trips ~ color*challenge + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod18s2,
                                 control = glmerControl(optimizer = "bobyqa"))
          
        #STAGE 1 Table
          # m_soc18_st1_t <- tab_model(m_mvis_18_s1b, m_fvis_18_s1b, m_trip_18_s1b, 
          #           pred.labels = c("Intercept", "Signal Dulled"),
          #           dv.labels = c("Daily Unique Male Visitors", "Daily Unique Female Visitors", "Daily Unique Trips to Other Boxes"))
          # saveRDS(m_soc18_st1_t, here::here("5_other_outputs/m_soc_18_st1_t.RDS"))
          
        #STAGE 2 Table
          # m_soc18_st2_t <- tab_model(m_mvis_18_s2b, m_fvis_18_s2, m_trip_18_s2b,
          #           pred.labels = c("Intercept", "Signal Dulled", "Predator", "Flight Reduction",
          #                           "Dulled * Predator", "Dulled * Flight", "Initial Brightness", "Dulled * Initial Bright", 
          #                           "Predator * Initial Bright", "Flight * Initial Bright", "Dulled * Predator * Bright",
          #                           "Dulled * Flight * Bright"),
          #           dv.labels = c("Daily Unique Male Visitors", "Daily Unique Female Visitors", "Daily Unique Trips to Other Boxes"))
          # saveRDS(m_soc18_st2_t, here::here("5_other_outputs/m_soc_18_st2_t.RDS"))
        
    ## Fit models for each response variable 2018   
        #STAGE 1
          m_fvis_19_s1 <- glmer(tot_f_vis ~ color*bbright_s + scale(offset) + (1|uby), family = "poisson", data = soc_mod19s1)
          m_fvis_19_s1b <- glmer(tot_f_vis ~ color + scale(offset) + (1|uby), family = "poisson", data = soc_mod19s1)
          m_trip_19_s1 <- glmer(tot_trips ~ color*bbright_s + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod19s1)
          m_trip_19_s1b <- glmer(tot_trips ~ color + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod19s1)
        #STAGE 2
          m_fvis_19_s2 <- glmer(tot_f_vis ~ color*challenge*bbright_s + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod19s2,
                                control = glmerControl(optimizer = "bobyqa"))
          m_fvis_19_s2b <- glmer(tot_f_vis ~ color*challenge + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod19s2,
                                 control = glmerControl(optimizer = "bobyqa"))
          m_trip_19_s2 <- glmer(tot_trips ~ color*challenge*bbright_s + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod19s2,
                                control = glmerControl(optimizer = "bobyqa"))
          m_trip_19_s2b <- glmer(tot_trips ~ color*challenge + scale(offset) + (1|uby) + (1|doy), family = "poisson", data = soc_mod19s2,
                                 control = glmerControl(optimizer = "bobyqa"))
          
        #STAGE 1 Table
          # m_soc19_st1_t <- tab_model(m_mvis_19_s1, m_fvis_19_s1b, m_trip_19_s1b, 
          #           pred.labels = c("Intercept", "Predator", "Initial Brightness", "Predator * Brightness"),
          #           dv.labels = c("Daily Unique Male Visitors", "Daily Unique Female Visitors", "Daily Unique Trips to Other Boxes"))
          # saveRDS(m_soc19_st1_t, here::here("5_other_outputs/m_soc_19_st1_t.RDS"))
          
        #STAGE 2 Table
          # m_soc19_st2_t <- tab_model(m_mvis_19_s2b, m_fvis_19_s2b, m_trip_19_s2b,
          #           pred.labels = c("Intercept", "Predator", "Signal Dulled", "Predator * Dulled"),
          #           dv.labels = c("Daily Unique Male Visitors", "Daily Unique Female Visitors", "Daily Unique Trips to Other Boxes"))
          # saveRDS(m_soc19_st2_t, here::here("5_other_outputs/m_soc_19_st2_t.RDS"))
        
## Social Interaction plots ----        
    # Making a plot for visits for experiment 1
        post1 <- as.data.frame(mvrnorm(n = 1e5, mu = fixef(m_fvis_18_s2), Sigma = vcov(m_fvis_18_s2))) 
        colnames(post1)[1] <- "Intercept"
          
        r <- seq(-2.5, 2.5, 0.1)
       
      # get maximum likelihood estimate for each group   
        mu_cc <- sapply(r, function(z)mean(post1$Intercept + post1$bbright_s*z))
        mu_cp <- sapply(r, function(z)mean(post1$Intercept + post1$bbright_s*z + post1$challengePredator +
                                             post1$`challengePredator:bbright_s`*z))
        mu_ch <- sapply(r, function(z)mean(post1$Intercept + post1$bbright_s*z + post1$challengeStress +
                                             post1$`challengeStress:bbright_s`*z))
        mu_dc <- sapply(r, function(z)mean(post1$Intercept + post1$bbright_s*z + post1$colorDull +
                                             post1$`colorDull:bbright_s`*z))
        mu_dp <- sapply(r, function(z)mean(post1$Intercept + post1$bbright_s*z + post1$colorDull +
                                             post1$`colorDull:bbright_s`*z +
                                             post1$`colorDull:challengePredator` +
                                             post1$`colorDull:challengePredator:bbright_s`*z))
        mu_dh <- sapply(r, function(z)mean(post1$Intercept + post1$bbright_s*z + post1$colorDull +
                                             post1$`colorDull:bbright_s`*z +
                                             post1$`colorDull:challengeStress` +
                                             post1$`colorDull:challengeStress:bbright_s`*z))
      
      # get confidence interval for each group across range of brightness 
        ci_cc <- sapply(r, function(z)rethinking::HPDI(post1$Intercept + post1$bbright_s*z))
        ci_cp <- sapply(r, function(z)rethinking::HPDI(post1$Intercept + post1$bbright_s*z + post1$challengePredator +
                                             post1$`challengePredator:bbright_s`*z))
        ci_ch <- sapply(r, function(z)rethinking::HPDI(post1$Intercept + post1$bbright_s*z + post1$challengeStress +
                                             post1$`challengeStress:bbright_s`*z))
        ci_dc <- sapply(r, function(z)rethinking::HPDI(post1$Intercept + post1$bbright_s*z + post1$colorDull +
                                             post1$`colorDull:bbright_s`*z))
        ci_dp <- sapply(r, function(z)rethinking::HPDI(post1$Intercept + post1$bbright_s*z + post1$colorDull +
                                             post1$`colorDull:bbright_s`*z +
                                             post1$`colorDull:challengePredator` +
                                             post1$`colorDull:challengePredator:bbright_s`*z))
        ci_dh <- sapply(r, function(z)rethinking::HPDI(post1$Intercept + post1$bbright_s*z + post1$colorDull +
                                             post1$`colorDull:bbright_s`*z +
                                             post1$`colorDull:challengeStress` +
                                             post1$`colorDull:challengeStress:bbright_s`*z))
        
      # get 95% ci for slope for each group to tell which ones have 'significant' slopes
        z <- c(0, 1)
        sl_cc <- HPDI((post1$Intercept + post1$bbright_s*z[2]) - (post1$Intercept + post1$bbright_s*z[1]))
        sl_cp <- HPDI((post1$Intercept + post1$bbright_s*z[2] + post1$challengePredator +
                                                         post1$`challengePredator:bbright_s`*z[2]) -
                        (post1$Intercept + post1$bbright_s*z[1] + post1$challengePredator +
                           post1$`challengePredator:bbright_s`*z[1]))
        sl_ch <- HPDI((post1$Intercept + post1$bbright_s*z[2] + post1$challengeStress +
                                                         post1$`challengeStress:bbright_s`*z[2]) -
                        (post1$Intercept + post1$bbright_s*z[1] + post1$challengeStress +
                           post1$`challengeStress:bbright_s`*z[1]))
        sl_dc <- HPDI((post1$Intercept + post1$bbright_s*z[2] + post1$colorDull +
                                                         post1$`colorDull:bbright_s`*z[2]) -
                        (post1$Intercept + post1$bbright_s*z[1] + post1$colorDull +
                           post1$`colorDull:bbright_s`*z[1]))
        sl_dp <- HPDI((post1$Intercept + post1$bbright_s*z[2] + post1$colorDull +
                                                         post1$`colorDull:bbright_s`*z[2] +
                                                         post1$`colorDull:challengePredator` +
                                                         post1$`colorDull:challengePredator:bbright_s`*z[2]) -
                        (post1$Intercept + post1$bbright_s*z[1] + post1$colorDull +
                           post1$`colorDull:bbright_s`*z[1] +
                           post1$`colorDull:challengePredator` +
                           post1$`colorDull:challengePredator:bbright_s`*z[1]))
        sl_dh <- HPDI((post1$Intercept + post1$bbright_s*z[2] + post1$colorDull +
                                                         post1$`colorDull:bbright_s`*z[2] +
                                                         post1$`colorDull:challengeStress` +
                                                         post1$`colorDull:challengeStress:bbright_s`*z[2]) -
                        (post1$Intercept + post1$bbright_s*z[1] + post1$colorDull +
                           post1$`colorDull:bbright_s`*z[1] +
                           post1$`colorDull:challengeStress` +
                           post1$`colorDull:challengeStress:bbright_s`*z[1]))
        
        
        p1d <- data.frame(color = c(rep("Control", length(r)*3), rep("Dulled", length(r)*3)),
                          challenge = rep(c(rep("Control", length(r)), rep("Predator", length(r)), rep("Handicap", length(r))), 2),
                          lo = c(ci_cc[1,], ci_cp[1,], ci_ch[1,], ci_dc[1,], ci_dp[1,], ci_dh[1,]),
                          hi = c(ci_cc[2,], ci_cp[2,], ci_ch[2,], ci_dc[2,], ci_dp[2,], ci_dh[2,]),
                          mu = c(mu_cc, mu_cp, mu_ch, mu_dc, mu_dp, mu_dh),
                          r = rep(r, 6))
        p1d$full <- paste(p1d$color, p1d$challenge, sep = "_")
        
        p1d_plot <- ggplot(p1d, mapping = aes(x = r, y = mu, color = color, linetype = full, fill = color)) +
          geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
          geom_line() +
          facet_grid(challenge ~ color) +
          guides(color = "none", fill = "none", linetype = "none") +
          scale_color_manual(values = c("orange", "slateblue")) +
          scale_fill_manual(values = c("orange", "slateblue")) +
          xlab("Initial Brightness (SD)") +
          ylab("Visitors to Nest \n (Model Predicted)") +
          theme_bw() +
          theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                  axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
          ggtitle("Experiment One") +
          scale_linetype_manual(values = c(1, 2, 2, 2, 2, 1))
          
        saveRDS(p1d_plot, here::here("5_other_outputs/social1_plot.rds"))
                                            
        
          
    # Making a plot for trips for experiment 2
        post2 <- as.data.frame(mvrnorm(n = 1e5, mu = fixef(m_trip_19_s2), Sigma = vcov(m_trip_19_s2))) 
        colnames(post2)[1] <- "Intercept"
        
        r <- seq(-2.5, 2.5, 0.1)
        
      # calculate maximum likelihood estimates  
        mu_cc2 <- sapply(r, function(z)mean(post2$Intercept + post2$bbright_s*z))
        mu_cp2 <- sapply(r, function(z)mean(post2$Intercept + post2$bbright_s*z + post2$colorPredator +
                                             post2$`colorPredator:bbright_s`*z))
        mu_dc2 <- sapply(r, function(z)mean(post2$Intercept + post2$bbright_s*z + post2$challengeDull +
                                             post2$`challengeDull:bbright_s`*z))
        mu_dp2 <- sapply(r, function(z)mean(post2$Intercept + post2$bbright_s*z + post2$challengeDull +
                                             post2$`challengeDull:bbright_s`*z +
                                             post2$`colorPredator:challengeDull` +
                                             post2$`colorPredator:challengeDull:bbright_s`*z))
       
      # calculate confidence intervals   
        ci_cc2 <- sapply(r, function(z)rethinking::HPDI(post2$Intercept + post2$bbright_s*z))
        ci_cp2 <- sapply(r, function(z)rethinking::HPDI(post2$Intercept + post2$bbright_s*z + post2$colorPredator +
                                             post2$`colorPredator:bbright_s`*z))
        ci_dc2 <- sapply(r, function(z)rethinking::HPDI(post2$Intercept + post2$bbright_s*z + post2$challengeDull +
                                             post2$`challengeDull:bbright_s`*z))
        ci_dp2 <- sapply(r, function(z)rethinking::HPDI(post2$Intercept + post2$bbright_s*z + post2$challengeDull +
                                             post2$`challengeDull:bbright_s`*z +
                                             post2$`colorPredator:challengeDull` +
                                             post2$`colorPredator:challengeDull:bbright_s`*z))
        
      # calculate cis of slopes
        sl_cc2 <- HPDI((post2$Intercept + post2$bbright_s*z[2]) -
                         (post2$Intercept + post2$bbright_s*z[1]))
        sl_cp2 <- HPDI((post2$Intercept + post2$bbright_s*z[2] + post2$colorPredator +
                                                          post2$`colorPredator:bbright_s`*z[2]) -
                         (post2$Intercept + post2$bbright_s*z[1] + post2$colorPredator +
                            post2$`colorPredator:bbright_s`*z[1]))
        sl_dc2 <- HPDI((post2$Intercept + post2$bbright_s*z[2] + post2$challengeDull +
                                                          post2$`challengeDull:bbright_s`*z[2]) -
                         (post2$Intercept + post2$bbright_s*z[1] + post2$challengeDull +
                            post2$`challengeDull:bbright_s`*z[1]))
        sl_dp2 <- HPDI((post2$Intercept + post2$bbright_s*z[2] + post2$challengeDull +
                                                          post2$`challengeDull:bbright_s`*z[2] +
                                                          post2$`colorPredator:challengeDull` +
                                                          post2$`colorPredator:challengeDull:bbright_s`*z[2]) -
                         (post2$Intercept + post2$bbright_s*z[1] + post2$challengeDull +
                            post2$`challengeDull:bbright_s`*z[1] +
                            post2$`colorPredator:challengeDull` +
                            post2$`colorPredator:challengeDull:bbright_s`*z[1]))

        p2d <- data.frame(color = c(rep("Control", length(r)*2), rep("Dulled", length(r)*2)),
                          challenge = rep(c(rep("Control", length(r)), rep("Predator", length(r))), 2),
                          lo = c(ci_cc2[1,], ci_cp2[1,], ci_dc2[1,], ci_dp2[1,]),
                          hi = c(ci_cc2[2,], ci_cp2[2,], ci_dc2[2,], ci_dp2[2,]),
                          mu = c(mu_cc2, mu_cp2, mu_dc2, mu_dp2),
                          r = rep(r, 4))
        p2d$full <- paste(p2d$color, p2d$challenge, sep = "_")
        
        p2d_plot <- ggplot(p2d, mapping = aes(x = r, y = mu, color = color, fill = color, linetype = full)) +
          geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
          geom_line() +
          facet_grid(challenge ~ color) +
          guides(color = "none", fill = "none", linetype = "none") +
          scale_color_manual(values = c("orange", "slateblue")) +
          scale_fill_manual(values = c("orange", "slateblue")) +
          xlab("Initial Brightness (SD)") +
          ylab("Trips to Other Nests \n (Model Predicted)") +
          theme_bw() +
          theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
          ggtitle("Experiment Two") +
          scale_linetype_manual(values = c(2, 2, 2, 1))
        
        saveRDS(p2d_plot, here::here("5_other_outputs/social2_plot.rds"))
        
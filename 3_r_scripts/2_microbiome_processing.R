## This script processes the raw microbiome sequence data into an object for analysis
  # It only needs to be run once to make the phyloseq object and then the object is saved
  # and can be re-loaded. Takes a very long time to run.


## Load packages ----
    # I think some of these are no longer actually required
    pacman::p_load(BiocStyle, ggplot2, gridExtra, dada2, phyloseq, phangorn, magrittr, sjPlot,
               plyr, dplyr, picante, decontam, GUniFrac, vegan, here, reshape2, DECIPHER)

## Settings to apply throughout ----
    ## Set random seed
      set.seed(47)

## Locate files ----
  ## Specify where fastq files are stored
    # I'm not syncing these to github
        miseq_path <- here::here("/1_raw_data/16s_sequences/") ## Set path to downloaded sequences
        head(list.files(miseq_path))

## Begin to process the read files ----
    # Sort ensures forward/reverse reads are in same order
        fnFs <- sort(list.files(miseq_path, pattern = "_R1.fastq"))
        fnRs <- sort(list.files(miseq_path, pattern = "_R2.fastq"))
    # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
        sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 7)
    # Specify the full path to the fnFs and fnRs
        fnFs <- file.path(miseq_path, fnFs)
        fnRs <- file.path(miseq_path, fnRs)

## Trim and filter, plot the read qualities ----
        plotQualityProfile(fnFs[c(1:3)])
        plotQualityProfile(fnRs[c(1:3)])

## Filter based on plots above to 180 bp long
      ## Define filenames for the new filtered files
          filt_path <- file.path(miseq_path, "filtered")
          if(!file_test("-d", filt_path)) dir.create(filt_path)
          filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
          filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))  

          
## Further filtering ----        
    ## Filter out reads based on the numbers that we decided on and some other
      ## expected error rates. Both forward and reverse reads from a pair have
      ## to pass in order to be included.
        out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(180, 181), trimLeft = c(19, 20),
                       maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
        # change trunclen back to 180,181, change trimleft to 19,20, make longer and 20, 20 for coi

## Main DADA2 functions ----
    ## Running DADA2
      ## DEREPLICATION (this is the section using dada2 to assign ASVs)
          derepFs <- derepFastq(filtFs, verbose = TRUE)
          derepRs <- derepFastq(filtRs, verbose = TRUE)
      ## Name the derep-class objects by the sample names
          names(derepFs) <- sampleNames
          names(derepRs) <- sampleNames 
      ## Unsupervised learning to build the models that characterize error rates
          errF <- learnErrors(filtFs, multithread = TRUE)
          errR <- learnErrors(filtRs, multithread = TRUE)  
      ## Look at the results of the error learning runs
          plotErrors(errF)
          plotErrors(errR)  

## Calling ASVs ----
    ## Inference of ASVs using the error models built above          
        dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
        dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
    
    ## Inspect the outcome
        dadaFs[[1]]
        dadaRs[[1]] 

## Merge forward and reverse reads ----
    ## Construct the sequence table with the cleaned dada data
        mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
        seqtabAll <- makeSequenceTable(mergers)
        seqtab <- seqtabAll[, nchar(colnames(seqtabAll)) %in% seq(247, 270)]  
        table(nchar(getSequences(seqtab)))

## Remove chimeras ----
    seqtab.nochim <- removeBimeraDenovo(seqtab)

## Final filtering table
    ## Produce table showing results of all filtering steps so far
    ## Check on table for all the filtering steps taken so far      
        getN <- function(x) sum(getUniques(x))
        track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
        colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
        rownames(track) <- sampleNames
        head(track)  

## Assign taxonomy ----
    ## Assign taxonomy. The silva file is downloaded and placed in the folder ahead of time.
      fastaRef <- here("/9_supporting_files/silva_nr_v132_train_set.fa.gz")
      taxTab <- assignTaxonomy(seqtab.nochim, refFasta = fastaRef, multithread = TRUE)          
      unname(head(taxTab))

## Build phylogeny ----
    ## Build phylogenetic tree based on samples identified
      seq <- getSequences(seqtab.nochim)
      names(seq) <- seq # This propagates to the tip labels of the tree
      alignment <- AlignSeqs(DNAStringSet(seq),anchor=NA,verbose=FALSE)

      phangAlign <- phyDat(as(alignment, "matrix"), type = "DNA")
      dm <- dist.ml(phangAlign)
      treeNJ <- NJ(dm) #note, tip order != sequence order
      fit = pml(treeNJ, data = phangAlign)
      fitGTR <- update(fit, k = 4, inv = 0.2)
      fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, ptGamma = TRUE,
                        rearrangement = "stochastic", control = pml.control(trace = 0))
      
    ## save the objects
      saveRDS(fitGTR, here("2_modified_data/fit_GTR.rds"))
      saveRDS(taxTab, here("2_modified_data/tax_tab.rds"))
      saveRDS(seqtab.nochim, here("2_modified_data/seqtab.rds"))


## Build phyloseq object ----
  ## Combine otu table, sample data, taxonomy table, and tree into phyloseq object
    samdf <- read.delim(here("1_raw_data/data_by_micro_sample.txt"))
    samdf$sampleID <- as.character(samdf$micro_number)
    rownames(samdf) <- samdf$sampleID
    ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                 sample_data(samdf),
                 tax_table(taxTab),
                 phy_tree(fitGTR$tree))
    ps

## Save original to disk
  ## Save the original unmodified phyloseq objec to disk
    #saveRDS(ps, here("2_modified_data/tres18_19_original_ps.rds"))
    ps <- readRDS(here("2_modified_data/tres18_19_original_ps.rds"))

## Remove non-bacteria ----
    ## Remove eukaryotes, archaea, chloroplasts, and mitochondria
    ps2<-subset_taxa(ps,
                     Kingdom == "Bacteria" &
                       Family != "mitochondria" &
                       Class != "Chloroplast"
    )
    ps2
    
## Remove contaminants with decontam ----
    ## Identify contaminants based on combination of prevalence and frequency then remove them
        sample_data(ps2)$is.neg <- sample_data(ps2)$sample_or_control == "neg_control"
        contamdf.prev <- isContaminant(ps2, conc = "quant_reading", method = "combined", neg = "is.neg")
        table(contamdf.prev$contaminant)
        ps.nc <- prune_taxa(!contamdf.prev$contaminant, ps2)
 
## Remove singletons ----
    ps.no.sing <- prune_taxa(taxa_sums(ps.nc) > 1,ps.nc)

## Save this phyloseq object to disk ----
    saveRDS(ps.no.sing, here("2_modified_data/tres_ps_no_sings.rds"))
    tres_ps <- readRDS(here("2_modified_data/tres_ps_no_sings.rds"))
    
## Plot reads ----
    ## Plot distribution of sample read depths
      sample_sum_df <- data.frame(sum = sample_sums(tres_ps))
    
      ggplot(sample_sum_df, aes(x = sum)) +
        geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
        ggtitle("Distribution of sample sequencing depth") +
        xlab("Read counts") +
        theme(axis.title.y = element_blank())

## Rarify for alpha diversity ----
  ## Rarify and then estimate diversity in phyloseq and picante
      pr <- prune_samples(sample_sums(tres_ps) >= 1000, tres_ps)
      pr1 <- rarefy_even_depth(pr, sample.size = 1000, rngseed = 121, trimOTUs = TRUE)
      
    # Agglomerate to taxonomic level before calculating alpha diversity
        prx <- tax_glom(pr1, "Genus")
        prx <- pr1  # If using pr1 then no agglomeration
      
    # richness and diversity metrics for each sample
        pr_rich <- estimate_richness(prx)
        OTU <- phyloseq::otu_table(prx)
        if(taxa_are_rows(OTU)){
          OTU <- t(OTU)
        }
        otutable <- as(OTU, "matrix")
        tree <- phyloseq::phy_tree(prx)
        pdt <- picante::pd(otutable, tree, include.root = FALSE)
        pr_rich$sampleID <- rownames(pr_rich)
        pdt$sampleID <- rownames(pdt)
        rich <- join(pdt, pr_rich, "sampleID")  
        out.rich <- join(samdf, rich, "sampleID")
        
    # Cast into a wide version that has each measure for captures 1, 2, 3 in the same row for models
        out.rich2 <- subset(out.rich, out.rich$sample_or_control == "sample")
        rich_wide <- tidyr::pivot_wider(out.rich2, id_cols = c(uby, band, unit, nest, year, color, challenge, full_treatment, age),
                                 names_from = cap_num, 
                                 values_from = c(PD, SR, Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher))
        rich_wide_df <- as.data.frame(rich_wide)
     
    # Make a separate dataframe for the 2018 and 2019 experiments   
        rwdf_18 <- subset(rich_wide_df, rich_wide_df$year == 2018)
        rwdf_19 <- subset(rich_wide_df, rich_wide_df$year == 2019)
     
    # Fit models for 2018 and make an html output table to add to appendix
        rwdf_18$initial <- rwdf_18$PD_1
        m18_PD3 <- lm(PD_3 ~ initial + color*challenge, data = rwdf_18)
        rwdf_18$initial <- rwdf_18$Shannon_1
        m18_Sh <- lm(Shannon_3 ~ initial + color*challenge, data = rwdf_18)
        rwdf_18$initial <- rwdf_18$InvSimpson_1
        m18_InvS <- lm(InvSimpson_3 ~ initial + color*challenge, data = rwdf_18)
        
        m_micro_18_t <- tab_model(m18_PD3, m18_Sh, m18_InvS, 
                  pred.labels = c("Intercept", "Pre-treatment Measure", "Signal Dulled", "Predator", "Flight Reduction",
                                  "Dulled * Predator", "Dulled * Flight"),
                  dv.labels = c("Faith's D 3rd Capture", "Shannon 3rd Capture", "Inverse Simpson 3rd Capture"))
        
        saveRDS(m_micro_18_t, here::here("5_other_outputs/m_micro_18_t.RDS"))
        
     # Fit models for 2019 and make an html output table to add to appendix
        rwdf_19$initial <- rwdf_19$PD_1
        m19_PD2 <- lm(PD_2 ~ initial + color*challenge, data = rwdf_19)
        m19_PD3 <- lm(PD_3 ~ initial + color*challenge, data = rwdf_19)
        rwdf_19$initial <- rwdf_19$Shannon_1
        m19_Sh2 <- lm(Shannon_2 ~ initial + color*challenge, data = rwdf_19)
        m19_Sh3 <- lm(Shannon_3 ~ initial + color*challenge, data = rwdf_19)
        rwdf_19$initial <- rwdf_19$InvSimpson_1
        m19_InvS2 <- lm(InvSimpson_2 ~ initial + color*challenge, data = rwdf_19)
        m19_InvS3 <- lm(InvSimpson_3 ~ initial + color*challenge, data = rwdf_19)
        
        m_micro_19_t <- tab_model(m19_PD2, m19_Sh2, m19_InvS2, m19_PD3, m19_Sh3, m19_InvS3,
                  pred.labels = c("Intercept", "Pre-treatment Measure", "Predator", "Signal Dulled", "Dulled * Predator"),
                  dv.labels = c("Faith's D 2nd Capture", "Shannon 2nd Capture", "Inverse Simpson 2nd Capture",
                                "Faith's D 3rd Capture", "Shannon 3rd Capture", "Inverse Simpson 3rd Capture"))
        
        saveRDS(m_micro_19_t, here::here("5_other_outputs/m_micro_19_t.RDS"))
      
## Group to phylum ----
  ## Create relative abundance table at phylum level
    phy_glom <- tres_ps %>%
      tax_glom(taxrank = "Phylum") %>%
      transform_sample_counts(function(x) {x / sum(x)}) %>%
      psmelt() %>%
      filter(Abundance > 0.001) %>%
      arrange(Phylum)
    
    phy_glom$sID_cap <- paste(phy_glom$sampleID, phy_glom$cap_num)
    library(reshape2)
    phy_glom2 <- dcast(phy_glom, sID_cap ~ Phylum, value.var = "Abundance")
    phy_glom2[is.na(phy_glom2)] <- 0
    phy_glom_x <- phy_glom[,c("sID_cap", "Sample", "band", "ubox", "year", "uby", "color", "challenge", "full_treatment",
                              "cap_doy", "cap_num", "age", "sampleID")]
    phy_glom3 <- join(phy_glom2, phy_glom_x, "sID_cap", match = "first")
    
## Group to order ----
    ## Create relative abundance table at order level
      ord_glom <- tres_ps %>%
        tax_glom(taxrank = "Order") %>%
        transform_sample_counts(function(x) {x / sum(x)}) %>%
        psmelt() %>%
        filter(Abundance > 0.001) %>%
        arrange(Order)
      
      ord_glom$sID_cap <- paste(ord_glom$sampleID, ord_glom$cap_num)
      library(reshape2)
      ord_glom2 <- dcast(ord_glom, sID_cap ~ Order, value.var = "Abundance")
      ord_glom2[is.na(ord_glom2)] <- 0
      ord_glom_x <- ord_glom[,c("sID_cap", "Sample", "band", "ubox", "year", "uby", "color", "challenge", "full_treatment",
                                "cap_doy", "cap_num", "age", "sampleID")]
      ord_glom3 <- join(ord_glom2, ord_glom_x, "sID_cap", match = "first")

## Plot Order ----
  # These have too many groups to be useful. Not included.
    order18 <- subset(ord_glom, ord_glom$year == 2018 & ord_glom$full_treatment != "")
    ggplot(order18, aes(x = full_treatment, y = Abundance, fill = Order)) +
      facet_grid(cap_num ~ .) +
      geom_bar(position = "fill", stat = "identity")
    
    order19 <- subset(ord_glom, ord_glom$year == 2019 & ord_glom$full_treatment != "")
    ggplot(order19, aes(x = full_treatment, y = Abundance, fill = Order)) +
      facet_grid(cap_num ~ .) +
      geom_bar(position = "fill", stat = "identity")
  
## Plot Phylum ----
  # Stacked barplot of phyla by treatment group for 2018
    phylum18 <- subset(phy_glom, phy_glom$year == 2018 & phy_glom$full_treatment != "")
    xlabs <- c("Control", "Predator", "Handicap", "Dulled", "Dull + Predator", "Dull + Handicap")
    bar18 <- ggplot(phylum18, aes(x = full_treatment, y = Abundance, fill = Phylum)) +
      facet_grid(cap_num ~ .) + theme(axis.text.x = element_text(angle = 90)) +
      geom_bar(position = "fill", stat = "identity") +
      scale_x_discrete(labels = xlabs) + theme(axis.title.x = element_blank())
    ggsave(here::here("3_r_scripts/stack_2018.pdf"), plot = bar18, width = 5.8, height = 6.2, units = "in", device = "pdf")
   
  # Stacked barplot of phyla by treatment group for 2019
    phylum19 <- subset(phy_glom, phy_glom$year == 2019 & phy_glom$full_treatment != "")
    xlabs <- c("Control", "Dulled", "Predator", "Predator + Dull")
    bar19 <- ggplot(phylum19, aes(x = full_treatment, y = Abundance, fill = Phylum)) +
      facet_grid(cap_num ~ .) +
      geom_bar(position = "fill", stat = "identity") +
      theme(axis.title.x = element_blank()) +
      scale_x_discrete(labels = xlabs)
    ggsave(here::here("3_r_scripts/stack_2019.pdf"), plot = bar19, width = 6, height = 7, units = "in", device = "pdf")

## Plot dissimilarity ----
  # Agglomerate to whatever level you want to ordinate at
    pr_glom <- tax_glom(pr1, "Genus")
    pr_use <- pr1   # No agglomeration if pr1 is used, agglomeration if pr_glom is used
  
  # Make separate phyloseq objects for each year
    ps_18 <- subset_samples(pr_use, year == 2018)
    ps_18 <- subset_samples(ps_18, full_treatment != "")
    ps_19 <- subset_samples(pr_use, year == 2019)
    ps_19 <- subset_samples(ps_19, full_treatment != "")

  # Ordination for 2018
    ps_18_prop <- transform_sample_counts(ps_18, function(otu) otu / sum(otu))
    nmds_bc_18 <- ordinate(ps_18_prop, method = "PCoA", distance = "bray")
    p18 <- plot_ordination(ps_18_prop, nmds_bc_18, color = "color", shape = "challenge", title = "Bray-Curtis PCoA") + 
      facet_wrap(cap_num ~ .) + theme_bw() + labs(color = "Color", shape = "Challenge", linetype = "Challenge") +
      scale_color_manual(values = c("orange", "slateblue")) + 
      scale_shape(labels = c("Control", "Predator", "Handicap")) + 
      stat_ellipse(aes(group = full_treatment, lty = challenge), level = 0.9) +
      scale_linetype_manual(values = c(1, 2, 3), labels = c("Control", "Predator", "Handicap"))
    ggsave(here::here("3_r_scripts/ordinate_2018.pdf"), plot = p18, width = 8, height = 4, units = "in", device = "pdf")
    
  # Ordination for 2019
    ps_19_prop <- transform_sample_counts(ps_19, function(otu) otu / sum(otu))
    nmds_bc_19 <- ordinate(ps_19_prop, method = "PCoA", distance = "bray")
    p19 <- plot_ordination(ps_19_prop, nmds_bc_19, color = "challenge", shape = "color", title = "Bray-Curtis PCoA") + 
      facet_wrap(cap_num ~ .) + theme_bw() + labs(color = "Color", shape = "Challenge", linetype = "Challenge") +
      scale_color_manual(values = c("orange", "slateblue")) + 
      stat_ellipse(aes(group = full_treatment, lty = color), level = 0.9)
    ggsave(here::here("3_r_scripts/ordinate_2019.pdf"), plot = p19, width = 9, height = 4, units = "in", device = "pdf")
    
## Permanova tests ----
    # Run separately for each capture and each year
        ps_18_1 <- subset_samples(ps_18, cap_num == 1)
        ps_18_3 <- subset_samples(ps_18, cap_num == 3)
        ps_19_1 <- subset_samples(ps_19, cap_num == 1)
        ps_19_2 <- subset_samples(ps_19, cap_num == 2)
        ps_19_3 <- subset_samples(ps_19, cap_num == 3)
    
        ps18_1_bray <- phyloseq::distance(ps_18_1, method = "bray")
        dat_temp <- data.frame(sample_data(ps_18_1))
        adonis2(ps18_1_bray ~ color*challenge, data = dat_temp)
        
        ps18_3_bray <- phyloseq::distance(ps_18_3, method = "bray")
        dat_temp <- data.frame(sample_data(ps_18_3))
        adonis2(ps18_3_bray ~ color*challenge, data = dat_temp)
        
        ps19_1_bray <- phyloseq::distance(ps_19_1, method = "bray")
        dat_temp <- data.frame(sample_data(ps_19_1))
        adonis2(ps19_1_bray ~ color*challenge, data = dat_temp)
        
        ps19_2_bray <- phyloseq::distance(ps_19_2, method = "bray")
        dat_temp <- data.frame(sample_data(ps_19_2))
        adonis2(ps19_2_bray ~ color*challenge, data = dat_temp)
        
        ps19_3_bray <- phyloseq::distance(ps_19_3, method = "bray")
        dat_temp <- data.frame(sample_data(ps_19_3))
        adonis2(ps19_3_bray ~ color*challenge, data = dat_temp)
    
      # The ps19_3 permanova is significant, so test for dispersion differences.
        beta <- betadisper(ps19_3_bray, dat_temp$full_treatment)
        permutest(beta)
        

        
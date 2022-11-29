
#Download required packages
install.packages("BiocManager")

BiocManager::install(c("Biobase", "AnnotationFilter", "SummarizedExperiment",
"tibble", "dbplyr", "BiocGenerics", "AnnotationDbi",
"S4Vectors", "robustbase", "BiocStyle", "RMariaDB"), force = TRUE)

library(Biobase)
library(AnnotationFilter)
library(SummarizedExperiment)
library(BiocGenerics)
library(AnnotationDbi)
library(S4Vectors)
library(robustbase)
library(BiocStyle)
library(RMariaDB)
library(usethis)
library(git2r)
library(devtools)


remotes::install_git("https://gitlab.gm.eurac.edu/metabolomics/BioCHRIStes", force=T, git = 'external')
#, auth_token = "glpat-UZ6XJp25AxEH97RDz-2E")

remotes::install_git("https://gitlab.gm.eurac.edu/metabolomics/BioCHRIStes",
                     credentials = git2r::cred_user_pass("dghasemisemeskandeh", "azm000n_amari.IR"))

library(BioCHRIStes)
library(biochristes7500)
browseVignettes("biochristes7500")

library(pander)
data(biochristes7500)
do.call(cbind, metadata(biochristes7500))
#-----------------------------------------------------#

#Metabolites data
rowData(biochristes7500)

#Custimized theme
blank_theme <- theme_minimal()+
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

#frequency table and Pie chart
rowData(biochristes7500) %>% View()
  as_tibble() %>% 
  #count(analyte_name) # 175 metabolites
  count(analyte_class) %>% #   6 categories
  #arrange(desc(n))
  ggplot(aes(x="", y = n, fill = fct_inorder(analyte_class))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  scale_fill_brewer("Blues") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5)) + 
  #geom_text(aes(y = n/3 + c(0, cumsum(n)[-length(n)]),
  #              label = scales::percent(n/100)), size=5) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.key.size = unit(1.3, 'cm'),
          legend.key.width = unit(0.9, 'cm'),
          legend.text  = element_text(size = 14))

ggsave("18-11-2022_Metabolites Categories.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-----------------------------------------------------#

# individuals characteristics
colData(biochristes7500)

#individual level data
concentrations(biochristes7500)[1:4, 1:3]
concentrations(biochristes7500, blessing = "none")[1:4, 1:3]

# Box plot
rowData(biochristes7500)["C14", ]
boxplot(split(log2(c14), biochristes7500$plate_name), main = "C14",
        ylab = expression(log[2]~signal), las = 2)

pts <- concentrations(biochristes7500, blessing = "none")["Putrescine", ]
plot(hist(log2(pts), breaks = 128), xlab = expression(log[2]~signal),
     main = "Putrescine")
#-----------------------------------------------------#

#preparing metabolites data for jointing with genotypes
chrisMass <- 
  concentrations(biochristes7500,
                 blessing = "none") %>%
  as_tibble() %>%
  mutate(Metabolite = row.names(concentrations(biochristes7500))) %>% 
  pivot_longer(cols = -Metabolite,
               names_to = "AID") %>% 
  pivot_wider(names_from = Metabolite) %>%
  mutate(across(!AID, function(x) log(x))) %>% 
  janitor::clean_names() %>% #janitor::row_to_names(1)
  rename(AID = aid)

#-----------------------------------------------------#
#Saving metabolites names for mediation analysis
metabolites <- chrisMass %>% select(-AID) %>% colnames() #%>% View
#targets    <- vcfmod    %>% select(-AID) %>% colnames()

#-----------------------------------------------------#

#Merge metabolites7500 with CHRIS baseline
vcfMass <-
  chris[c("AID", "Age", "Sex", "eGFRw.log.Res")] %>%
  inner_join(PCs_13K,   by = "AID") %>% #dim()
  inner_join(chrisMass, by = "AID") %>% #dim()
  inner_join(vcfmod,    by = "AID") #%>% dim()
  
#-----------------------------------------------------#
#--------------- Step 3: metabolite as outcome -------
#-----------------------------------------------------#

#Mediation Analysis -> Step 3: Metabolites as the outcome

#function for iterating retrieving 
#the regression coefficients

step3_to_table <- 
  function(mytrait, 
           mytarget, 
           myformula,
           data){
    res3 <- map2_dfr(
      .x = mytrait,
      .y = myformula,
      .f = function(trait, formula){
        res1 <- map_dfr(
          .x = mytarget,
          .f = function(SNP)
          {
            m <- lm(as.formula(formula), data = data)
            p <- summary(m)$coefficients[2, c(1,2,4)]
            return(p)
          }
        )
        res2 <- as.data.frame(res1)
        colnames(res2) <- c("Estimate",
                            "SE",
                            "Pvalue")
        Nsnp       <- length(mytarget)
        res2$SNPid <- rep(colnames(mytarget)[1:Nsnp], 1)
        return(res2)
      }
    )
    res3$pheno <- rep(colnames(mytrait), each = length(unique(res3$SNPid)))
    res4 <- 
      res3 %>%
      pivot_wider(id_cols     = SNPid,
                  names_from  = pheno,
                  values_from = c(Estimate, SE , Pvalue),
                  names_glue  = "{pheno}_{.value}")
    return(res3)
  }
#---------#

MMA_Step3_raw <-
  step3_to_table(vcfMass[metabolites],
                 vcfMass[targets],
                 paste("trait ~ SNP + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
                 vcfMass)

#Appending Locus column to the results
results_Step3_long <-
  repSNPs %>%                    
  select(SNPid, Locus) %>% 
  inner_join(MMA_Step3_raw,
             by = c("SNPid")) %>% 
  #change the order of the columns
  select(contains(c("SNPid",
                    "Locus",
                    "pheno",
                    "Estimate",
                    "SE",
                    "Pvalue"))) %>%
  #Significant association between variants and metabolites
  mutate(related = ifelse(Pvalue <= 0.05/1750,
                          "Yes",
                          "No")) #%>% filter(related == "Yes") %>% View

write.csv(results_Step3_long, "29-Nov-2022_MMD_Step3_SNPs associated with metabolites.csv", row.names = F, quote = F)

#heatmap of step 3
results_Step3_long %>%
  mutate(SNP = factor(SNPid,
                      levels = str_sort(unique(SNPid),
                                        numeric = TRUE,
                                        decreasing = TRUE))) %>%
  ggplot(aes(x = pheno, y = SNP, fill = related)) +
  geom_tile(color = "grey80") +
  #scale_fill_gradient(low = "red", high = "white")+
  coord_equal(ratio = 1) + 
  #viridis::scale_fill_viridis(name="Estimate")#(option="magma") +
  #ggthemes::theme_tufte(base_family="Helvetica") +
  #theme_classic()+
  theme(legend.title= element_text(size=7, face="bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size=6, face="bold", angle=60, vjust=1.05, hjust=1),
        axis.text.y = element_text(size=6, face="bold"),
        axis.title  = element_blank())

ggsave("29-Nov-22_MMD_Step3_SNPs associated with metabolites_Pvalue.png", 
       last_plot(), width = 16, height = 12, pointsize = 4, dpi = 300, units = "in")

#pretty heatmap
pdf('24-Nov-22_Heatmap of MMD_Step3_SNPs associated with metabolites.pdf', width=8, height = 6)
library(pheatmap)
pheatmap(results_Step3_long$Estimate,
         cluster_cols = T, 
         cluster_rows = T, 
         show_rownames = T, 
         labels_row = results_Step3_long$SNPid,
         border_color = NA, #'Black',
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = "45")

dev.off()


#----------------- Step 2: trait as covariate -----------------

step2_to_table <- 
  function(mytrait,
           mytarget,
           data){
    results1 <- lapply(mytrait,
                       function(trait){
                         map_df(mytarget,
                                function(SNP){
                                  myformula <- as.formula(
    eGFRw.log.Res ~ SNP + trait + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
                                  m <- lm(myformula, data = data)
                                  s <- coef(m)[2]
                                  p <- summary(m)$coefficients[2, c(1,2,4)] #Beta, se, Pvalue
                                  #t <- tidy(m)[2, c(1,2,4)]#broom
                                  return(p)
                                })
                       })
    results2 <- as.data.frame(do.call(cbind, results1)) #%>% clean_names() #library(janitor)
    names(results2) <- str_replace_all(names(results2),
                                       c("Estimate"   = "Estimate",
                                         "Std. Error" = "SE",
                                         "Pr\\([^\\(]*\\)" = "Pvalue"))
    #colnames(results2) <- names(mytrait)
    nloc <- length(mytarget)
    results3 <- cbind(SNPid = names(mytarget),
                      Locus = repSNPs$Locus[1:nloc],
                      results2)
    return(results3)
  }

MMA_Step2_raw <-
  step2_to_table(vcfMass[metabolites],
                 vcfMass[targets],
                 vcfMass)

#Turning the results to longer format for merging with step3
results_Step2_long <- 
  MMA_Step2_raw %>% 
  pivot_longer(cols = -c(SNPid, Locus),
               names_to = c("trait", "value"),
               names_pattern = "(.+).(Estimate|SE|Pvalue)$",
               values_to = c("score")) %>%
  pivot_wider(names_from = "value",
              values_from = "score") %>% 
  group_by(Locus, SNPid) %>%
  mutate(outlier = ifelse(is_outlier(Estimate),
                          "Yes",
                          "No")) %>% 
  ungroup()

write.csv(results_Step2_long,
          "29-Nov-2022_MMD_Step2_SNPs adjusted for metabolites_long format.csv",
          row.names = F,
          quote = F)

#---------#
#Heatmap for step 2

library(pheatmap)
png("28-Nov-2022_Heatmap_MMD_Step2_SNPs adjusted for metabolites.png", 
    units="in", res = 300, width=10, height=12)
#pdf('pheatmap2.pdf', width=18, height = 18)


results_Step2_long <-
  repSNPs %>%                    
  select(SNPid, Locus, BETA) %>% 
  inner_join(MMA_Step2_raw,
             by = "SNPid") %>% 
  #filter(SNPid %in% leadingSNPs) %>%
  select(SNPid, Locus, BETA, contains("Estimate")) %>%
  mutate(across(contains("Estimate"), function(x) x / BETA))

pheatmap(results_Step2_long[,-c(1:3)],
         cluster_cols = F,
         cluster_rows = F, 
         show_rownames = T, 
         labels_row = MMA_Step2_raw$SNPid, 
         border_color = NA, 
         fontsize = 4, 
         angle_col = "45")

dev.off()

#----------------- Merge 3 Steps -----------------

results_Step2_long %>%
  inner_join(results_Step3_long,
             by = c("SNPid" = "SNPid",
                    "Locus" = "Locus",
                    "trait" = "pheno"),
             suffix = c("_Step2", "_Step3")) %>%
  inner_join(repSNPs[c("SNPid",
                       "Locus",
                       "MARKER_ID",
                       "BETA",
                       "SEBETA",
                       "PVALUE")],
             by = c("SNPid", "Locus")) %>%
  mutate(
    REF = str_split(MARKER_ID, "\\W", simplify = TRUE)[,5],
    ALT = str_split(MARKER_ID, "\\W", simplify = TRUE)[,6],
    outlierRelated = case_when(
    outlier == "Yes" & related == "Yes" ~ "Yes/Yes",
    outlier == "Yes" & related != "Yes" ~ "Yes/No",
    outlier != "Yes" & related == "Yes" ~ "No/Yes",
    outlier != "Yes" & related != "Yes" ~ "No/No"),
    mediator = ifelse(outlier == "Yes" & related == "Yes", "Yes", "No"),
    SNP = factor(SNPid,
                 levels = str_sort(unique(SNPid),
                                   numeric = TRUE,
                                   decreasing = TRUE))) %>%
  rename(Estimate_GWAS = BETA, 
         SE_GWAS       = SEBETA, 
         Pvalue_GWAS   = PVALUE) %>%
  select(SNPid, Locus, REF, ALT, 
         Estimate_GWAS, SE_GWAS, Pvalue_GWAS, everything()) %>% #View()
  #filter(mediator == "Yes") %>% 
  write.csv(., "29-Nov-22_Heatmap_MMD_Step 1&2&3_outlierRelated traits.csv", row.names = FALSE)
  ggplot(aes(x = trait, y = SNP, fill = outlierRelated)) +
    geom_tile() +
    theme_classic() +
  scale_fill_manual(values=c('white', "grey50", "#FF6666"))+ #'#999999', "#E69F00"
  labs(x    = "",
       y    = "")+
  theme(legend.title= element_text(size=7, face="bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text  = element_text(size = 8),
        axis.text.x = element_text(size=4, face="bold", angle=90, vjust=1.05, hjust=1),
        axis.text.y = element_text(size=4, face="bold"),
        axis.title  = element_text(size=12,face="bold"))

ggsave("29-Nov-22_Heatmap_MMD_Step 1&2&3_outlierRelated traits.png", 
       last_plot(), width = 10, height = 9, pointsize = 4, dpi = 300, units = "in")

 




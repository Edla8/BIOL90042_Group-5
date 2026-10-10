install.packages("tidyverse")

library(tidyverse)
load("group5_data/prs.Rda")
load("group5_data/rnaseq.Rda")
load("group5_data/rnaseq_gene_info.Rda")
load("group5_data/drugscreen.Rda")
load("group5_data/gait.Rda")

# First question - calculating PRS 
# First you need to data wrangle using pivot_longer,
# merge ca_co_genotypes dataset with gwas_top_hits via the locus
# this will then bring beta value close to each locus for each patient
# summarise() to allow calculation for each patient

################################################################
# Population Genetics


# genotype_beta includes merged data-set
 genotype_beta <- ca_co_genotypes |> 
   pivot_longer(
   cols= contains("locus"), 
   names_to = "locus", 
   values_to = "presence")  |> 
   left_join(gwas_top_hits, 
             by = "locus") 
 # collapses columns that have the word "locus" into long tables
 
 
 # PRS calculation on merged dataset (presence * beta), 
 # leading to aggregate table showing PRS value for each patient
 genotype_beta |>  summarise(PRS = sum(presence * beta), .by = donorid) 
 
 # grabs patient with the highest PRS value
 genotype_beta |>  
   summarise(PRS = sum(presence * beta), .by = donorid) |> 
   slice_max(PRS)
 
 # grabs patient with the lowest PRS value
 genotype_beta |>  
   summarise(PRS = sum(presence * beta), .by = donorid) |> 
   slice_min(PRS)
 
 
 # making the patient PRS data frame a variable
 patient_prs <- genotype_beta |>  summarise(PRS = sum(presence * beta), .by = c(donorid,status)) 
 
 # Making boxplot
 ggplot(patient_prs, aes(x = status, y = PRS)) + geom_boxplot() + geom_point(position = position_jitter(width = 0.2)) + 
   labs(
     title = "Case vs Control PRS",
     x = "experimental status",
     y = "Count",
   )
 
 
 # cheking normality
 patient_prs |> na.omit() |>  ggplot(aes(sample = PRS)) +
   geom_qq() +
   geom_qq_line() + facet_wrap(~status)
 
 # follows the qqline therefore it's linear
 
 # Observations should be independent (check how they did experiment)
 # equal variance?
 
 patient_prs |> group_by(status) |> summarise(variance = var(PRS))
 leveneTest(PRS ~ status, data = patient_prs)
 
 
 #t-test
 t.test(PRS ~ status, data = patient_prs, var.equal = FALSE)
 # p-value gives 7.235e-09 -> <0.05 -> reject null hypothesis
  
 
 # do manual calculation to double check - Done
 # do linear regression on this for statistical significance

#######################################################################


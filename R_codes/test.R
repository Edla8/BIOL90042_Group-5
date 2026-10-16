library(tidyverse)
load("data/prs.Rda")
load("data/rnaseq.Rda")
load("data/rnaseq_gene_info.Rda")
load("data/drugscreen.Rda")
load("data/gait.Rda")

install.packages("car")
library(car)


dge

# First question - calculating PRS 
# First you need to data wrangle using pivot_longer,
# merge ca_co_genotypes dataset with gwas_top_hits via the locus
# this will then bring beta value close to each locus for each patient
# summarise() to allow calculation for each patient


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
 
 # Now we do statistical testing
 # we need to check assumptions of t-test
 # Observations should be independent (check how they did experiment)
 # Normality = ?
 # Equal variance = ?

 
 
 # checking normality
 # Pipeline: Using qqplot, then do quantative Shapiro test 
 
 patient_prs |> na.omit() |>  ggplot(aes(sample = PRS)) +
   geom_qq() +
   geom_qq_line() + facet_wrap(~status) # plots doesn't seem to fit well- double chekc with quantitative test
 
 patient_prs |> na.omit() |> 
   filter(status == "case") |> 
   pull(PRS) |> 
   shapiro.test() # returns p-value = 0.2709 -> we do not reject null
 
 
 patient_prs |> na.omit() |> 
   filter(status == "ctrl") |> 
   pull(PRS) |> 
   shapiro.test() # returns p-value = 0.08706 -> we do not reject null
 
# Overall - PRS in both case and control samples are not normally 
# distributed

 # equal variance?
 # Used levenTest (imported from Car package)
 
 patient_prs |> group_by(status) |> summarise(variance = var(PRS))
 leveneTest(PRS ~ status, data = patient_prs)
 
# Instead we do wilcoxon - 
 # assumptions: outcome variable is ordinal?
 # the groups should have similarly shaped distributions?
wilcox.test(PRS ~ status, data=patient_prs) #p-value not significant p-value = 2.279e-08
 

  

#######################################################################


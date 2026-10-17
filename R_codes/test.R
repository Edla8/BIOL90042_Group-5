library(tidyverse)
load("data/prs.Rda")
load("data/rnaseq.Rda")
load("data/rnaseq_gene_info.Rda")
load("data/drugscreen.Rda")
load("data/gait.Rda")


# Step 1: Data Wrangling & Polygenic Risk Score (PRS) Calculation
# Reshape genotype data from wide to long format using pivot_longer()
# and join with GWAS top hits to map allele effect sizes (beta values) by locus


# genotype_beta stores the merged tidy data frame
genotype_beta <- ca_co_genotypes |> 
  pivot_longer(
    cols = contains("locus"), 
    names_to = "locus", 
    values_to = "presence"
  ) |> 
  left_join(gwas_top_hits, 
            by = "locus") 
# Collapses wide locus columns into tidy long-form observations


# Calculate aggregate PRS per donor using .by argument in summarise()
genotype_beta |> summarise(PRS = sum(presence * beta), .by = donorid) 

# Identify the patient with the highest PRS value using slice_max()
genotype_beta |>  
  summarise(PRS = sum(presence * beta), .by = donorid) |>  
  slice_max(PRS)

# Identify the patient with the lowest PRS value using slice_min()
genotype_beta |>  
  summarise(PRS = sum(presence * beta), .by = donorid) |>  
  slice_min(PRS)


# Create the patient-level PRS data frame retaining metadata (status)
patient_prs <- genotype_beta |> summarise(PRS = sum(presence * beta), .by = c(donorid, status)) 

# Visualize PRS distributions across treatment groups with a boxplot and jittered points
ggplot(patient_prs, aes(x = status, y = PRS)) + geom_boxplot() + geom_point(position = position_jitter(width = 0.2)) + 
  labs(
    title = "Case vs Control PRS",
    x = "Treatment Group",
    y = "PRS",
  )

# Extract summary statistics for the case group pipeline
patient_prs |> filter(status == "case") |> summarise(
  mean = mean(PRS),
  median = median(PRS),
  sd = sd(PRS),
  variance = var(PRS),
  iqr = IQR(PRS),
  n = n()
)

# Extract summary statistics for the control group pipeline
patient_prs |> filter(status == "ctrl") |> summarise(
  mean = mean(PRS),
  median = median(PRS),
  sd = sd(PRS),
  variance = var(PRS),
  iqr = IQR(PRS),
  n = n()
)

# Step 2: Statistical Testing & Assumption Checks
# Verify core statistical test assumptions: independence, normality, and equal variance


# Check normality visually using faceted Q-Q plots
patient_prs |> na.omit() |> ggplot(aes(sample = PRS)) +
  geom_qq() +
  geom_qq_line() + facet_wrap(~status) 

# Perform Shapiro-Wilk normality test for the case group
patient_prs |> na.omit() |>  
  filter(status == "case") |>  
  pull(PRS) |>  
  shapiro.test() 

# Perform Shapiro-Wilk normality test for the control group
patient_prs |> na.omit() |>  
  filter(status == "ctrl") |>  
  pull(PRS) |>  
  shapiro.test() 

# Note: PRS distributions in both case and control groups deviate from normality


# Step 3: Non-Parametric Hypothesis Testing
# Apply the Wilcoxon rank-sum test due to non-normal distributions
wilcox.test(PRS ~ status, data = patient_prs) 

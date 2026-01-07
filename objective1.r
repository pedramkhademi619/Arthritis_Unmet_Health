# ==============================================================================
# FINAL GOLDEN SCRIPT — TABLE OUTPUT VERSION
# CCHS 2017–2018 | Bootstrap-weighted prevalence estimates
# Author: Melika Zarei
# Purpose:
#   - Validate sample logic
#   - Produce StatsCan-compliant estimates with readable tables
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. LIBRARIES
# ------------------------------------------------------------------------------
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("knitr")) install.packages("knitr")

library(tidyverse)
library(knitr)

cat("\n============================================================\n")
cat("CCHS ANALYSIS PIPELINE — INITIALIZING\n")
cat("============================================================\n")

# ------------------------------------------------------------------------------
# 1. LOAD MAIN DATA
# ------------------------------------------------------------------------------
cat("\n[1/6] Loading main CCHS data file\n")
cat("------------------------------------------------------------\n")

main_data <- read_csv(
  "data.csv",
  col_types = cols(.default = "n"),
  show_col_types = FALSE
)

cat("Main file loaded\n")
cat("Total records:", format(nrow(main_data), big.mark=","), "\n")

# ------------------------------------------------------------------------------
# 2. LOAD BOOTSTRAP WEIGHTS (FIXED WIDTH)
# ------------------------------------------------------------------------------
cat("\n[2/6] Loading bootstrap weights (fixed-width format)\n")
cat("------------------------------------------------------------\n")

bsw_widths <- c(20, 7, rep(7, 1000))
bsw_names  <- c("ADM_RNO", "FWGT_BSW", paste0("BSW", 1:1000))

bsw_data <- read_fwf(
  "bsw.txt",
  col_positions = fwf_widths(bsw_widths, bsw_names),
  col_types = cols(.default = "n"),
  show_col_types = FALSE
)

cat("Bootstrap weights loaded\n")
cat("Rows:", format(nrow(bsw_data), big.mark=","), "\n")

# ------------------------------------------------------------------------------
# 3. MERGE DATA
# ------------------------------------------------------------------------------
cat("\n[3/6] Merging survey data with bootstrap weights\n")
cat("------------------------------------------------------------\n")

main_data <- main_data %>% mutate(ADM_RNO = as.numeric(ADM_RNO))
bsw_data  <- bsw_data  %>% mutate(ADM_RNO = as.numeric(ADM_RNO))

data <- inner_join(main_data, bsw_data, by = "ADM_RNO")

cat("Merge complete\n")
cat("Combined records:", format(nrow(data), big.mark=","), "\n")

rm(main_data, bsw_data)
gc(verbose = FALSE)

# ------------------------------------------------------------------------------
# 4. DEFINE STUDY POPULATION
# ------------------------------------------------------------------------------
cat("\n[4/6] Defining analytical study population\n")
cat("------------------------------------------------------------\n")

bsw_cols <- paste0("BSW", 1:1000)

valid_provinces <- data %>%
  group_by(GEO_PRV) %>%
  summarise(Has_Data = mean(UCN_005 %in% c(1, 2), na.rm = TRUE) > 0) %>%
  filter(Has_Data) %>%
  pull(GEO_PRV)

cat("Participating provinces:",
    paste(valid_provinces, collapse = ", "), "\n")

analysis_data <- data %>%
  filter(DHHGAGE >= 13) %>%        # Age 65+
  filter(CCC_050 == 1) %>%         # Arthritis
  filter(GEO_PRV %in% valid_provinces) %>%
  mutate(
    Has_Unmet_Need = if_else(UCN_005 == 1, 1, 0, missing = 0),
    Age_Group = if_else(DHHGAGE == 16, "80+", "65-79"),
    Sex_Label = if_else(DHH_SEX == 1, "Male", "Female")
  )

cat("Final study population size:",
    format(nrow(analysis_data), big.mark=","), "\n")

# ------------------------------------------------------------------------------
# 5. BOOTSTRAP ESTIMATION FUNCTION
# ------------------------------------------------------------------------------
calculate_stats_with_bootstrap <- function(df, group_var, target_var,
                                           bsw_cols, weight_var = "WTS_M") {

  groups <- unique(df[[group_var]])
  out <- list()

  for (g in groups) {

    gd <- df[df[[group_var]] == g, ]
    if (nrow(gd) == 0) next

    # Main estimate
    theta <- sum(gd[[weight_var]] * gd[[target_var]]) /
             sum(gd[[weight_var]]) * 100

    # Bootstrap replicates
    W <- as.matrix(gd[, bsw_cols])
    y <- gd[[target_var]]

    theta_b <- (t(W) %*% y) / colSums(W) * 100
    var_b   <- mean((theta_b - theta)^2)
    se      <- sqrt(var_b)
    cv      <- (se / theta) * 100

    out[[g]] <- tibble(
      Group = g,
      Unweighted_n = nrow(gd),
      Estimate_percent = round(theta, 2),
      CI_lower = round(max(0, theta - 1.96 * se), 2),
      CI_upper = round(min(100, theta + 1.96 * se), 2),
      CV_percent = round(cv, 1)
    )
  }

  bind_rows(out)
}

# ------------------------------------------------------------------------------
# 6. TABLE PRINTER
# ------------------------------------------------------------------------------
print_table <- function(df, title) {
  cat("\n", title, "\n", sep = "")
  cat(rep("-", nchar(title)), "\n")
  cat(kable(df, align = "c"), sep = "\n")
  cat("\n")
}

# ------------------------------------------------------------------------------
# 7. ANALYSIS & TABLE OUTPUT
# ------------------------------------------------------------------------------
cat("\n[6/6] Running final analyses\n")
cat("------------------------------------------------------------\n")

# ---- Part A: Prevalence -------------------------------------------------------
age_tab <- calculate_stats_with_bootstrap(
  analysis_data, "Age_Group", "Has_Unmet_Need", bsw_cols
)

sex_tab <- calculate_stats_with_bootstrap(
  analysis_data, "Sex_Label", "Has_Unmet_Need", bsw_cols
)

print_table(
  age_tab,
  "Table 1. Prevalence of Unmet Health Care Needs by Age Group"
)

print_table(
  sex_tab,
  "Table 2. Prevalence of Unmet Health Care Needs by Sex"
)

# ---- Part B: Reasons (Conditional Subsample) ---------------------------------
unmet_pop <- analysis_data %>% filter(Has_Unmet_Need == 1)
unmet_pop$All <- "Unmet_Pop"

cat("\nSubpopulation with unmet needs:",
    format(nrow(unmet_pop), big.mark=","), "\n")

reasons <- list(
  Cost = "UCN_010G",
  Transportation = "UCN_010J",
  Dr_Unnecessary = "UCN_010I",
  Choice = "UCN_010H",
  Inadequate = "UCN_010F",
  Other = "UCN_010K"
)

reason_tables <- list()

for (r in names(reasons)) {
  var <- reasons[[r]]
  unmet_pop$Temp <- ifelse(unmet_pop[[var]] == 1, 1, 0)

  res <- calculate_stats_with_bootstrap(
    unmet_pop, "All", "Temp", bsw_cols
  ) %>% mutate(Reason = r)

  reason_tables[[r]] <- res
}

final_reason_table <- bind_rows(reason_tables) %>%
  select(Reason, everything(), -Group)

print_table(
  final_reason_table,
  "Table 3. Reasons for Unmet Health Care Needs (Conditional Sample)"
)

cat("\n============================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")




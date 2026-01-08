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
cat("Total records:", format(nrow(main_data), big.mark = ","), "\n")

# ------------------------------------------------------------------------------
# 2. LOAD BOOTSTRAP WEIGHTS (FIXED WIDTH)
# ------------------------------------------------------------------------------
cat("\n[2/6] Loading bootstrap weights (fixed-width format)\n")
cat("------------------------------------------------------------\n")

bsw_widths <- c(20, 7, rep(7, 1000))
bsw_names <- c("ADM_RNO", "FWGT_BSW", paste0("BSW", 1:1000))

bsw_data <- read_fwf(
  "bsw.txt",
  col_positions = fwf_widths(bsw_widths, bsw_names),
  col_types = cols(.default = "n"),
  show_col_types = FALSE
)

cat("Bootstrap weights loaded\n")
cat("Rows:", format(nrow(bsw_data), big.mark = ","), "\n")

# ------------------------------------------------------------------------------
# 3. MERGE DATA
# ------------------------------------------------------------------------------
cat("\n[3/6] Merging survey data with bootstrap weights\n")
cat("------------------------------------------------------------\n")

main_data <- main_data %>% mutate(ADM_RNO = as.numeric(ADM_RNO))
bsw_data <- bsw_data %>% mutate(ADM_RNO = as.numeric(ADM_RNO))

data <- inner_join(main_data, bsw_data, by = "ADM_RNO")

cat("Merge complete\n")
cat("Combined records:", format(nrow(data), big.mark = ","), "\n")

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

cat(
  "Participating provinces:",
  paste(valid_provinces, collapse = ", "), "\n"
)

analysis_data <- data %>%
  filter(DHHGAGE >= 13) %>% # Age 65+
  filter(CCC_050 == 1) %>% # Arthritis
  filter(GEO_PRV %in% valid_provinces) %>%
  mutate(
    # --- Core variables ---
    Has_Unmet_Need = if_else(UCN_005 == 1, 1, 0, missing = 0),
    Age_Group = if_else(DHHGAGE == 16, "80+", "65-79"),
    Sex_Label = if_else(DHH_SEX == 1, "Male", "Female"),

    # --- START OF NEW VARIABLES ---

    # 1. Province Labels
    Province_Label = case_when(
      GEO_PRV == 10 ~ "NL",
      GEO_PRV == 11 ~ "PE",
      GEO_PRV == 12 ~ "NS",
      GEO_PRV == 13 ~ "NB",
      GEO_PRV == 24 ~ "QC",
      GEO_PRV == 35 ~ "ON",
      GEO_PRV == 46 ~ "MB",
      GEO_PRV == 47 ~ "SK",
      GEO_PRV == 48 ~ "AB",
      GEO_PRV == 59 ~ "BC",
      TRUE ~ "Territories/Other"
    ),

    # 2. Functional Limitations (WDM Variables)

    # Mobility limitation
    Mobility_Label = case_when(
      WDM_015 == 1 ~ "No Mobility Limit",
      WDM_015 %in% c(2, 3, 4) ~ "Has Mobility Limit",
      TRUE ~ NA_character_
    ),

    # Cognitive limitation
    Cognition_Label = case_when(
      WDM_020 == 1 ~ "No Cognition Limit",
      WDM_020 %in% c(2, 3, 4) ~ "Has Cognition Limit",
      TRUE ~ NA_character_
    )

    # --- END OF NEW VARIABLES ---
  )

cat(
  "Final study population size:",
  format(nrow(analysis_data), big.mark = ","), "\n"
)


# ------------------------------------------------------------------------------
# 5. BOOTSTRAP ESTIMATION FUNCTION (ROBUST VERSION)
# ------------------------------------------------------------------------------
# این بخش را بین مرحله ۴ (Define Study Population) و مرحله ۶ (Table 1 Helper) قرار دهید

calculate_stats_with_bootstrap <- function(df, group_var, target_var,
                                           bsw_cols, weight_var = "WTS_M") {
  # 1. Get valid groups (removing NA categories from the loop list)
  groups <- unique(na.omit(df[[group_var]]))
  out <- list()

  for (g in groups) {
    # 2. ROBUST SUBSETTING: Use which() to ignore NAs.
    # This ensures we don't get NA rows when the grouping variable has NAs
    gd <- df[which(df[[group_var]] == g), ]

    if (nrow(gd) == 0) next

    # Main estimate
    theta <- sum(gd[[weight_var]] * gd[[target_var]]) /
      sum(gd[[weight_var]]) * 100

    # Bootstrap replicates
    W <- as.matrix(gd[, bsw_cols])
    y <- gd[[target_var]]

    # Matrix multiplication for speed
    theta_b <- (t(W) %*% y) / colSums(W) * 100
    var_b <- mean((theta_b - theta)^2)
    se <- sqrt(var_b)
    cv <- (se / theta) * 100

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

calc_col_pct <- function(df, target_col) {
  # Create a dummy variable for each level of the target column
  levels <- unique(na.omit(df[[target_col]]))
  res_list <- list()

  for (l in levels) {
    # --- START OF FIX ---
    # Use dplyr::if_else which handles NA values gracefully via the `missing` argument.
    # This ensures NAs in the target column become 0 in the flag, not NA.
    df$temp_flag <- if_else(df[[target_col]] == l, 1, 0, missing = 0)
    # --- END OF FIX ---

    # Use existing function with a dummy group
    df$Total_Group <- "Total Population"

    res <- calculate_stats_with_bootstrap(df, "Total_Group", "temp_flag", bsw_cols)
    res$Characteristic <- paste(target_col, ":", l)
    res_list[[l]] <- res %>% select(Characteristic, Estimate_percent, CI_lower, CI_upper)
  }
  bind_rows(res_list)
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

cat(
  "\nSubpopulation with unmet needs:",
  format(nrow(unmet_pop), big.mark = ","), "\n"
)

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


# ------------------------------------------------------------------------------
# EXTRA: STRATIFIED TABLES (REQUESTED BY JESSICA)
# ------------------------------------------------------------------------------

# ---- Part C: Prevalence by Province ------------------------------------------
prov_tab <- calculate_stats_with_bootstrap(
  analysis_data, "Province_Label", "Has_Unmet_Need", bsw_cols
)

print_table(
  prov_tab,
  "Table 4. Prevalence of Unmet Health Care Needs by Province"
)

# ---- Part D: Prevalence by Functional Limitations ----------------------------

# 1. Mobility
mob_tab <- calculate_stats_with_bootstrap(
  analysis_data, "Mobility_Label", "Has_Unmet_Need", bsw_cols
)
print_table(
  mob_tab,
  "Table 5. Prevalence by Mobility Limitation (Walking/Climbing)"
)

# 2. Cognition
cog_tab <- calculate_stats_with_bootstrap(
  analysis_data, "Cognition_Label", "Has_Unmet_Need", bsw_cols
)
print_table(
  cog_tab,
  "Table 6. Prevalence by Cognition Limitation (Memory/Concentration)"
)

# ---- Part E: Table 1 - Study Sample Description (Estimating Population %) ----
# Note: To create Table 1 (Descriptive), we treat the demographic variable
# as the outcome to see its distribution in the population.

cat("\n--- TABLE 1 GENERATION (Descriptive Characteristics) ---\n")

# Helper to calculate column percentage
calc_col_pct <- function(df, target_col) {
  levels <- unique(na.omit(df[[target_col]]))
  res_list <- list()

  for (l in levels) {
    # تغییر مهم: افزودن missing = 0
    df$temp_flag <- if_else(df[[target_col]] == l, 1, 0, missing = 0)

    df$Total_Group <- "Total Population"

    # فراخوانی تابع اصلاح شده بالا
    res <- calculate_stats_with_bootstrap(df, "Total_Group", "temp_flag", bsw_cols)

    res$Characteristic <- paste(target_col, ":", l)
    res_list[[l]] <- res %>% select(Characteristic, Estimate_percent, CI_lower, CI_upper)
  }
  bind_rows(res_list)
}

# Run for key demographics
t1_age <- calc_col_pct(analysis_data, "Age_Group")
t1_sex <- calc_col_pct(analysis_data, "Sex_Label")
t1_mob <- calc_col_pct(analysis_data, "Mobility_Label")
t1_cog <- calc_col_pct(analysis_data, "Cognition_Label") # <-- Add this line

# Combine all characteristics into the final table
table1_final <- bind_rows(t1_age, t1_sex, t1_mob, t1_cog) # <-- And add it here

print_table(
  table1_final,
  "Table 1. Weighted Sample Characteristics (Descriptive)"
)




cat("\n============================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

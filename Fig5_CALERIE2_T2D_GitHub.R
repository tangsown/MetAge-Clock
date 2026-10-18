# Run from the repository root. CALERIE2 has both scores; T2D has Light-MetAge only.
# The public T2D CSV omits individual covariates and 1-year HbA1c outcomes.
# Full analysis requires an ethics-approved T2D CSV with the same file name and
# Age_baseline, Sex, BMI_baseline, FPG_baseline_mmol_L, and
# HbA1c_target_lt6_5_1y columns. The regression code below is retained.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(nlme)
  library(emmeans)
  library(lmerTest)
  library(broom)
  library(gghalves)
})

data_dir <- Sys.getenv("FIG5_DATA_DIR", unset = ".")
if (!file.exists(file.path(data_dir, "CALERIE2_Fig5_participants.csv"))) {
  data_dir <- file.path("outputs", "Fig5_GitHub_source_data")
}
cr_file <- file.path(data_dir, "CALERIE2_Fig5_participants.csv")
t2d_file <- file.path(data_dir, "T2D373_Fig5_participants.csv")
if (!all(file.exists(c(cr_file, t2d_file)))) stop("Fig5 source CSV files not found.")

cr_wide <- utils::read.csv(cr_file, check.names = FALSE,
                           na.strings = c("", "NA"))
t2d_wide <- utils::read.csv(t2d_file, check.names = FALSE,
                            na.strings = c("", "NA"))
restricted_t2d_fields <- c(
  "Age_baseline", "Sex", "BMI_baseline", "FPG_baseline_mmol_L",
  "HbA1c_target_lt6_5_1y"
)
missing_t2d_fields <- setdiff(restricted_t2d_fields, names(t2d_wide))
if (length(missing_t2d_fields)) {
  stop("Full Fig 5 analysis requires restricted T2D fields: ",
       paste(missing_t2d_fields, collapse = ", "),
       ". See README.md for the access process.")
}
cr_visits <- c("Baseline", "1m", "3m", "6m", "9m", "12m", "18m", "24m")
t2d_visits <- c("Baseline", "2w", "3m", "6m", "9m", "12m")
treatment_levels <- paste("Group", LETTERS[1:4])

stopifnot(
  nrow(cr_wide) == 220L, nrow(t2d_wide) == 373L,
  !anyDuplicated(cr_wide$Participant_ID),
  !anyDuplicated(t2d_wide$Participant_ID),
  all(cr_wide$Group %in% c("Control", "Caloric Restriction")),
  all(t2d_wide$Group %in% LETTERS[1:4])
)

cr_long <- cr_wide %>%
  pivot_longer(
    cols = matches("^(Light_MetAge_Accel|Light_MetAge|MetAge_Accel|MetAge)_"),
    names_to = c(".value", "Visit"),
    names_pattern = paste0(
      "^(Light_MetAge_Accel|Light_MetAge|MetAge_Accel|MetAge)_",
      "(Baseline|1m|3m|6m|9m|12m|18m|24m)$"
    )
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Control", "Caloric Restriction")),
    Visit = factor(Visit, levels = cr_visits),
    Sex = factor(Sex, levels = c("Male", "Female")),
    Delta_MetAge = MetAge - MetAge[Visit == "Baseline"][match(
      Participant_ID, Participant_ID[Visit == "Baseline"]
    )],
    Delta_Light_MetAge = Light_MetAge - Light_MetAge[Visit == "Baseline"][match(
      Participant_ID, Participant_ID[Visit == "Baseline"]
    )]
  ) %>%
  arrange(Participant_ID, Visit)

t2d_long <- t2d_wide %>%
  pivot_longer(
    cols = matches("^Light_MetAge(_Accel)?_"),
    names_to = c(".value", "Visit"),
    names_pattern =
      "^(Light_MetAge_Accel|Light_MetAge)_(Baseline|2w|3m|6m|9m|12m)$"
  ) %>%
  mutate(
    Group = factor(Group, levels = LETTERS[1:4]),
    Treatment_group = factor(paste("Group", Group), levels = treatment_levels),
    Visit = factor(Visit, levels = t2d_visits),
    Sex = factor(Sex, levels = c("Male", "Female")),
    Delta_Light_MetAge = Light_MetAge - Light_MetAge[Visit == "Baseline"][match(
      Participant_ID, Participant_ID[Visit == "Baseline"]
    )]
  ) %>%
  arrange(Participant_ID, Visit)

stopifnot(
  nrow(cr_long) == 220L * length(cr_visits),
  nrow(t2d_long) == 373L * length(t2d_visits),
  !anyDuplicated(cr_long[c("Participant_ID", "Visit")]),
  !anyDuplicated(t2d_long[c("Participant_ID", "Visit")])
)

summarise_scores <- function(data, group, score) {
  data %>%
    filter(is.finite(.data[[score]])) %>%
    group_by(across(all_of(c(group, "Visit")))) %>%
    summarise(
      n = n_distinct(Participant_ID),
      mean = mean(.data[[score]]), sd = sd(.data[[score]]),
      se = sd / sqrt(n),
      ci_low = mean - qt(0.975, pmax(n - 1L, 1L)) * se,
      ci_high = mean + qt(0.975, pmax(n - 1L, 1L)) * se,
      .groups = "drop"
    )
}

# CALERIE2: randomized-group by visit and differences in change from baseline.
cr_model_data <- cr_long %>%
  mutate(Participant_ID = factor(Participant_ID))
fit_cr <- function(outcome) {
  nlme::lme(
    fixed = stats::as.formula(paste0(
      outcome, " ~ Group * Visit + Age_baseline + Sex"
    )),
    random = ~ 1 | Participant_ID,
    data = cr_model_data, na.action = na.exclude, method = "REML",
    control = nlme::lmeControl(opt = "optim")
  )
}
cr_metage_model <- fit_cr("MetAge")
cr_light_model <- fit_cr("Light_MetAge")
cr_global_p <- function(model) {
  as.numeric(anova(model, type = "marginal")["Group:Visit", "p-value"])
}
cr_interactions <- function(model) {
  as.data.frame(summary(model)$tTable) %>%
    tibble::rownames_to_column("term") %>%
    filter(grepl("Group.*:Visit", term)) %>%
    transmute(
      Visit = sub("^.*Visit", "", term), estimate_years = Value,
      se_years = Std.Error, df = DF,
      ci_low_years = estimate_years + qt(0.025, df) * se_years,
      ci_high_years = estimate_years + qt(0.975, df) * se_years,
      p_raw = `p-value`
    ) %>%
    mutate(
      p_holm_7_followups = p.adjust(p_raw, method = "holm"),
      Visit = factor(Visit, levels = cr_visits)
    ) %>%
    arrange(Visit) %>%
    mutate(Visit = as.character(Visit))
}
cr_metage_interactions <- cr_interactions(cr_metage_model)
cr_light_interactions <- cr_interactions(cr_light_model)

cr_summary <- list(
  MetAge = summarise_scores(cr_long, "Group", "MetAge"),
  Change = summarise_scores(cr_long, "Group", "Delta_MetAge"),
  Acceleration = summarise_scores(cr_long, "Group", "MetAge_Accel"),
  Light_MetAge = summarise_scores(cr_long, "Group", "Light_MetAge"),
  Light_Change = summarise_scores(cr_long, "Group", "Delta_Light_MetAge"),
  Light_Acceleration = summarise_scores(cr_long, "Group", "Light_MetAge_Accel")
)

s_table1_cr <- bind_rows(lapply(names(cr_summary)[1:3], function(name) {
  cr_summary[[name]] %>% mutate(Outcome = name, .before = 1)
}))
s_table1_cr_light <- bind_rows(lapply(names(cr_summary)[4:6], function(name) {
  cr_summary[[name]] %>% mutate(Outcome = name, .before = 1)
}))
extract_cr_global <- function(model, outcome) {
  test <- as.data.frame(anova(model, type = "marginal")) %>%
    tibble::rownames_to_column("Term") %>%
    filter(Term == "Group:Visit")
  tibble(
    Outcome = outcome, Test = "Group x visit", Visit = "All visits",
    numerator_df = test$numDF, denominator_df = test$denDF,
    F_statistic = test$`F-value`, p_value = test$`p-value`
  )
}
s_table2_cr_global <- bind_rows(
  extract_cr_global(cr_metage_model, "MetAge"),
  extract_cr_global(cr_light_model, "Light-MetAge")
)
s_table2_cr_followups <- bind_rows(
  cr_metage_interactions %>% mutate(Outcome = "MetAge", .before = 1),
  cr_light_interactions %>% mutate(Outcome = "Light-MetAge", .before = 1)
) %>%
  as_tibble()

# T2D: baseline-adjusted maintenance-phase mixed models (3m to 12m).
t2d_baseline <- t2d_long %>%
  filter(Visit == "Baseline") %>%
  transmute(
    Participant_ID,
    Baseline_Light_MetAge = Light_MetAge,
    Baseline_Light_MetAge_Accel = Light_MetAge_Accel
  )
t2d_model_data <- t2d_long %>%
  filter(Visit %in% c("3m", "6m", "9m", "12m")) %>%
  select(Participant_ID, Treatment_group, Visit,
         Light_MetAge, Light_MetAge_Accel) %>%
  left_join(t2d_baseline, by = "Participant_ID") %>%
  filter(if_all(c(Light_MetAge, Light_MetAge_Accel,
                  Baseline_Light_MetAge, Baseline_Light_MetAge_Accel),
                is.finite)) %>%
  mutate(
    Visit = factor(as.character(Visit), levels = c("3m", "6m", "9m", "12m")),
    Participant_ID = factor(Participant_ID)
  ) %>%
  as.data.frame()
stopifnot(nrow(t2d_model_data) == 936L,
          n_distinct(t2d_model_data$Participant_ID) == 330L)

fit_t2d <- function(outcome, baseline) {
  base <- paste(outcome, "~", baseline)
  no_group <- lmerTest::lmer(
    as.formula(paste(base, "+ Visit + (1 | Participant_ID)")),
    data = t2d_model_data, REML = FALSE
  )
  parallel <- lmerTest::lmer(
    as.formula(paste(base, "+ Treatment_group + Visit + (1 | Participant_ID)")),
    data = t2d_model_data, REML = FALSE
  )
  full <- lmerTest::lmer(
    as.formula(paste(base, "+ Treatment_group * Visit + (1 | Participant_ID)")),
    data = t2d_model_data, REML = FALSE
  )
  list(
    full = full,
    any_group = stats::anova(no_group, full),
    interaction = stats::anova(parallel, full)
  )
}
t2d_light_fit <- fit_t2d("Light_MetAge", "Baseline_Light_MetAge")
t2d_accel_fit <- fit_t2d("Light_MetAge_Accel", "Baseline_Light_MetAge_Accel")

extract_global <- function(fit, outcome, test) {
  row <- fit[[test]]
  tibble(
    Outcome = outcome, Test = test,
    chi_square = as.numeric(row[["Chisq"]][2]),
    df = as.numeric(row[["Df"]][2]),
    p_value = as.numeric(row[["Pr(>Chisq)"]][2]),
    model_participants = n_distinct(t2d_model_data$Participant_ID),
    model_observations = nrow(t2d_model_data)
  )
}
s_table3_t2d_global <- bind_rows(
  extract_global(t2d_light_fit, "Light-MetAge", "any_group"),
  extract_global(t2d_light_fit, "Light-MetAge", "interaction"),
  extract_global(t2d_accel_fit, "Light-MetAge acceleration", "any_group"),
  extract_global(t2d_accel_fit, "Light-MetAge acceleration", "interaction")
)

metformin_by_visit <- function(fit, outcome) {
  em <- emmeans::emmeans(
    fit$full, specs = ~ Treatment_group | Visit,
    lmer.df = "satterthwaite"
  )
  contrast <- emmeans::contrast(
    em,
    method = list("(A+C)/2 - (B+D)/2" = c(0.5, -0.5, 0.5, -0.5)),
    adjust = "none"
  )
  as_tibble(summary(contrast, infer = c(TRUE, TRUE), adjust = "none")) %>%
    mutate(Outcome = outcome, p_holm_4_visits = p.adjust(p.value, "holm"),
           .before = 1)
}
t2d_metformin_contrasts <- bind_rows(
  metformin_by_visit(t2d_light_fit, "Light-MetAge"),
  metformin_by_visit(t2d_accel_fit, "Light-MetAge acceleration")
) %>%
  mutate(p_holm_8_outcome_visits = p.adjust(p.value, "holm"))

t2d_summary <- list(
  Light_MetAge = summarise_scores(t2d_long, "Treatment_group", "Light_MetAge"),
  Change = summarise_scores(t2d_long, "Treatment_group", "Delta_Light_MetAge"),
  Acceleration = summarise_scores(t2d_long, "Treatment_group", "Light_MetAge_Accel")
)
s_table3_t2d_descriptive <- bind_rows(lapply(names(t2d_summary), function(name) {
  t2d_summary[[name]] %>% mutate(Outcome = name, .before = 1)
}))
s_table4_t2d_contrasts <- t2d_metformin_contrasts %>%
  transmute(
    Outcome, Visit = as.character(Visit), Contrast = contrast,
    estimate_years = estimate, se_years = SE, df,
    ci_low_years = lower.CL, ci_high_years = upper.CL,
    p_raw = p.value, p_holm_4_visits, p_holm_8_outcome_visits
  )

# Fig5G: HbA1c <6.5% at one year.
lr_data <- t2d_wide %>%
  transmute(
    Participant_ID,
    Outcome = as.integer(HbA1c_target_lt6_5_1y),
    Sex = factor(Sex, levels = c("Male", "Female")),
    BMI = BMI_baseline,
    Baseline_FPG_mmol_L = FPG_baseline_mmol_L,
    Light_MetAge = Light_MetAge_Baseline,
    Treatment_group = relevel(
      factor(paste("Group", Group), levels = treatment_levels), ref = "Group D"
    )
  ) %>%
  filter(complete.cases(Outcome, Sex, BMI, Baseline_FPG_mmol_L,
                        Light_MetAge, Treatment_group))
lr_sd <- sd(lr_data$Light_MetAge)
lr_data <- lr_data %>%
  mutate(Light_MetAge_perSD = (Light_MetAge - mean(Light_MetAge)) / lr_sd)
lr_primary <- glm(
  Outcome ~ Light_MetAge_perSD + Sex + BMI + Baseline_FPG_mmol_L,
  data = lr_data, family = binomial()
)
lr_group_adjusted <- glm(
  Outcome ~ Light_MetAge_perSD + Sex + BMI + Baseline_FPG_mmol_L +
    Treatment_group,
  data = lr_data, family = binomial()
)
tidy_or <- function(model) {
  broom::tidy(model) %>%
    mutate(
      OR = exp(estimate),
      CI_low = exp(estimate - qnorm(0.975) * std.error),
      CI_high = exp(estimate + qnorm(0.975) * std.error),
      Predictor = recode(term,
        Light_MetAge_perSD = "Light-MetAge", SexFemale = "Female vs male",
        BMI = "BMI", Baseline_FPG_mmol_L = "Baseline FPG",
        .default = term
      )
    )
}
lr_primary_results <- tidy_or(lr_primary)
lr_group_results <- tidy_or(lr_group_adjusted)
s_table5_t2d_logistic <- bind_rows(
  lr_primary_results %>% filter(term != "(Intercept)") %>%
    mutate(Model = "Primary", .before = 1),
  lr_group_results %>% filter(term == "Light_MetAge_perSD") %>%
    mutate(Model = "Plus randomized group", .before = 1)
) %>%
  transmute(
    Model, Predictor, OR, ci_low = CI_low, ci_high = CI_high,
    p_value = p.value, n = nrow(lr_data), events = sum(lr_data$Outcome == 1L),
    light_metage_sd_years = lr_sd
  )

# Fig5A-C and supplementary CALERIE2 panels.
pt_to_mm <- 25.4 / 72.27
cr_colors <- c("Control" = "#377EB8", "Caloric Restriction" = "#E64B35")
make_cr_plot <- function(data, y_label, title_text, zero_line = FALSE) {
  dodge <- position_dodge(width = 0.3)
  p <- ggplot(data, aes(Visit, mean, color = Group, group = Group))
  if (zero_line) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.5 * pt_to_mm)
  }
  p +
    geom_line(linewidth = 0.75 * pt_to_mm, position = dodge) +
    geom_point(size = 1, position = dodge) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12,
                  linewidth = 0.5 * pt_to_mm, position = dodge) +
    scale_color_manual(values = cr_colors) +
    labs(title = title_text,
         subtitle = "Points show group means; error bars show 95% CIs",
         x = "Visit", y = y_label, color = "Randomized group") +
    theme_bw(base_family = "Arial", base_size = 7) +
    theme(
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, face = "plain"),
      axis.ticks = element_line(linewidth = 0.5 * pt_to_mm),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.5 * pt_to_mm,
                                      color = "grey90"),
      panel.border = element_rect(fill = NA, linewidth = 0.75 * pt_to_mm),
      legend.position = "bottom",
      legend.title = element_text(size = 7, face = "plain"),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "plain"),
      plot.subtitle = element_text(size = 6, face = "plain"),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}
format_p <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.4f", p)
add_global_p <- function(plot, p) {
  plot + annotate(
    "text", x = Inf, y = Inf,
    label = paste0("Group \u00d7 visit: ", format_p(p)),
    family = "Arial", size = 6 * pt_to_mm, hjust = 1.08, vjust = 1.25
  )
}
significance_symbol <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"
  )
}
add_cr_interactions <- function(plot, summary_data, tests) {
  span <- diff(range(c(summary_data$ci_low, summary_data$ci_high), na.rm = TRUE))
  if (!is.finite(span) || span == 0) span <- 1
  annotation <- summary_data %>%
    group_by(Visit) %>%
    summarise(y = max(ci_high), .groups = "drop") %>%
    mutate(Visit = as.character(Visit)) %>%
    inner_join(tests %>% mutate(symbol = significance_symbol(p_holm_7_followups)),
               by = "Visit") %>%
    mutate(x = match(Visit, cr_visits),
           bracket_y = y + 0.055 * span,
           tip_y = bracket_y - 0.025 * span,
           label_y = bracket_y + 0.025 * span)
  plot +
    geom_segment(data = annotation,
      aes(x = x - 0.11, xend = x + 0.11, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE, linewidth = 0.5 * pt_to_mm) +
    geom_segment(data = annotation,
      aes(x = x - 0.11, xend = x - 0.11, y = bracket_y, yend = tip_y),
      inherit.aes = FALSE, linewidth = 0.5 * pt_to_mm) +
    geom_segment(data = annotation,
      aes(x = x + 0.11, xend = x + 0.11, y = bracket_y, yend = tip_y),
      inherit.aes = FALSE, linewidth = 0.5 * pt_to_mm) +
    geom_text(data = annotation, aes(x = x, y = label_y, label = symbol),
      inherit.aes = FALSE, family = "Arial", size = 7 * pt_to_mm) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))
}

p_cr_metage <- add_global_p(make_cr_plot(
  cr_summary$MetAge, "MetAge (years)", "MetAge across CALERIE2 follow-up"
), cr_global_p(cr_metage_model))
p_cr_delta <- add_cr_interactions(make_cr_plot(
  cr_summary$Change, "Change from baseline (years)",
  "Within-person change in MetAge from baseline", zero_line = TRUE
), cr_summary$Change, cr_metage_interactions)
p_cr_accel <- make_cr_plot(
  cr_summary$Acceleration, "MetAge-Accel (years)",
  "MetAge-Accel across CALERIE2 follow-up"
)
p_cr_column <- (p_cr_metage / p_cr_delta / p_cr_accel) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

p_cr_light <- add_global_p(make_cr_plot(
  cr_summary$Light_MetAge, "Light-MetAge (years)",
  "Light-MetAge across CALERIE2 follow-up"
), cr_global_p(cr_light_model))
p_cr_light_delta <- add_cr_interactions(make_cr_plot(
  cr_summary$Light_Change, "Change from baseline (years)",
  "Within-person change in Light-MetAge from baseline", zero_line = TRUE
), cr_summary$Light_Change, cr_light_interactions)
p_cr_light_accel <- make_cr_plot(
  cr_summary$Light_Acceleration, "Light-MetAge-Accel (years)",
  "Light-MetAge-Accel across CALERIE2 follow-up"
)
p_cr_light_column <- (p_cr_light / p_cr_light_delta / p_cr_light_accel) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Fig5D-F: four randomized T2D groups.
t2d_colors <- c(
  "Group A" = "#BC3C29FF", "Group B" = "#0072B5FF",
  "Group C" = "#E18727FF", "Group D" = "#20854EFF"
)
t2d_shapes <- c("Group A" = 16, "Group B" = 17,
                "Group C" = 15, "Group D" = 18)
make_t2d_plot <- function(data, y_label, title_text,
                          zero_line = FALSE, corner_label = NULL) {
  dodge <- position_dodge(width = 0.36)
  p <- ggplot(data, aes(Visit, mean, color = Treatment_group,
                        shape = Treatment_group, group = Treatment_group))
  if (zero_line) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.75 * pt_to_mm)
  }
  p <- p +
    geom_line(linewidth = 1 * pt_to_mm, position = dodge, na.rm = TRUE) +
    geom_point(size = 1.6, position = dodge, na.rm = TRUE) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.30,
                  linewidth = 0.75 * pt_to_mm, position = dodge,
                  na.rm = TRUE) +
    scale_color_manual(name = "Randomized group", values = t2d_colors,
                       breaks = names(t2d_colors), drop = FALSE) +
    scale_shape_manual(name = "Randomized group", values = t2d_shapes,
                       breaks = names(t2d_shapes), drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    labs(title = title_text,
         subtitle = "Points show group means; error bars show 95% CIs",
         x = "Visit", y = y_label) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
           shape = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_bw(base_family = "Arial", base_size = 7) +
    theme(
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, face = "plain"),
      axis.ticks = element_line(linewidth = 0.75 * pt_to_mm),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.75 * pt_to_mm,
                                      color = "grey90"),
      panel.border = element_rect(fill = NA, linewidth = 1 * pt_to_mm),
      legend.position = "bottom",
      legend.title = element_text(size = 7, face = "plain"),
      legend.text = element_text(size = 6),
      legend.key.width = grid::unit(8, "mm"),
      plot.title = element_text(size = 8, face = "plain"),
      plot.subtitle = element_text(size = 6, face = "plain"),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
  if (!is.null(corner_label)) {
    p <- p + annotate("text", x = Inf, y = Inf, label = corner_label,
      hjust = 1.05, vjust = 1.25, family = "Arial", size = 6 * pt_to_mm)
  }
  p
}
p_t2d_light <- make_t2d_plot(
  t2d_summary$Light_MetAge, "Light-MetAge (years)",
  "Light-MetAge across T2D parent-trial follow-up",
  corner_label = paste0("Maintenance group \u00d7 visit: ",
                        format_p(s_table3_t2d_global$p_value[2]))
)
p_t2d_delta <- make_t2d_plot(
  t2d_summary$Change, "Change from baseline (years)",
  "Within-person change in Light-MetAge from baseline", zero_line = TRUE
)
p_t2d_accel <- make_t2d_plot(
  t2d_summary$Acceleration, "Light-MetAge acceleration (years)",
  "Light-MetAge acceleration across T2D parent-trial follow-up",
  corner_label = paste0("Maintenance group \u00d7 visit: ",
                        format_p(s_table3_t2d_global$p_value[4]))
)
p_t2d_column <- (p_t2d_light / p_t2d_delta / p_t2d_accel) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Fig5G: adjusted logistic-regression table and forest plot.
make_lr_forest <- function(results) {
  data <- results %>%
    filter(term != "(Intercept)") %>%
    mutate(
      row_id = rev(seq_len(n())),
      or_text = sprintf("%.2f (%.2f-%.2f)", OR, CI_low, CI_high),
      p_text = if_else(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  stripes <- data %>% filter(row_number() %% 2L == 0L)
  y_header <- max(data$row_id) + 0.70
  y_limits <- c(0.45, y_header + 0.25)
  table_panel <- ggplot(data, aes(y = row_id)) +
    geom_rect(data = stripes, aes(ymin = row_id - 0.47,
               ymax = row_id + 0.47), xmin = -Inf, xmax = Inf,
               inherit.aes = FALSE, fill = "#F5F5F5") +
    geom_text(aes(x = 0.01, label = Predictor), hjust = 0,
              size = 7 * pt_to_mm, family = "Arial") +
    geom_text(aes(x = 0.54, label = or_text), hjust = 0.5,
              size = 7 * pt_to_mm, family = "Arial") +
    geom_text(aes(x = 0.96, label = p_text), hjust = 1,
              size = 7 * pt_to_mm, family = "Arial") +
    annotate("text", x = c(0.01, 0.54, 0.96), y = y_header,
             label = c("Predictor", "Adjusted OR (95% CI)", "P value"),
             hjust = c(0, 0.5, 1), fontface = "bold",
             size = 7 * pt_to_mm, family = "Arial") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(plot.margin = margin(2, 4, 2, 4))
  x_min <- max(0.10, min(data$CI_low) * 0.78)
  x_max <- max(1.10, max(data$CI_high) * 1.22)
  x_breaks <- c(0.25, 0.50, 0.75, 1, 1.25, 1.5, 2, 3)
  x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  forest_panel <- ggplot(data, aes(y = row_id, x = OR)) +
    geom_rect(data = stripes, aes(ymin = row_id - 0.47,
               ymax = row_id + 0.47), xmin = -Inf, xmax = Inf,
               inherit.aes = FALSE, fill = "#F5F5F5") +
    geom_vline(xintercept = 1, linetype = "longdash",
               linewidth = 0.5 * pt_to_mm, color = "#777777") +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y",
                  width = 0.16, linewidth = 0.75 * pt_to_mm,
                  color = "#B11235") +
    geom_point(shape = 15, size = 1.2, color = "#B11235") +
    scale_x_log10(
      limits = c(x_min, x_max), breaks = x_breaks,
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    labs(x = "Adjusted odds ratio", y = NULL) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line.y = element_blank(),
      axis.line.x = element_line(linewidth = 0.75 * pt_to_mm),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.5 * pt_to_mm),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 6, color = "black"),
      axis.title.x = element_text(size = 7),
      plot.margin = margin(2, 4, 2, 1)
    )
  (table_panel + forest_panel + plot_layout(widths = c(1.65, 1))) +
    plot_annotation(
      title = "Baseline Light-MetAge and HbA1c <6.5% at 1 year",
      subtitle = "Outcome: attained versus not attained; Light-MetAge per 1 SD",
      theme = theme(
        plot.title = element_text(family = "Arial", size = 8,
                                  face = "plain", hjust = 0.5),
        plot.subtitle = element_text(family = "Arial", size = 6, hjust = 0.5)
      )
    )
}
p_t2d_lr <- make_lr_forest(lr_primary_results)

# Fig5H: complete Baseline/12m pairs, with model-based 12m annotations.
endpoint_data <- t2d_long %>%
  filter(Visit %in% c("Baseline", "12m")) %>%
  mutate(
    Comparison_group = factor(
      if_else(Group %in% c("A", "C"), "A+C", "B+D"),
      levels = c("A+C", "B+D")
    ),
    Endpoint_visit = factor(as.character(Visit), levels = c("Baseline", "12m"))
  )
prepare_pairs <- function(data, outcome) {
  paired <- data %>%
    transmute(Participant_ID, Comparison_group, Endpoint_visit,
              Value = .data[[outcome]]) %>%
    filter(is.finite(Value)) %>%
    group_by(Participant_ID, Comparison_group) %>%
    filter(n() == 2L, n_distinct(Endpoint_visit) == 2L) %>%
    ungroup()
  order <- sort(unique(paired$Participant_ID))
  paired %>%
    mutate(
      group_number = if_else(Comparison_group == "A+C", 1, 2),
      visit_offset = if_else(Endpoint_visit == "Baseline", -0.12, 0.12),
      pair_index = match(Participant_ID, order),
      pair_jitter = (((pair_index * 37L) %% 101L) / 100 - 0.5) * 0.05,
      x_center = group_number + visit_offset,
      x_point = x_center + pair_jitter
    )
}
paired_light <- prepare_pairs(endpoint_data, "Light_MetAge")
paired_accel <- prepare_pairs(endpoint_data, "Light_MetAge_Accel")

summarise_pairs <- function(data, outcome) {
  data %>%
    select(Participant_ID, Comparison_group, Endpoint_visit, Value) %>%
    pivot_wider(names_from = Endpoint_visit, values_from = Value) %>%
    mutate(change = `12m` - Baseline) %>%
    group_by(Comparison_group) %>%
    summarise(
      Outcome = outcome, n_pairs = n(),
      baseline_mean = mean(Baseline), baseline_sd = sd(Baseline),
      baseline_median = median(Baseline),
      baseline_q1 = unname(quantile(Baseline, 0.25)),
      baseline_q3 = unname(quantile(Baseline, 0.75)),
      month12_mean = mean(`12m`), month12_sd = sd(`12m`),
      month12_median = median(`12m`),
      month12_q1 = unname(quantile(`12m`, 0.25)),
      month12_q3 = unname(quantile(`12m`, 0.75)),
      change_mean = mean(change), change_sd = sd(change),
      change_ci_low = change_mean - qt(0.975, n_pairs - 1L) *
        change_sd / sqrt(n_pairs),
      change_ci_high = change_mean + qt(0.975, n_pairs - 1L) *
        change_sd / sqrt(n_pairs),
      .groups = "drop"
    ) %>%
    mutate(Comparison_group = as.character(Comparison_group)) %>%
    relocate(Outcome)
}
s_table4_t2d_paired <- bind_rows(
  summarise_pairs(paired_light, "Light-MetAge"),
  summarise_pairs(paired_accel, "Light-MetAge acceleration")
)

endpoint_p_label <- function(outcome) {
  row <- t2d_metformin_contrasts %>%
    filter(Outcome == outcome, as.character(Visit) == "12m")
  stopifnot(nrow(row) == 1L)
  paste0("Adjusted 12m P = ", format_p(row$p.value),
         "\nHolm-8 P = ", format_p(row$p_holm_8_outcome_visits))
}
make_paired_violin <- function(data, y_label, title_text, p_label) {
  range_y <- range(data$Value)
  span <- diff(range_y)
  bracket_y <- range_y[2] + 0.05 * span
  tip <- 0.025 * span
  ggplot(data, aes(y = Value)) +
    gghalves::geom_half_violin(
      data = data %>% filter(Endpoint_visit == "Baseline"),
      aes(x = x_center, fill = Endpoint_visit,
          group = interaction(Endpoint_visit, Comparison_group)),
      side = "l", position = "identity", width = 0.34,
      trim = TRUE, alpha = 0.4, color = "grey35",
      linewidth = 0.75 * pt_to_mm
    ) +
    gghalves::geom_half_violin(
      data = data %>% filter(Endpoint_visit == "12m"),
      aes(x = x_center, fill = Endpoint_visit,
          group = interaction(Endpoint_visit, Comparison_group)),
      side = "r", position = "identity", width = 0.34,
      trim = TRUE, alpha = 0.4, color = "grey35",
      linewidth = 0.75 * pt_to_mm
    ) +
    geom_line(aes(x = x_point,
                  group = interaction(Comparison_group, Participant_ID)),
              color = "grey65", alpha = 0.3, linewidth = 0.5 * pt_to_mm) +
    geom_boxplot(aes(x = x_center,
                     group = interaction(Endpoint_visit, Comparison_group)),
                 width = 0.15, outlier.shape = NA, fill = "white",
                 alpha = 0.82, color = "grey25", linewidth = 0.75 * pt_to_mm) +
    geom_point(aes(x = x_point, fill = Endpoint_visit), shape = 21,
               size = 0.7, stroke = 0.5 * pt_to_mm,
               color = "grey15", alpha = 0.5) +
    annotate("segment",
      x = c(1.12, 1.12, 2.12), xend = c(1.12, 2.12, 2.12),
      y = c(bracket_y - tip, bracket_y, bracket_y),
      yend = c(bracket_y, bracket_y, bracket_y - tip),
      linewidth = 0.75 * pt_to_mm, color = "black") +
    annotate("text", x = 1.62, y = bracket_y + 0.01 * span,
             label = p_label, hjust = 0.5, vjust = 0,
             family = "Arial", size = 6 * pt_to_mm) +
    scale_fill_manual(name = "Visit",
      values = c("Baseline" = "#BC3C29FF", "12m" = "#0072B5FF")) +
    scale_x_continuous(breaks = c(1, 2),
      labels = c("A+C (Met)", "B+D (no Met)"), limits = c(0.62, 2.38),
      expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.03))) +
    labs(title = title_text, x = NULL, y = y_label) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(shape = 21, size = 2,
                                                  alpha = 0.75))) +
    theme_bw(base_family = "Arial", base_size = 7) +
    theme(
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, face = "plain"),
      axis.ticks = element_line(linewidth = 0.75 * pt_to_mm),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.5 * pt_to_mm,
                                        color = "grey90"),
      panel.border = element_rect(fill = NA, linewidth = 1 * pt_to_mm),
      legend.position = "bottom",
      legend.title = element_text(size = 7, face = "plain"),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "plain"),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}
p_t2d_paired_light <- make_paired_violin(
  paired_light, "Light-MetAge (years)", "Light-MetAge",
  endpoint_p_label("Light-MetAge")
)
p_t2d_paired_accel <- make_paired_violin(
  paired_accel, "Light-MetAge acceleration (years)",
  "Light-MetAge acceleration", endpoint_p_label("Light-MetAge acceleration")
)
p_t2d_paired_row <- (p_t2d_paired_light | p_t2d_paired_accel) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

stopifnot(
  nrow(s_table1_cr) == 48L,
  nrow(s_table2_cr_followups) == 14L,
  nrow(s_table3_t2d_descriptive) == 72L,
  nrow(s_table3_t2d_global) == 4L,
  nrow(s_table4_t2d_paired) == 4L,
  nrow(s_table4_t2d_contrasts) == 8L,
  nrow(s_table5_t2d_logistic) == 5L
)

cat("CALERIE2 randomized:", nrow(cr_wide),
    "; T2D randomized:", nrow(t2d_wide), "\n")
cat("CALERIE2 group x visit P:", cr_global_p(cr_metage_model), "\n")
cat("T2D Light-MetAge group x visit P:",
    s_table3_t2d_global$p_value[2], "\n")
print(s_table2_cr_global, n = Inf, width = Inf)
print(s_table2_cr_followups, n = Inf, width = Inf)
print(s_table3_t2d_global, n = Inf, width = Inf)
print(s_table4_t2d_paired, n = Inf, width = Inf)
print(s_table4_t2d_contrasts, n = Inf, width = Inf)
print(s_table5_t2d_logistic, n = Inf, width = Inf)

if (interactive() && !identical(Sys.getenv("FIG5_SKIP_PLOTS"), "1")) {
  print(p_cr_column)
  print(p_cr_light_column)
  print(p_t2d_column)
  print(p_t2d_lr)
  print(p_t2d_paired_row)
}

# Optional exports. Uncomment only after reviewing the objects in RStudio.
table_output_dir <- file.path("outputs", "Fig5_tables")
figure_output_dir <- file.path("outputs", "Fig5_figures")
# dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)
# writexl::write_xlsx(
#   list(MetAge = s_table1_cr, Light_MetAge = s_table1_cr_light),
#   file.path(table_output_dir, "Table_S1_CALERIE2_descriptive.xlsx")
# )
# writexl::write_xlsx(
#   list(Global_tests = s_table2_cr_global,
#        Followup_contrasts = s_table2_cr_followups),
#   file.path(table_output_dir, "Table_S2_CALERIE2_LMM.xlsx")
# )
# writexl::write_xlsx(
#   list(Descriptive = s_table3_t2d_descriptive,
#        Global_tests = s_table3_t2d_global),
#   file.path(table_output_dir, "Table_S3_T2D_trajectories.xlsx")
# )
# writexl::write_xlsx(
#   list(Paired_descriptive = s_table4_t2d_paired,
#        Model_contrasts = s_table4_t2d_contrasts),
#   file.path(table_output_dir, "Table_S4_T2D_AC_vs_BD.xlsx")
# )
# writexl::write_xlsx(
#   list(Logistic_regression = s_table5_t2d_logistic),
#   file.path(table_output_dir, "Table_S5_T2D_logistic.xlsx")
# )

# dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
# ggsave(file.path(figure_output_dir, "Fig5_CALERIE2.pdf"), p_cr_column,
#        device = grDevices::cairo_pdf, width = 3.6, height = 6.2,
#        units = "in", scale = 1, family = "Arial")
# ggsave(file.path(figure_output_dir, "Fig5_T2D_trajectories.pdf"), p_t2d_column,
#        device = grDevices::cairo_pdf, width = 88, height = 155,
#        units = "mm", scale = 1, family = "Arial")
# ggsave(file.path(figure_output_dir, "Fig5_T2D_logistic.pdf"), p_t2d_lr,
#        device = grDevices::cairo_pdf, width = 3.7, height = 1.65,
#        units = "in", scale = 1, family = "Arial")
# ggsave(file.path(figure_output_dir, "Fig5_T2D_paired.pdf"), p_t2d_paired_row,
#        device = grDevices::cairo_pdf, width = 3.5, height = 2.2,
#        units = "in", scale = 1, family = "Arial")

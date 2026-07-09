# Load Library
library(tidyverse)
library(afex)
library(sjPlot)
library(emmeans)
options(scipen = "999")

# Read the data 
df <- read.csv('demo.csv')
data <- read.csv('data.csv')

# Between group analyses for baseline information 
var.test(Steps ~ PD, data = df) # variance of heterogenity
t.test(Steps ~ PD, data = df, alternative = "two.sided",var.equal = TRUE) # t-test

# Create the Line graph point with the error bar 
df <- read.csv("data_for_graph.csv")
df$Timepoints <- as.factor(df$Timepoints)

# Collapse to AA and TRB 
df <- df %>%
  mutate(
    Group_collapsed = fct_collapse(
      as.factor(Session),
      AA  = c("A", "B"),
      TRB = c("C", "D")
    ),
    # order levels
    Group_collapsed = fct_relevel(Group_collapsed, "AA", "TRB")
  )

# Collapsed to HOA-AA, HOA-TRB, PD-AA, PD-TRB
df_sub <- df %>%
  filter(Group %in% c("HOA", "PD"),
         Group_collapsed %in% c("AA", "TRB")) %>%
  mutate(Panel = factor(paste0(Group, "-", Group_collapsed),
                        levels = c("HOA-AA", "HOA-TRB", "PD-AA", "PD-TRB")))

sum_df <- df_sub %>%
  group_by(Panel, Timepoints) %>%
  summarise(
    mean_vo2 = mean(VO2_mL_kg_min_mean, na.rm = TRUE),
    sd_vo2   = sd(VO2_mL_kg_min_mean, na.rm = TRUE),
    n        = dplyr::n(),
    se_vo2   = sd_vo2 / sqrt(n),
    .groups  = "drop"
  )

# Create the graph
ggplot(
  sum_df,
  aes(x = Timepoints, y = mean_vo2, color = Panel, linetype = Panel, group = Panel)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean_vo2 - se_vo2, ymax = mean_vo2 + se_vo2),
                width = 0.12, linewidth = 0.6) +
  scale_color_manual(values = c(
    "HOA-AA"  = "#1338BE", 
    "HOA-TRB" = "#1338BE",  
    "PD-AA"   = "#EC9706",  
    "PD-TRB"  = "#EC9706"   
  ), name = "Group") +
  scale_linetype_manual(values = c(
    "HOA-AA"  = "longdash",   
    "PD-AA"   = "longdash",
    "HOA-TRB" = "solid",    
    "PD-TRB"  = "solid"
  ), name = "Group") +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 1),
  ) +
  labs(
    x = "Timepoint",
    y = "VO2 (mL/kg/min)",
    title = "VO2 Across Timepoints",
    subtitle = "Same color for HOA & PD within group; AA dotted, TRB solid (mean ± SE)"
  ) +
  theme_classic(base_size = 12)

# Subset to HR or RPE session
#dfs <- subset(df, df$Session == "B" | df$Session == "D") # RPE session
dfs <- subset(df, df$Session == "A" | df$Session == "C") # HR 

# Collapsed to HOA-AA, HOA-TRB, PD-AA, PD-TRB
df_sub <- dfs %>%
  filter(Group %in% c("HOA", "PD"),
         Group_collapsed %in% c("AA", "TRB")) %>%
  mutate(Panel = factor(paste0(Group, "-", Group_collapsed),
                        levels = c("HOA-AA", "HOA-TRB", "PD-AA", "PD-TRB")))

sum_df <- df_sub %>%
  group_by(Panel, Timepoints) %>%
  summarise(
    mean_vo2 = mean(VO2_mL_kg_min_mean, na.rm = TRUE),
    sd_vo2   = sd(VO2_mL_kg_min_mean, na.rm = TRUE),
    n        = dplyr::n(),
    se_vo2   = sd_vo2 / sqrt(n),
    .groups  = "drop"
  )

# Create the graph
ggplot(
  sum_df,
  aes(x = Timepoints, y = mean_vo2, color = Panel, linetype = Panel, group = Panel)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean_vo2 - se_vo2, ymax = mean_vo2 + se_vo2),
                width = 0.12, linewidth = 0.6) +
  scale_color_manual(values = c(
    "HOA-AA"  = "#1338BE", 
    "HOA-TRB" = "#1338BE",  
    "PD-AA"   = "#EC9706",  
    "PD-TRB"  = "#EC9706"   
  ), name = "Group") +
  scale_linetype_manual(values = c(
    "HOA-AA"  = "longdash",   
    "PD-AA"   = "longdash",
    "HOA-TRB" = "solid",    
    "PD-TRB"  = "solid"
  ), name = "Group") +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 1),
  ) +
  ylim(6,16) +
  labs(
    x = "Timepoint",
    y = "VO2 (mL/kg/min)",
    title = "VO2 Across Timepoints",
    subtitle = "Same color for HOA & PD within group; AA dotted, TRB solid (mean ± SE)"
  ) +
  theme_classic(base_size = 12)

# Convert relevant columns to factors
data <- data %>%
  mutate(
    PID = as.factor(PID),
    Group = as.factor(Group),
    Session = as.factor(Session),
    Timepoints = as.factor(Timepoints)
  )

# Mixed-design 
mod <- aov_car(
  FS ~ Group * Session * Timepoints + 
    Error(PID / (Session * Timepoints)),
  data = data,
  factorize = FALSE
)

mod_summary <- mod$anova_table

## Post hoc for Session: Timepoints
# Get estimated marginal means
session_emm <- emmeans(mod, ~ Session | Timepoints)

# Compute sample size per cell from raw data
cell_ns <- data %>%
  group_by(Session, Timepoints) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Session = as.character(Session),
         Timepoints = as.character(Timepoints))

# Create mapping for Timepoints
time_map <- data.frame(Timepoints_emm = c("X1", "X2", "X3"),
                       Timepoints_raw = c("1", "2", "3"))

# Prepare EMM summary and apply mapping
session_emm_summary <- summary(session_emm) %>%
  mutate(Session = as.character(Session),
         Timepoints = as.character(Timepoints)) %>%
  left_join(time_map, by = c("Timepoints" = "Timepoints_emm")) %>%
  mutate(Timepoints = Timepoints_raw) %>%
  select(-Timepoints_raw)

# Join with cell_ns and compute SD
session_emm_summary <- session_emm_summary %>%
  left_join(cell_ns, by = c("Session", "Timepoints")) %>%
  mutate(SD = SE * sqrt(as.numeric(n)))

# Pairwise comparisons
session_pairs <- pairs(session_emm, adjust = "tukey") %>%
  summary(infer = TRUE) %>%
  mutate(Session1 = sub(" - .*", "", contrast),
         Session2 = sub(".* - ", "", contrast),
         Session1 = as.character(Session1),
         Session2 = as.character(Session2),
         Timepoints = as.character(Timepoints)) %>%
  left_join(time_map, by = c("Timepoints" = "Timepoints_emm")) %>%
  mutate(Timepoints = Timepoints_raw) %>%
  select(-Timepoints_raw)

# Merge means and SDs for both groups
session_pairs_with_means <- session_pairs %>%
  left_join(session_emm_summary %>% rename(Session1 = Session), by = c("Session1", "Timepoints")) %>%
  rename(Mean1 = emmean, SD1 = SD) %>%
  left_join(session_emm_summary %>% rename(Session2 = Session), by = c("Session2", "Timepoints")) %>%
  rename(Mean2 = emmean, SD2 = SD)

# Extract only the requested columns
final_df <- session_pairs_with_means %>%
  select(contrast, Timepoints, df, t.ratio, p.value,
         Session1, Session2, Mean1, SD1, Mean2, SD2)


## Post hoc for intervention

# Get EMMs for Group within each Session × Timepoints
group_emm <- emmeans(mod, ~ Group | Session * Timepoints)

# Summarize EMMs
group_emm_summary <- summary(group_emm) %>%
  mutate(Group = as.character(Group),
         Session = as.character(Session),
         Timepoints = as.character(Timepoints))

# Compute sample size per cell from raw data
cell_ns_group <- data %>%
  group_by(Group, Session, Timepoints) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Group = as.character(Group),
         Session = as.character(Session),
         Timepoints = as.character(Timepoints))

# Create mapping for Timepoints
time_map <- data.frame(Timepoints_emm = c("X1", "X2", "X3"),
                       Timepoints_raw = c("1", "2", "3"))

# Apply mapping to EMM summary
group_emm_summary <- group_emm_summary %>%
  left_join(time_map, by = c("Timepoints" = "Timepoints_emm")) %>%
  mutate(Timepoints = Timepoints_raw) %>%
  select(-Timepoints_raw)

# Join and compute SD
group_emm_summary <- group_emm_summary %>%
  left_join(cell_ns_group, by = c("Group", "Session", "Timepoints")) %>%
  mutate(SD = SE * sqrt(as.numeric(n)))

# Pairwise comparisons for Group
group_pairs <- pairs(group_emm, adjust = "tukey") %>%
  summary(infer = TRUE) %>%
  mutate(Group1 = sub(" - .*", "", contrast),
         Group2 = sub(".* - ", "", contrast),
         Group1 = as.character(Group1),
         Group2 = as.character(Group2),
         Session = as.character(Session),
         Timepoints = as.character(Timepoints)) %>%
  left_join(time_map, by = c("Timepoints" = "Timepoints_emm")) %>%
  mutate(Timepoints = Timepoints_raw) %>%
  select(-Timepoints_raw)

# Merge means and SDs for both groups
group_pairs_with_means <- group_pairs %>%
  left_join(group_emm_summary %>% rename(Group1 = Group), by = c("Group1", "Session", "Timepoints")) %>%
  rename(Mean1 = emmean, SD1 = SD) %>%
  left_join(group_emm_summary %>% rename(Group2 = Group), by = c("Group2", "Session", "Timepoints")) %>%
  rename(Mean2 = emmean, SD2 = SD)

# Extract final columns
final_group_df <- group_pairs_with_means %>%
  select(contrast, Session, Timepoints, df, t.ratio, p.value,
         Group1, Group2, Mean1, SD1, Mean2, SD2)

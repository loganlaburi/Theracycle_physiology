# Load Library
library(tidyverse)
library(afex)
library(emmeans)
options(scipen = "999")

# Read the data
data <- read.csv('data.csv')

# Convert relevant columns to factors
data <- data %>%
  mutate(
    PID = as.factor(PID),
    Group = as.factor(Group),
    Session = as.factor(Session),
    Timepoints = as.factor(Timepoints)
  )

# Average sessions A and B to AA
avg_ab <- data %>%
  filter(Session %in% c("A", "B")) %>%
  group_by(PID, Group, Timepoints) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(Session = "AA")


# Average sessions C and D to TRB
avg_cd <- data %>%
  filter(Session %in% c("C", "D")) %>%
  group_by(PID,  Group, Timepoints) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(Session = "RB")

# Average sessions HR
avg_hr <- data %>%
  filter(Session %in% c("A", "C")) %>%
  group_by(PID, Group, Timepoints) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(Session = "HR")

# Average Sessions RPE
avg_rpe <- data %>%
  filter(Session %in% c("B", "D")) %>%
  group_by(PID, Group, Timepoints) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(Session = "RPE")

# Combine the result
# df <- bind_rows(avg_ab, avg_cd)
df <- bind_rows(avg_hr, avg_rpe)

# Mixed-design 
mod <- aov_car(
  SBP_mmHg ~ Group * Session * Timepoints + 
    Error(PID / (Session * Timepoints)),
  data = df,
  factorize = FALSE
)

mod_summary <- mod$anova_table

# ---------------- post hoc for Session: Timepoints
# Get estimated marginal means
session_emm <- emmeans(mod, ~ Session | Timepoints)

# Compute sample size per cell from raw data
cell_ns <- df %>%
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
cell_ns_group <- df %>%
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

## post hoc for Group: session

# Get estimated marginal means
session_emm <- emmeans(mod, ~ Group | Session)

# Summarize EMMs
session_emm_summary <- summary(session_emm) %>%
  mutate(Group = as.character(Group),
         Session = as.character(Session))

# Compute sample size per cell from raw data
cell_ns <- df %>%
  group_by(Group, Session) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Group = as.character(Group),
         Session = as.character(Session))

# Join and compute SD
session_emm_summary <- session_emm_summary %>%
  left_join(cell_ns, by = c("Group", "Session")) %>%
  mutate(SD = SE * sqrt(as.numeric(n)))

# Pairwise comparisons for Group within Session
session_pairs <- pairs(session_emm, adjust = "tukey") %>%
  summary(infer = TRUE) %>%
  mutate(Group1 = sub(" - .*", "", contrast),
         Group2 = sub(".* - ", "", contrast),
         Group1 = as.character(Group1),
         Group2 = as.character(Group2),
         Session = as.character(Session))

# Merge means and SDs for both groups
session_pairs_with_means <- session_pairs %>%
  left_join(session_emm_summary %>% rename(Group1 = Group), by = c("Group1", "Session")) %>%
  rename(Mean1 = emmean, SD1 = SD) %>%
  left_join(session_emm_summary %>% rename(Group2 = Group), by = c("Group2", "Session")) %>%
  rename(Mean2 = emmean, SD2 = SD)

# Extract only requested columns
final_df <- session_pairs_with_means %>%
  select(contrast, Session, df, t.ratio, p.value,
         Group1, Group2, Mean1, SD1, Mean2, SD2)

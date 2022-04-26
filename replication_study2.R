# ---
# title: 'Replication code for Study 2 in _The Affective Style of Politics: 
#    Evidence from Surveys and Laboratory Experiments_'
# author: "Julie Hassing Nielsen & Dan M&oslash;nster"
# date: "April 2022"
# ---

# The groundhog package is used to ensure reproducibility.

packages <- c("tidyr",
              "dplyr",
              "ggplot2",
              "scales",
              "directlabels",
              "table1",
              "psych",
              "lme4",
              "lmerTest",
              "emmeans",
              "kableExtra",
              "effectsize",
              "texreg")
library(groundhog)
groundhog_day <- "2021-11-10"
groundhog.library(packages, groundhog_day)


# Set variables to control figure dimensions
fig.width <- 12
fig.height <- 10

# load experiment data
data <- read.csv("study2.csv", header=TRUE)


# Descriptive statistics reported in the article text
n.sessions <- length(unique(data$session.id))
n.sessions.4.subjs <- sum(table(data$session.id) == 4)
n.sessions.6.subjs <- sum(table(data$session.id) == 6)
age.mean <- round(mean(data$age), 1)
age.sd <- round(sd(data$age), 1)
n.female <- sum(data$sex == "Female")
n.male <- sum(data$sex == "Male")
female.pct <- round(table(data$sex)["Female"] / (sum(table(data$sex))) * 100, 1)
lr.mean <- round(mean(data$left_right, na.rm = TRUE), 1)
lr.sd <- round(sd(data$left_right, na.rm = TRUE), 1)
trust.mean <- round(mean(data$trust_index, na.rm = TRUE), 1)
trust.sd <- round(sd(data$trust_index, na.rm = TRUE), 1)
pol.active.pct <- round(table(data$not_active)["Inactive"] /
                          sum(table(data$not_active)) * 100, 1)
payoff.mean <- round(mean(data$public_goods.1.player.payoff, na.rm = TRUE), 1)
payoff.sd <- round(sd(data$public_goods.1.player.payoff, na.rm = TRUE), 1)

## Table A.11

table_a11 <- 
  table1(~ aff_1 + aff_2 + aff_3 + aff_4 + aff_5 + aff_6 + aff_7 + aff_8  +
           aff_9  + aff_10  + aff_11  + aff_12  + aff_13 + aff_14 + aff_15  +
           aff_16  + aff_17 + aff_18  + aff_19  + aff_20, 
                 data = data,
                 digits = 3)

markdown::markdownToHTML(text = table1::t1kable(table_a11),
                         output = "table_a11.html")


# Reorder treatment factor levels so the order is:
# "Neutral", "Joy", "Disgust", "Fear"
data$treatment <- factor(data$treatment,
                         levels = c("Neutral", "Joy", "Disgust", "Fear"))

# Convert valence into long format
data.valence.long <- gather(data, SAM.valence.round, SAM.valence.value,
                            matches("iaps.([1-6]).player.valence_response"))

data.valence.long$SAM.valence.round <- as.numeric(
  gsub("iaps.([1-6]).player.valence_response",
       "\\1",
       data.valence.long$SAM.valence.round))

data.valence.long$round <- data.valence.long$SAM.valence.round

# Convert to long format for SAM arousal
data.arousal.long <- gather(data, SAM.arousal.round, SAM.arousal.value,
                            matches("iaps.([1-6]).player.arousal_response"))
data.arousal.long$SAM.arousal.round <- as.numeric(
  gsub("iaps.([1-6]).player.arousal_response", 
       "\\1", data.arousal.long$SAM.arousal.round))
data.arousal.long$round <- data.arousal.long$SAM.arousal.round


# Average the valence and arousal over all six rounds.

sam.arousal.by.round <- data.arousal.long %>% 
  select(pid, round, arousal = SAM.arousal.value)

# Drop rows with missing responses
sam.arousal.by.round <-
  sam.arousal.by.round[complete.cases(sam.arousal.by.round), ]

# Alternative
sam.valence.by.round <- data.valence.long %>% 
  select(pid, round, valence = SAM.valence.value)

sam.valence.by.round <-
  sam.valence.by.round[complete.cases(sam.valence.by.round), ]

# Next step: merge data frames
sam.response.by.round <- merge(sam.arousal.by.round, sam.valence.by.round,
                               by=c("pid", "round"))


# Rename left_right to ideology and level_of_educ to education for
# convenience of presentation.

trust_affective_data <- data %>%
  select(pid, Contribution = contribution, treatment,
         Social = soc_trust, Local = trust_loc,
         National = trust_nat, EU = trust_eu, trust_index,
         Tolerating = tolerating.normalized,
         Concealing = concealing.normalized,
         Adjusting = adjusting.normalized,
         age, sex, ideology = left_right, salience)


# Take the average SAM response over all rounds.
sam_average_response <- sam.response.by.round %>% 
  group_by(pid) %>% 
  summarise(valence = mean(valence),
            arousal = mean(arousal),
            .groups = "drop")

sam_average <- merge(sam_average_response, trust_affective_data,
                     by = "pid")


## Table A.13

table_a13_valence <- table1(~ valence | treatment, data = sam_average)
table_a13_arousal <- table1(~ arousal | treatment, data = sam_average)

markdown::markdownToHTML(text = table1::t1kable(table_a13_valence),
                         output = "table_a13_valence.html")
markdown::markdownToHTML(text = table1::t1kable(table_a13_arousal),
                         output = "table_a13_arousal.html")

## Table A.14

image.by.round <- data %>%
  select(pid, treatment, matches('iaps.[1-6].player.iaps_image')) %>%
  gather(iaps.stim, image.number, matches('iaps.[1-6].player.iaps_image'))

# Replace iaps.stim by the number of the stimulus
image.by.round$round <- as.numeric(gsub("iaps.([1-6]).player.iaps_image", "\\1",
                                        image.by.round$iaps.stim))
image.by.round <- image.by.round[complete.cases(image.by.round), ]

sam.response.by.image.round <- merge(sam.response.by.round, image.by.round,
                                     by=c("pid", "round"))
sam.response.by.image <- sam.response.by.image.round %>% 
  group_by(image.number, treatment) %>% 
  summarise( Arousal = mean(arousal),
             Arousal.SD   = sd(arousal),
             Valence = mean(valence),
             Valence.SD   = sd(valence),
             image = unique(image.number))
             
#Add an identifier
sam.response.by.image$Dataset <- "This study"
# Convert image to factor
sam.response.by.image$image <- as.factor(sam.response.by.image$image)

# Keep only the same columns as in the db.ratings loaded below
sam.response.by.image.unscaled <- sam.response.by.image %>% 
  ungroup() %>% 
  select(Image = image, Treatment = treatment,
         Valence, Valence.SD,
         Arousal, Arousal.SD, Dataset)

sam.response.by.image.scaled <- sam.response.by.image.unscaled
sam.response.by.image.scaled$Dataset <- "This study (rescaled)"
# Scale to go from 1-5 Likert to 1 <- 9 Likert scale
sam.response.by.image.scaled$Valence <-
  sam.response.by.image.scaled$Valence * 9 / 5
sam.response.by.image.scaled$Valence.SD <-
  sam.response.by.image.scaled$Valence.SD * 9 / 5
sam.response.by.image.scaled$Arousal <-
  sam.response.by.image.scaled$Arousal * 9 / 5
sam.response.by.image.scaled$Arousal.SD <-
  sam.response.by.image.scaled$Arousal.SD * 9 / 5

# Get the ratings from IAPS.TechManual.1-20.2008.pdf
db.ratings <- read.csv("iaps_db.csv")
db.ratings$Image <- as.factor(db.ratings$Image)
db.ratings$Dataset <- "IAPS TechManual"

# Sort Image factor levels according to treatment
db.ratings$Treatment <- factor(db.ratings$Treatment,
                               levels = c("Neutral", "Joy", "Disgust", "Fear"))
# Print this to test that the result is correct
# db.ratings[order(db.ratings$Treatment), c("Image", "Treatment")]
image.levels <- as.character(db.ratings[order(db.ratings$Treatment),
                                        c("Image", "Treatment")]$Image)


compare.ratings <- rbind(sam.response.by.image.unscaled, db.ratings)

# Order levels of Treatment and Image
compare.ratings$Treatment <- 
  factor(compare.ratings$Treatment,
         levels = c("Neutral", "Joy", "Disgust", "Fear"))
compare.ratings$Image <- 
  factor(compare.ratings$Image, levels = image.levels)

markdown::markdownToHTML(text = kbl(compare.ratings, digits = 2),
                         output = "table_a14.html")


## Figure A.2
compare.ratings <- rbind(sam.response.by.image.scaled, db.ratings)

ggplot(compare.ratings,
       aes(x = Image, y = Arousal, linetype = Dataset, colour = Treatment)) +
  geom_point(position=position_dodge(width = 0.6)) +
  scale_linetype_manual(values=c("solid", "longdash")) +
  geom_errorbar(aes(ymin = Arousal - Arousal.SD,
                    ymax = Arousal + Arousal.SD),
                width=.1,
                position=position_dodge(width = 0.6)) +
  ylim(0, 10) +
  coord_flip() +
  theme_classic()
ggsave("figure_a2.pdf", width = fig.width, height = fig.height, units = "cm")

## Figure A.3
ggplot(compare.ratings,
       aes(x = Image, y = Valence, linetype = Dataset, colour = Treatment)) +
  geom_point(position=position_dodge(width = 0.6)) +
  scale_linetype_manual(values=c("solid", "longdash")) +
  geom_errorbar(aes(ymin = Valence - Valence.SD,
                    ymax = Valence + Valence.SD),
                width=.1,
                position=position_dodge(width = 0.6)) +
  ylim(0, 10) +
  coord_flip() +
  theme_classic()
ggsave("figure_a3.pdf", width = fig.width, height = fig.height, units = "cm")

## Figure A.4

# Summarise SAM valence by treatment for plotting
SAM.valence.by.treatment <- data.valence.long %>%
  select(Treatment = treatment, Round = round, SAM.valence.value)  %>%
  group_by(Treatment, Round) %>% 
  summarise(N = sum(!is.na(SAM.valence.value)),
            mean = mean(SAM.valence.value, na.rm = TRUE),
            sd = sd(SAM.valence.value, na.rm = TRUE),
            se = sd / sqrt(N),
            .groups = "drop")

# Summarise arousal by treatment for plots
SAM.arousal.by.treatment <- data.arousal.long %>%
  select(Treatment = treatment, Round = round, SAM.arousal.value)  %>%
  group_by(Treatment, Round) %>% 
  summarise(N = sum(!is.na(SAM.arousal.value)),
            mean = mean(SAM.arousal.value, na.rm = TRUE),
            sd = sd(SAM.arousal.value, na.rm = TRUE),
            se = sd / sqrt(N),
            .groups = "drop")

SAM.summmary.by.treatment <- SAM.valence.by.treatment
colnames(SAM.summmary.by.treatment) <-
  c("Treatment", "Round", "V.N", "V.mean", "V.sd", "V.se" )
SAM.summmary.by.treatment <- 
  merge(SAM.summmary.by.treatment, SAM.arousal.by.treatment,
        by = c("Treatment", "Round"))
colnames(SAM.summmary.by.treatment) <- c("Treatment", "Round",
                                         "V.N", "V.mean", "V.sd", "V.se" ,
                                         "A.N", "A.mean", "A.sd", "A.se" )
ggplot(data = SAM.summmary.by.treatment,
       aes(x = V.mean, y = A.mean, color = Treatment)) +
  geom_point() +
  geom_errorbar(aes(ymin = A.mean - A.se, ymax = A.mean + A.se)) +
  geom_errorbarh(aes(xmin = V.mean - V.se, xmax = V.mean + V.se)) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(breaks = 1:5) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE,
              clip = "on") +
  xlab("Valence") +
  ylab("Arousal") +
  theme_classic()
ggsave("figure_a4.pdf", width = fig.width, height = fig.height, units = "cm")


## Table A.15
tolerating_items <- c("aff_3", "aff_6", "aff_11", "aff_14", "aff_17")

tolerating_cor <- round(cor(data[ , tolerating_items],
                            use = "pairwise.complete.obs"), 2)

tol_rel <- psych::alpha(data[ , tolerating_items])

tol_note <- paste("Note: Scale reliability coefficient:",
                  round(tol_rel$total$std.alpha, 2))

markdown::markdownToHTML(text = kbl(replace(tolerating_cor,
                                            upper.tri(tolerating_cor),
                                            "")) %>% 
                           add_footnote(tol_note, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a15_tolerating.html")




concealing_items <- c("aff_1", "aff_5", "aff_9", "aff_10", "aff_13",
                       "aff_15", "aff_18", "aff_20")

concealing_cor <- round(cor(data[ , concealing_items],
                               use = "pairwise.complete.obs"), 2)

con_rel <- psych::alpha(data[ , concealing_items])

con_note <- paste("Note: Scale reliability coefficient:",
                     round(con_rel$total$std.alpha, 2))

markdown::markdownToHTML(text = kbl(replace(concealing_cor,
                                            upper.tri(concealing_cor),
                                            "")) %>% 
                           add_footnote(con_note, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a15_concealing.html")


adjusting_items <- c("aff_2", "aff_4", "aff_7", "aff_8",
                     "aff_12", "aff_16", "aff_19")

adjusting_without8 <-  c("aff_2", "aff_4", "aff_7",
                         "aff_12", "aff_16", "aff_19")

adjusting_cor <- round(cor(data[ , adjusting_items],
                            use = "pairwise.complete.obs"), 2)

adj_rel <- psych::alpha(data[ , adjusting_items])
adj_rel_no8 <- psych::alpha(data[ , adjusting_without8])

adj_note <- paste0("Note: Scale reliability coefficient with item 8: ",
                  round(adj_rel$total$std.alpha, 2),
                  "; without item 8: ",
                  round(adj_rel_no8$total$std.alpha, 2))

markdown::markdownToHTML(text = kbl(replace(adjusting_cor,
                                            upper.tri(adjusting_cor),
                                            "")) %>% 
                           add_footnote(adj_note, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a15_adjusting.html")


## Figure A.5

# Perform a factor analysis on the data from the affective styles survey from
# study 2 (lab experiment) to get factor loadings for the affective styles
# dimensions.

# Construct a vector with the column names for the affective style questions
aff_cols <- as.vector(mapply(paste0, rep("aff_", 20), 1:20))

# Extract just the affective styles data from the questionnaire
aff_data <- data[, aff_cols]

# Use only complete cases
aff_data <- aff_data[complete.cases(aff_data), ]


#
# Factor analysis with three factors
#
fit_fa <- factanal(aff_data, 3, rotation="varimax", scores = "Bartlett")

pdf("figure_a5.pdf", width = 7, height = 5) 
scree(fit_fa$correlation, pc = FALSE, main = "")
dev.off()


## Table A.16

# Extract the loadings
fa_loadings <- fit_fa$loadings
fa_loadings <- as.data.frame(fa_loadings[1:20, 1:3])

# Get the loadings with a cutoff of 0.3
fa_loadings_cut <- round(fa_loadings * (abs(fa_loadings) > 0.3), digits = 2) 

# Name the rows and columns
fa_loadings_cut <- fa_loadings_cut %>% 
  rename(Concealing = Factor1,
         Adjusting = Factor2,
         Tolerating = Factor3)

rownames(fa_loadings_cut) <- aff_cols

factor_note <- paste("Factor analysis with varimax rotation with",
                     "absolute loadings > 0.3. Blank spaces indicate",
                     "loadings below 0.3.")

markdown::markdownToHTML(text = kbl(replace(fa_loadings_cut,
                                            fa_loadings_cut == 0, "")) %>% 
                           add_footnote(factor_note, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a16.html")


## Figure 2a: Manipulation check for SAM valence

ggplot(data = SAM.valence.by.treatment,
       aes(x = Round, y = mean,  group = Treatment)) +
  geom_line(color = "grey") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                color = "black", size = 0.5, width = 0.1) +
  geom_dl(aes(label = Treatment), method =  "last.bumpup",
          position = position_nudge(0.1)) +
  expand_limits(x = 6.4) +
  xlab("Stimulus number") +
  ylab("Mean SAM valence") +
  theme_classic() +
  scale_shape_manual(values = c(0, 1, 2, 4)) +
  scale_x_continuous(breaks = 1:6, limits = c(1, 6.8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  theme(legend.position = "none")

ggsave("figure_2a.pdf",
       width = fig.width, height = fig.height, units = "cm")


ggplot(data = SAM.arousal.by.treatment,
       aes(x = Round, y = mean, group = Treatment)) +
  geom_line(color = "grey") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                color = "black", size = 0.5, width = 0.1) +
  geom_dl(aes(label = Treatment), method =  "last.bumpup",
          position = position_nudge(0.1)) +
  expand_limits(x = 6.4) +
  xlab("Stimulus number") +
  ylab("Mean SAM arousal") +
  theme_classic() +
  scale_shape_manual(values = c(0, 1, 2, 4)) +
  scale_x_continuous(breaks = 1:6, limits = c(1, 6.8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  theme(legend.position = "none")

ggsave("figure_2b.pdf",
       width = fig.width, height = fig.height, units = "cm")



## Table 3: Manipulation test for SAM valence and SAM arousal

# A linear mixed effects model takes into account repeated measures
# per subject.

valence.test.lmm  <- lmer(SAM.valence.value ~ treatment  + (1 | round),
                          data = data.valence.long)

arousal.test.lmm  <- lmer(SAM.arousal.value ~ treatment + (1 | round),
                          data = data.arousal.long)

htmlreg(list(valence.test.lmm, arousal.test.lmm),
        digits = 2,
        custom.coef.names = c("(Intercept)", "Joy",
                              "Disgust", "Fear"),
        custom.model.names = c("SAM valence", "SAM arousal"),
        caption = "Manipulation test for SAM valence and SAM arousal (LMM)",
        caption.above = TRUE,
        include.aic = FALSE, include.bic = FALSE,
        include.loglik = FALSE, include.rsquared = TRUE,
        doctype = TRUE, html.tag = TRUE,
        head.tag = TRUE, body.tag = TRUE,
        file = "table_3.html")



## Table 4: pairwise differences for valence and arousal

valence_pairwise <- pairs(emmeans(valence.test.lmm, "treatment"))

arousal_pairwise <- pairs(emmeans(arousal.test.lmm, "treatment"))

v_a_pairwise <- cbind(as.data.frame(valence_pairwise),
                      as.data.frame(arousal_pairwise)[, -1])
# Format the p-values
v_a_pairwise[ , "p.value" == colnames(v_a_pairwise)][
  v_a_pairwise[ ,  "p.value" == colnames(v_a_pairwise)] < 0.0001] <- "< 0.0001"


pairwise_table <- kable(v_a_pairwise, "html", booktabs = TRUE,
                        caption = "Pairwise differences between treatments",
                        label = "man_diff",
                        digits = c(0, 3, 3, 1, 2, 3, 3, 3, 1, 2, 3)) %>% 
  add_header_above(c(" ", "Valence" = 5, "Arousal" = 5))

write(pairwise_table, file = "table_4.html", sep = "\n")



## Table 5

# Models

model_1_beh <- lm(data = sam_average,
              Contribution ~ Tolerating + Concealing + Adjusting)

model_1_beh_std <- standardize(model_1_beh, method = "refit", two_sd = FALSE,
                                     include_response = TRUE)


model_2_beh <- lm(data = sam_average,
              Contribution ~ Tolerating + Concealing + Adjusting +
                valence + arousal)

model_2_beh_std <- standardize(model_2_beh, method = "refit", two_sd = FALSE,
                                     include_response = TRUE)


model_3_beh <- lm(data = sam_average,
              Contribution ~ Tolerating + Concealing + Adjusting +
                valence + arousal + age + sex + ideology + salience)

model_3_beh_std <- standardize(model_3_beh, method = "refit", two_sd = FALSE,
                                     include_response = TRUE)


model_4_soc <- lm(data = sam_average,
              Social ~ Tolerating + Concealing + Adjusting +
                age + sex + ideology + salience)

model_4_soc_std <- standardize(model_4_soc, method = "refit", two_sd = FALSE,
                                     include_response = TRUE)

model_5_pol <- lm(data = sam_average,
              trust_index ~ Tolerating + Concealing + Adjusting +
                + age + sex + ideology + salience)

model_5_pol_std <- standardize(model_5_pol, method = "refit", two_sd = FALSE,
                                     include_response = TRUE)

htmlreg(list(model_1_beh_std,
             model_2_beh_std,
             model_3_beh_std,
             model_4_soc_std,
             model_5_pol_std),
       digits = 2,
       custom.header = list("Behavioral trust" = 1:3,
                            "Social trust" = 4,
                            "Political trust" = 5),
       custom.coef.names = c("(Intercept)", "Tolerating", "Concealing",
                             "Adjusting", "Valence", "Arousal",
                             "Age", "Sex (male)",
                             "Ideology", "Salience"),
       caption = "The impact of affective styles on behavioral trust",
       include.aic = FALSE, include.bic = FALSE,
       include.loglik = FALSE, include.rsquared = TRUE,
       file = "table_5.html")

## Table 6: The impact of affective styles on SAM valence and arousal

# This tests whether affective styles moderate the effect of the IAPS stimuli
# on the self-reported SAM valence and arousal.

# Valence

valence_simple <- lmer(data = data.valence.long,
                       SAM.valence.value ~ tolerating.normalized +
                         concealing.normalized +
                         adjusting.normalized +
                         treatment + (1 | round))

valence_simple_std <- standardize(valence_simple,
                                  method = "refit", two_sd = FALSE,
                                  include_response = TRUE)

valence_moderation <- lmer(data = data.valence.long,
                           SAM.valence.value ~ (tolerating.normalized +
                                                  concealing.normalized +
                                                  adjusting.normalized) *
                             treatment + (1 | round))

valence_moderation_std <- standardize(valence_moderation,
                                      method = "refit", two_sd = FALSE,
                                      include_response = TRUE)



# Arousal
arousal_simple <- lmer(data = data.arousal.long,
                       SAM.arousal.value ~ tolerating.normalized +
                         concealing.normalized +
                         adjusting.normalized +
                         treatment + (1 | round))

arousal_simple_std <- standardize(arousal_simple,
                                  method = "refit", two_sd = FALSE,
                                  include_response = TRUE)

arousal_moderation <- lmer(data = data.arousal.long,
                           SAM.arousal.value ~ (tolerating.normalized +
                                                  concealing.normalized +
                                                  adjusting.normalized) *
                             treatment + (1 | round))

arousal_moderation_std <- standardize(arousal_moderation,
                                      method = "refit", two_sd = FALSE,
                                      include_response = TRUE)

htmlreg(list(valence_simple_std, valence_moderation_std,
             arousal_simple_std, arousal_moderation_std),
        digits = 2,
        custom.header = list("Valence" = 1:2, "Arousal" = 3:4),
        custom.coef.names = c("(Intercept)", "Tolerating", "Concealing",
                              "Adjusting",
                              "Joy", "Disgust", "Fear",
                              "Tolerating:Joy", "Tolerating:Disgust",
                              "Tolerating:Fear",
                              "Concealing:Joy", "Concealing:Disgust",
                              "Concealing:Fear",
                              "Adjusting:Joy", "Adjusting:Disgust",
                              "Adjusting:Fear "),
        caption = "The impact of affective styles on SAM valence and arousal",
        include.aic = FALSE, include.bic = FALSE,
        include.loglik = FALSE, include.rsquared = TRUE,
        file = "table_6.html")



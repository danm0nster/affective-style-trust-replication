# ---
# title: 'Replication code for Study 1 in _The Affective Style of Politics:
#   Evidence from Surveys and Laboratory Experiments_'
# author: "Julie Hassing Nielsen & Dan M&oslash;nster"
# date: "January 2022"
# ---

# The groundhog package is used to ensure reproducibility, by loading
# the version of all R-packages used originally.
library(groundhog)
groundhog_day <- "2021-11-10"
packages <- c("dplyr",
              "table1",
              "markdown",
              "psych",
              "kableExtra",
              "lme4",
              "lmerTest",
              "texreg",
              "effectsize")
groundhog.library(packages, groundhog_day)

# load survey data
DK <- read.csv("study1_dk.csv")
US <- read.csv("study1_us.csv")

# Affective styles variables are created. Because there are different numbers of
# items in each of the affective styles and to enable comparison, we normalize
# the three affective styles, so that they range from 0-100. For adjusting we
# use the normalized adjusting without item 8 unless otherwise stated.

# Function to construct affective style variables.
construct_affective_styles <- function(data) {
  af_columns <- as.vector(mapply(paste0, rep("af", 20), 1:20))
  for (af_var in af_columns) {
    if (!(af_var %in% colnames(data))) {
      stop("Some columns are missing! Wrong name?")
    }
  }
  data$tolerating <- data$af3 + data$af6 + data$af11 + data$af14 + data$af17
  data$tolerating_normalized <- (data$tolerating - 5) / (5*5 - 5) * 100
  data$concealing <- data$af1 + data$af5 + data$af9 + data$af10 + 
    data$af13 + data$af15 + data$af18 + data$af20
  data$concealing_normalized <- (data$concealing - 8) / (8*5 - 8) * 100
  # item 8 is excluded from adjusting
  data$adjusting_without8 <- data$af2 + data$af4 + data$af7 + 
    data$af12 + data$af16 + data$af19
  data$adjusting_without8_normalized <- 
    (data$adjusting_without8 - 6) / (6*5 - 6) * 100
  # For comparison, create a variable where item 8 is included in adjusting
  data$adjusting <- data$af2 + data$af4 + data$af7 +
    data$af8 + data$af12 + data$af16 + data$af19
  data$adjusting_normalized <- (data$adjusting - 7) / (7*5 - 7) * 100
  return(data)
}

DK <- construct_affective_styles(DK)
US <- construct_affective_styles(US)


# Descriptive statistics

## Online Appendix: Table A1 

#Descriptive statistics of included variables in the study, excluding dependent
#variables, affective style and big five personality traits : Denmark and the US

# As two separate tables
DK_descriptive <-
  table1::table1(~age + education + ideology + sex + salience + region + urban,
               data = DK)
markdown::markdownToHTML(text = table1::t1kable(DK_descriptive),
                         output = "descriptive_dk.html")


US_descriptive <-
  table1::table1(~age + education + ideology + sex + salience + region + urban,
               data = US)
markdown::markdownToHTML(text = table1::t1kable(US_descriptive),
                         output = "descriptive_us.html")


## Online Appendix: Table A2

# Descriptive statistics of dependent variables: political trust index and
# social trust: Denmark and the US.

# As two separate tables
trust_table_data <- rbind(DK %>% select(soc_trust, trust_index) %>% 
                            mutate(Social_trust = soc_trust,
                                   Political_trust_index = trust_index,
                                   Country = "DK"),
                          US %>% select(soc_trust, trust_index) %>% 
                            mutate(Social_trust = soc_trust,
                                   Political_trust_index = trust_index,
                                   Country = "US"))

table_a2_trust <- 
  table1::table1(~Social_trust + Political_trust_index | Country, 
                 data = trust_table_data,
                 overall = FALSE)

markdown::markdownToHTML(text = table1::t1kable(table_a2_trust),
                         output = "table_a2.html")

## Online Appendix: Table A3

# Correlation matrix between political trust index variables: Denmark and the
# US.

# Select variables for the data frame alpha test;
# do the analysis only with the selected variables.
trust_items_dk <- c("trust_eu", "trust_nat", "trust_loc", "soc_trust")
trust_items_us <- c("trust_fed", "trust_nat", "trust_loc", "soc_trust")

cor_dk <- round(cor(DK[ , trust_items_dk], use = "pairwise.complete.obs"), 2)
cor_us <- round(cor(US[ , trust_items_us], use = "pairwise.complete.obs"), 2)

a_dk <- psych::alpha(DK[ , trust_items_dk])
a_us <- psych::alpha(US[ , trust_items_us])

cor_note_dk <- paste("Note: Scale reliability coefficient:",
                     round(a_dk$total$raw_alpha, 2))
cor_note_us <- paste("Note: Scale reliability coefficient:",
                     round(a_us$total$raw_alpha, 2))

markdown::markdownToHTML(text = kbl(replace(cor_dk,
                                            upper.tri(cor_dk),
                                            "")) %>% 
                           add_footnote(cor_note_dk, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a3_dk.html")
markdown::markdownToHTML(text = kbl(replace(cor_us,
                                            upper.tri(cor_us),
                                            "")) %>%
                           add_footnote(cor_note_us,notation = "none") %>% 
                           kable_styling(),
                         output = "table_a3_us.html")

## Online Appendix: Table A4

# Descriptive statistics of personality (B5) survey items: Denmark and the US.

B5_table_data <- rbind(DK %>% select(extraverted, reserved, critical, 
                                     sympathetic, disorganized, dependable, 
                                     anxious, calm, conventional, open) %>% 
                         mutate(Country = "DK"),
                       US %>% select(extraverted, reserved, critical, 
                                     sympathetic, disorganized, dependable, 
                                     anxious, calm, conventional, open) %>% 
                         mutate(Country = "US"))

B5_table <- table1::table1(~extraverted + reserved + critical +
                             sympathetic + disorganized + dependable +
                             anxious + calm + conventional + open | Country,
                           data = B5_table_data,
                           overall = FALSE)

markdown::markdownToHTML(text = table1::t1kable(B5_table),
                         output = "table_a4.html")


## Online Appendix: Table A5

# Correlation matrix between the Big Five Personality Traits (B5): Denmark and
# US.

cor_B5_dk <- round(cor(DK[ , c("extraversion", "agreeableness",
                               "conscientiousness",
                               "neuroticism", "openness")],
                       use = "pairwise.complete.obs"),
                   2)

cor_B5_us <- round(cor(US[ , c("extraversion", "agreeableness",
                               "conscientiousness",
                               "neuroticism", "openness")],
                       use = "pairwise.complete.obs"),
                   2)

markdown::markdownToHTML(text = kbl(replace(cor_B5_dk,
                                            upper.tri(cor_B5_dk),
                                            "")) %>% 
                           kable_styling(),
                         output = "table_a5_dk.html")
markdown::markdownToHTML(text = kbl(replace(cor_B5_us,
                                            upper.tri(cor_B5_us),
                                            "")) %>%
                           kable_styling(),
                         output = "table_a5_us.html")


## Online Appendix: Table A6

# Correlation matrix between affective style items: Tolerating.

tolerating_items <- c("af3", "af6", "af11", "af14", "af17")

tolerating_cor_dk <- round(cor(DK[ , tolerating_items],
                               use = "pairwise.complete.obs"), 2)
tolerating_cor_us <- round(cor(US[ , tolerating_items],
                               use = "pairwise.complete.obs"), 2)

tol_rel_dk <- psych::alpha(DK[ , tolerating_items])
tol_rel_us <- psych::alpha(US[ , tolerating_items])

tol_note_dk <- paste("Note: Scale reliability coefficient:",
                     round(tol_rel_dk$total$std.alpha, 2))
tol_note_us <- paste("Note: Scale reliability coefficient:",
                     round(tol_rel_us$total$std.alpha, 2))

markdown::markdownToHTML(text = kbl(replace(tolerating_cor_dk,
                                            upper.tri(tolerating_cor_dk),
                                            "")) %>% 
                           add_footnote(tol_note_dk, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a6_dk.html")

markdown::markdownToHTML(text = kbl(replace(tolerating_cor_us,
                                            upper.tri(tolerating_cor_us),
                                            "")) %>%
                           add_footnote(tol_note_us, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a6_us.html")


## Online Appendix: Table A7

# Correlation matrix between affective style items: Concealing.

concealing_styles <- c("af1", "af5", "af9", "af10",
                       "af13", "af15", "af18", "af20")

concealing_cor_dk <- round(cor(DK[ , concealing_styles],
                               use = "pairwise.complete.obs"), 2)
concealing_cor_us <- round(cor(US[ , concealing_styles],
                               use = "pairwise.complete.obs"), 2)

con_rel_dk <- psych::alpha(DK[ , concealing_styles])
con_rel_us <- psych::alpha(US[ , concealing_styles])

con_note_dk <- paste("Note: Scale reliability coefficient:",
                     round(con_rel_dk$total$std.alpha, 2))
con_note_us <- paste("Note: Scale reliability coefficient:",
                     round(con_rel_us$total$std.alpha, 2))

markdown::markdownToHTML(text = kbl(replace(concealing_cor_dk,
                                            upper.tri(concealing_cor_dk),
                                            "")) %>% 
                           add_footnote(con_note_dk, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a7_dk.html")
markdown::markdownToHTML(text = kbl(replace(concealing_cor_us,
                                            upper.tri(concealing_cor_us),
                                            "")) %>%
                           add_footnote(con_note_us, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a7_us.html")




## Online Appendix: Table A8

# Correlation matrix between affective style items: Adjusting.

# The analysis is conducted without item 8 as a part of the affective style:
# Adjusting, because item 8 performs poorly in the correlation between the
# affective style items. Hence, below we make both an adjusting variable, which
# includes item 8 `adjusting_styles` and one that excludes item 8
# `adjusting_without8`.

adjusting_styles <- c("af2", "af4", "af7", "af8", "af12", "af16", "af19")
adjusting_without8 <- c("af2", "af4", "af7", "af12", "af16", "af19")

adjusting_cor_dk <- round(cor(DK[ , adjusting_styles],
                               use = "pairwise.complete.obs"), 2)
adjusting_cor_us <- round(cor(US[ , adjusting_styles],
                               use = "pairwise.complete.obs"), 2)

adj_rel_dk <- psych::alpha(DK[ , adjusting_styles])
adj_rel_us <- psych::alpha(US[ , adjusting_styles])

adj_no8_rel_dk <- psych::alpha(DK[ , adjusting_without8])
adj_no8_rel_us <- psych::alpha(US[ , adjusting_without8])

adj_note_dk <- paste("Note: Scale reliability coefficient in DK with item 8:",
                     round(adj_rel_dk$total$std.alpha, 2),
                     "; without item 8:",
                     round(adj_no8_rel_dk$total$std.alpha, 2))
adj_note_us <- paste("Note: Scale reliability coefficient in US with item 8:",
                     round(adj_rel_us$total$std.alpha, 2),
                     "; without item 8:",
                     round(adj_no8_rel_us$total$std.alpha, 2))

markdown::markdownToHTML(text =  kbl(replace(adjusting_cor_dk,
                                             upper.tri(adjusting_cor_dk),
                                             "")) %>% 
                           add_footnote(adj_note_dk, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a8_dk.html")
markdown::markdownToHTML(text = kbl(replace(adjusting_cor_us,
                                            upper.tri(adjusting_cor_us),
                                            "")) %>%
                           add_footnote(adj_note_us, notation = "none") %>% 
                           kable_styling(),
                         output = "table_a8_us.html")



## Online Appendix: Table A9

# Results of factor analysis with affective style items: Denmark and the US. 

# Exploratory factor analysis (Denmark and US) 

factor_analysis <- function(country_data) {
  # Construct a vector with the column names for the affective style questions
  aff_cols <- as.vector(mapply(paste0, rep("af", 20), 1:20))
  # Extract just the affective styles data from the questionnaire
  aff_data <- country_data[, aff_cols]
  # Use only complete cases
  aff_data <- aff_data[complete.cases(aff_data), ]
  #
  # Factor analysis with three factors
  #
  fit_fa <- factanal(aff_data, 3, rotation = "varimax", scores = "Bartlett")
  # Extract the loadings
  fa_loadings <- fit_fa$loadings
  fa_loadings <- as.data.frame(fa_loadings[1:20, 1:3])
  # Get the loadings with a cutoff of 0.3
  fa_loadings <- round(fa_loadings * (abs(fa_loadings) > 0.3), digits = 2) 
  # Name the rows and columns
  # Concealing has the highest loading on af1
  con_col <- which(fa_loadings[1, ] == max(fa_loadings[1, ]))
  # Tolerating has the highest loading on af6
  tol_col <- which(fa_loadings[6, ] == max(fa_loadings[6, ]))
  # Adjusting has the highest loading on af19
  adj_col <- which(fa_loadings[19, ] == max(fa_loadings[19, ]))
  colnames(fa_loadings)[con_col] <- "Concealing"
  colnames(fa_loadings)[tol_col] <- "Tolerating"
  colnames(fa_loadings)[adj_col] <- "Adjusting"
  rownames(fa_loadings) <- aff_cols
  # Make sure factors are in the same order every time
  fa_loadings <- fa_loadings[, c("Concealing", "Adjusting", "Tolerating")]
  return(fa_loadings)
}

factor_note <- paste("Factor analysis with varimax rotation with",
                     "absolute loadings > 0.3. Blank spaces indicate",
                     "loadings below 0.3.")

DK_fa_loadings <- factor_analysis(DK)

DK_fa_table <- kbl(replace(DK_fa_loadings,
                            DK_fa_loadings == 0,
                            ""),
                    format = "html", row.names = TRUE) %>% 
  add_footnote(factor_note, notation = "none")

US_fa_loadings <- factor_analysis(US)

US_fa_table <- kbl(replace(US_fa_loadings,
                           US_fa_loadings == 0,
                           ""),
                   format = "html", row.names = TRUE) %>% 
  add_footnote(factor_note, notation = "none")


markdown::markdownToHTML(text = DK_fa_table,
                         output = "table_a9_dk.html")

markdown::markdownToHTML(text = US_fa_table,
                         output = "table_a9_us.html")



## Table 1: Affective style items

# Construct a vector with the column names for the affective style questions

af_columns <- as.vector(mapply(paste0, rep("af", 20), 1:20))
af_descriptions <- 
  c("1: People usually can’t tell how I am feeling inside",
    "2: I have my emotions well under control",
    "3: I can tolerate having strong emotions",
    "4: I can avoid getting upset by taking a different perspective on things",
    "5: I often suppress my emotional reactions to things",
    "6: It’s ok if people see me being upset",
    "7: I can calm down very quickly",
    "8: I am able to let go of my feelings",
    "9: I am good at hiding my feelings",
    "10: People usually cannot tell when I am upset",
    "11: It’s ok to feel negative emotions at times",
    "12: I know exactly what to do to get myself into a better mood",
    "13: People usually cannot tell when I am sad",
    "14: I can tolerate being upset",
    "15: I can act in a way such that people do not see me being upset",
    "16: I can get into a better mood quite easily",
    "17: There is nothing wrong with feeling very emotional",
    "18: I could easily fake emotions",
    "19: I can get out of a bad mood very quickly",
    "20: I can hide my anger well if I have to")
tolerating_items <- c(3, 6, 11, 14, 17)
concealing_items <- c(1, 5, 9, 10, 13, 15, 18, 20)
adjusting_items <- c(2, 4, 7, 8, 12, 16, 19)

item_descriptions <- data.frame(data_col = af_columns, 
                                description = af_descriptions)

affective <- DK[ , af_columns]
affective$country <- "DK"
affective_US <- US[ , af_columns]
affective_US$country <- "US"
affective <- rbind(affective, affective_US)
rm(affective_US)


for (n in 1:nrow(item_descriptions)) {
  table1::label(affective[[as.character(item_descriptions$data_col[n])]]) <- 
    item_descriptions$description[n]
}

table1::label(affective$af1) <- af_descriptions[1]
table1::label(affective$af2) <- af_descriptions[2]

# Generate the variable in the right order for the table1() call below
table_1 <- 
  table1::table1(~ af3 + af6 + af11 + af14 + af17 +
                   af1  + af5  + af9  + af10  + af13  + af15  + af18  + af20 +
                   af2  + af4  + af7  + af8  + af12  + af16  + af19 | country, 
                 data = affective,
                 overall = FALSE,
                 digits = 2)

markdown::markdownToHTML(text = table1::t1kable(table_1),
                         output = "table_1.html")


## Table 2

# Linear mixed effects model: DK, Model 1-4

# Model 1 where region, urban and education are included as random effects
# with a random intercept
lmm_trust_1_dk <- lmer(soc_trust ~ ideology + salience + age_mean + sex +
                         extraversion + agreeableness + conscientiousness +
                         neuroticism + openness +
                         (1 | region) + (1 | urban) + (1 | education),
                       data = DK)

# Model 2 where region, urban, and education are included as random effects
# with a random intercept 
lmm_trust_2_dk <- lmer(soc_trust ~ ideology + salience + age_mean + sex +
                         extraversion + agreeableness + conscientiousness +
                         neuroticism + openness + concealing_normalized +
                         tolerating_normalized + adjusting_without8_normalized +
                         (1 | region) + (1 | urban) + (1 | education),
                       data = DK)

# Model 3 where region, urban, and education are included as random effects
# with a random intercept 
lmm_trust_3_dk <- lmer(trust_index ~ ideology + salience + age_mean + sex +
                         extraversion + agreeableness + conscientiousness +
                         neuroticism + openness +
                         (1 | region) + (1 | urban) + (1 | education),
                       data = DK)

# Model 4 where region, urban, and education are included as random effects
# with a random intercept 
lmm_trust_4_dk <- lmer(trust_index ~ ideology + salience + age_mean + sex +
                         extraversion + agreeableness + conscientiousness +
                         neuroticism + openness + concealing_normalized +
                         tolerating_normalized + adjusting_without8_normalized +
                         (1 | region) + (1 | urban) + (1 | education),
                       data = DK)


# Linear mixed effects model: US Model 5-8

# Model 1 where region, urban and education are included as random effects
# with a random intercept
lmm_trust_1_us <- lmer(soc_trust ~ ideology + salience + age_mean + sex +
                         extraversion + agreeableness + conscientiousness +
                         neuroticism + openness +
                         (1 | region) + (1 | urban) + (1 | education),
                       data = US)


# Model 2 where region, urban, and education are included as random effects
# with a random intercept 
lmm_trust_2_us <- lmer(soc_trust ~ ideology + salience + age_mean + sex +
                         extraversion + agreeableness + conscientiousness +
                         neuroticism + openness + concealing_normalized +
                         tolerating_normalized + adjusting_without8_normalized +
                         (1 | region) + (1 | urban) + (1 | education),
                       data = US)


# Model 3 where region, urban, and education are included as random effects
# with a random intercept 
lmm_trust_3_us <- lmer(trust_index ~ ideology + salience + age_mean + sex +
                         extraversion + agreeableness + conscientiousness +
                         neuroticism + openness +
                         (1 | region) + (1 | urban) + (1 | education),
                       data = US)


# Model 4 where region, urban, and education are included as random effects
# with a random intercept 
lmm_trust_4_us <- lmer(trust_index ~ ideology + salience + age_mean + sex +
                         extraversion + agreeableness + conscientiousness +
                         neuroticism + openness + concealing_normalized +
                         tolerating_normalized + adjusting_without8_normalized +
                         (1 | region) + (1 | urban) + (1 | education),
                       data = US)

# standardize models, using the effectsize package
lmm_trust_1_dk_std <- standardize(lmm_trust_1_dk, method = "refit",
                                  include_response = TRUE)
lmm_trust_2_dk_std <- standardize(lmm_trust_2_dk, method = "refit",
                                  include_response = TRUE)
lmm_trust_3_dk_std <- standardize(lmm_trust_3_dk, method = "refit",
                                  include_response = TRUE)
lmm_trust_4_dk_std <- standardize(lmm_trust_4_dk, method = "refit",
                                  include_response = TRUE)
lmm_trust_1_us_std <- standardize(lmm_trust_1_us, method = "refit",
                                  include_response = TRUE)
lmm_trust_2_us_std <- standardize(lmm_trust_2_us, method = "refit",
                                  include_response = TRUE)
lmm_trust_3_us_std <- standardize(lmm_trust_3_us, method = "refit",
                                  include_response = TRUE)
lmm_trust_4_us_std <- standardize(lmm_trust_4_us, method = "refit",
                                  include_response = TRUE)


trust_models <- list(lmm_trust_1_dk_std, lmm_trust_2_dk_std,
                     lmm_trust_3_dk_std, lmm_trust_4_dk_std,
                     lmm_trust_1_us_std, lmm_trust_2_us_std,
                     lmm_trust_3_us_std, lmm_trust_4_us_std)

# Compute conditional R-squared, to be included in table
trs <- list()
n <- 1
for (m in trust_models) {
  # Extract results from model
  tr <- texreg::extract(m, include.aic = FALSE, include.bic = FALSE,)
  # Compute the R-squared of the model
  rsq <- MuMIn::r.squaredGLMM(m)
  # Add the R-squared to the GOF-part
  tr@gof.names <- c(tr@gof.names, "Conditional R2")
  tr@gof <- c(tr@gof, rsq[, "R2c"])
  tr@gof.decimal <- c(tr@gof.decimal, TRUE)
  trs[n] <- tr
  n <- n + 1
}

htmlreg(trs,
        custom.header = list("Denmark" = 1:4, "US" = 5:8),
        custom.model.names = c("Social trust", "Social trust",
                               "Political trust", "Political trust",
                               "Social trust", "Social trust",
                               "Political trust", "Political trust"),
        custom.coef.names = c("(Intercept)", "Ideology", "Salience", "Age",
                              "Sex (Male)", "Extraversion", "Agreeableness",
                              "Conscientiousness", "Neuroticism", "Openness",
                              "Concealing", "Tolerating", "Adjusting"),
        digits = 2,
        caption.above = TRUE,
        caption = "Determinants of Political Trust",
        doctype = TRUE,
        html.tag = TRUE,
        head.tag = TRUE,
        body.tag = TRUE,
        file = "table_2.html")

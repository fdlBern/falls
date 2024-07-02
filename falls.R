################################################################################
######### Risk factors for in-hospital falls among medical patients ############
################### A systematic review and meta-analysis ######################
################################################################################
# Heinzmann J, Rossen ML et al.

library(tidyverse)
library(meta)

img_save <- FALSE # Should figures be saved?

# Read dataextraction file from SPSS -------------------------------------------
df <- read_rds("df.RDS")

# Study characteristics --------------------------------------------------------
tab_char_all <- df %>%
  filter(is.na(model) | model  == "*") %>%
  mutate(id = floor(id)) %>%
  distinct(id, .keep_all = TRUE)

tab_char_all %>%  
  arrange(year) %>%
  select(c(model, author_year, country, design, cases_n, control_n, hospsize, ward)) %>%
  rename(c("Author, publication year" = 'author_year',
           "Country" = 'country',
           "Study design" = 'design',
           "Hospital size" = 'hospsize',
           "Sample size (cases)" = 'cases_n',
           "Sample size (controls)" = 'control_n',
           "Ward type" = 'ward'))

tab_char_all %>%
  filter(falls == "yes") %>%
  filter(is.na(model)| model  == "*") %>%
  count(design) %>%
  rename("Study design" = 'design') %>%
  arrange(desc(n))

tab_char_all %>%
  filter(falls == "yes") %>%
  filter(is.na(model)| model  == "*") %>%
  arrange(year) %>%
  select(author_year, QUIPS1, QUIPS2, QUIPS3, QUIPS4, QUIPS5, QUIPS6, quips) %>%
  rename("Author, publication year" = author_year) %>% 
  rename("Item 1" = QUIPS1) %>%
  rename("Item 2" = QUIPS2) %>%
  rename("Item 3" = QUIPS3) %>%
  rename("Item 4" = QUIPS4) %>%
  rename("Item 5" = QUIPS5) %>%
  rename("Item 6" = QUIPS6) %>%
  rename("Conclusion" = quips)

tab_char_all %>%
  filter(falls == "yes") %>%
  filter(is.na(model)| model  == "*") %>%
  count(quips) %>%
  rename("Overall QUIPS appraisal" = quips)

tab_char_all %>%
  group_by(quips, falls, injfalls) %>%
  summarise(n = n())

# Filter and calculate s.e. ----------------------------------------------------
# outcome = falls, effect emasure = OR
data_orf <- df %>%
  filter(outcome == "Falls" & !is.na(oddsratio)) # for main analysis 
  #filter(outcome == "Falls" & excl_age == "60 to 75 years or less" & !is.na(oddsratio)) # sensitivity analysis
  #filter(outcome == "Falls" & excl_age == "No" & !is.na(oddsratio)) # sensitivity analysis

data_orf$se <- ifelse(is.na(data_orf$se) & data_orf$oddsratio > 0,
                      ((log(data_orf$ci_up) - log(data_orf$ci_low))/3.92),
                      data_orf$se)

data_orf <- data_orf %>%
  filter(se > 0) %>%
  filter(ci_low > 0)


# Categorisation and count tables ----------------------------------------------
# All studies
ev_subgroup <- data_orf %>%
  group_by(rf) %>%
  summarise(n = n()) %>%
  filter(n > 4) 

# Studies with adjusted effect measures
ev_subgroup_adj <- data_orf %>%
  filter(adj_age == "yes - model" | adj_age == "matched controls") %>%
  filter(adj_sex == "yes - model" | adj_sex == "matched controls") %>%
  group_by(rf) %>%
  summarise(n = n()) 

ev_subgroup <- left_join(ev_subgroup, ev_subgroup_adj, 
                         by = "rf",
                         suffix = c("_all", "_adj"))

# Studies excl. QUIPS with high-risk of bias
ev_subgroup_Q <- data_orf %>%
  filter(quips == "low or moderate risk") %>%
  group_by(rf) %>%
  summarise(n = n())

ev_subgroup <- left_join(ev_subgroup, ev_subgroup_Q, 
                         by = "rf")

#. Studies with adjusted effect measures excl. QUIPS with high-risk of bias
ev_subgroup_Qadj <- data_orf %>%
  filter(adj_age == "yes - model" | adj_age == "matched controls") %>%
  filter(adj_sex == "yes - model" | adj_sex == "matched controls") %>%
  filter(quips == "low or moderate risk") %>%
  group_by(rf) %>%
  summarise(n = n())

ev_subgroup <- left_join(ev_subgroup, ev_subgroup_Qadj, 
                         by = "rf",
                         suffix = c("_Q", "_Qadj"))



# Meta-analysis (falls, or) ====================================================

## Filter / define data sets ---------------------------------------------------

data_all <- data_orf

data_adj <- data_orf %>%
  filter(adj_age == "yes - model" | adj_age == "matched controls") %>%
  filter(adj_sex == "yes - model" | adj_sex == "matched controls")

data_q <- data_orf %>%
  filter(quips == "low or moderate risk")

data_qadj <- data_orf %>%
  filter(quips == "low or moderate risk") %>%
  filter(adj_age == "yes - model" | adj_age == "matched controls") %>%
  filter(adj_sex == "yes - model" | adj_sex == "matched controls")


## create a tribble ------------------------------------------------------------
ma_chart <- tribble(
  ~ev_name, ~quips, ~n, 
  ~effect, ~lower_or, ~upper_or, 
  ~I2, ~lower_I2, ~upper_I2, 
  ~lower_predict, ~upper_predict)

# define loop vectors
ev_subgroup <- ev_subgroup %>%
  #mutate_all(~ replace(., is.na(.), 0)) %>% 
  filter(!is.na(rf))

ev_all <- ev_subgroup$rf[ev_subgroup$n_all >= 5]
ev_adj <- ev_subgroup$rf[ev_subgroup$n_adj >= 5]
ev_q <- ev_subgroup$rf[ev_subgroup$n_Q >= 5]
ev_qadj <- ev_subgroup$rf[ev_subgroup$n_Qadj >= 5]


## Calculate for unadjusted ORs ------------------------------------------------
ma_chart_all <- ma_chart

for (i in ev_all) {
  metaanalyse <- data_all %>%
    filter(rf == i) %>%
    metagen(TE = log(oddsratio),
            seTE = se,
            studlab = ref,
            data = ., 
            sm = "OR",
            fixed = FALSE,
            random = TRUE,
            method.tau = "REML",
            hakn = TRUE,
            title = "Falls risk factors",
            subgroup = quips,
            tau.common = TRUE, 
            prediction = TRUE,
            prediction.subgroup = TRUE,
            control=list(stepadj=0.5, maxiter=10000))
  
  # Fill tibble with results
  ma_chart_all <- ma_chart_all %>%
    add_row(ev_name = i,
            quips = "pooled",
            n = metaanalyse[["k.study"]],
            effect = exp(metaanalyse[["TE.random"]]),
            lower_or = exp(metaanalyse[["lower.random"]]),
            upper_or = exp(metaanalyse[["upper.random"]]),
            I2 = metaanalyse[["I2"]],
            lower_I2 = metaanalyse[["lower.I2"]],
            upper_I2 = metaanalyse[["upper.I2"]],
            lower_predict = exp(metaanalyse[["lower.predict"]]),
            upper_predict = exp(metaanalyse[["upper.predict"]]))
  
  if(!is.na(metaanalyse$k.study.w["low or moderate risk"])) {
    ma_chart_all <- ma_chart_all %>%
      add_row(ev_name = i,
              quips = "low or moderate risk",
              n = metaanalyse[["k.study.w"]][["low or moderate risk"]],
              effect = exp(metaanalyse[["TE.random.w"]][["low or moderate risk"]]),
              lower_or = exp(metaanalyse[["lower.random.w"]][["low or moderate risk"]]),
              upper_or = exp(metaanalyse[["upper.random.w"]][["low or moderate risk"]]),
              I2 = metaanalyse[["I2.w"]][["low or moderate risk"]],
              lower_I2 = metaanalyse[["lower.I2.w"]][["low or moderate risk"]],
              upper_I2 = metaanalyse[["upper.I2.w"]][["low or moderate risk"]],
              lower_predict = exp(metaanalyse[["lower.predict.w"]][["low or moderate risk"]]), 
              upper_predict = exp(metaanalyse[["upper.predict.w"]][["low or moderate risk"]]))
  }
  
  # Create and save Forest plots for individual effect measures
  if(img_save)  {
    png(file = paste0(i, "1_orf_all_forestplot.png"), 
        width = 1600, height = 1400, res = 150)  }
    
    forest.meta(metaanalyse, #all
                sortvar = studlab, 
                header.line = "both",
                text.addline1 = " ",
                text.addline2 = paste0("Risk factor: ", i, 
                                       ", all studies"),
                fs.addline = 12,
                ff.addline = "bold",
                subgroup = k.w > 4,
                prediction = case_when(
                  metaanalyse$k >= 5 ~ TRUE,
                  metaanalyse$k < 5 ~ FALSE), 
                prediction.subgroup = k.w > 4,
                leftcols = c("studlab", "TE", "seTE"),
                leftlabs = c("Publication year, author", "logOR", "SE"),
                rightlabs = c("OR", "95% CI", "Weight"),
                xlim = case_when(exp(metaanalyse[["upper.predict"]]) < 7 
                                 ~ c(0.2, 5),
                                 exp(metaanalyse[["upper.predict"]]) >= 7 
                                 ~ c(0.2, 20)),
                col.square = "navy",
                col.square.lines = "navy",
                col.diamond = case_when(
                  metaanalyse[["k.w"]][["low or moderate risk"]] >= 5 ~ "grey",
                  metaanalyse[["k.w"]][["low or moderate risk"]] < 5 ~ "transparent"), 
                col.diamond.lines = case_when(
                  metaanalyse[["k.w"]][["low or moderate risk"]] >= 5 ~ "navy",
                  metaanalyse[["k.w"]][["low or moderate risk"]] < 5 ~ "transparent"),
                pooled.totals = FALSE,
                comb.fixed = FALSE,
                print.tau2 = TRUE,
                print.Q = TRUE,
                print.pval.Q = TRUE,
                print.I2 = TRUE,
                digits.tau2 = 2,
                digits.TE = 2,
                digits.se = 2,
                digits.Q = 1,
                fs.hetstat = 9,
                hetlab = "Het.: ",
                resid.hetlab = "Res. het.: ",
                label.test.subgroup.random = "Subgroup diff.: ")
    if(img_save)  {  dev.off()  }
  
  if(img_save) {
    png(file = paste0(i, "2_orf_all_funnelplot.png"), 
        width = 1800, height = 800, res = 150) }
    
    # Create funnel plots for individual effect measures
    col.contour = c("gray75", "gray85", "gray95")
    funnel.meta(metaanalyse,
                studlab = FALSE,
                contour = c(0.9, 0.95, 0.99),
                col.contour = col.contour)
    # Add a legend
    legend("topright", 
           legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
           fill = col.contour)
    # Add a title
    title(paste0(i, " (all studies)"))
    if(img_save)  {  dev.off()  }
} 
ma_chart_all


# Calculate for adjusted -------------------------------------------------------

ma_chart_adj <- ma_chart

for (i in ev_adj) {
  metaanalyse <- data_adj %>%
    filter(rf == i) %>%
    metagen(TE = log(oddsratio),
            seTE = se,
            studlab = ref,
            data = ., 
            sm = "OR",
            fixed = FALSE,
            random = TRUE,
            method.tau = "REML",
            hakn = TRUE,
            title = "Falls risk factors",
            subgroup = quips,
            tau.common = TRUE, 
            prediction = TRUE,
            prediction.subgroup = TRUE,
            control=list(stepadj=0.5, maxiter=10000))
  
  # Fill tibble with results
  ma_chart_adj <- ma_chart_adj %>%
    add_row(ev_name = i,
            quips = "pooled",
            n = metaanalyse[["k.study"]],
            effect = exp(metaanalyse[["TE.random"]]),
            lower_or = exp(metaanalyse[["lower.random"]]),
            upper_or = exp(metaanalyse[["upper.random"]]),
            I2 = metaanalyse[["I2"]],
            lower_I2 = metaanalyse[["lower.I2"]],
            upper_I2 = metaanalyse[["upper.I2"]],
            lower_predict = exp(metaanalyse[["lower.predict"]]), 
            upper_predict = exp(metaanalyse[["upper.predict"]]))
  
  if(!is.na(metaanalyse$k.study.w["low or moderate risk"])) {
    ma_chart_adj <- ma_chart_adj %>%
      add_row(ev_name = i,
              quips = "low or moderate risk",
              n = metaanalyse[["k.study.w"]][["low or moderate risk"]],
              effect = exp(metaanalyse[["TE.random.w"]][["low or moderate risk"]]),
              lower_or = exp(metaanalyse[["lower.random.w"]][["low or moderate risk"]]),
              upper_or = exp(metaanalyse[["upper.random.w"]][["low or moderate risk"]]),
              I2 = metaanalyse[["I2.w"]][["low or moderate risk"]],
              lower_I2 = metaanalyse[["lower.I2.w"]][["low or moderate risk"]],
              upper_I2 = metaanalyse[["upper.I2.w"]][["low or moderate risk"]],
              lower_predict = exp(metaanalyse[["lower.predict.w"]][["low or moderate risk"]]), 
              upper_predict = exp(metaanalyse[["upper.predict.w"]][["low or moderate risk"]]))
  }
  
  # Create and save Forest plots for individual effect measures
  if(img_save)  {
    png(file = paste0(i, "3_orf_adj_forestplot.png"), 
        width = 1600, height = 1200, res = 150)  }
    
    forest.meta(metaanalyse, #adjusted
                header.line = "both",
                text.addline1 = " ",
                text.addline2 = paste0("Risk factor: ", i,
                                       ", studies adjusted for sex and age"),
                fs.addline = 12,
                ff.addline = "bold",
                subgroup = k.w > 4,
                sortvar = studlab, 
                prediction = case_when(
                  metaanalyse$k >= 5 ~ TRUE,
                  metaanalyse$k < 5 ~ FALSE), 
                prediction.subgroup = k.w > 4,
                leftcols = c("studlab", "TE", "seTE"),
                leftlabs = c("Publication year, author", " logOR ", " SE "),
                rightlabs = c("OR (adj.)", "95% CI", "Weight"),
                xlim = case_when(exp(metaanalyse[["upper.predict"]]) < 7 
                                 ~ c(0.2, 5),
                                 exp(metaanalyse[["upper.predict"]]) >= 7 
                                 ~ c(0.2, 20)),
                col.square = "navy",
                col.square.lines = "navy",
                col.diamond = case_when(
                  metaanalyse[["k.w"]][["low or moderate risk"]] >= 5 ~ "grey",
                  metaanalyse[["k.w"]][["low or moderate risk"]] < 5 ~ "transparent"), 
                col.diamond.lines = case_when(
                  metaanalyse[["k.w"]][["low or moderate risk"]] >= 5 ~ "navy",
                  metaanalyse[["k.w"]][["low or moderate risk"]] < 5 ~ "transparent"),
                pooled.totals = FALSE,
                comb.fixed = FALSE,
                print.tau2 = TRUE,
                print.Q = TRUE,
                print.pval.Q = TRUE,
                print.I2 = TRUE,
                digits.TE = 2,
                digits.se = 2,
                digits.tau2 = 2,
                digits.Q = 1,
                fs.hetstat = 9,
                hetlab = "Het.: ",
                resid.hetlab = "Res. het.: ",
                label.test.subgroup.random = "Subgroup diff.: ")
    if(img_save)  {  dev.off()  }
  
  if(img_save) {
    png(file = paste0(i, "4_orf_adj_funnelplot.png"), 
        width = 1800, height = 800, res = 150)  }
    
    # Create funnel plots for individual effect measures
    col.contour = c("gray75", "gray85", "gray95")
    funnel.meta(metaanalyse,
                studlab = FALSE,
                contour = c(0.9, 0.95, 0.99),
                col.contour = col.contour)
    # Add a legend
    legend("topright", 
           legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
           fill = col.contour)
    # Add a title
    title(paste0(i, " (studies adjusted for sex and age)"))
    if(img_save)  { dev.off() }
}
ma_chart_adj


# Calculate Egger's p-values ---------------------------------------------------
ma_chart_egger <- tribble(
  ~ev_name, ~subgroup, ~n, ~egger_p)

ev_all_egger <- ev_subgroup$rf[ev_subgroup$n_all >= 10]
ev_adj_egger <- ev_subgroup$rf[ev_subgroup$n_adj >= 10]
ev_q_egger <- ev_subgroup$rf[ev_subgroup$n_Q >= 10]
ev_qadj_egger <- ev_subgroup$rf[ev_subgroup$n_Qadj >= 10]

# all
for (i in ev_all_egger) {
  metaanalyse <- data_all %>%
    filter(rf == i) %>%
    metagen(TE = log(oddsratio),
            seTE = se,
            studlab = ref,
            data = ., 
            sm = "OR",
            fixed = FALSE,
            random = TRUE,
            method.tau = "REML",
            hakn = TRUE,
            title = "Falls risk factors")
  
  ma_chart_egger <- ma_chart_egger %>%
    add_row(ev_name = i,
            subgroup = "non-adjusted, incl. high-risk bias",
            n = metaanalyse$k.study,
            egger_p = metabias(metaanalyse, method.bias = "Egger")$pval)
}

# adjusted
for (i in ev_adj_egger) {
  metaanalyse <- data_adj %>%
    filter(rf == i) %>%
    metagen(TE = log(oddsratio),
            seTE = se,
            studlab = ref,
            data = ., 
            sm = "OR",
            fixed = FALSE,
            random = TRUE,
            method.tau = "REML",
            hakn = TRUE,
            title = "Falls risk factors")
  
  ma_chart_egger <- ma_chart_egger %>%
    add_row(ev_name = i,
            subgroup = "adjusted, incl. high-risk bias",
            n = metaanalyse$k.study,
            egger_p = metabias(metaanalyse, method.bias = "Egger")$pval)
}

# non-adjusted, QUIPS
for (i in ev_q_egger) {
  metaanalyse <- data_q %>%
    filter(rf == i) %>%
    metagen(TE = log(oddsratio),
            seTE = se,
            studlab = ref,
            data = ., 
            sm = "OR",
            fixed = FALSE,
            random = TRUE,
            method.tau = "REML",
            hakn = TRUE,
            title = "Falls risk factors")
  
  ma_chart_egger <- ma_chart_egger %>%
    add_row(ev_name = i,
            subgroup = "non-adjusted, QUIPS",
            n = metaanalyse$k.study,
            egger_p = metabias(metaanalyse, method.bias = "Egger")$pval)
}

# adjusted, QUIPS
for (i in ev_qadj_egger) {
  metaanalyse <- data_qadj %>%
    filter(rf == i) %>%
    metagen(TE = log(oddsratio),
            seTE = se,
            studlab = ref,
            data = ., 
            sm = "OR",
            fixed = FALSE,
            random = TRUE,
            method.tau = "REML",
            hakn = TRUE,
            title = "Falls risk factors")
  
  ma_chart_egger <- ma_chart_egger %>%
    add_row(ev_name = i,
            subgroup = "adjusted, QUIPS",
            n = metaanalyse$k.study,
            egger_p = metabias(metaanalyse, method.bias = "Egger")$pval)
}
ma_chart_egger


## plot grouped meta-analyses results ------------------------------------------

# for all studies
if(img_save)  {
  tiff("ggplot_orf_grouped.tiff", unit = "mm",
       width = 2*107, height = 3*80, res = 300)  }
ma_chart_all %>%
  filter(quips == "pooled" & n >= 5) %>%
  ggplot(aes(x = reorder(ev_name, effect),
             middle = effect,
             lower = lower_or, ymin = lower_predict, 
             upper = upper_or, ymax = upper_predict, 
             fill = I2)) +
  geom_boxplot(stat = "identity") +
  ggtitle("") +
  ylab("Odds ratio") +
  xlab("Risk factor") + 
  labs(fill = "I\u00b2") +
  scale_fill_continuous(limits = c(0,1), low = "white", high = "#AD002AFF") +
  coord_flip() +
  geom_hline(yintercept=1) +
  scale_y_continuous(trans='log2',
                     breaks = c(0.5, 4, 32),
                     labels = c("0.5", "4", "32"))  +
  theme(axis.text.x = element_text(size = 17), 
        axis.text.y = element_text(size = 13), 
        axis.title = element_text(size = 19), 
        legend.title = element_text(size = 17), 
        legend.text = element_text(size = 13)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 12))  +
  theme(legend.position = "right") 
if(img_save)  { dev.off()  }

# for adjusted studies
if(img_save)  {
  tiff("ggplot_orf_adj_grouped.tiff", unit = "mm",
       width = 2*107, height = 3*80, res = 300)  }
ma_chart_adj %>%
  filter(quips == "pooled" & n >= 5) %>%
  ggplot(aes(x = reorder(ev_name, effect), 
             middle = effect, 
             lower = lower_or, ymin = lower_predict, 
             upper = upper_or, ymax = upper_predict, 
             fill = I2
  )) +
  geom_boxplot(stat = "identity") +
  ggtitle("") +
  ylab("Odds ratio") +
  xlab("Risk factor") + 
  labs(fill = "I\u00b2") +
  scale_fill_continuous(limits = c(0,1), low = "white", high = "#AD002AFF") +
  coord_flip()+
  #theme(legend.position="none") +
  geom_hline(yintercept=1) +
  scale_y_continuous(trans='log2',
                     breaks = c(0.125, 1, 32),
                     labels = c("0.125", "1", "32"))  +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 12))  +
  theme(legend.position = "right") 
if(img_save)  { dev.off()  }


# for QUIPS excl. high-risk
if(img_save)  {
  #tiff("ggplot_orf_q_grouped.tiff", unit = "mm", width = 2*107, height = 3*80, res = 300)  
  pdf("Figure_2.pdf", width = 10, height = 15)
  }

ma_chart_all %>%
  filter(quips == "low or moderate risk" & n >= 5) %>%
  ggplot(aes(x = reorder(ev_name, effect), 
             middle = effect, 
             lower = lower_or, ymin = lower_predict, 
             upper = upper_or, ymax = upper_predict, 
             fill = I2
  )) +
  geom_boxplot(stat = "identity") +
  ggtitle("") +
  ylab("Odds ratio") +
  xlab("Risk factor") + 
  labs(fill = "I\u00b2") +
  scale_fill_continuous(limits = c(0,1), low = "white", high = "#AD002AFF") +
  coord_flip()+
  #theme(legend.position="none") +
  geom_hline(yintercept=1) +
  scale_y_continuous(trans='log2',
                     breaks = c(0.125, 1, 4, 32),
                     labels = c("0.125", "1", "4", "32")) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 12))  +
  theme(legend.position = "right") 
if(img_save)  { dev.off()  }


# for adjusted & QUIPS excl. high-risk
if(img_save)  {
  #tiff("ggplot_orf_qadj_grouped.tiff", unit = "mm", width = 2*107, height = 2*80, res = 300)  
  pdf("Figure_3.pdf", width = 10, height = 10)
}
ma_chart_adj %>%
  filter(quips == "low or moderate risk" & n >= 5) %>%
  ggplot(aes(x = reorder(ev_name, effect), 
             middle = effect, 
             lower = lower_or, ymin = lower_predict, 
             upper = upper_or, ymax = upper_predict, 
             fill = I2
  )) +
  geom_boxplot(stat = "identity") +
  ggtitle("") +
  ylab("Odds ratio") +
  xlab("Risk factor") + 
  labs(fill = "I\u00b2") +
  scale_fill_continuous(limits = c(0,1), low = "white", high = "#AD002AFF") +
  coord_flip() +
  geom_hline(yintercept=1) +
  scale_y_continuous(trans='log2',
                     breaks = c(0.25, 1, 4, 16),
                     labels = c("0.25", "1", "4", "16")) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 12))  +
  theme(legend.position = "right") 
if(img_save)  { dev.off()  }


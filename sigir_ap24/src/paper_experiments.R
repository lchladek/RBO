source("src/common.R")
source("rbo/low.R")
source("src/rbo_vectorized.R")

library(dplyr)
library(latex2exp)
library(glue)
library(ggplot2)
theme_set(theme_bw() +
            theme(plot.subtitle = element_text(hjust = 0.5),
                  plot.margin = margin(0,5,0,5)))
.W <- 7500
.RES <- 1000
path <- "paper/"


# 3 permutations in TREC data ######################################################################

d <- lapply(.WEB_YEAR, function(year) {
  read_inputs(glue("scratch/01-trec-download/web{year}"),"random") |>
    group_by(sys,topic) |>
    summarise(nperms = prod(factorial(table(score)))) |>
    mutate(year = year) |>
    ungroup()
}) |> bind_rows()

nrow(d)
# number of permutations greater than...
d |> count(nperms>1)
d |> count(nperms>1e6)
d |> count(nperms>1e80)

# 4.1 TREC Data ####################################################################################

d <- rio::import_list(list.files("output/trec-stats/", full.names = TRUE), rbind = TRUE) |>
  select(-`_file`)

# number of runs
d |> count(year,sys) |> nrow()
# number of groups
d |> count(year,group) |> nrow()

d <- lapply(.WEB_YEAR, function(year) {
  d <- rio::import(glue("output/rbo-trec/web{year}.csv"))
  d_stats <- rio::import(glue("output/trec-stats/web{year}.csv"))

  left_join(d, d_stats, by = c("sys1"="sys", "topic"="topic")) |>
    left_join(d_stats, by = c("sys2"="sys","topic"="topic"), suffix = c("1","2"))
}) |> bind_rows() |>
  filter(n_t1 > 0 & n_t2 > 0) # only rankings with ties

# number of pairs
d |> count(p)
# ranking lengths
d |> summarize(mean((n1+n2)/2),
               max(n1),max(n2))

# analyzing bounds

bounds_trec <- d |> mutate(res_su = (rbo_hig_max - rbo_low_min),
                                 res_s = (rbo_hig_ext - rbo_low_ext),
                                 res_u = (rbo_a_max - rbo_a_min),
                                 dif_a = abs(rbo_a_ext - rbo_avg_ext),
                                 dif_min = abs(rbo_a_min - rbo_avg_min),
                                 dif_max = abs(rbo_a_max - rbo_avg_max))

# mean, mid, and large diff in residuals TREC

prop_trec <- bounds_trec |>
  group_by(p) |>
  summarize(mean(res_u), max(res_u), sum(between(res_u,.01,.1)) / n(), sum(res_u>.1) / n(),
            mean(res_s), max(res_s), sum(between(res_s,.01,.1)) / n(), sum(res_s>.1) / n(),
            mean(res_su), max(res_su), sum(between(res_su,.01,.1)) / n(), sum(res_su>.1) / n()) |>
  round(2) |> as.data.frame()

xtable::xtable(prop_trec, digits = 2)

# 4.1 USE CASE #####################################################################################

# Read runs
d <- read_inputs("scratch/01-trec-download/web2009/", "random")
d_before <- d |> filter(sys == "UMHOOsd")
d_after <- d |> filter(sys == "UMHOOsdp")

# For every topic, compute RBO scores
r <- future_lapply(unique(d_before$topic), function(top) {
  # per-topic rankings
  dd_before <- d_before |> filter(topic == top)
  dd_after <- d_after |> filter(topic == top)
  # transform to our format
  before <- extract_ranking(dd_before$doc, dd_before$score)
  after <- extract_ranking(dd_after$doc, dd_after$score)
  # calculate %ties
  ties_before <- lengths(before)
  ties_before <- sum(ties_before[ties_before>1])/sum(ties_before)
  ties_after <- lengths(after)
  ties_after <- sum(ties_after[ties_after>1])/sum(ties_after)

  m <- low_rankings(before, after)
  data.frame(ties_before = ties_before,
             ties_after = ties_after,
             rbo_low_min = rbo(m$x, m$y, .95, score = "min")$min,
             rbo_avg_min = rbo_avg_min(before, after, .95))
}) |> bind_rows() |> # append AP scores (copy-pasted from TREC summaries)
  mutate(
    ap_before = c(0.539903, 0.441785, 0.000727, 0.115018, 0.175231, 0.101039,
                  0.034696, 0.05222, 0.10495, 0.336969, 0.105364, 0.211183, 0.0591,
                  0.017638, 0.134902, 0.364686, 0.082242, 0.106601, 0, NA, 0.414696,
                  0.506828, 0.007539, 0.530083, 0.290742, 0.063088, 0.298102, 0.398537,
                  0.031245, 0.188065, 0.264901, 0.03317, 0.47328, 0.034604, 0.405827,
                  0.043304, 0.062112, 0.115024, 0.101431, 0.141821, 0.55667, 0.108844,
                  0.443683, 0.035068, 0.226496, 0.754691, 0.530057, 0.090825, 0.295222,
                  0.065158),
    ap_after=c(0.539903, 0.441785, 0.000727, 0.115018, 0.175231, 0.101039,
               0.034696, 0.05222, 0.10495, 0.336969, 0.105364, 0.211183, 0.0591,
               0.017638, 0.134902, 0.364686, 0.082242, 0.106601, 0, NA, 0.414696,
               0.506828, 0.007539, 0.530083, 0.290742, 0.063088, 0.298102, 0.398537,
               0.010618, 0.188065, 0.264901, 0.03317, 0.47328, 0.034604, 0.405827,
               0.043304, 0.062112, 0.115024, 0.101431, 0.141821, 0.55667, 0.108844,
               0.443683, 0.035068, 0.226496, 0.754691, 0.530057, 0.090825, 0.295222,
               0.065158)
  )

# % of ties
mean(r$ties_before)
mean(r$ties_after)
# RBO above threshold?
sum(r$rbo_avg_min>.99)
sum(r$rbo_low_min>.99)

png(glue("{path}/UMHOOsd_ap.png"), width = .W/3, height = .W/3*.98, res = .RES)
ggplot(r, aes(ap_before, ap_after)) +
  geom_abline(intercept = 0, slope = 1, linewidth = .25, color = "red") +
  geom_point(size = 1, stroke = NA) +
  labs(x = TeX("UMHOOsd"), y = TeX("UMHOOsdp"), subtitle = "AP")
dev.off()

png(glue("{path}/UMHOOsd_rbo.png"), width = .W/3*1.1, height = .W/3, res = .RES)
ggplot(r, aes(rbo_avg_min, rbo_low_min)) +
  geom_abline(intercept = 0, slope = 1, linewidth = .25, color = "red") +
  geom_point(size = 1, stroke = NA) +
  labs(x = TeX("$RBO^{avg}_{MIN}$"), y = TeX("$RBO^{low}_{MIN}$"), subtitle = "RBO") +
  coord_cartesian(xlim=c(.98,1),ylim=c(.98,1))
dev.off()



# SYNTETHIC DATA

synt_data <- read.csv("output/rbo-synthetic/rbo.csv")
dim(synt_data)
head(synt_data)

# general stats about generated rankings and number of ties

synt_data |>
  mutate(tot_n = (n_x + n_y), frac_ties = ((n_t_x/n_x + n_t_y/n_y)/2)) |>
  summarize(mean_n = mean((n_x + n_y)/2),
            mean_t = mean(frac_ties))


# analyzing bounds

bounds <- synt_data |> mutate(res_su = (rbo_hig_max - rbo_low_min),
                               res_s = (rbo_hig_ext - rbo_low_ext),
                               res_u = (rbo_a_max - rbo_a_min),
                               dif_a = abs(rbo_a_ext - rbo_avg_ext),
                               dif_min = abs(rbo_a_min - rbo_avg_min),
                               dif_max = abs(rbo_a_max - rbo_avg_max))

#checking differences on average

bounds |> select(dif_a, dif_min, dif_max) |> colMeans()

# mid and large diff in residuals SYNT by s

synt_s_1 <- bounds |>
  group_by(p) |>
  filter(n_x <= 25) |>
  summarize(mean(res_u), max(res_u), sum(between(res_u,.01,.1)) / n(), sum(res_u>.1) / n(),
            mean(res_s), max(res_s), sum(between(res_s,.01,.1)) / n(), sum(res_s>.1) / n(),
            mean(res_su), max(res_su), sum(between(res_su,.01,.1)) / n(), sum(res_su>.1) / n()) |>
  round(2) |> as.data.frame()

xtable::xtable(synt_s_1, digits = 2)

synt_s_2 <- bounds |>
  group_by(p) |> filter(n_x > 25 & n_x <=50) |>
  summarize(mean(res_u), max(res_u), sum(between(res_u,.01,.1)) / n(), sum(res_u>.1) / n(),
            mean(res_s), max(res_s), sum(between(res_s,.01,.1)) / n(), sum(res_s>.1) / n(),
            mean(res_su), max(res_su), sum(between(res_su,.01,.1)) / n(), sum(res_su>.1) / n()) |>
  round(2) |> as.data.frame()

xtable::xtable(synt_s_2, digits = 2)

synt_s_3 <- bounds |>
  group_by(p) |> filter(n_x > 50) |>
  summarize(mean(res_u), max(res_u), sum(between(res_u,.01,.1)) / n(), sum(res_u>.1) / n(),
            mean(res_s), max(res_s), sum(between(res_s,.01,.1)) / n(), sum(res_s>.1) / n(),
            mean(res_su), max(res_su), sum(between(res_su,.01,.1)) / n(), sum(res_su>.1) / n()) |>
  round(2) |> as.data.frame()

xtable::xtable(synt_s_3, digits = 2)

# mid and large diff in residuals SYNT by tie prop

synt_t_1 <- bounds |>
  group_by(p) |>
  mutate(ties_frac = ( (n_t_x/n_x) + (n_t_y/n_y) )/2 ) |>
  filter(ties_frac <= .4) |>
  summarize(mean(res_u), max(res_u), sum(between(res_u,.01,.1)) / n(), sum(res_u>.1) / n(),
            mean(res_s), max(res_s), sum(between(res_s,.01,.1)) / n(), sum(res_s>.1) / n(),
            mean(res_su), max(res_su), sum(between(res_su,.01,.1)) / n(), sum(res_su>.1) / n()) |>
  round(2) |> as.data.frame()

xtable::xtable(synt_t_1, digits = 2)

synt_t_2 <- bounds |>
  group_by(p) |>
  mutate(ties_frac = ( (n_t_x/n_x) + (n_t_y/n_y) )/2 ) |>
  filter(ties_frac > .4 & ties_frac <= .6) |>
  summarize(mean(res_u), max(res_u), sum(between(res_u,.01,.1)) / n(), sum(res_u>.1) / n(),
            mean(res_s), max(res_s), sum(between(res_s,.01,.1)) / n(), sum(res_s>.1) / n(),
            mean(res_su), max(res_su), sum(between(res_su,.01,.1)) / n(), sum(res_su>.1) / n()) |>
  round(2) |> as.data.frame()

xtable::xtable(synt_t_2, digits = 2)

synt_t_3 <- bounds |>
  group_by(p) |>
  mutate(ties_frac = ( (n_t_x/n_x) + (n_t_y/n_y) )/2 ) |>
  filter(ties_frac > .6) |>  summarize(mean(res_u), max(res_u), sum(between(res_u,.01,.1)) / n(), sum(res_u>.1) / n(),
                                        mean(res_s), max(res_s), sum(between(res_s,.01,.1)) / n(), sum(res_s>.1) / n(),
                                        mean(res_su), max(res_su), sum(between(res_su,.01,.1)) / n(), sum(res_su>.1) / n()) |>
  round(2) |> as.data.frame()

xtable::xtable(synt_t_3, digits = 2)


## Scatter plots for syntetic data to be used in the paper

#p=0.95
png(glue("{path}/syn_s_u_95.png"), width = .W/3, height = .W/3*.9, res = .RES)
bounds |> filter(p == .95) |>
  ggplot(aes(res_s, res_u, alpha = abs(res_s-res_u))) +
  geom_point(size = .25, shape = 19) +
  geom_abline(intercept = 0, slope = 1, linewidth = .25, color = "red") +
  labs(x = TeX("$RES_S$"), y = TeX("$RES_U$")) +
  scale_alpha_continuous(range = c(.02,.3)) + # for better display of large differences
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1)
dev.off()

#p=0.9
png(glue("{path}/syn_s_u_9.png"), width = .W/3, height = .W/3*.9, res = .RES)
bounds |> filter(p == .9) |>
  ggplot(aes(res_s, res_u, alpha = abs(res_s-res_u))) +
  geom_point(size = .25, shape = 19) +
  geom_abline(intercept = 0, slope = 1, linewidth = .25, color = "red") +
  labs(x = TeX("$RES_S$"), y = TeX("$RES_U$")) +
  scale_alpha_continuous(range = c(.02,.3)) + # for better display of large differences
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1)
dev.off()

#p=0.8
png(glue("{path}/syn_s_u_8.png"), width = .W/3, height = .W/3*.9, res = .RES)
bounds |> filter(p == .8) |>
  ggplot(aes(res_s, res_u, alpha = abs(res_s-res_u))) +
  geom_point(size = .25, shape = 19) +
  geom_abline(intercept = 0, slope = 1, linewidth = .25, color = "red") +
  labs(x = TeX("$RES_S$"), y = TeX("$RES_U$")) +
  scale_alpha_continuous(range = c(.02,.3)) + # for better display of large differences
  scale_y_continuous(limits = c(0,1)) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "none")
dev.off()



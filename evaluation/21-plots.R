# Loading the library
suppressPackageStartupMessages(source("src/lib.R"))

library(latex2exp)
library(plyr)
library(viridis)
library(ggrastr)

theme_set( theme_light(base_family = "Linux Biolinum O"))
theme_update(legend.title.position = "left",
             legend.position = "inside",
             legend.position.inside = c(0.9, 0.7),
             legend.title = element_text(angle = 90, hjust = 0.5))


# Plot of K_n with varying p-values
dp <- 0.001
ps <- seq(0.74, 0.91, dp)
select_ps <- c(0.8, 0.85, 0.9)
n <- 7
continuous_k_n <- ldply(ps, function(p) {
  ms <- seq(1, n, 0.01)
  y <- (1 - p) / p * sapply(ms - 1, function(m) {
    sum(sapply(seq(50), function(d) { p^(d + m) / (d + m)} ))
  })
  intm = ms
  intm[which(ms != as.integer(ms))] <- NA
  selp <- rep(as.character(select_ps[which(select_ps==p)][1]), length(ms))
  data.frame(p, selp, m=ms, intm, y)
})

k_n_select <- ldply(select_ps, function(p) {
  ms <- seq(1, n, 0.01)
  y <- (1 - p) / p * sapply(ms - 1, function(m) {
    sum(sapply(seq(50), function(d) { p^(d + m) / (d + m)} ))
  })
  intm = ms
  intm[which(ms != as.integer(ms))] <- NA
  lab <- as.character(p)
  data.frame(p, m=ms, intm, y, lab)
})

summary_p <- k_n_select[which(!is.na(k_n_select$intm)), c("p", "m", "y")]
rownames(summary_p) <- NULL
summary_p$y <- round(summary_p$y, 3)
print(summary_p)

ggplot(continuous_k_n, aes(x=intm, y=y, color=p)) +
  scale_y_continuous(limits=c(0, 0.47)) +
  scale_x_continuous(breaks=seq(n), minor_breaks = c()) +
  scale_shape_manual(values=c(1, 3, 6)) +
  scale_color_distiller() +
  rasterise(geom_point(aes(x=m), size=0.1), dpi=300) +
  geom_point(data=k_n_select, aes(shape=lab), size=1.5, color="black") +
  geom_line(data=k_n_select, aes(x=m, shape=lab), linewidth=0.4, color="black") +
  labs(x=TeX("n"), y=TeX("$K_n$"), color="p", shape="select p") +
  theme(axis.title.y = element_text(angle = 0, vjust=0.5))

ggsave("plots/Kn_for_varying_p.pdf", device=cairo_pdf,
       width = 12, height = 9, units= "cm")


# The example in the paper: comparison of basic and culling convolutoin
pair <- lapply(list(x="(a b c)", y="(a b) c"), from_string)

est_pmf_basic <- estimate_pmf(pair, p=0.9, convolution_fn=basic_convolution)
est_pmf_cull <- estimate_pmf(pair, p=0.9, convolution_fn=cull_all_convolution)
real_pmf <- complete_pmf(pair, p=0.9)

ggplot(mapping=aes(x=rbo, xend=rbo, y=prob, yend=0)) +
  scale_y_continuous(limits=c(0, 0.65)) +
  scale_x_continuous(limits=c(.3, .7)) +
  scale_color_brewer(type="qual", direction=-1, palette=3) +
  geom_segment(data=real_pmf, aes(color="Actual"), linewidth=1.5) +
  geom_segment(data=est_pmf_basic, aes(color="Basic convolution"), linewidth=1.5, position=position_nudge(x=0.005)) + 
  geom_segment(data=est_pmf_cull, aes(color="Culling convolution"), linewidth=1.5, position=position_nudge(x=0.01)) +
  geom_point(data=real_pmf, aes(color="Actual"), size=3) +
  geom_point(data=est_pmf_basic, aes(color="Basic convolution"), size=3, position=position_nudge(x=0.005)) +
  geom_point(data=est_pmf_cull, aes(color="Culling convolution"), size=3, position=position_nudge(x=0.01)) +
  labs(x="RBO value", y="Probability", colour="Estimation method")

ggsave("plots/counterexample.pdf", device=cairo_pdf,
       width = 12, height = 6, units= "cm")

emdist::emdw(real_pmf$rbo, real_pmf$prob, est_pmf_basic$rbo, est_pmf_basic$prob)
emdist::emdw(real_pmf$rbo, real_pmf$prob, est_pmf_cull$rbo, est_pmf_cull$prob)

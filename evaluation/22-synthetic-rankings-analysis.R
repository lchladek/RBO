
# Loading the library
suppressPackageStartupMessages(source("src/lib.R"))

library(ggrastr)
library(ggridges)
library(latex2exp)

theme_set( theme_light(base_family = "Linux Biolinum O"))
theme_update(legend.title.position = "left",
             legend.position = "inside",
             legend.position.inside = c(0.9, 0.3),
             legend.title = element_text(angle = 90, hjust = 0.5))


# Adjust this to the data directory
setwd("data/experiments")

numeric_dirs <- na.omit(suppressWarnings(as.integer(dir())))
sort(numeric_dirs, decreasing=TRUE)

latest_experiment <- numeric_dirs[1]
if (is.na(latest_experiment)) {
  stop("Experiment not found in the specified folder")
}
message(paste("Using the latest experiment: ", latest_experiment))
setwd(as.character(latest_experiment))
setwd("0")

# Import outputs from Gourd, in CSV format with filename 'stdout'
run_dirs <- dir()
files <- lapply(run_dirs, function(d) {
  fname <- paste(d, "stdout", sep="/")
  if (file.exists(fname)) {
    return (fname)
  } else {
    message(paste("No file at", fname))
    return (NULL)
  }
})

message(paste("Found", length(files), "simulation output files"))
all_data <- data.frame()
for (f in files) {
  contents <- read.csv(f)
  nentries <- nrow(contents)
  
  for (col in colnames(contents)) {
    if (length(which(col==c("ranking_x", "ranking_y", "kl_divergence",
                            "earth_movers_distance")))>0) {
      next
    } else {
      contents[[col]] <- suppressWarnings(as.numeric(contents[[col]]))
      nas <- which(is.na(contents[[col]]))
      if (length(nas)>0) {
        message(paste("invalid values for", col, ":", length(nas)))
        contents <- contents[-nas,]
      }
    }
  }
  not_strings <- which(!is.na(suppressWarnings(as.numeric(contents$ranking_x))))
  if (length(not_strings)>0) {
    contents <- contents[-not_strings,]
    message(paste("invalid values for ranking_x :", length(not_strings)))
  }
  not_strings <- which(!is.na(suppressWarnings(as.numeric(contents$ranking_y))))
  if (length(not_strings)>0) {
    contents <- contents[-not_strings,]
    message(paste("invalid values for ranking_y :", length(not_strings)))
  }
  contents$length <- suppressWarnings(as.integer(contents$length))
  nas <- which(is.na(contents$length))
  if (length(nas)>0) {
    message(paste("invalid values for length :", length(nas)))
    contents <- contents[-nas,]
  }
  contents$earth_movers_distance <- as.numeric(contents$earth_movers_distance)
  contents$kl_divergence <- as.numeric(contents$kl_divergence)
  
  message(paste(f, ":", nrow(contents), "entries  (", nentries - nrow(contents), "bad rows )"))
  all_data <- rbind(all_data, contents)
}

lo_hi_off <- c(which(all_data$pmf.q0 != all_data$rbo_low), which(all_data$pmf.q1 != all_data$rbo_high))
if (length(lo_hi_off) > 0) {
  message("RBO low/high failed!")
  print(lo_hi_off)
}

item_count_category <- function(item_count) {
  if (item_count < 12) "XS"
  else if (item_count < 18) "S"
  else if (item_count < 24) "M"
  else if (item_count < 30) "L"
  else if (item_count >= 30) "XL"
}

all_data$size <- factor(sapply(all_data$item_count, item_count_category), c("XS", "S", "M", "L", "XL"))

emd_data <- all_data[!is.na(all_data$earth_movers_distance),]
emd_data <- emd_data[emd_data$p == 0.9,]
message(paste(nrow(emd_data), "rankings with Earth Mover's Distance"))

# data coverage summary
ggplot(emd_data, aes(x=length)) +
  geom_bar()
ggplot(emd_data, aes(x=size, fill=size)) +
  geom_bar()
ggplot(emd_data, aes(x=p, fill=as.factor(p))) +
  geom_bar()

print(
  emd_data %>%
    group_by(size) %>%
    count()
)

ggplot(emd_data, aes(x=earth_movers_distance)) + 
  geom_histogram()

ggplot(emd_data, aes(x=size, fill=size, y=earth_movers_distance)) + 
  geom_violin()

ggplot(emd_data, aes(x=permutations, y=earth_movers_distance, color=size)) +
  scale_x_log10() +
  geom_point(size=0.5)

ggplot(emd_data, aes(x=pmf.mean, y=est.mean, color=size)) +
  geom_point(size=0.1)

ggplot(emd_data, aes(x=pmf.q0.025, y=est.q0.025)) +
  scale_x_continuous(limits=c(0.5, 1)) +
  #scale_alpha_distiller(palette = "Spectral", limits=c(0, 1)) +
  geom_bin2d(aes(alpha=after_stat(ndensity)),fill='black', bins = 200) #+
  #geom_abline(aes(slope=1, intercept = 0))

ggplot(emd_data, aes(x=pmf.q0.5, y=est.q0.5)) +
  scale_fill_distiller(direction=1, trans = "log10") +
  geom_bin2d(aes(fill=after_stat(ndensity)), bins = 200)

ggplot(emd_data, aes(x=pmf.q0.5, y=est.q0.5-pmf.q0.5)) +
  scale_x_binned(n.breaks=10) +
  geom_violin()

ggplot(emd_data, aes(x=pmf.q0.1, y=est.q0.1, color=permutations)) +
 scale_x_continuous(limits=c(0.45, 0.55)) +
  scale_y_continuous(limits=c(0.45, 0.55)) +
  geom_point(size=0.1, alpha=1)

ggplot(emd_data, aes(x=est.q0.1-pmf.q0.1, y=earth_movers_distance, colour = size)) +
  geom_point(size=0.1, alpha=0.5)

ggplot(emd_data, aes(x=est.q0.9-pmf.q0.9, y=est.q0.1-pmf.q0.1, colour = size)) +
  geom_point(size=0.1, alpha=0.5)

ggplot(emd_data, aes(x=pmf.mean, y=(pmf.mean-est.mean)/pmf.mean, colour = size)) +
  geom_point(size=0.5)

ggplot(emd_data, aes(x=rbo_a, y=(pmf.mean-est.mean)/pmf.mean, colour = size)) +
  geom_point(size=0.1, alpha=0.5)

ggplot(emd_data, aes(x=rbo_a, y=(pmf.mean-est.mean)/pmf.mean, colour = size)) +
  geom_point(size=0.5)

# always overestimates the uncertainty
ggplot(emd_data, aes(x=pmf.q0.975-pmf.q0.025, y=est.q0.975-est.q0.025, colour = size)) +
  geom_point(size=0.1, alpha=0.5) +
  geom_abline(aes(slope=1, intercept = 0), alpha=0.5, linewidth=0.3)

ggplot(emd_data, aes(x=pmf.q1, y=est.q1)) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(0, 1)) +
  scale_fill_distiller(direction=1, trans = "log10") +
  geom_bin2d(aes(fill=after_stat(density)), bins = 200) +
  geom_abline(aes(slope=1, intercept = 0), alpha=0.5, linewidth=0.3)

ggsave("~/RP/lchladek/RBO/plots/range.pdf", device=cairo_pdf,
       width = 12, height = 12, units= "cm")

ggplot(emd_data, aes(x=pmf.q1-pmf.q0, y=est.q1-est.q0-(pmf.q1-pmf.q0))) +
  scale_x_continuous(limits=c(-0.01, .5)) +
  scale_y_continuous(limits=c(-0.01, .2)) +
  scale_fill_fermenter(direction=1, trans = "log10") +
  geom_bin2d(aes(fill=after_stat(ndensity)), bins = 200)

ggsave("~/RP/lchladek/RBO/plots/dist_q1.pdf", device=cairo_pdf,
       width = 12, height = 12, units= "cm")

ggplot(emd_data, aes(x=pmf.mean, y=est.mean)) +
  scale_fill_distiller(direction=1, trans = "log10") +
  geom_bin2d(aes(fill=after_stat(ndensity)), bins = 200) +
  geom_abline(aes(slope=1, intercept = 0), alpha=0.5, linewidth=0.3) +
  labs(x=TeX("$RBO_a$"), y="Estimated PMF mean", color="Absolute error", fill="Normalised density")

ggsave("~/RP/lchladek/RBO/plots/RBO_a.pdf", device=cairo_pdf,
       width = 12, height = 12, units= "cm")

ggplot(emd_data, aes(x=pmf.mean, y=est.mean)) +
  scale_x_continuous(limits=c(0.46, .54)) +
  scale_y_continuous(limits=c(0.46, .54)) +
  scale_fill_distiller(direction=1) +
  scale_color_grey() +
  geom_bin2d(aes(fill=after_stat(ndensity)), bins = 100) +
  geom_abline(aes(slope=1, intercept = 0, , color="0"), alpha=0.5, linewidth=0.5) +
  geom_abline(aes(slope=1, intercept = 0.01, color="0.01"), alpha=0.5,linewidth=0.5) +
  geom_abline(aes(slope=1, intercept = -0.01, color="0.01"), alpha=0.5,linewidth=0.5) +
  geom_abline(aes(slope=1, intercept = 0.02, color="0.02"), alpha=0.5,linewidth=0.5) +
  geom_abline(aes(slope=1, intercept = -0.02, color="0.02"), alpha=0.5,linewidth=0.5) +
  labs(x=TeX("$RBO_a$"), y="Estimated PMF mean", color="Absolute error", fill="Normalised density")

ggsave("~/RP/lchladek/RBO/plots/RBO_a_zoomed.pdf", device=cairo_pdf,
       width = 6, height = 6, units= "cm")

ggplot(emd_data, aes(x=pmf.mean, y=est.mean, color=size)) +
  scale_x_continuous(limits=c(0.47, .53)) +
  scale_y_continuous(limits=c(0.47, .53)) +
  geom_point( size=0.5)

sizes_counts <- sapply(levels(emd_data$size), function(s) {length(which(emd_data$size==s))})
emd_data$size_prop <- 1/sizes_counts[emd_data$size]

ggplot(emd_data, aes(weight=size_prop, fill=size, y=after_stat(density))) +
  scale_x_continuous(limits=c(-0.025, 0.05)) +
  geom_histogram(aes(x=(est.q0.975-pmf.q0.975)))

ggplot(emd_data, aes(weight=size_prop/4, fill=size, y=after_stat(count))) +
  scale_y_continuous(limits=c(0, 0.73), labels=scales::percent) +
  scale_fill_brewer(palette=9, direction=-1) +
  scale_x_continuous(limits=c(-0.08, 0.022), labels=scales::percent, n.breaks=8) +
  geom_histogram(aes(x=(est.q0.025-pmf.q0.025)/pmf.q0.025), binwidth = 0.01) +
  labs(x="Relative error in the 2.5th percentile", y="Proportion", fill="Ranking size")

ggsave("~/RP/lchladek/RBO/plots/dist_2.5.pdf", device=cairo_pdf,
       width = 6, height = 12, units= "cm")

ggplot(emd_data, aes(weight=size_prop, fill=size, y=after_stat(count/4))) +
  scale_y_continuous(limits=c(0, 0.73), labels=scales::percent) +
  scale_fill_brewer(palette=9, direction=-1) +
  scale_x_continuous(limits=c(-0.022, 0.1), labels=scales::percent, n.breaks=8) +
  geom_histogram(aes(x=(est.q0.975-pmf.q0.975)/pmf.q0.975), binwidth=0.01) +
  labs(x="Relative error in the 97.5th percentile")


ggsave("~/RP/lchladek/RBO/plots/dist_97.5.pdf", device=cairo_pdf,
       width = 6, height = 12, units= "cm")

p_sequence=c(0.75,0.8,0.85,0.9,0.95)
p_seq_indices <- function(first_index) {seq_along(p_sequence) - 1 + first_index}

p_first_indices <- which(all_data$p == p_sequence[1])

# select only those which are valid for a p sequence
valid_low_indices <- sapply(low_indices, function(first_index) {
  seq_indices <- p_seq_indices(first_index)
  is_p_sequence <- as.logical(prod(
    all_data$p[seq_indices] == p_sequence
  ))
  if (!is_p_sequence) {
    return (FALSE)
  }
  f <- as.factor(all_data$ranking_x[seq_indices])
  return (length(levels(f)) == 1)
})

message(paste("Found", sum(valid_low_indices), "valid p-sequences out of", length(low_indices)))

p_data <- all_data[as.vector(
  sapply(low_indices[valid_low_indices], function(index) {
    seq_along(p_sequence) - 1 + index
  })
),]

p_data$p_pair_index <- as.factor(
  as.vector(sapply(rank(low_indices[valid_low_indices]), rep, length(p_sequence)))
)

ggplot(p_data, aes(x=pmf.mean, y=est.mean)) +
  scale_fill_distiller(direction=1, trans = "log") +
  geom_bin2d(aes(fill=after_stat(density)), bins = 200) +
  geom_abline(aes(slope=1, intercept = 0), alpha=0.5, linewidth=0.3)

ggplot(p_data, aes(x=size, fill=as.factor(size))) +
  geom_bar()

p_data$range_mse <- with(p_data, ((pmf.q0-pmf.q1)-(est.q0-est.q1))^2)

q0.025_diff <- with(emd_data, pmf.q0.025 - pmf.q0)
q0.975_diff <- with(emd_data, pmf.q1 - pmf.q0.975)

max_q0.975_diff<- max(q0.025_diff)
max_q0.025_diff<- max(q0.975_diff)

max9row <- emd_data[which(q0.975_diff == max(q0.975_diff)),]
max2row <- emd_data[which(q0.025_diff == max(q0.025_diff)),]

row_df <- function(row) {
  qs <- seq(0, 1, 0.01) 
  pmf.q <- sapply(qs, function(q) {row[[paste("pmf.q", q, sep='')]]})
  est.q <- sapply(qs, function(q) {row[[paste("est.q", q, sep='')]]})
  data.frame(q=qs, pmf.q, est.q)
}


ggplot(row_df(max9row), aes(x=q)) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(0, 1)) +
  geom_line(aes(y=pmf.q, color="real")) +
  geom_line(aes(y=est.q, color="estimated"))

ggplot(row_df(max2row)) +
  scale_y_continuous(limits=c(0, 1)) +
  geom_histogram(aes(x=pmf.q, y=after_stat(ndensity), fill="real"), alpha=0.5) +
  geom_histogram(aes(x=est.q, y=after_stat(ndensity), fill="est"), alpha=0.5)

ggplot(row_df(max2row), aes(x=q)) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(0, 1)) +
  geom_line(aes(y=pmf.q, color="real")) +
  geom_line(aes(y=est.q, color="estimated"))


ggplot(max9row_df, aes(x=q)) +
  geom_line(aes(y=pmf.q))


mse_comp <- emd_data %>% group_by(size) %>%
  summarise(mean_emd=mean(earth_movers_distance),
            mse_q0=mean((est.q0-pmf.q0)^2),
            mse_mean=mean((est.mean-pmf.mean)^2),
            mse_var=mean((est.var-pmf.var)^2),
            mse_q0.025=mean((est.q0.025-pmf.q0.025)^2),
            mse_q0.05=mean((est.q0.05-pmf.q0.05)^2),
            mse_q0.95=mean((est.q0.95-pmf.q0.95)^2),
            mse_q0.975=mean((est.q0.975-pmf.q0.975)^2),
            mse_q1=mean((est.q1-pmf.q1)^2))

mse_comp <- rbind(mse_comp, 
                  mse_comp <- emd_data %>%
                    summarise(size="ALL",
                              mean_emd=mean(earth_movers_distance),
                              mse_mean=mean((est.mean-pmf.mean)^2),
                              mse_var=mean((est.var-pmf.var)^2),
                              mse_q0=mean((est.q0-pmf.q0)^2),
                              mse_q0.025=mean((est.q0.025-pmf.q0.025)^2),
                              mse_q0.05=mean((est.q0.05-pmf.q0.05)^2),
                              mse_q0.95=mean((est.q0.95-pmf.q0.95)^2),
                              mse_q0.975=mean((est.q0.975-pmf.q0.975)^2),
                              mse_q1=mean((est.q1-pmf.q1)^2)))

as.data.frame(mse_comp)
write.csv(mse_comp, "mse_comp2.csv")

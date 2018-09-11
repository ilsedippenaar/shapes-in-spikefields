jaws_color <- 'orange'
zorin_color <- 'purple'

load_mat <- function(name, folder_name, convert=T) {
  if (!exists("data_dir", mode = "character")) source("config.R")
  if (missing(folder_name)) folder_name <- sys.frame(0)$folder_name
  d <- readMat(file.path(data_dir, folder_name, paste(name, '.mat', sep = '')))
  if (convert)  d[['x']] <- as.matrix(as.numeric(d[['x']])) # make x a column vector
  d
}

std_dev_plot_from_mat <- function(name, ...) {
  d <- load_mat(name)
  x <-  d[['x']]
  y <- d[['y']]
  std_err <- d[['std.err']]
  std_dev_plot(x, y, std_err, ...)
}

std_dev_plot_from_mat_monkeys <- function(name, ...) {
  jaws <- load_mat(paste0(name, '_jaws'))
  zorin <- load_mat(paste0(name, '_zorin'))
  x <- rbind(zorin$x, jaws$x)
  y <- rbind(zorin$y, jaws$y)
  std_err <- rbind(zorin$std.err, jaws$std.err)
  name <- t(c(rep('Z', length(zorin$x)), 
               rep('J', length(jaws$x))))
  std_dev_plot_monkeys(x, y, std_err, name, ...)
}

std_dev_plot <- function(x, y, std_err, color = 'blue', legend_labels=NULL) {
  n <- dim(y)[2]
  if (n > 1 & length(color) != n) {
    colors <- matlab::fliplr(matlab::jet.colors(n))
  } else {
    colors <- color
  }
  means <- data.frame(x=x, y=y)
  colnames(means) <- c('x', 1:n)
  stds <- data.frame(x=x, std_err=std_err)
  colnames(stds) <- c('x', 1:n)
  p <- means %>% 
    gather('name', y, -x) %>%
    left_join(stds %>% gather('name', std_err, -x), by = c('x', 'name')) %>%
    ggplot(aes(x = x)) +
    #geom_line(aes(y = y, col = name), size=1) +
    geom_ribbon(aes(ymin = y-std_err, ymax = y+std_err, fill = name), alpha=0.2) +
    scale_color_manual(labels=legend_labels, name = NULL, values = colors) +
    scale_fill_manual(values = colors) +
    guides(fill=F) +
    theme_classic() +
    theme(text = element_text(size=16))
  if (is.null(legend_labels)) p <- p + guides(color = F)
  p
}

std_dev_plot_monkeys <- function(x, y, std_err, name, colors = c('Z'=zorin_color, 'J'=jaws_color)) {
  data.frame(x=x, y=y, std_err=std_err, name=as.factor(name)) %>%
    ggplot(aes(x=x)) +
    geom_line(aes(y = y, col = name), size=1) +
    geom_ribbon(aes(ymin = y-std_err, ymax = y+std_err, fill = name), alpha=0.2) +
    scale_color_manual(name = NULL, values = colors) +
    scale_fill_manual(values = colors) +
    guides(fill=F) +
    theme_classic() +
    theme(text = element_text(size=16))
}

plot_low_high_phase_hist <- function(name, n_bins=60) {
  breaks <- seq(0, 360, length.out = n_bins+1)
  phases <- lapply(load_mat(name, convert = F)$x, function(x) x[[1]])
  phases <- lapply(phases, function(x) rad2deg(x) %% 360)
  phases <- vapply(phases, function(x) table(cut(x, breaks = breaks)) / length(x), numeric(n_bins))
  phases <- data.frame(phases)
  colnames(phases) <- 1:4
  freq_names <- rep(c("Low (7.5-25 Hz)", "High (31-100 Hz)"), each=2)
  freq_names <- factor(freq_names, levels=freq_names[c(1,3)])
  phases <- phases %>%
    mutate(x = rollmean(seq(0, 360, length.out = n()+1), 2)) %>%
    gather(num, y, -x) %>%
    mutate(num = as.numeric(num),
           success = ifelse(num %% 2 == 1, 'yes', 'no'),
           f = freq_names[num]) %>%
    select(-num) %>%
    group_by(f) %>%
    spread(success, y) %>%
    mutate(rate = yes / (yes+no))
  phase_means <- phases %>%
    group_by(f) %>%
    summarize(v = mean(rate*exp(1i*deg2rad(x)))) %>%
    mutate(x = rad2deg(Arg(v)) %% 360,
           xend = x,
           y = 0,
           yend = Mod(v))
  phases %>%
    ggplot(aes(x = x, y=rate)) +
    geom_polygon(fill=NA, col='#0000FF') +
    geom_segment(data = phase_means, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = F, size=0.5, col='darkgrey') +
    coord_polar(start = -pi/2, direction = -1) +
    facet_wrap(~ f) +
    theme(panel.grid.major.y = element_line(size = 0.5, color = 'darkgrey'),
          axis.ticks.y = element_blank(),
          panel.spacing = unit(0.2, 'cm')) +
    scale_x_continuous(name = NULL) +
    scale_y_continuous(name = NULL, labels = NULL)
}

plot_phase_success <- function(name, n_bins=60) {
  breaks <- seq(0, 360, length.out = n_bins+1)
  phases <- lapply(load_mat(name, convert = F)$x, function(x) x[[1]])
  phases <- lapply(phases, function(x) rad2deg(x) %% 360)
  phases <- vapply(phases, function(x) table(cut(x, breaks = breaks)) / length(x), numeric(n_bins))
  phases <- data.frame(phases)
  colnames(phases) <- 1:2
  phases <- phases %>%
    mutate(x = rollmean(seq(0, 360, length.out = n()+1), 2)) %>%
    gather(num, y, -x) %>%
    mutate(num = as.numeric(num),
           success = ifelse(num == 1, 'yes', 'no')) %>%
    select(-num) %>%
    spread(success, y) %>%
    mutate(rate = yes / (yes+no))
  phase_means <- phases %>%
    summarize(v = mean(rate*exp(1i*deg2rad(x)))) %>%
    mutate(x = rad2deg(Arg(v)) %% 360,
           xend = x,
           y = 0,
           yend = Mod(v))
  phases %>%
    ggplot(aes(x = x, y=rate)) +
    geom_polygon(fill=NA, col='#0000FF') +
    geom_segment(data = phase_means, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = F, size=0.5, col='darkgrey') +
    coord_polar(start = -pi/2, direction = -1) +
    theme(panel.grid.major.y = element_line(size = 0.5, color = 'darkgrey'),
          axis.ticks.y = element_blank(),
          panel.spacing = unit(0.2, 'cm')) +
    scale_x_continuous(name = NULL) +
    scale_y_continuous(name = NULL, labels = NULL)
}

plot_rxn_time_power_bin <- function(name, monkey_name) {
  mat <- load_mat(name, convert = F)
  pow_lab <- as.factor(mat$x)
  freq_lab <- as.factor(c("Low (7.5-25 Hz)", "High (31-100 Hz)"))
  levels(freq_lab) <- freq_lab
  freqs <- lapply(mat$y, function(x) 
    lapply(x[[1]], function(y) y[[1]]))
  d <- data.frame(y = numeric(0), pow_lab = numeric(0), freq_lab = numeric(0))
  for (i in 1:length(freqs)) {
    for (j in 1:length(freqs[[i]])) {
      e <- data.frame(y = freqs[[i]][[j]], pow_lab = pow_lab[j], freq_lab = freq_lab[i])
      d <- rbind(d, e)
    }
  }
  d %>% 
    ggplot(aes(x=pow_lab, y=y)) +
    geom_boxplot(aes(col=monkey_name)) +
    scale_color_manual(values = c('Z'=zorin_color, 'J'=jaws_color)) +
    guides(col = F) +
    facet_wrap(~freq_lab) +
    labs(y = "Reaction Time (ms)", x = "Power bin")
}

plot_binned_coh <- function(name) {
  d <- load_mat(name)
  dist_edges <- d[['dist.edges']]
  legend_labels <- c()
  for (i in 1:(length(dist_edges)-1)) {
    legend_labels[i] <- TeX(sprintf('%.0f - %.0f $\\mu m$', dist_edges[i], dist_edges[i+1]))
  }
  idxs <- d$x < 100
  x <- d$x[idxs]
  y <- d$y[idxs,]
  std_err <- d$std.err[idxs,]
  
  std_dev_plot(x, y, std_err, legend_labels = legend_labels) +
    xlab('Frequency (Hz)') + ylab('Coherence') + 
    theme(legend.text = element_text(size=12), 
          legend.position = c(0.98,0.98),
          legend.justification = c('right', 'top'),
          legend.key.size = unit(0.1, units = 'in'))
}

plot_heatmap <- function(name, freq_cutoff=100, convert=F) {
  d <- load_mat(name, convert = convert)
  if (!convert) {
    days <- as.Date(unlist(lapply(d$x, function(x) x[[1]])))
    trial = 1:length(days)
  } else {
    trial <- d$x
  }
  freqs <- d$freqs
  vals <- d$y
  expand.grid(f=freqs, trial=trial) %>%
    mutate(z = as.numeric(vals),
           z = ifelse(is.infinite(z), NaN, z)) %>%
    filter(f <= freq_cutoff) %>%
    ggplot(aes(trial, f)) +
    geom_raster(aes(fill=z)) +
    scale_fill_gradientn(colors = inferno(1000), name = NULL) +
    labs(x = "Day number", y = "Frequency")
}
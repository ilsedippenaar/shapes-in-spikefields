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
    geom_line(aes(y = y, col = name), size=1) +
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
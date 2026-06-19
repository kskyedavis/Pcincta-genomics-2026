###########################
# ADJUSTED EEMS FUNCTIONS #
###########################

make_eems_plots2 <- function(mcmcpath, longlat = TRUE, dpi = 250,
                            add_grid = FALSE, col_grid = "#BBBBBB",
                            add_demes = FALSE, col_demes = "#000000",
                            add_outline = FALSE, col_outline = "#FFFFFF",
                            eems_colors = NULL, prob_level = 0.9,
                            m_colscale = NULL, q_colscale = NULL,
                            add_abline = FALSE) {
  check_mcmcpath_contents2(mcmcpath)
  func_params <- list(add_grid = add_grid, add_demes = add_demes,
                      add_outline = add_outline,
                      eems_colors = eems_colors, prob_level = prob_level,
                      col_grid = col_grid, col_demes = col_demes,
                      col_outline = col_outline,
                      m_colscale = m_colscale, q_colscale = q_colscale,
                      add_abline = add_abline)
  plot_params <- check_plot_params2(func_params)
  dimns <- read_dimns2(mcmcpath[1], longlat, dpi)
  plots <- list()
  
  p <- eems_contours2(mcmcpath, dimns, longlat, plot_params, is_mrates = TRUE)
  plots$mrates01 <- p[[1]]
  plots$mrates02 <- p[[2]]
  p <- eems_contours2(mcmcpath, dimns, longlat, plot_params, is_mrates = FALSE)
  plots$qrates01 <- p[[1]]
  plots$qrates02 <- p[[2]]
  
  dissimilarities <- pairwise_dist2(mcmcpath, longlat, plot_params)
  p <- plot_pairwise_dissimilarities_2(dissimilarities, add_abline)
  plots$rdist01 <- p[[1]]
  plots$rdist02 <- p[[2]]
  plots$rdist03 <- p[[3]]
  
  plots$pilogl01 <- plot_log_posterior2(mcmcpath)
  plots
}

eems_contours2 <- function(mcmcpath, dimns, longlat, plot_params, is_mrates) {
  if (is_mrates)
    message("Generate effective migration surface ",
            "(posterior mean of m rates). ",
            "See plots$mrates01 and plots$mrates02.")
  else
    message("Generate effective diversity surface ",
            "(posterior mean of q rates). ",
            "See plots$qrates01 and plots$qrates02.")
  zrates <- rep(0, dimns$nmrks)
  pr_gt0 <- rep(0, dimns$nmrks)
  pr_lt0 <- rep(0, dimns$nmrks)
  niters <- 0
  # Loop over each directory in mcmcpath to average the contour plots
  for (path in mcmcpath) {
    voronoi <- read_voronoi2(path, longlat, is_mrates)
    rslt <- tiles2contours(voronoi$tiles, voronoi$rates,
                           cbind(voronoi$xseed, voronoi$yseed),
                           dimns$marks, dimns$dist_metric)
    zrates <- zrates + rslt$zrates
    niters <- niters + rslt$niters
    pr_gt0 <- pr_gt0 + rslt$pr_gt0
    pr_lt0 <- pr_lt0 + rslt$pr_lt0
  }
  zrates <- zrates / niters
  pr_gt0 <- pr_gt0 / niters
  pr_lt0 <- pr_lt0 / niters
  p1 <- filled_eems_contour2(dimns, zrates, plot_params, is_mrates)
  p2 <- filled_prob_contour2(dimns, pr_gt0 - pr_lt0, plot_params, is_mrates)
  list(p1, p2)
}

plot_log_posterior2 <- function(mcmcpath) {
  message("Generate posterior probability trace. ",
          "See plots$pilog01.")
  rleid <- function(x) {
    r <- rle(x)
    rep(seq_along(r$lengths), r$lengths)
  }
  pl_df <- NULL
  for (path in mcmcpath) {
    pl <- read_matrix2(file.path(path, "mcmcpilogl.txt"))
    pl_df <- bind_rows(pl_df, as_data_frame(pl) %>% mutate(path))
  }
  pl_df <- pl_df %>%
    setNames(c("pi", "logl", "path")) %>%
    mutate(mcmcpath = factor(rleid(path))) %>%
    group_by(mcmcpath) %>%
    mutate(iter = row_number(), pilogl = pi + logl)
  ggplot(pl_df, aes(x = iter, y = pilogl, color = mcmcpath)) +
    geom_path() +
    labs(x = "MCMC iteration  (after burn-in and thinning)",
         y = "log posterior",
         title = "Have the MCMC chains converged?",
         subtitle = "If not, restart EEMS and/or increase numMCMCIter") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
}

plot_pairwise_dissimilarities_2 <- function(dissimilarities, add_abline) {
  message("Generate average dissimilarities within and between demes. ",
          "See plots$rdist01, plots$rdist02 and plots$rdist03.")
  p1 <- ggplot(dissimilarities$between %>% filter(size > 1),
               aes(fitted, obsrvd)) +
    geom_point(shape = 1) +
    theme_minimal() +
    labs(x = expression(paste("Fitted dissimilarity between demes  ",
                              Delta[alpha * beta], " - (",
                              Delta[alpha * alpha], " + ",
                              Delta[beta * beta], ") / 2")),
         y = expression(paste("Observed dissimilarity between demes  ",
                              D[alpha * beta], " - (",
                              D[alpha * alpha], " + ",
                              D[beta * beta], ") / 2")),
         title = expression(paste("Dissimilarities between pairs of ",
                                  "sampled demes (", alpha, ", ", beta, ")")),
         subtitle = paste("Singleton demes, if any, are excluded from this",
                          "plot (but not from EEMS)"))
  p2 <- ggplot(dissimilarities$within %>% filter(size > 1),
               aes(fitted, obsrvd)) +
    geom_point(shape = 1) +
    theme_minimal() +
    labs(x = expression(paste("Fitted dissimilarity within demes  ",
                              Delta[alpha * alpha])),
         y = expression(paste("Observed dissimilarity within demes ",
                              D[alpha * alpha])),
         title = expression(paste("Dissimilarities within sampled ",
                                  "demes ", alpha)),
         subtitle = paste("Singleton demes, if any, are excluded from ",
                          "this plot (but not from EEMS)"))
  p3 <- ggplot(dissimilarities$ibd %>% filter(size > 1),
               aes(fitted, obsrvd)) +
    geom_point(shape = 1) +
    theme_minimal() +
    labs(x = "Great circle distance between demes (km)",
         y = expression(paste("Observed dissimilarity between demes  ",
                              D[alpha * beta], " - (",
                              D[alpha * alpha], " + ",
                              D[beta * beta], ") / 2")),
         title = expression(paste("Dissimilarities between pairs of ",
                                  "sampled demes (", alpha, ", ", beta, ")")),
         subtitle = paste("Singleton demes, if any, are excluded from this",
                          "plot (but not from EEMS)"))
  if (add_abline) {
    p1 <- p1 + geom_smooth(method = "lm", se = FALSE)
    p2 <- p2 + geom_smooth(method = "lm", se = FALSE)
  }
  list(p1, p2, p3)
}
  
  
check_plot_params2 <- function(pars) {
    
    if (is.logical(pars$add_grid)) pars$add_grid <- pars$add_grid[1]
    else pars$add_grid <- FALSE
    if (is_color2(pars$col_grid)) pars$col_grid <- pars$col_grid[1]
    else pars$col_grid <- "#BBBBBB"
    
    if (is.logical(pars$add_outline)) pars$add_outline <- pars$add_outline[1]
    else pars$add_outline <- FALSE
    if (is_color2(pars$col_outline)) pars$col_outline <- pars$col_outline[1]
    else pars$col_outline <- "#EEEEEE"
    
    if (is.logical(pars$add_demes)) pars$add_demes <- pars$add_demes[1]
    else pars$add_demes <- FALSE
    if (is_color2(pars$col_demes)) pars$col_demes <- pars$col_demes[1]
    else pars$col_demes <- "#000000"
    if (is.logical(pars$add_seeds)) pars$add_seeds <- pars$add_seeds[1]
    else pars$add_seeds <- TRUE
    
    if (is.numeric(pars$m_colscale)) pars$m_colscale <- pars$m_colscale
    else pars$m_colscale <- c(-2.5, 2.5)
    if (is.numeric(pars$q_colscale)) pars$q_colscale <- pars$q_colscale
    else pars$q_colscale <- c(-0.1, 0.1)
    
    if (length(pars$eems_colors) < 2 || any(!is_color2(pars$eems_colors)))
      pars$eems_colors <- default_eems_colors2()
    
    if (is.null(pars$prob_level)) prob_level <- 0.9
    else prob_level <- pars$prob_level
    prob_level <- prob_level[prob_level > 0.5 & prob_level < 1]
    if (length(prob_level) != 1) prob_level <- 0.9
    pars$prob_level <- prob_level
    pars
  }

default_eems_colors2 <- function() {
    # Use the default DarkOrange to Blue color scheme, which combines
    # two color schemes from the `dichromat` package. These are based
    # on a collection of color schemes for scientific graphics:
    # See http://geog.uoregon.edu/datagraphics/color_scales.htm
    # To reproduce the default eems colors,
    # let oranges be dichromat::colorschemes$BluetoDarkOrange.12[12:7]
    # and blues be dichromat::colorschemes$BrowntoBlue.12[7:12]
    c("#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC",
      "#FBFBFB",
      "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99")
  }

is_color2 <- function(x) {
    if (is.null(x)) return(FALSE)
    sapply(x, function(x) {
      tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)
    })
  }
  
read_matrix2 <- function(file, ncol = 2) {
    stopifnot(file.exists(file))
    matrix(scan(file, what = numeric(), quiet = TRUE),
           ncol = ncol, byrow = TRUE)
  }
  
read_vector2 <- function(file) {
    stopifnot(file.exists(file))
    scan(file, what = numeric(), quiet = TRUE)
  }

get_dist_metric2 <- function(mcmcpath) {
    dist_metric <- "euclidean"
    lines <- tolower(readLines(file.path(mcmcpath, "eemsrun.txt")))
    for (line in lines) {
      if (grepl("\\s*distance\\s*=\\s*", line))
        dist_metric <- gsub("\\s*distance\\s*=\\s*(\\w+)", "\\1", line)
    }
    if (dist_metric != "euclidean" && dist_metric != "greatcirc")
      stop("eemsrun.txt should specify `euclidean` or `greatcirc` distance.")
    dist_metric
  }

read_dimns2 <- function(mcmcpath, longlat, nmrks = 100) {
    outer <- read_matrix2(file.path(mcmcpath, "outer.txt"))
    ipmap <- read_vector2(file.path(mcmcpath, "ipmap.txt"))
    demes <- read_matrix2(file.path(mcmcpath, "demes.txt"))
    edges <- read_matrix2(file.path(mcmcpath, "edges.txt"))
    dist_metric <- get_dist_metric2(mcmcpath)
    if (!longlat) {
      outer <- outer[, c(2, 1)]
      demes <- demes[, c(2, 1)]
    }
    xlim <- range(outer[, 1])
    ylim <- range(outer[, 2])
    aspect <- (diff(ylim) / diff(xlim)) / cos(mean(ylim) * pi / 180)
    aspect <- abs(aspect)
    if (aspect > 1) {
      nxmrks <- nmrks
      nymrks <- round(nxmrks * aspect)
    } else {
      nymrks <- nmrks
      nxmrks <- round(nymrks / aspect)
    }
    # Construct a rectangular "raster" of equally spaced pixels/marks
    xmrks <- seq(from = xlim[1], to = xlim[2], length = nxmrks)
    ymrks <- seq(from = ylim[1], to = ylim[2], length = nymrks)
    marks <- cbind(rep(xmrks, times = nymrks), rep(ymrks, each = nxmrks))
    # Exclude pixels that fall outside the habitat outline
    outer_poly <-
      sp::SpatialPolygons(list(Polygons(list(Polygon(outer, hole = FALSE)), "1")))
    marks <- sp::SpatialPoints(marks)[outer_poly, ]
    marks <- marks@coords
    outer <- as_data_frame(outer) %>% setNames(c("x", "y"))
    ipmap <- data_frame(id = ipmap) %>% count(id)
    demes <- as_data_frame(demes) %>% setNames(c("x", "y")) %>%
      mutate(id = row_number()) %>%
      left_join(ipmap) %>% arrange(id) %>%
      mutate(n = if_else(is.na(n), 0L, n))
    # edges <- bind_cols(demes[edges[, 1], ] %>% select(x, y),
    #                     demes[edges[, 2], ] %>% select(x, y)) %>%
    #    setNames(c("x", "y", "xend", "yend"))
    edges <- cbind(demes[edges[, 1], c("x", "y")],
                   demes[edges[, 2], c("x", "y")])
    colnames(edges) <- c("x", "y", "xend", "yend")
    edges <- as_tibble(edges)
    list(marks = marks, nmrks = nrow(marks), xlim = xlim, ylim = ylim,
         outer = outer, demes = demes, edges = edges,
         dist_metric = dist_metric)
}
  
filled_contour_graph2 <- function(p, dimns, plot_params) {
    if (plot_params$add_grid) {
      p <- p + geom_segment(data = dimns$edges,
                            aes(x = x, y = y, xend = xend, yend = yend),
                            color = plot_params$col_grid) +
        coord_fixed()
      
    }
    if (plot_params$add_demes) {
      p <- p + geom_point(data = dimns$demes %>% filter(n > 0),
                          aes(x = x, y = y, size = n), shape = 1,
                          color = plot_params$col_demes) +
        scale_size_continuous(guide = FALSE)
    }
    if (plot_params$add_outline) {
      p <- p + geom_path(data = dimns$outer, aes(x = x, y = y),
                         color = plot_params$col_outline)
    }
    p
}

check_mcmcpath_contents2 <- function(mcmcpath) {
  for (path in mcmcpath) {
    for (file in c("rdistJtDobsJ.txt", "rdistJtDhatJ.txt", "rdistoDemes.txt",
                   "mcmcmtiles.txt", "mcmcmrates.txt", "mcmcxcoord.txt",
                   "mcmcycoord.txt", "mcmcqtiles.txt", "mcmcqrates.txt",
                   "mcmcwcoord.txt", "mcmczcoord.txt", "mcmcpilogl.txt",
                   "outer.txt", "demes.txt", "edges.txt", "ipmap.txt",
                   "eemsrun.txt")) {
      if (!file.exists(file.path(path, file)))
        stop("Each EEMS output folder should include ", file)
    }
  }
}

read_voronoi2 <- function(mcmcpath, longlat, is_mrates) {
  if (is_mrates) {
    rates <- read_vector2(file.path(mcmcpath, "mcmcmrates.txt"))
    tiles <- read_vector2(file.path(mcmcpath, "mcmcmtiles.txt"))
    xseed <- read_vector2(file.path(mcmcpath, "mcmcxcoord.txt"))
    yseed <- read_vector2(file.path(mcmcpath, "mcmcycoord.txt"))
  } else {
    rates <- read_vector2(file.path(mcmcpath, "mcmcqrates.txt"))
    tiles <- read_vector2(file.path(mcmcpath, "mcmcqtiles.txt"))
    xseed <- read_vector2(file.path(mcmcpath, "mcmcwcoord.txt"))
    yseed <- read_vector2(file.path(mcmcpath, "mcmczcoord.txt"))
  }
  if (!longlat) {
    tempi <- xseed
    xseed <- yseed
    yseed <- tempi
  }
  list(rates = log10(rates), tiles = tiles, xseed = xseed, yseed = yseed)
}

get_dist_metric2 <- function(mcmcpath) {
  dist_metric <- "euclidean"
  lines <- tolower(readLines(file.path(mcmcpath, "eemsrun.txt")))
  for (line in lines) {
    if (grepl("\\s*distance\\s*=\\s*", line))
      dist_metric <- gsub("\\s*distance\\s*=\\s*(\\w+)", "\\1", line)
  }
  if (dist_metric != "euclidean" && dist_metric != "greatcirc")
    stop("eemsrun.txt should specify `euclidean` or `greatcirc` distance.")
  dist_metric
}

theme_void2 <- function() {
  theme(line = element_blank(), rect = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size = rel(0.8)),
        legend.title = element_text(hjust = 0),
        legend.text.align = 1)
}

filled_contour_rates2 <- function(z, dimns) {
  w <- cbind(dimns$marks, z)
  colnames(w) <- c("x", "y", "z")
  ggplot(as_data_frame(w), aes(x = x, y = y)) +
    geom_tile(aes(fill = z)) +
    scale_x_continuous(limits = dimns$xlim) +
    scale_y_continuous(limits = dimns$ylim) +
    coord_quickmap() +
    theme_void2() +
    theme(legend.text.align = 1)
}

filled_eems_contour2 <- function(dimns, zmean, plot_params, is_mrates) {
  if (is_mrates) {
    title <- "log(m)"
    colscale <- plot_params$m_colscale
  } else {
    title <- "log(q)"
    colscale <- plot_params$q_colscale
  }
  limits <- range(zmean, colscale, na.rm = TRUE, finite = TRUE)
  p <- filled_contour_rates2(zmean, dimns)
  p <- filled_contour_graph2(p, dimns, plot_params) +
    scale_fill_gradientn(colors = plot_params$eems_colors,
                         limits = limits, name = title)
  p
}

filled_prob_contour2 <- function(dimns, probs, plot_params, is_mrates) {
  probs <- (probs + 1) / 2
  probs[probs < 0] <- 0
  probs[probs > 1] <- 1
  if (is_mrates) r <- "m" else r <- "q"
  breaks <- c(1 - plot_params$prob_level, plot_params$prob_level)
  labels <- c(paste0("P{log(", r, ") < 0} = ", plot_params$prob_level),
              paste0("P{log(", r, ") > 0} = ", plot_params$prob_level))
  p <- filled_contour_rates2(probs, dimns)
  p <- filled_contour_graph2(p, dimns, plot_params) +
    geom_contour(aes(z = z), breaks = breaks, color = "white") +
    scale_fill_gradientn(colors = default_eems_colors2(),
                         limits = c(0, 1), name = "",
                         breaks = breaks, labels = labels)
  p
}

decompose_distances2 <- function(diffs, sizes = NULL) {
  # Diffs can have NAs on the main diagonal; these elements correspond to demes
  # with a single observation. For such deme a, no dissimilarities between
  # two distinct individuals are observed. I approximate diffs(a,a) with the
  # average diffs(b,b) computed across demes b with multiple samples.
  if (!is.null(sizes))
    diag(diffs)[sizes < 2] <- mean(diag(diffs)[sizes >= 2])
  within <- diag(diffs)
  selfsim <- matrix(within, nrow(diffs), ncol(diffs))
  between <- diffs - (selfsim + t(selfsim)) / 2
  between <- between[upper.tri(between, diag = FALSE)]
  list(within = within, between = between)
}

geo_distm2 <- function(coord, longlat, plot_params) {
  if (!longlat) coord <- coord[, c(2, 1)]
  dist <- sp::spDists(coord, longlat = TRUE)
  dist <- dist[upper.tri(dist, diag = FALSE)]
  dist
}

pairwise_dist2 <- function(mcmcpath, longlat, plot_params) {
  # List of observed demes, with number of samples taken collected Each row
  # specifies: x coordinate, y coordinate, n samples
  obs_demes <- read_matrix2(file.path(mcmcpath[1], "rdistoDemes.txt"), ncol = 3)
  sizes <- obs_demes[, 3]
  npops <- nrow(obs_demes)
  demes <- seq(npops)
  diffs_obs <- matrix(0, npops, npops)
  diffs_hat <- matrix(0, npops, npops)
  for (path in mcmcpath) {
    tempi <- read_matrix2(file.path(path, "rdistoDemes.txt"), ncol = 3)
    if (sum(dim(obs_demes) != dim(tempi)) || sum(obs_demes != tempi)) {
      message("EEMS results for at least two different population grids. ",
              "Plot pairwise dissimilarity for each grid separately.")
      return(list(between = data_frame(), within = data_frame(),
                  ibd = data_frame()))
    }
    diffs_obs <- diffs_obs +
      as.matrix(read.table(file.path(path, "rdistJtDobsJ.txt")))
    diffs_hat <- diffs_hat +
      as.matrix(read.table(file.path(path, "rdistJtDhatJ.txt")))
  }
  diffs_obs <- diffs_obs / length(mcmcpath)
  diffs_hat <- diffs_hat / length(mcmcpath)
  alpha <- matrix(demes, nrow = npops, ncol = npops)
  beta <- t(alpha)
  tempi <- matrix(sizes, npops, npops)
  smaller_deme <- pmin(tempi, t(tempi))
  smaller_deme <- smaller_deme[upper.tri(smaller_deme, diag = FALSE)]
  alpha <- alpha[upper.tri(alpha, diag = FALSE)]
  beta <- beta[upper.tri(beta, diag = FALSE)]
  # Under pure isolation by distance, we expect the genetic dissimilarities
  # between demes increase with the geographic distance separating them
  dist <- geo_distm2(obs_demes[, 1:2], longlat, plot_params)
  if (sum(sizes > 1) < 2) {
    message("There should be at least two observed demes ",
            "to plot pairwise dissimilarities")
    return(NULL)
  }
  bw_obs <- decompose_distances2(diffs_obs, sizes)
  bw_hat <- decompose_distances2(diffs_hat)
  b_component <- data_frame(alpha_x = obs_demes[, 1][alpha],
                            alpha_y = obs_demes[, 2][alpha],
                            beta_x = obs_demes[, 1][beta],
                            beta_y = obs_demes[, 2][beta],
                            fitted = bw_hat$between,
                            obsrvd = bw_obs$between,
                            size = smaller_deme)
  w_component <- data_frame(alpha_x = obs_demes[, 1][demes],
                            alpha_y = obs_demes[, 2][demes],
                            fitted = bw_hat$within,
                            obsrvd = bw_obs$within,
                            size = sizes)
  g_component <- data_frame(alpha_x = obs_demes[, 1][alpha],
                            alpha_y = obs_demes[, 2][alpha],
                            beta_x = obs_demes[, 1][beta],
                            beta_y = obs_demes[, 2][beta],
                            fitted = dist,
                            obsrvd = bw_obs$between,
                            size = smaller_deme)
  list(between = b_component, within = w_component, ibd = g_component)
}

load_required_packages2 <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE))
      stop("The ", package, " package is required. ",
           "Please install it first.")
    else
      message("Loading ", package, ".")
  }
}





plot_voronoi_tiles2 <- function(mcmcpath, longlat, num_draws = 1,
                               add_seeds = TRUE, eems_colors = NULL,
                               m_colscale = NULL, q_colscale = NULL) {
  check_mcmcpath_contents2(mcmcpath)
  load_required_packages2("deldir")
  func_params <- list(add_seeds = add_seeds, eems_colors = eems_colors,
                      m_colscale = m_colscale, q_colscale = q_colscale)
  plot_params <- check_plot_params2(func_params)
  
  plots <- list()
  for (draw in seq(num_draws)) {
    p <- random_voronoi_diagram_2(sample(mcmcpath, 1), longlat,
                                 plot_params, num_draws, is_mrates = TRUE)
    plots[[paste0("mtiles", draw)]] <- p
  }
  for (draw in seq(num_draws)) {
    p <- random_voronoi_diagram_2(sample(mcmcpath, 1), longlat,
                                 plot_params, num_draws, is_mrates = FALSE)
    plots[[paste0("qtiles", draw)]] <- p
  }
  plots
}

random_voronoi_diagram_2 <- function(mcmcpath, longlat, plot_params,
                                    num_draws = 1, is_mrates = TRUE) {
  message("Plotting Voronoi tessellation of estimated effective rates")
  if (is_mrates) {
    colscale <- plot_params$m_colscale
    title <- "log(m)"
  } else {
    colscale <- plot_params$q_colscale
    title <- "log(q)"
  }
  dimns <- read_dimns2(mcmcpath, longlat)
  voronoi <- read_voronoi2(mcmcpath, longlat, is_mrates)
  # Choose one saved posterior draw at random
  iter <- sample(seq_along(voronoi$tiles), 1)
  message(paste0("Draw ", iter))
  # Jump over stored parameters for draws 1 to (iter - 1)
  skip <- sum(voronoi$tiles[iter:1][-1])
  now_tiles <- voronoi$tiles[iter]
  now_rates <- voronoi$rates[(skip + 1):(skip + now_tiles)]
  now_xseed <- voronoi$xseed[(skip + 1):(skip + now_tiles)]
  now_yseed <- voronoi$yseed[(skip + 1):(skip + now_tiles)]
  outer <- as_data_frame(dimns$outer) %>% setNames(c("x", "y"))
  seeds <- data_frame(x = now_xseed, y = now_yseed)
  poly <- NULL
  if (now_tiles == 1) {
    poly <- data_frame(x = rep(dimns$xlim, each = 2),
                       y = rep(dimns$ylim, times = 2),
                       rate = now_rates)
  } else {
    tile_list <-
      deldir::tile.list(deldir::deldir(now_xseed, now_yseed,
                                       rw = c(dimns$xlim, dimns$ylim)))
    for (t in seq(now_tiles))
      poly <- bind_rows(poly, as_data_frame(tile_list[[t]][c("x", "y")]) %>%
                          mutate(id = t, rate = now_rates[t]))
  }
  limits <- range(poly$rate, colscale, na.rm = TRUE, finite = TRUE)
  p <- ggplot(poly, aes(x = x, y = y)) +
    geom_polygon(aes(fill = rate, group = id), color = "white") +
    geom_path(data = outer) +
    scale_fill_gradientn(colors = plot_params$eems_colors,
                         limits = limits, name = title) +
    theme_void()
  if (plot_params$add_seeds) p <- p + geom_point(data = seeds, shape = 1)
  p
}

plot_resid_heatmap2 <- function(datapath, mcmcpath,
                               hm_colors = NULL, hm_scale = NULL) {
  load_required_packages2(c("Matrix", "scales"))
  check_mcmcpath_contents2(mcmcpath)
  p <- plot_resid_heatmap_2(datapath, mcmcpath, hm_colors, hm_scale)
  list(residhm = p)
}

plot_resid_heatmap_2 <- function(datapath, mcmcpath,
                                hm_colors, hm_scale) {
  message("Generate heatmap of n-by-n matrix of residuals ",
          "(observed - fitted). See plots$residhm.")
  stopifnot(file.exists(paste0(datapath, ".diffs")))
  diffs <- as.matrix(read.table(paste0(datapath, ".diffs")))
  n <- nrow(diffs)
  m <- length(mcmcpath)
  delta <- matrix(0, n, n)
  for (path in mcmcpath) {
    ipmap <- read_vector2(file.path(path, "ipmap.txt"))
    indicator <- as.matrix(Matrix::spMatrix(n, n_distinct(ipmap),
                                            i = seq(n), j = ipmap,
                                            x = rep(1, n)))
    d_hat <- as.matrix(read.table(file.path(path, "rdistJtDhatJ.txt")))
    delta <- delta + indicator %*% d_hat %*% t(indicator)
  }
  resid <- abs(diffs - (delta / m))
  diag(resid) <- NA
  hm_limits <- range(resid, hm_scale, na.rm = TRUE, finite = TRUE)
  colnames(resid) <- seq_len(ncol(resid))
  tiles <- resid %>%
    as.data.frame() %>%
    as_data_frame() %>%
    mutate(row = row_number()) %>%
    gather(col, value, - row) %>%
    mutate(col = as.integer(col))
  ggplot(tiles, aes(x = col, y = row, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = hm_colors, limits = hm_limits,
                         name = "|err|", na.value = "white") +
    coord_equal() +
    theme_void()
}




plot_population_grid2 <- function (mcmcpath, longlat, add_demes = TRUE, col_demes = "black",
                                    col_outline = "black", col_grid = "gray50"){
  func_params <- list(add_outline = TRUE, add_grid = TRUE, 
                      add_demes = add_demes, col_demes = col_demes, col_grid = col_grid, 
                      col_outline = col_outline)
  plot_params <- check_plot_params2(func_params)
  dimns <- read_dimns2(mcmcpath, longlat = T)
  p <- ggplot() + coord_quickmap() + theme_void2()
  p <- filled_contour_graph2(p, dimns, plot_params)
  out <- list(popgrid = p)
  return(out)
}




tiles2contours_standardize <- function(tiles, rates, seeds, marks, distm) {
  .Call('_reemsplots2_tiles2contours_standardize', PACKAGE = 'reemsplots2', tiles, rates, seeds, marks, distm)
}

tiles2contours <- function(tiles, rates, seeds, marks, distm) {
  .Call('_reemsplots2_tiles2contours', PACKAGE = 'reemsplots2', tiles, rates, seeds, marks, distm)
}

#' @title Plot distances to MAP
#' @description Plot distances to MAP for each prior type
#' @param results A list of results from the simulation containing the distances to MAP
#' @param prior_types A character vector of prior types
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot. By default, it is 600
#' @param width An integer of the width of the plot. By default, it is 10
#' @param height An integer of the height of the plot. By default, it is 10
#' @return A ggplot object
#' @export

plt_distances_to_MAP <- function(results, prior_types, path = NULL, dpi = 600, width = 10, height = 10) {
    all_distances <- data.frame()
    for (prior_type in prior_types) {
        distances_to_MAP <- lapply(results, function(x) x[[prior_type]]$distance_vector)
        distances_to_MAP <- do.call(rbind, distances_to_MAP)
        distances_to_MAP <- as.data.frame(distances_to_MAP)
        distances_to_MAP$iteration <- seq_len(nrow(distances_to_MAP))
        distances_to_MAP <- reshape2::melt(distances_to_MAP, id.vars = "iteration")
        distances_to_MAP$prior_type <- prior_type
        all_distances <- rbind(all_distances, distances_to_MAP)
    }

    p <- ggplot(all_distances) +
        stat_ecdf(aes(value, color = prior_type)) +
        facet_wrap(~variable)
    if (!is.null(path)) {
        ggsave(path, dpi = dpi, width = width, height = height)
    }
    print(p)
}

#' @title Mean distances to MAP and standard deviation of Gaussian approximation around the median
#' @description Gets the mean distances to MAP and the standard deviation of the Gaussian approximation around the median for each prior type
#' @param results A list of results from the simulation containing the distances to MAP
#' @param prior_types A character vector of prior types
#' @return A list of vectors containing the mean distances to MAP and the standard deviation of the Gaussian approximation around the median for each prior type
#' @export

mean_distance_to_MAP_and_std_dev_of_Gaussian_approximation <- function(results, prior_types) {
    mean_distances <- lapply(prior_types, function(prior_type) {
        mean_distances <- sapply(1:5, function(i) {
            all_distances <- sapply(seq_along(results), function(j) {
                results[[j]][[prior_type]]$distance_vector[i]
            })
            mean(all_distances)
        })
        names(mean_distances) <- parameter_names

        std_dev_estimates_Gaussian_median <- do.call(rbind, lapply(results, function(x) x[[prior_type]]$std_dev_estimates_Gaussian_median))
        mean_std_dev <- colMeans(std_dev_estimates_Gaussian_median)
        list(mean_distances = mean_distances, mean_std_dev = mean_std_dev)
    })
    mean_distances
}


#' @title Plot credible interval lengths and get mean lengths
#' @description Plot credible interval lengths and get mean lengths for each prior type and approximation type
#' @param results A list of results from the simulation containing the credible intervals
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param parameter_names A character vector of parameter names
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @return A list of data frames containing the mean lengths for each prior type and approximation type
#' @export

plt_CI_lengths_and_get_mean_lengths <- function(results, prior_types, approximation_types, parameter_names, path = NULL, dpi = 600, width = 10, height = 10) {
    lengths_df <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types, function(approximation_type) {
            lengths <- lapply(parameter_names, function(parameter_name) {
                all_lengths <- sapply(seq_along(results), function(j) {
                    length <- diff(results[[j]][[prior_type]]$credible_intervals[[approximation_type]][parameter_name, ])
                    length
                })
                all_lengths
            })
            lengths <- do.call(cbind, lengths)
            colnames(lengths) <- parameter_names
            lengths
        })
    })
    mean_lengths_df <- lapply(lengths_df, function(prior_type) {
        lapply(prior_type, function(approximation_type) {
            colMeans(approximation_type)
        })
    })
    all_lengths <- data.frame()
    for (prior_type in prior_types) {
        for (approximation_type in approximation_types) {
            lengths <- lengths_df[[prior_type]][[approximation_type]]
            lengths <- as.data.frame(lengths)
            lengths$iteration <- seq_len(nrow(lengths))
            lengths <- reshape2::melt(lengths, id.vars = "iteration") # Necessary to use ggplot as it expects a data frame in long format
            lengths$prior_type <- prior_type
            lengths$approximation_type <- approximation_type
            all_lengths <- rbind(all_lengths, lengths)
        }
    }

    p <- ggplot(all_lengths) +
        stat_ecdf(aes(value, color = prior_type, linetype = approximation_type)) +
        facet_wrap(~variable)
    if (!is.null(path)) {
        ggsave(path, dpi = dpi, width = width, height = height)
    }
    print(p)
    mean_lengths_df
}

#' @title Plot the frequency of the true parameter being in the credible interval
#' @description Plot the frequency of the true parameter being in the credible interval for each prior type and approximation type
#' @param results A list of results from the simulation containing the credible intervals
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param parameter_names A character vector of parameter names
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @export

plt_frequency_true_parameter_in_CI <- function(results = results, prior_types = prior_types, approximation_types = approximation_types, parameter_names = parameter_names, path = NULL, dpi = 600, width = 10, height = 10) {
    within_ci <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types, function(approximation_type) {
            rowMeans(sapply(results, function(x) x[[prior_type]][["true_parameter_within_c_interval"]][[approximation_type]]))
        })
    })
    # Convert the list to a data frame
    within_ci <- do.call(rbind, lapply(names(within_ci), function(prior_type) {
        do.call(rbind, lapply(names(within_ci[[prior_type]]), function(approximation_type) {
            data.frame(
                prior_type = prior_type,
                approximation_type = approximation_type,
                parameter = factor(names(within_ci[[prior_type]][[approximation_type]]), levels = c("log_kappa", "v1", "v2", "log_sigma_u", "log_sigma_epsilon")),
                value = c(within_ci[[prior_type]][[approximation_type]])
            )
        }))
    }))

    # Plot the points using ggplot
    p <- ggplot(within_ci, aes(x = approximation_type, y = value, color = prior_type)) +
        geom_point() +
        geom_text(aes(label = round(value, 2)), vjust = -0.5) +
        facet_wrap(~parameter) +
        labs(x = "Approximation Type", y = "Value") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    if (!is.null(path)) {
        ggsave(path, dpi = dpi, width = width, height = height)
    }
    print(p)
}

#' @title Plot KL divergences
#' @description Plot KL divergences for each prior type and approximation type
#' @param results A list of results from the simulation containing the KL divergences
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @export

plt_KL_and_get_mean_KL <- function(results, prior_types, approximation_types, path = NULL, dpi = 600, width = 10, height = 10) {
    KL <- lapply(prior_types, function(prior_type) {
        lapply(KL_approx_types, function(approximation_type) {
            sapply(results, function(x) {
                kl_values <- x[[prior_type]]$importance[[paste0("KL_divergence_", approximation_type, "_Gaussian_median")]]
                kl_values <- replace(kl_values, is.na(kl_values), Inf) # Replace NA values with Inf
                kl_values
            })
        })
    })

    all_KL <- data.frame()
    for (prior_type in prior_types) {
        for (approximation_type in approximation_types) {
            KL_divergences <- KL[[prior_type]][[approximation_type]]
            KL_divergences <- as.data.frame(KL_divergences)
            KL_divergences$iteration <- seq_len(nrow(KL_divergences))
            KL_divergences <- reshape2::melt(KL_divergences, id.vars = "iteration") # Necessary to use ggplot as it expects a data frame in long format
            KL_divergences$prior_type <- prior_type
            KL_divergences$approximation_type <- approximation_type
            all_KL <- rbind(all_KL, KL_divergences)
        }
    }

    p <- ggplot(all_KL) +
        stat_ecdf(aes(value, color = prior_type, linetype = approximation_type)) +
        facet_wrap(~variable)

    if (!is.null(path)) {
        ggsave(path, dpi = dpi, width = width, height = height)
    }
    print(p)
    KL_Gaussian_median_mean <- lapply(prior_types, function(prior_type) {
        lapply(KL_approx_types, function(approximation_type) {
            mean(KL[[prior_type]][[approximation_type]])
        })
    })
    KL_Gaussian_median_mean
}

#' @title Plot probabilities
#' @description Plot probabilities for each prior type and approximation type that the marginal posterior is smaller than the true parameter value
#' @param results A list of results from the simulation containing the probabilities
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param parameter_names A character vector of parameter names
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @export

plt_probabilities <- function(results, prior_types, approximation_types, parameter_names, path = NULL, dpi = 600, width = 10, height = 10) {
    all_probabilities <- data.frame()
    for (prior_type in prior_types) {
        for (approximation_type in approximation_types) {
            probabilities <- sapply(results, function(x) x[[prior_type]][["importance"]][[paste0("probabilities_", approximation_type)]])
            for (i in seq_along(parameter_names)) {
                df <- data.frame(prob = unlist(probabilities[i, ]), parameter = parameter_names[[i]], prior = prior_type, approximation = approximation_type)
                all_probabilities <- rbind(all_probabilities, df)
            }
        }
    }

    p <- ggplot(all_probabilities) +
        stat_ecdf(aes(prob, col = prior, linetype = approximation)) +
        geom_abline(slope = 1, intercept = 0, color = "red") + # Add this line
        facet_wrap(~parameter)

    if (!is.null(path)) {
        ggsave(path, dpi = dpi, width = width, height = height)
    }
    print(p)
}

#' @title Plot KS statistics and p-values
#' @description Plot KS statistics and p-values for each prior type and approximation type
#' @param results A list of results from the simulation containing the probabilities
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param parameter_names A character vector of parameter names
#' @param path1 A character string of the path to save the plot for the KS statistics
#' @param path2 A character string of the path to save the plot for the p-values
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @export

plt_KS <- function(results, prior_types, approximation_types, parameter_names, path1 = NULL, path2 = NULL, dpi = 600, width = 10, height = 10) {
    all_probabilities <- data.frame()
    for (prior_type in prior_types) {
        for (approximation_type in approximation_types) {
            probabilities <- sapply(results, function(x) x[[prior_type]][["importance"]][[paste0("probabilities_", approximation_type)]])
            for (i in seq_along(parameter_names)) {
                df <- data.frame(prob = unlist(probabilities[i, ]), parameter = parameter_names[[i]], prior = prior_type, approximation = approximation_type)
                all_probabilities <- rbind(all_probabilities, df)
            }
        }
    }
    KS_results <- data.frame()

    for (i in seq_along(parameter_names)) {
        # Get the probabilities for the current parameter
        probabilities <- all_probabilities[all_probabilities$parameter == parameter_names[[i]], ]
        # Calculate the KS statistic for each prior and approximation type
        for (prior_type in prior_types) {
            for (approximation_type in approximation_types) {
                KS_result <- ks.test(probabilities[probabilities$prior == prior_type & probabilities$approximation == approximation_type, ]$prob, "punif")
                # Add the KS statistic and p-value for the current parameter to the data frame
                KS_results <- rbind(KS_results, data.frame(parameter = parameter_names[[i]], prior = prior_type, approximation = approximation_type, statistic = KS_result$statistic, p_value = KS_result$p.value))
            }
        }
    }

    p1 <- ggplot(KS_results) +
        geom_point(aes(x = parameter, y = statistic, color = prior, shape = approximation))

    p2 <- ggplot(KS_results) +
        geom_point(aes(x = parameter, y = p_value, color = prior, shape = approximation))

    if (!is.null(path1)) {
        ggsave(path1, dpi = dpi, width = width, height = height)
    }
    if (!is.null(path2)) {
        ggsave(path2, dpi = dpi, width = width, height = height)
    }
    print(p1)
    print(p2)
}

#' @title Plot complexity and get mean complexity
#' @description Plot complexity and get mean complexity for each prior type and approximation type
#' @param results A list of results from the simulation containing the complexity
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @return A list of vectors containing the mean complexity for each prior type and approximation type
#' @export

plt_complexity_and_get_mean_complexity <- function(results, prior_types, approximation_types, path = NULL, dpi = 600, width = 10, height = 10) {
    complexity <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types[2:3], function(approximation_type) {
            sapply(results, function(x) {
                complexity <- x[[prior_type]]$importance[[paste0("complexity_", approximation_type)]]
                complexity
            })
        })
    })

    all_complexity <- data.frame()
    for (prior_type in prior_types) {
        for (approximation_type in approximation_types[2:3]) {
            complexity_values <- complexity[[prior_type]][[approximation_type]]
            complexity_values <- as.data.frame(complexity_values)
            complexity_values$iteration <- seq_len(nrow(complexity_values))
            complexity_values <- reshape2::melt(complexity_values, id.vars = "iteration") # Necessary to use ggplot as it expects a data frame in long format
            complexity_values$prior_type <- prior_type
            complexity_values$approximation_type <- approximation_type
            all_complexity <- rbind(all_complexity, complexity_values)
        }
    }

    p <- ggplot(all_complexity) +
        stat_ecdf(aes(value, color = prior_type, linetype = approximation_type)) +
        facet_wrap(~variable)

    if (!is.null(path)) {
        ggsave(path, dpi = dpi, width = width, height = height)
    }
    print(p)

    complexity_mean <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types[2:3], function(approximation_type) {
            mean(complexity[[prior_type]][[approximation_type]])
        })
    })
    complexity_mean
}

#' @title Plot k diagnostics
#' @description Plot k diagnostics for each prior type
#' @param results A list of results from the simulation containing the k diagnostics
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @export

plt_k_diagnostics <- function(results, prior_types, path = NULL, dpi = 600, width = 10, height = 10) {
    all_k_diagnostics <- data.frame(matrix(ncol = length(prior_types), nrow = length(results)))
    colnames(all_k_diagnostics) <- prior_types
    for (prior_type in prior_types) {
        all_k_diagnostics[[prior_type]] <- sapply(results, function(x) x[[prior_type]]$importance$k_diagnostic)
    }
    all_k_diagnostics <- as.data.frame(all_k_diagnostics)
    all_k_diagnostics$iteration <- seq_len(nrow(all_k_diagnostics))
    all_k_diagnostics <- reshape2::melt(all_k_diagnostics, id.vars = "iteration")

    p <- ggplot(all_k_diagnostics) +
        stat_ecdf(aes(value, color = variable))

    if (!is.null(path)) {
        ggsave(path, dpi = dpi, width = width, height = height)
    }
    print(p)
}

plt_weights_cdf <- function(results, prior_types, path = NULL, dpi = 600, width = 10, height = 10) {
    all_weights <- list()
    for (prior_type in prior_types) {
        weights_unsmoothed <- c(unlist(sapply(results, function(x) x[[prior_type]]$importance$log_unnormalized_weights)))
        weights_smoothed <- c(unlist(sapply(results, function(x) x[[prior_type]]$importance$log_unnormalized_weights_smoothed)))
        weights <- data.frame(weights_unsmoothed = weights_unsmoothed, weights_smoothed = weights_smoothed)
    }
    p <- ggplot(weights) +
        stat_ecdf(aes(weights_unsmoothed, color = "Unsmoothed")) +
        stat_ecdf(aes(weights_smoothed, color = "Smoothed")) +
        labs(x = "Log weight", y = "Cumulative density") +
        theme(legend.position = "bottom")

    if (!is.null(path)) {
        ggsave(path, dpi = dpi, width = width, height = height)
    }
    print(p)
}


## Example for j=1
j <- 1
hist(results_not_pc[[j]]$pc$importance$log_unnormalized_weights_smoothed, main = "Log weights", xlab = "Log weight")
hist(results_not_pc[[j]]$pc$importance$log_unnormalized_weights, main = "Log weights", xlab = "Log weight")

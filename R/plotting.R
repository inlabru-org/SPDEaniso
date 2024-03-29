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

plt_distances_to_MAP <- function(results, prior_types, path = NULL, dpi = 100, width = 10, height = 10) {
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

plt_CI_lengths_and_get_mean_lengths <- function(results, prior_types, approximation_types, parameter_names, path = NULL, dpi = 100, width = 10, height = 10) {
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

plt_frequency_true_parameter_in_CI <- function(results = results, prior_types = prior_types, approximation_types = approximation_types, parameter_names = parameter_names, path = NULL, dpi = 100, width = 10, height = 10) {
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

plt_KL_and_get_mean_KL <- function(results, prior_types, approximation_types, path = NULL, dpi = 100, width = 10, height = 10) {
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

plt_probabilities <- function(results, prior_types, approximation_types, parameter_names, path = NULL, dpi = 100, width = 10, height = 10) {
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
    all_probabilities$parameter <- factor(all_probabilities$parameter, levels = parameter_names)
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

plt_KS <- function(results, prior_types, approximation_types, parameter_names, path1 = NULL, path2 = NULL, dpi = 100, width = 10, height = 10) {
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
    KS_results$parameter <- factor(KS_results$parameter, levels = parameter_names)

    p1 <- ggplot(KS_results) +
        geom_point(aes(x = parameter, y = statistic, color = prior, shape = approximation))
    if (!is.null(path1)) {
        ggsave(path1, dpi = dpi, width = width, height = height)
    }

    p2 <- ggplot(KS_results) +
        geom_point(aes(x = parameter, y = p_value, color = prior, shape = approximation))
    if (!is.null(path2)) {
        ggsave(path2, dpi = dpi, width = width, height = height)
    }
    print(p1)
    print(p2)
    # We only get the ones using last approximation type
    KS_results <- KS_results[KS_results$approximation == approximation_types[length(approximation_types)], ]
    # Since we only have one approximation type, we can eliminate it from the names
    KS_results$approximation <- NULL
    aa <- KS_results %>%
        mutate(result = paste0("statistic: ", round(statistic, 2), ", p_value: ", round(p_value, 2)))

    # Convert the table from long format to wide format
    wide_table <- aa %>%
        pivot_wider(names_from = "parameter", values_from = "result")


    # Convert the wide table to a LaTeX table
    latex_table <- xtable(KS_results)

    # Print the LaTeX table
    print(latex_table, type = "latex", include.rownames = FALSE)
    KS_results
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

plt_complexity_and_get_mean_complexity <- function(results, prior_types, approximation_types, path = NULL, dpi = 100, width = 10, height = 10) {
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

plt_k_diagnostics <- function(results, prior_types, path = NULL, dpi = 100, width = 10, height = 10) {
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

plt_weights_cdf <- function(results, prior_types, path = NULL, dpi = 100, width = 10, height = 10) {
    all_weights <- list()
    for (prior_type in prior_types) {
        weights_unsmoothed <- c(unlist(sapply(results, function(x) x[[prior_type]]$importance$log_unnormalized_weights)))
        weights_smoothed <- c(unlist(sapply(results, function(x) x[[prior_type]]$importance$log_unnormalized_weights_smoothed)))
        weights <- data.frame(prior_type = prior_type, weights_unsmoothed = weights_unsmoothed, weights_smoothed = weights_smoothed)
        all_weights[[prior_type]] <- weights
    }
    all_weights_df <- do.call(rbind, all_weights)

    # Reshape the data to long format
    all_weights_long <- pivot_longer(all_weights_df, c(weights_unsmoothed, weights_smoothed), names_to = "weight_type", values_to = "weight")

    p <- ggplot(all_weights_long) +
        stat_ecdf(aes(weight, color = prior_type, linetype = weight_type)) +
        labs(x = "Log weight", y = "Cumulative density") +
        theme(legend.position = "bottom") +
        xlim(c(-5, 0))

    if (!is.null(path)) {
        ggsave(path, dpi = dpi, width = width, height = height)
    }
    print(p)
}

prior_posterior_plotter <- function(theta_fixed = map_pc$par, log_priors = log_priors,
                                    log_posteriors = log_posteriors, l = 2, n_points = 10, n_parameters_to_plot = 3, path = NULL) {
    function_types <- list(prior = "prior", posterior = "posterior", path = NULL)
    ########## UNNORMALIZED Gaussian_median APPROXIMATION TO THE POSTERIOR############

    ### UNNORMALIZED LOG FUNCTION SO THEY ALL START AT 0###
    unnormalize_prior_and_posterior <- function(log_prior) {
        function(log_kappa, v1, v2, log_sigma_u, log_sigma_epsilon) {
            log_prior(
                log_kappa = log_kappa, v = c(v1, v2),
                log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
            ) - log_prior(
                log_kappa = theta_fixed[1], v = theta_fixed[2:3],
                log_sigma_u = theta_fixed[4], log_sigma_epsilon = theta_fixed[5]
            )
        }
    }

    unnormalized_priors <- lapply(log_priors, unnormalize_prior_and_posterior)
    unnormalized_posteriors <- lapply(log_posteriors, unnormalize_prior_and_posterior)
    unnormalized_priors_and_posteriors <- list(prior = unnormalized_priors, posterior = unnormalized_posteriors)

    # Restricting the functions to one parameter
    restricting_function_to_one_parameter <- function(f, x0) {
        f_list <- lapply(seq_along(x0)[1:n_parameters_to_plot], function(i) {
            function(x) {
                x0_copy <- x0
                x0_copy[i] <- x
                unname(do.call(f, as.list(x0_copy)))
            }
        })
        names(f_list) <- names(x0[1:n_parameters_to_plot])
        f_list
    }
    restricted_priors_and_posteriors <- lapply(unnormalized_priors_and_posteriors, function(f) {
        lapply(f, restricting_function_to_one_parameter, theta_fixed)
    })


    # Getting data for plotting
    partitions <- lapply(seq_along(theta_fixed)[1:n_parameters_to_plot], function(i) {
        seq(theta_fixed[i] - l, theta_fixed[i] + l, length.out = n_points)
    })
    names(partitions) <- names(theta_fixed[1:n_parameters_to_plot])

    plot_data <- do.call(rbind, lapply(names(restricted_priors_and_posteriors), function(function_type) {
        do.call(rbind, lapply(names(restricted_priors_and_posteriors[[function_type]]), function(prior_type) {
            do.call(rbind, lapply(names(restricted_priors_and_posteriors[[function_type]][[prior_type]]), function(parameter_name) {
                # Calculate the function values
                x_values <- partitions[[parameter_name]]
                y_values <- sapply(x_values, restricted_priors_and_posteriors[[function_type]][[prior_type]][[parameter_name]])
                # Normalize y values so max is 0
                y_values <- y_values - max(y_values)

                data.frame(
                    x = x_values,
                    Value = y_values,
                    Parameter = parameter_name,
                    FunctionType = function_type,
                    PriorType = prior_type,
                    stringsAsFactors = FALSE
                )
            }))
        }))
    }))


    # Create a single ggplot object instead of a list of plots
    p <- ggplot(plot_data, aes(x = x, y = exp(Value), color = PriorType, linetype = FunctionType)) +
        geom_line() +
        geom_vline(
            data = data.frame(Parameter = names(theta_fixed[1:n_parameters_to_plot]), xintercept = unlist(theta_fixed[1:n_parameters_to_plot])),
            aes(xintercept = xintercept), color = "blue"
        ) +
        facet_wrap(~Parameter, ncol = 2) +
        labs(y = "Density")

    # If you want to save the plot
    if (!is.null(path)) {
        ggsave(path, p, device = "pdf", dpi = 100, width = 10, height = 5)
    }
    p
}

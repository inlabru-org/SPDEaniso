#' @title Plot distances to MAP
#' @description Plot distances to MAP for each prior type
#' @param results A list of results from the simulation containing the distances to MAP
#' @param prior_types A character vector of prior types
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @return A ggplot object
#' @export

plt_distances_to_MAP <- function(results, prior_types, path = NULL, dpi = 600) {
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

    ggplot(all_distances) +
        stat_ecdf(aes(value, color = prior_type)) +
        facet_wrap(~variable)
    if (!is.null(path)) {
        ggsave(path, dpi = dpi)
    }
}

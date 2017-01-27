fit_exp_decay <- function(metricsTable, metric = c("median", "mean")) {
    metric <- match.arg(metric)
    
    metricsTable <- dplyr::select_(metricsTable, "bin", "window", metric) %>%
        dplyr::rename_(metric = metric)
    
    tt <- dplyr::filter(metricsTable, window == "top_window")
    bt <- dplyr::filter(metricsTable, window != "top_window")
    # diff is the diference between actual value and each randomisation
    bt$diff <- vapply(
        seq_len(nrow(bt)),
        function(i) {
            unlist(tt[which(tt$bin == unlist(bt[i, "bin"])), "metric"]) - unlist(bt[i, "metric"])
        },
        0
    )
    
    # we fit a exponential decay, started values guessed from example datasets 
    return(nls(
        diff ~ x1 + x2 * exp(x3 * bin),
        data = bt,
        start = list(x1 = 0.1, x2 = 0.5, x3 = -0.5)
    ))
}


determine_bin_cutoff <- function(exp_decay_output, threshold = 0.1) {
    mm <- coef(exp_decay_output)
    mf <- function(x) mm[1] + mm[2] * exp(mm[3] * x)
    y_th <- mm[1] + (mf(1) - mm[1]) * threshold
    return(
        log((y_th - mm[1]) / mm[2]) / mm[3]
    )
}

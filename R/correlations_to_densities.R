
.applyDensity <- function(x, n = 64, absolute_cc = TRUE) {
    if (absolute_cc) {
        myDens <- density(abs(x), n = n, from = 0, to = 1, cut = 0)
    } else {
        myDens <- density(x, n = n, from = -1, to = 1, cut = 0)
    }
    return(list(
        cor_coef = myDens$x,
        density = myDens$y
    ))
}

correlations_to_densities <- function(df, absolute_cc = TRUE) {
    group_by(df, bin, window) %>%
        do(do.call(data_frame, .applyDensity(.$cor_coef, absolute_cc = absolute_cc))) %>%
        ungroup
}

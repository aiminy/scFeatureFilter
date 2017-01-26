get_mean_median <- function(df) {
    dplyr::group_by(df, bin, window) %>%
        dplyr::summarise(mean = mean(abs(cor_coef)), median = median(abs(cor_coef))) %>%
        dplyr::ungroup()
}

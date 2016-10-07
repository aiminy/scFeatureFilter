get_mean_median <- function(df) {
    dplyr::group_by(df, bin, window) %>%
        dplyr::summarise_if(is.numeric, funs(mean, median)) %>%
        dplyr::ungroup
}

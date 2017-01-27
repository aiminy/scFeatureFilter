

getExpressionTable <- function(binedTable, bin_cutoff = NULL, as.matrix = FALSE) {
    
    if (!is.null(bin_cutoff)) {
        binedTable <- dplyr::filter(binedTable, bin < bin_cutoff)
    }
    
    binedTable <- select(binedTable, -mean, -sd, -cv, -bin)
    
    if (as.matrix) {
        geneName <- binedTable$geneName
        binedTable <- select(binedTable, -geneName) %>%
            as.matrix
    }
    
    return(binedTable)
}


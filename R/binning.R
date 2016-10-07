#' Select top highly expressed genes.
#'
#' Selects the group of genes in the dataset that will be considered more highly expressed,
#' the top window, using the filtering method chosen by the user.
#'
#' There are three selection methods available:
#'
#' \itemize{
#'
#'      \item \code{"window_size"}: genes are ranked by mean expression, and a subset of the size
#'      indicated in \code{parameter} is selected from the top.
#'
#'      \item \code{"min_expression"}: genes where all expression values are above a minimum
#'      expression threshold indicated in \code{parameter} are selected.
#'
#'      \item \code{"mean_expression"}: the \code{mean} column is checked, and all genes with mean
#'      expression above the threshold indicated in \code{parameter} are selected.
#' }
#'
#' There are no restrictions to the \code{parameter} argument, however, the value should
#' be coherent with the characteristics of the data set provided and the chosen method.
#'
#' In general, it is adviseable to avoid generating top windows much larger than 250 genes,
#' to prevent excessively long computation time as well as to preserve the quality of the
#' analysis, as the top window should only include a subset of reliable values. As a rule,
#' the bigger the top window is, the more likely is that the reliability of the values is
#' compromised, given the characteristics of single cell RNA sequencing data.
#'
#' @param dataset A data frame, containing genes as rows and cells as columns, and where
#' the mean expression value for each gene has been added as a column.
#'
#' @param method A string indicating the method to use when creating the top window. If no
#' method is indicated, \code{"window_size"} will be used.
#'
#' @param parameter An integer. Indicates the numeric parameter to use in the previously
#' chosen method.
#'
#' @return A list with two elements, both data frames: the generated top window, and
#' the rest of the genes.

extract_top_genes <- function(dataset,
                              window_size = NULL,
                              mean_expression = NULL,
                              min_expression = NULL
                              ){
    divided_data <- list()
    expr_values <- data.frame()
    sorted_values <- data.frame()
 
    if (is.null(window_size) & is.null(mean_expression) & is.null(min_expression)) stop("Need to provide one of the following parameter: window_size, mean_expression or min_expression")

    if (is.numeric(window_size)){

        # sort data by mean expression - original row names are lost
        sorted_values <- arrange(dataset, desc(mean))
        # select the top x genes (x=window size selected)
        divided_data$topgenes <- sorted_values[1:window_size, ]
        # assign bin 0 to the top window
        divided_data[[1]]$bin <- 1
        # display mean FPKM value of the last gene in the top window
        message(paste("Mean expression of last top gene:",
                      sorted_values[window_size, ]$mean))
        # store the rest of the genes a the second element of the list
        divided_data$restofgenes <- sorted_values[-(1:nrow(divided_data$topgenes)), ]

    } else if (is.numeric(mean_expression)){
        
        # select top genes and store in first position of the list
        divided_data$topgenes <- subset(dataset, dataset$mean > mean_expression)
        # assign bin 0 to the top window
        divided_data[[1]]$bin <- 1
        # display number of genes in the window
        message(paste("Number of genes in top window:",
                      nrow(divided_data$topgenes)))
        # store rest of genes in second position
        divided_data$restofgenes <- subset(dataset, dataset$mean < mean_expression)
        
    } else if (is.numeric(min_expression)){

        # extract expression values
        expr_values <- dplyr::select(dataset, -geneName, -mean, -sd, -cv)
        # select top genes and store in first position of the list
        divided_data$topgenes <- subset(dataset,
                                        rowSums(expr_values > min_expression) == ncol(expr_values))
        # assign bin 0 to the top window
        divided_data[[1]]$bin <- 1
        # display number of genes in the window
        message(paste("Number of genes in top window:",
                      nrow(divided_data$topgenes)))
        # store rest of genes in second position
        divided_data$restofgenes <- subset(dataset,
                                           rowSums(expr_values > min_expression) != ncol(expr_values))


        
    } else {
        stop("The second parameter should be numeric.")
    }

    return(divided_data)
}

#' Bin genes by mean expression.
#'
#' Divides the genes that were not included in the top window in windows of the same size,
#' by their mean expression.
#'
#' There are two binning methods available:
#'
#' \itemize{
#'      \item \code{"window_number"}: Divides the genes into the number of windows specified in
#'      \code{parameter}, regardless of their size.
#'
#'      \item \code{"window_size"}: Divides the genes into windows of the size specified in
#'      \code{parameter}, regardless of the number of windows generated.
#' }
#'
#' This function uses the \code{ntile} function, in the \code{dplyr} package to assign a bin
#' number to each gene based on the value contained in the \code{mean} column, corresponding
#' to its mean expression. These are then added as a the column \code{bin} using the \code{mutate}
#' function, also in the \code{dplyr} package.
#'
#' \strong{Important note:} This function is designed to take the list output by the
#' \code{extract_top_window} function as an argument, operating only on the second element
#' of it.
#'
#' Once the genes in it have been binned, both elements of the list are bound
#' together in a data frame and returned. The output is similar, but a new column \code{bin}
#' is added, which indicates the window number assigned to each gene.
#'
#' @param dataset A list, containing the top window generated by \code{extract_top_genes}
#' as the first element, and the rest of undivided genes as the second.
#'
#' @param method A string, indicating the method to be used to bin the genes by mean
#' expression.
#'
#' @param parameter An integer. Indicates the numeric parameter to use in the previously
#' chosen method. Values are not restricted, but should be coherent with the method of choice.
#'
#' @return A data frame containing the binned genes.

bin_scdata <- function(dataset, window_number = NULL, window_size = NULL){

    if (is.null(window_number) & is.null(window_size)) stop("Need to provide window_number or window_size.")

    if (is.numeric(window_number)){
        # bin into the selected number of windows
        windows <- ntile(desc(dataset[[2]]$mean), window_number) + 1
    } else if (is.numeric(window_size)){
        # calculate number of windows of selected window size possible and bin
        windows <- ntile(desc(dataset[[2]]$mean), trunc(nrow(dataset$restofgenes)/window_size)) + 1
    } else {
        stop("The second parameter should be numeric.")
    }
    dataset$restofgenes <- mutate(dataset$restofgenes, bin = windows)

    if (is.numeric(window_number)){
        # print size of the desired number of windows created
        message(paste("Window size:", length(which(dataset[[2]]$bin == 2))))
    } else if (is.numeric(window_size)){
        # print number of windows of desired size created
        message(paste("Number of windows:", max(dataset[[2]]$bin)))
    }

    # bind all the data and correct bin number
    dataset <- dplyr::bind_rows(dataset) %>% 
        dplyr::select(geneName, mean, sd, cv, bin, everything())
    return(dataset)
}

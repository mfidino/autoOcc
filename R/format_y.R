#' @title Format detection data into a 3D array
#'
#' @description Take a data.frame with detection data and convert it to a
#' three dimensional array to be used in \code{\link{auto_occ}}.
#'
#'
#' @param x A data.frame in mostly long format that should have one columns that denote
#' sites sampled, one columns that denotes the primary sampling period, and multiple columns
#' for the detection history (one column per secondary sampling period within a primary sampling
#' period). See details for more information.
#'
#' @param site_column Either a character that is the column name of the site column or a
#' numeric or integer the represents the site column's location.
#'
#' @param time_column Either a character that is the column name of the primary sampling period column or a
#' numeric or integer the represents the primary sampling period column's location.
#'
#' @param history_columns Either a character object of length 1 that contains a regular expression to
#' identify the columns that house a species detection history or a numeric or an integer vector that
#' represents the column locations that house a species detection history.
#'
#' @param report Some QA/QC information that the function shares about how the data
#' get ordered within the detection array and also the names of the history_columns
#' that were identified by either the regex or numeric/integer vector. Defaults to \code{TRUE}.
#'
#' @details
#'
#' This function is just a quick way to ensure that more standard representations
#' of a species detection history get converted into what would be accepted by
#' \code{\link{auto_occ}}. The classes of the associated columns in \code{x} may vary somewhat
#'  \code{format_y}. Mor specifically:
#'
#'  \itemize{
#'    \item{The site column}{Should preferably be either of class \code{character} or
#'      \code{factor}. If this column is either a \code{numeric} or \code{integer},
#'      \code{format_y} will provide a warning that it is converting the column
#'      to a character class with padded zeroes on the left of the numbers to ensure
#'      correct sorting (e.g., \code{"001"} instead of \code{"1"}).
#'    }
#'    \item{The primary sampling period column}{Can either be of class \code{factor},
#'    \code{character}, \code{numeric}, or \code{integer}. As a \code{factor}, the rows
#'    of \code{x} will be sorted based on the levels of the primary sampling period column.
#'    As a \code{character}, then their order of appearance from the top of \code{x} downwards
#'    will be used to sort them. As such, when this column is a \code{character}, it
#'    essentially assumes that the rows are already sorted in the way they should be
#'    for analysis, with the first primary sampling period on top of \code{x}, followed by the
#'    second primary sampling period, etc. If this column is either a \code{numeric} or
#'    \code{integer}, then their numeric order will be used to sort the rows of \code{x}.
#'    }
#'    \item{The detection history columns}{should either be of class \code{numeric} or
#'    \code{integer}. Each element of these columns must either be a 1 if a species was
#'    detected on a given secondary sampling period within a primary sampling period, a 0 if
#'    the species was not detected, or an NA if no sampling occurred. This function will
#'    return an error if any element within the detection history columns are not a \code{0},
#'    \code{1}, or \code{NA}.
#'    }
#'  }
#'
#' @examples
#  load in data
#' data("opossum_det_hist")
#'
#'
#' # a quick look at the top six rows of this data.frame
#' #  so that you can see the ordering of the data.
#' #       Site Season Week_1 Week_2 Week_3 Week_4
#' # 1 D02-BMT1   JA19     NA      0     NA      0
#' # 2 D02-HUP0   JA19     NA      0      0      0
#' # 3 D02-JBS1   JA19     NA      0      1      0
#' # 4 D02-MDU1   JA19     NA     NA     NA     NA
#' # 5 D02-MOP1   JA19     NA      0      0      0
#' # 6 D02-RCP1   JA19     NA      0      1      1
#' #
#'
#' opossum_y <- format_y(
#'   x = opossum_det_hist,
#'   site_column = "Site",
#'   time_column = "Season",
#'   history_columns = "^Week", # the detection history columns all start with the word Week
#'   report = FALSE #only setting to FALSE for example, defaults to TRUE
#' )
#'
#'
#'
#' @export

format_y <- function(x, site_column, time_column, history_columns, report = TRUE){

  if(!class(site_column) %in% c("character","numeric","integer")){
    stop("site_column input must be either a character, numeric, or integer")
  }
  if(is.integer(site_column)|is.numeric(site_column)){
    site_column <- colnames(x)[site_column]
  }
  if(!class(time_column) %in% c("character","numeric","integer")){
    stop("time_column input must be either a character, numeric, or integer")
  }
  if(is.integer(time_column)|is.numeric(time_column)){
    time_column <- colnames(x)[time_column]
  }
  if(
    class(history_columns) %in% c("integer","numeric") &
    length(history_columns) == 1
  ){
    stop("An integer or numeric of length 1 was supplied to history_columns. When specifying an integer or numeric history_columns must be a vector that identifies which columns in x represent a species detection history.")
  }
  if(is.character(history_columns) & length(history_columns)>1){
    stop("A character object with length > 1 supplied to history_columns. When this argument is a character object regex is used to identify the associated columns, and so the character object must be a scalar.")
  }
  if(!class(history_columns) %in% c("character","numeric","integer")){
    stop("history_columns input must be either a character, numeric, or integer")
  }
  # check for any duplicates in history column input,
  #  error out if so.
  if(any(duplicated(history_columns))){
    stop("history_columns cannot have duplicates")
  }
  if(!class(x[,site_column]) %in% c("factor", "character")){
    my_warning <- paste0(
      "site_column was of class ",
      paste0(class(x[,site_column]),collpase = ", "),"converted to character with padded zeroes."
    )
    warning(my_warning)
    tmp <- x[,site_column]
    max_pad <- max(nchar(as.character(tmp)))
    x[,site_column] <- sprintf(
      paste0("%0",max_pad,"d"),
      as.numeric(tmp)
    )
  }

  # get indices, number of sites
  nsite <- length(
    unique(
      x[,site_column]
      )
  )
  # number of primary sampling periods
  nseason <- length(
    unique(
      x[,time_column]
      )
  )
  if(report){
    cat("\n\nTEMPORAL ORDERING\n-----------------\n\n")
  }
  if(is.factor(x[,time_column])){
    if(report){
      cat("Primary sampling period column is a factor, factor levels to order temporally.\n")
      cat(paste0("Ordering: ", paste0(levels(x[,time_column]), collapse = ", ") ))
    }
    x <- x[order(x[,time_column]),]
  }
  if(is.character(x[,time_column])){
    if(report){
      cat("Primary sampling period column is a character vector, using their order of appearance from top of x to order temporally.\n")
      cat(paste0("Ordering: ", paste0(unique(x[,time_column]), collapse = ", ")))
      tmp <- x[,time_column]
      tmp <- factor(tmp, levels = unique(tmp))
      x <- x[order(tmp),]
      rm(tmp)
    }
  }
  if(is.numeric(x[,time_column])|is.integer(x[,time_column])){
    if(report){
      cat("Primary sampling period column is a numeric or integer vector, sorting numerically to order temporally.\n")
      tmp <- unique(x[,time_column])
      cat(paste0("Ordering: ", paste0(tmp[order(tmp)], collapse = ", ")))
      rm(tmp)
    }
    x <- x[order(x[,time_column]),]
  }
  # if character, do your regex
  if(is.character(history_columns)){
    history_columns <- grep(history_columns, colnames(x))
    # check class of each column
    history_classes <- sapply(
      x[,history_columns],
      class
    )
    if(any(!history_classes %in% c("numeric", "integer"))){
      stop(
        paste0(
          "One of the columns identified by inputted value to 'history_columns'",
          " was not a numeric or integer.",
          "\n\nColumns queried: ", paste0(
            names(history_classes),
            collapse = ", "
          )
        )
      )
    }
  }
  history_columns <- colnames(x)[history_columns]
  if(length(history_columns) == 1){
    stop("Value input to 'history_columns' only returned a single column. Detection histories must have > 1 secondary sample.")
  }
  if(report){
    cat("\n\nDETECTION HISTORIES\n-------------------\n\n")
    cat(
      paste0(
        length(history_columns),
        " detection history columns found.\nColumn names: ",
        paste0(
          history_columns,
          collapse = ", "
        ),"\n"
      )
    )
  }

  # number of secondary samples
  nrep <- length(
    unique(
      history_columns
    )
  )

  # create the y array
  y <- array(
    NA,
    dim = c(nsite, nseason, nrep),
    dimnames = list(
      unique(x[,site_column]),
      unique(x[,time_column]),
      history_columns
    )
  )
  # split data.frame up by primary sampling period
  split_x <- split(
    x,
    factor(x[,time_column], levels = unique(x[,time_column]))
  )
  for(i in 1:length(split_x)){
    tmp_df <- data.frame(
      site = unlist(dimnames(y)[1]),
      season = names(split_x)[i]
    )
    colnames(tmp_df) <- c(site_column, time_column)
    tmp_df <- merge(
      tmp_df,
      split_x[[i]],
      by = c(site_column, time_column)
    )
    y[,i,] <- as.matrix(tmp_df[,history_columns])
  }
  # check to make sure there are only 0, 1, or NA present
  y_vals <- suppressWarnings(unique(as.numeric(y)))
  if(any(!y_vals %in% c(NA, 0, 1))){
    stop("Elements of the y array must be either 0, 1, or NA. Check your detection history columns for elements that do not take those values.")
  }
  return(y)
}

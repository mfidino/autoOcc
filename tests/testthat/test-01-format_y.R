test_that(
  "format_y",{
    f <- function(
      x,
      site_column,
      time_column,
      history_columns,
      report = FALSE
    ){
      to_return <- format_y(
        x = x,
        site_column = site_column,
        time_column = time_column,
        history_columns = history_columns,
        report = report
      )
      return(to_return)
    }
    data("opossum_det_hist")
    expect_silent(
      f(
        x = opossum_det_hist,
        site_column = "Site",
        time_column = "Season",
        history_columns = "^Week"
      )
    )
    # Error because it is 'Site' not 'site'
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = "site",
        time_column = "Season",
        history_columns = "^Week"
      )
    )
    # regex too specific, only returns one column.
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = "Site",
        time_column = "Season",
        history_columns = "^Week_1"
      )
    )
    # integers input
    expect_silent(
      f(
        x = opossum_det_hist,
        site_column = 1,
        time_column = 2,
        history_columns = 3:6
      )
    )
    # same thing, but only one column for history_columns
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = 1,
        time_column = 2,
        history_columns = 3
      )
    )
    expect_output(
     ay <-  f(
        x = opossum_det_hist,
        site_column = 1,
        time_column = 2,
        history_columns = 3:6,
        report = TRUE
      )
    )
    opo <- opossum_det_hist
    opo$Season <- factor(
      opo$Season,
      levels = unique(opo$Season)
    )
    expect_output(
      ay <-  f(
        x = opo,
        site_column = 1,
        time_column = 2,
        history_columns = 3:6,
        report = TRUE
      )
    )
    opo$Season <- as.numeric(opo$Season)
    expect_output(
      ay <-  f(
        x = opo,
        site_column = 1,
        time_column = 2,
        history_columns = 3:6,
        report = TRUE
      )
    )
    # way more locations then needed
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = 1,
        time_column = 2,
        history_columns = 3:10
      )
    )
    # sloppy regex
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = 1,
        time_column = 2,
        history_columns = "Site|Week"
      )
    )
    # put an incorrect number in detection history
    my_data <- opossum_det_hist
    my_data$Week_3[4] <- 5
    expect_error(
      f(
        x = my_data,
        site_column = 1,
        time_column = 2,
        history_columns = "Week"
      )
    )
    # site column identifier cant be a factor
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = factor("yo"),
        time_column = 2,
        history_columns = "Site|Week"
      )
    )
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = 1,
        time_column = factor("yo"),
        history_columns = "Site|Week"
      )
    )
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = 1,
        time_column = 2,
        history_columns = factor("yo")
      )
    )
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = 1,
        time_column = 2,
        history_columns = c("Week_1", "Week_2")
      )
    )
    expect_error(
      f(
        x = opossum_det_hist,
        site_column = 1,
        time_column = 2,
        history_columns = c(3,3,3,3)
      )
    )
    my_data <- opossum_det_hist
    my_data$site <- as.numeric(factor(my_data$Site))
    expect_warning(
      f(
        x = my_data,
        site_column = "site",
        time_column = 2,
        history_columns = "^Week"
      )
    )
    my_data$site <- factor(my_data$Site)
    expect_silent(
      f(
        x = my_data,
        site_column = "site",
        time_column = 2,
        history_columns = "^Week"
      )
    )
  }
)

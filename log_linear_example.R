# simulating a simple loglinear model

# example. We are interested in the average number of minutes
#  between raccoon images at sites that vary along a gradient
#  of forest cover. Our hypothesis is that in areas with lots
#  of forest cover, raccoon abundance is greater, and so
#  the amount of time between raccoon images should be lower.

set.seed(111)

n <- 300

# a covariate
forest_cover <- runif(n, 0, 1)

# scale the covariate for the regression
forest_cover_scaled <- scale(forest_cover)

# known intercept and slope term
#  these are on the log scale
b0 <- 4
b1 <- -0.75

# the mean response for each data point
response_mean <- b0 + b1 * as.numeric(forest_cover_scaled)

# add variation
response_logscale <- response_mean + rnorm(n, mean = 0, sd = 0.5)

# convert to real scale
response <- exp(response_logscale)


model_df <- data.frame(
  raccoon_time = response,
  forest = forest_cover_scaled
)

# model with the covariate
m1 <- lm(
  log(raccoon_time+1) ~ forest,
  data = model_df
)

# null model for comparison
m2 <- lm(
  log(raccoon_time+1) ~ 1,
  data = model_df
)

aic_tab <- AIC(m1, m2)
aic_tab$delta <- aic_tab$AIC - min(aic_tab$AIC)

# aic weight, or the relative likelihood of each model
#  in the model set.
aic_tab$aicwt <- exp(-0.5 *aic_tab$delta) /
  sum(exp(-0.5 *aic_tab$delta) )

aic_tab <- aic_tab[order(aic_tab$delta),]

# cumulative weight of each model
aic_tab$cumltvWt <- cumsum(aic_tab$aicwt)

aic_tab

# looks like model 1 fit better!. Let's see how well
#  it did in returning the model parameters we
#  used to simulate these data.
summary(m1)
b0
b1
# using forest cover, not forest cover scaled,
#  because we dont care about the scaled covariate aside
#  from fitting the model
plot(
  model_df$raccoon_time ~ forest_cover,
  ylab = "Minutes between images",
  xlab = "Forest cover (proportion)",
  las = 1,
  bty = "l"
)

# make a model prediction. Start with a sequence of values
#  that are relevant to the 'real world.' In other words
#  the non-scaled covariate values. We know that the
#  proportional data ranges from 0 to 1, but when
#  analyzing real data it may not run that full range.
#  so it is always a good idea to check and see
#  what realm you can reasonably make predictions
#  in vs. extrapolation.
(forest_range <- range(forest_cover))


# here, it looks like it is just about 0 to 1, so we
#  will stick with that. The covariate was called
#  'forest' in our model, so we need to stick with
#  that.
new_data <- data.frame(
  forest = seq(0,1, length.out = 500)
)

# scale it by the data
new_data_scaled <- new_data
new_data_scaled$forest <-
  (new_data_scaled$forest - mean(forest_cover)) / sd(forest_cover)

my_preds <- predict(
  m1,
  newdata = new_data_scaled,
  interval = "confidence"
)

# exponentiate, which is the reverse of the log link.
my_preds <- data.frame(exp(my_preds))


plot(
  model_df$raccoon_time ~ forest_cover,
  ylab = "Minutes between images",
  xlab = "Forest cover (proportion)",
  las = 1,
  bty = "l"
)
lines(
  my_preds$fit ~ new_data$forest,
  lwd = 4,
  col = "darkgreen"
)
lines(
  my_preds$lwr ~ new_data$forest,
  lwd = 2,
  col = "darkgreen",
  lty = 2
)
lines(
  my_preds$upr ~ new_data$forest,
  lwd = 2,
  col = "darkgreen",
  lty = 2
)

# some little 'facts' we could get from our model.

# average amount of time between raccoon images at an 'average'
#  location. 54.52 minutes (95% CI = 51.63, 57.56). This only
#  works because we scaled the covariate (subtracted the mean).
exp(
  qnorm(
    c(0.025,0.5,0.975),
    3.99852,
    0.02773
  )
)
# little example of a join

install.packages("dplyr")
library(dplyr)

my_join <- left_join(
  df_1,
  df_2,
  by = "siteID"
)

my_join <- data.frame(my_join)


my_unique_sites <- data.frame(
  sites = unique(racdata$sites)
)
write.csv(my_unique_sites, "mason_unique.csv")


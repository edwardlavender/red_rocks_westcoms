################################
################################
#### analyse_eggs.R

#### This script: 
# 1) Analyses egg-area scaling relationships 

#### Steps preceding this script: 
# 1) Define global parameters     (define_global_param.R)


################################
################################
#### Set up

#### Wipe workspace and source essential packages and variables
source("./R/define_global_param.R")

#### Load data (egg counts)
counts      <- readxl::read_excel("./data-raw/eggs/egg_counts.xlsx", sheet = "rov")
counts$size <- counts$area

#### Essential packages 
library(magrittr)
library(prettyGraphics)

### Define global parameters 
# Define seed
set.seed(1)
# Define threshold number of eggs of interest 
th_eggs <- 80
# Define function to add threshold eggs and corresponding areas to plots 
add_ths <- function(.th_area, .th_eggs = th_eggs){
  lines(c(0, .th_area), c(.th_eggs, .th_eggs), lty = 3)
  arrows(.th_area, .th_eggs, .th_area, 0, length = 0.1)
}

################################
################################
#### Modelling 

#### Implement model 
mod <- lm(count ~ 0 + size, data = counts)

#### Visualise model 
pretty_predictions_1d(mod)
calc_size <- function(y, c, m) (y - c)/m
calc_size(th_eggs, 0, coef(mod)[1])

#### Implement a simulation approach 
# ... to estimate the area within which the threshold number of eggs would be found
n_sim <- 1000
sim_size_glm_ls <- 
  pbapply::pblapply(1:n_sim, function(i){
    ind <- sample(1:nrow(counts), nrow(counts), replace = TRUE)
    dat <- counts[ind, ]
    mod <- lm(count ~ 0 + size, data = dat)
    th <- calc_size(th_eggs, 0, coef(mod)[1])
    out <- list(ind = ind, dat = dat, mod = mod, th = th)
    return(out)
  })

#### Examine estimates for the scalinng parameter
sapply(sim_size_glm_ls, function(elm) coef(elm$mod)) %>% 
  utils.add::basic_stats(p = NULL, 
                         f = list(median = median, 
                                  q = function(x) quantile(x, p = 0.975)))

#### Visualise results
## Set up figure to save 
png("./fig/egg_scaling.png", 
    height = 5, width = 5, units = "in", res = 600)
pp <- par(oma = c(1, 1, 1, 1))
## Define suitable ylimits 
ylim <- 
  pbapply::pblapply(sim_size_glm_ls, function(elm) {
  dat <- elm$dat
  mod <- elm$mod
  x <- seq(min(dat$size), max(dat$size), length.out = 100)
  p <- predict(elm$mod, newdata = data.frame(size = x), type = "response")
  ylim <- range(c(p, dat$count))
  return(ylim)
}) %>% unlist() %>% range()
## Define blank plot 
pretty_plot(counts$area, counts$count, 
            pretty_axis_args = list(x = list(x = range(counts$size), y = range(ylim))), 
            xlab = "", ylab = "", 
            type = "n")
## Add predictions 
pbapply::pblapply(sim_size_glm_ls, function(elm) {
    dat <- elm$dat
    mod <- elm$mod
    x <- seq(min(dat$size), max(dat$size), length.out = 100)
    p <- predict(elm$mod, newdata = data.frame(size = x), type = "response")
    lines(x, p, lwd = 0.1, col = scales::alpha("grey", 0.5))
    return(invisible())
  }) %>% invisible()
## Add observations 
points(counts$area, counts$count,
       pch = 21, 
       col = scales::alpha("dimgrey", 0.9), 
       bg = scales::alpha("dimgrey", 0.9))
## Add estimates area thresholds
sim_size_glm <- sapply(sim_size_glm_ls, function(elm) elm$th)
utils.add::basic_stats(sim_size_glm)
th_areas <- quantile(sim_size_glm, c(0.025, 0.975))
th_areas
add_ths(th_areas[1])
add_ths(th_areas[2])
mtext(side = 1, expression("Area (" * m^2 * ")"), line = 2.5)
mtext(side = 2, "Count", line = 2.5)
par(pp)
dev.off()


#### End of code.
################################
################################
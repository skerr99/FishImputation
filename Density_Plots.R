library(mvtnorm)
library(tidyverse)
library(mice)
library(plotly)


# create dataset
mu <- c(0, 0)
sigma <- matrix(c(1, .5, .5, 1), nrow = 2)
dat <- as.data.frame(mvtnorm::rmvnorm(n=500, mean=mu, sigma = sigma))

# duplicate dataset and add missingness
dat_missing <- dat
missing_type <- 'mar'

if(missing_type == 'mar'){
  med <- median(dat$V2)
  n1 <- sum(dat$V2 < med)
  dat_missing$V1[dat_missing$V2 < med & runif(n1) < .8] <- NA
  n2 <- sum(dat_missing$V2 >= med)
  dat_missing$V1[dat_missing$V2 >= med & runif(n2) < .2] <- NA
} 
if(missing_type == 'mcar'){
  missing_amt <-  nrow(dat) * .5
  dat_missing$V1[1:missing_amt] <- NA
}


# number of imputations (also the amount of iterations in the loop below)
m <- 10

# impute data
imputed <- mice(dat_missing, print=FALSE, m=m)

# note: all densities calculated have parameters from=-4, to=4, n=512 so they can be compared easier later 
#   - in geom_density, it's bounds=c(-4,4), n=512

imp_stacked <- data.frame() # will have raw imputation data binded to it (V1, V2)
full_imp_den <- data.frame() # will have imputation density curve data binded to it (x,y)
bws <- c() # will track nrd0 bandwidths of every imputed V1
plot <- ggplot() # will have individual imputation curves added to it
sds <- c()
iqrs <- c()
# loop from 1 to 10 (amount of imputations stored in `imputed`)
for(i in 1:m){
  
  # pull ith imputed data and bind to `imp_stacked`
  imp_data <- complete(imputed, i)
  imp_stacked <- bind_rows(imp_stacked, imp_data)
  
  # pull density of ith imputed data and bind to `full_imp_den`
  imp_den <- data.frame(x=density(imp_data$V1, from=-4, to=4, n=512)$x, y=density(imp_data$V1, from=-4, to=4, n=512)$y)
  full_imp_den <- bind_rows(full_imp_den, imp_den)
  
  # pull bandwidth of ith imputed data and append to `bws`
  bw <- bw.nrd0(imp_data$V1)
  bws <- c(bws, bw)
  
  # pull stds and iqrs
  sd <- sd(imp_data$V1)
  sds <- c(sds, sd)
  iqr <- IQR(imp_data$V1)
  iqrs <- c(iqrs, iqr)
  
  # create plot of ith imputed data and add to plot object
  imp_plot <- geom_density(data=imp_data, aes(x=V1, color='Individual Imputation'), size=.5, alpha=.5, bounds = c(-4,4), n=512)
  plot <- plot + imp_plot
  
}




# create "Rounded + Grouped" estimation of imputations combined
imp_den_grouped <- full_imp_den %>% 
  group_by(x) %>% summarize(y=mean(y))


# apply average bandwidth to stacked raw imputed data to create "Averaged BW" estimation of imputations combined
imp_avg_bw_den <- density(imp_stacked$V1, bw = mean(bws), from=-4, to=4, n=512)
imp_avg_bw <- data.frame(x= imp_avg_bw_den$x, y=imp_avg_bw_den$y)


# calculate the nrd0 bw for the stacked imputed data using n=500 and the average SDs and IQRs from individual imputations
stacked_bw <- 0.9 * min(mean(sds), mean(iqrs) / 1.34) * 500^-0.2


# create interactive `tot_plot`, which includes: 
#   "truth" -- blue
#   random sample -- red
#   random sample with missingness -- gold
#   all ten individual imputation curves (stored in `plot` object from loop) -- lightgray
#   all three estimations of the imputations combined -- dark gray, black, and purple

tot_plot <-
  plot + 
  geom_density(data=dat, aes(x=V1, color='Original')) + 
  theme_classic() +
  geom_density(data=dat_missing, aes(x=V1, color='Missing'), bounds=c(-4,4)) + 
  geom_density(data=imp_stacked, aes(x=V1, color='Imputed Combined: Stacked'), bw=stacked_bw, bounds=c(-4,4)) +
  geom_line(data=imp_avg_bw, aes(x=x, y=y, color='Imputed Combined: Averaged BW')) +
  geom_line(data=imp_den_grouped, aes(x=x, y=y, color='Imputed Combined: Grouped')) +
  xlim(-4, 4) +
  stat_function(aes(color='Truth'), fun = dnorm, args = list(mean = 0, sd = 1)) +
  scale_color_manual(
    values = c(
    'Truth' = 'blue',
    'Original' = 'red', 
    'Missing' = 'gold', 
    'Individual Imputation' = 'lightgray',
    'Imputed Combined: Stacked' = 'darkgray',
    'Imputed Combined: Rounded+Grouped' = 'black',
    'Imputed Combined: Averaged BW' = 'purple'))

ggplotly(tot_plot)



# pull y vectors of all density curves
orig_y <- density(dat$V1, from=-4, to=4, n=512)$y # red
missing_y <- density(na.omit(dat_missing$V1), from=-4, to=4, n=512)$y # gold
imp_stacked_y <- density(imp_stacked$V1, from=-4, to=4, n=512, bw=stacked_bw)$y # gray
imp_avg_bw_y <- imp_avg_bw$y # purple
imp_grouped_y <- imp_den_grouped$y 

orig_x <- density(dat$V1, from=-4, to=4, n=512)$x
truth_y <- dnorm(orig_x, 0, 1) # BLUE


# find a bunch of mean squared errors between these and truth (blue curve)
orig_mse <- sum((orig_y - truth_y)^2)
missing_mse <- sum((missing_y - truth_y)^2)
stacked_mse <- sum((imp_stacked_y - truth_y)^2)
avg_bw_mse <- sum((imp_avg_bw_y - truth_y)^2)
grouped_mse <- sum((imp_grouped_y - truth_y)^2)



# order resulting MSEs and display
results <- sort(c('Original'=orig_mse, 'Missing'=missing_mse, 'Stacked'=stacked_mse, 'Avg BW'=avg_bw_mse, 'Grouped'=grouped_mse))
results


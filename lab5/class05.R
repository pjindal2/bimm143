# class 5 Data Visualization 

#rnorm gives us x anmount of data points in a normal distribution 
x <- rnorm(1000)

# some summary stats
mean(x)
sd(x)

# gives several stats (box-whisker plot)
summary(x)

boxplot(x)
hist(x)
rug(x)

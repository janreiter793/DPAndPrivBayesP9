if(!("pwr" %in% rownames(installed.packages()))) {
  install.packages("pwr")
}
library(pwr)
library(magrittr)
library(tidyverse)

power_result <- pwr.r.test(n = 50, r = 0.67, sig.level = 0.05, alternative = "two.sided")
power_result$power

spearman_rank_coef <- 0
for(i in 1:10) {
   cat(i, "/", 10, "\n")
   x1 <- runif(10000)
   eps <- rnorm(10000, sd = 0.1, mean = 0)
   x2 <- x1 + eps
   X <- data.frame(x1, x2)
   while((X %>% 
          filter(x1 <= 1, x1 >= 0,
                 x2 <= 1, x2 >= 0) %>% 
          nrow) < 10000) {
     temp <-
       X %>% 
       filter(x1 > 1 | x1 < 0 | x2 > 1 | x2 < 0) %>% 
       nrow
     x1 <- runif(temp)
     eps <- rnorm(temp, sd = 0.1, mean = 0)
     x2 <- x1 + eps
     
     temp <- data.frame(x1, x2)
     X %<>% 
       filter(x1 <= 1, x1 >= 0,
              x2 <= 1, x2 >= 0) %>% 
       rbind(temp)
   }
   d <- rank(X$X1) - rank(X$X2)
   spearman_rank_coef <- spearman_rank_coef + 6 * sum(d^2) / (10000 * (10000^2 - 1))
}
spearman_rank_coef <- spearman_rank_coef / 10

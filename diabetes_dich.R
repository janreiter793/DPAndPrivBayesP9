library(tidyverse)
library(magrittr)

data <- read.csv("C:/Users/janre/Documents/uni/9. Semester/Projekt/Eksamen/Presentation/diabetes.csv")
data %>% glimpse

data %<>% 
  mutate(ageUnder45 = if_else(Age < 45, 1, 0),
         overweight = if_else(BMI >= 25, 1, 0),
         everPregnant = (Pregnancies > 0))

data %>% 
  select(Outcome, ageUnder45, overweight, everPregnant) %>% 
  write.csv("C:/Users/janre/Documents/uni/9. Semester/Projekt/Eksamen/Presentation/diabetes_dichotomised.csv",
                   row.names = FALSE)

read.csv("C:/Users/janre/Documents/uni/9. Semester/Projekt/Eksamen/Presentation/diabetes_dichotomised.csv")

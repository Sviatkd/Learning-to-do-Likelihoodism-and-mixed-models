library(tidyverse)
library(arm)
library(rstanarm)
library(readxl)
library(lmtest)
Maths_and_Stats_grades <- read_excel("C:/Users/alexa/Downloads/Maths and Stats grades.xlsx")
car::scatterplot(Adv.Stats ~ Adv.Maths | Gender, data = Maths_and_Stats_grades)
#Full Bayes we are sure that the relationship between maths and stats is near 1 to 1 positive but allow it to vary#
BayesReg <- stan_glm(Adv.Stats ~ Adv.Maths + factor(Gender),family=gaussian(),prior = normal(location = c(0.8,0),autoscale = TRUE),data=Maths_and_Stats_grades)
#Get the Credible interval#
posterior_interval(BayesReg,prob = 0.95)
#Describe Prior and Posterior#
describe_prior(BayesReg)
describe_posterior(BayesReg,rope_range = c(-0.1, 0.1))
#Predict Scores for the Sample#
pointpredictions <- posterior_epred(BayesReg,newdata=Maths_and_Stats_grades)
#Get MCMC SDs from the model#
sigma_draws <- as.data.frame(BayesReg)$sigma
risk_table <- data.frame(
  Adv.Maths = Maths_and_Stats_grades$Adv.Maths
)
#Set up the CDF for P(Score<40) and apply to predictions and sampled sds#
risk_grid_draws <- pnorm(40, mean = pointpredictions, sd = sigma_draws)
risk_table$Risk_Mean <- colMeans(risk_grid_draws)
risk_table$Risk_Lower <- apply(risk_grid_draws, 2, quantile, probs = 0.025)
risk_table$Risk_Upper <- apply(risk_grid_draws, 2, quantile, probs = 0.975)
ggplot(data =risk_table,mapping = aes(x=Adv.Maths,y=Risk_Mean,colour = Maths_and_Stats_grades$Gender,fill=Maths_and_Stats_grades$Gender))+geom_jitter()+geom_line(linewidth = 1)+ geom_ribbon(aes(ymin = Risk_Lower, ymax = Risk_Upper), alpha = 0.3, colour = NA) +
  scale_y_continuous(labels = scales::percent_format(),breaks = scales::breaks_pretty(n = 10)) +
  labs(
    x = "Advanced Maths Score",
    y = "Probability of Failing (< 40)",
    title = "Risk of Failing Advanced Stats",
    subtitle = "Shaded regions represent 95% posterior credible intervals"
  ) +
  theme_minimal()
library(brms)
my_priors <- set_prior("normal(0.8, 1)", class = "b", coef = "Adv.Maths")

# Fit the Bayesian regression model using brms
BayesReg_brms <- brm(
  Adv.Stats ~ Adv.Maths + Gender,
  data = Maths_and_Stats_grades,
  family = gaussian(),
  prior = my_priors,
  save_pars = save_pars(all = TRUE)
)

# Get the Credible Interval
posterior_interval(BayesReg_brms, prob = 0.95)

# Describe Prior and Posterior
describe_prior(BayesReg_brms)
describe_posterior(BayesReg_brms, rope_range = c(-0.1, 0.1))
plot(BayesReg_brms)
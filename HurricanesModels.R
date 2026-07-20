library(sandwich)
library(tidyverse)
library(lme4)
library(stargazer)
library(modelsummary)
library(mgcv)
library(gt)
panel <- pivot_longer(HurricanePanel,cols = H0:H5)
ggplot(panel,aes(x=Period,y=value,colour = factor(name)))+geom_line() +geom_smooth(method = "loess", aes(group = 1), colour = "darkblue", size = 1.2, se = FALSE)
ggplot(panel,aes(x=Period,y=c(0,diff(value)),colour = factor(name)))+geom_line() + geom_smooth(method = "loess", aes(group = 1), colour = "darkblue", size = 1.2, se = FALSE)
pooledols <- lm(value~poly(Period,degree = 2,raw = TRUE),data = panel)
fdmodel <- lm(c(0,diff(value))~poly(Period,degree = 2,raw = TRUE),data = panel)
fdmodelnull <- lm(c(0,diff(value))~poly(Period,degree = 1,raw = TRUE),data = panel)
femodel <- lm(value~poly(Period,degree = 2,raw = TRUE)+factor(name),data = panel)
ggplot(panel,aes(x=Period,y=residuals(femodel),colour = factor(name)))+geom_line() +geom_smooth(method = "loess", aes(group = 1), colour = "darkblue", size = 1.2, se = FALSE)
car::crPlots(femodel)
femodelinter <- lm(value~factor(name)*poly(Period,degree = 2,raw = TRUE),data = panel)
ggplot(panel,aes(x=Period,y=residuals(femodelinter),colour = factor(name)))+geom_line() +geom_smooth(method = "loess", aes(group = 1), colour = "darkblue", size = 1.2, se = FALSE)
remodelint <- lmer(value~poly(Period,degree = 2,raw = TRUE)+(1|name),data=panel)
remodelintslope <- lmer(value~poly(Period,degree = 2,raw = TRUE)+(Period|name),data=panel)
pooledols3 <- lm(value~poly(Period,degree = 3,raw = TRUE),data = panel)
fdmodel3 <- lm(c(0,diff(value))~poly(Period,degree = 3,raw = TRUE),data = panel)
fdmodelnull3 <- lm(c(0,diff(value))~poly(Period,degree = 3,raw = TRUE),data = panel)
femodel3 <- lm(value~poly(Period,degree = 3,raw = TRUE)+factor(name),data = panel)
ggplot(panel,aes(x=Period,y=residuals(femodel3),colour = factor(name)))+geom_line() +geom_smooth(method = "loess", aes(group = 1), colour = "darkblue", size = 1.2, se = FALSE)
car::crPlots(femodel)
femodelinter3 <- lm(value~factor(name)*poly(Period,degree = 3,raw = FALSE),data = panel)
ggplot(panel,aes(x=Period,y=residuals(femodelinter3),colour = factor(name)))+geom_line() +geom_smooth(method = "loess", aes(group = 1), colour = "darkblue", size = 1.2, se = FALSE)
remodelint3 <- lmer(value~poly(Period,degree = 3,raw = TRUE)+(1|name),data=panel)
remodelintslope3 <- lmer(value~poly(Period,degree = 3,raw = FALSE)+(Period|name),data=panel)
remodelintslope3quad <- lmer(value~poly(Period,degree = 3,raw = FALSE)+(poly(Period,degree = 2,raw = FALSE)|name),data=panel)
remodelintslope3full <- lmer(value~poly(Period,degree = 3,raw = FALSE)+(poly(Period,degree = 3,raw = FALSE)|name),data=panel)
gbicp = function(model) {as.numeric(-2*logLik(model)+log(nobs(model))*sum(diag(clubSandwich::vcovCR(model, cluster = panel$name,type = "CR2")%*%solve(vcov(model)))))}
ticp = function(model) {as.numeric(-2*logLik(model)+2*sum(diag(clubSandwich::vcovCR(model, cluster = panel$name,type = "CR2")%*%solve(vcov(model)))))}
rawtakeuchiscores<- t(tibble(ticp(pooledols),
ticp(femodel),
ticp(remodelint),
ticp(pooledols3),
ticp(femodel3),
ticp(remodelint3),
ticp(remodelintslope3),
ticp(femodelinter3)))
Deltatic <- rawtakeuchiscores - min(rawtakeuchiscores)
Rellikelihood <- exp(-0.5*Deltatic)
Weight <- Rellikelihood/sum(Rellikelihood)
Results <- data.frame(rawtakeuchiscores,Deltatic,round(Rellikelihood,4),round(Weight,4))
gt(Results,rownames_to_stub = TRUE)
ticpsandwich = function(model) {as.numeric(-2*logLik(model)+2*sum(diag(solve(vcov(model))%*%sandwich(model))))}
rawtakeuchiscores<- t(tibble(ticp(pooledols),
                             ticp(remodelint),
                             ticp(pooledols3),
                             ticp(remodelint3),
                             ticp(remodelintslope3)))
Deltatic <- rawtakeuchiscores - min(rawtakeuchiscores)
Rellikelihood <- exp(-0.5*Deltatic)
Weight <- Rellikelihood/sum(Rellikelihood)
Results <- data.frame(rawtakeuchiscores,Deltatic,round(Rellikelihood,4),round(Weight,4))
gt(Results,rownames_to_stub = TRUE)
# Convert the storm name column to a clean factor variable
panel$name <- as.factor(panel$name)
gamm_model <- gam(value ~ s(Period, bs = "cr") + 
                        s(name, bs = "re") + 
                        s(name, Period, bs = "re"), 
                      data = panel, 
                      method = "REML")

gamm_modelfull <- gam(value ~ s(Period, bs = "cr") + 
                    s(name, bs = "re") + 
                    s(name, Period, bs = "fs"), 
                  data = panel, 
                  method = "REML")
library(ggeffects)
# Mixed model individual storm paths using ggeffects
mixed_predictions <- ggpredict(remodelintslope3full, terms = c("Period [all]", "name"), type = "random")

plot(mixed_predictions) + theme_minimal()+ggtitle("BLUPs for CAT 5 storm windspeeds")
library(car)
library(lme4)

# 1. Extract the fixed effects and variance-covariance matrix
beta <- fixef(remodelintslope3)
V <- clubSandwich::vcovCR(remodelintslope3full,cluster = panel$name,type = "CR2")

# 2. Assign clean, simple names that R's parser won't choke on
names(beta) <- c("b0", "b1", "b2", "b3")
colnames(V) <- rownames(V) <- c("b0", "b1", "b2", "b3")

# 3. Define the simplified formulas using our clean names
formula_root1 <- "(-2*b2 + sqrt((2*b2)^2 - 12*b3*b1)) / (6*b3)"
formula_root2 <- "(-2*b2 - sqrt((2*b2)^2 - 12*b3*b1)) / (6*b3)"

# 4. Run deltaMethod using the cleaned vector and matrix
max_res <- deltaMethod(beta, formula_root1, vcov. = V)
min_res <- deltaMethod(beta, formula_root2, vcov. = V)

# View results
print(root1_res)
print(root2_res)
# 1. Define the formula for the inflection point using the clean names
formula_inflection <- "-b2 / (3 * b3)"

# 2. Run deltaMethod
inflection_res <- deltaMethod(beta, formula_inflection, vcov. = V)

# View the estimate, SE, and Wald confidence intervals
print(inflection_res)
library(ggplot2)

# Create a lookup table using your exact outputs
hurricane_landmarks <- data.frame(
  Feature = c("Peak Intensity (RMW)", "Max Weakening Rate", "Dissipation Floor"),
  Estimate = c(19.4327, 38.3602, 57.2876),
  Lower_CI = c(15.5155, 29.7703, 44.0006),
  Upper_CI = c(23.3500, 46.9500, 70.5750),
  Color = c("#e41a1c", "#377eb8", "#4daf4a") # Red, Blue, Green
)
mixed_predictions <- ggpredict(remodelintslope3, terms = c("Period [all]", "name"), type = "random")

plot(mixed_predictions) + theme_minimal()+ggtitle("Best Linear Unbiased Predictions for CAT 5 storm windspeeds+sample mean")+ 
  stat_summary(
    aes(x = Period, y = value),
    fun = mean, 
    geom = "line", 
    color = "darkblue", 
    linewidth = 1.1,
  data=panel,inherit.aes = FALSE)+geom_hline(yintercept =0)+labs(y="Windspeed(kt)",x="6 Hour observation Period",caption ="SD(intercept)=99.2,SD(Period)=0.308,Corr(Intercept~SlopePeriod)=0.3" )+theme_light()
ggplot(panel,aes(x=Period,y=value,colour = factor(name)))+geom_line()+ 
  stat_summary(
    aes(x = Period, y = value),
    fun = mean, 
    geom = "line", 
    color = "darkblue", 
    linewidth = 1.1
  )+
  
  # 1. Break the FILL scale memory for the background bands
  new_scale_fill() +
  geom_rect(data = hurricane_landmarks, 
            aes(xmin = Lower_CI, xmax = Upper_CI, ymin = -Inf, ymax = Inf, fill = Feature), 
            alpha = 0.15, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Peak Intensity (RMW)" = "#e41a1c", 
                               "Max Weakening Rate" = "#377eb8", 
                               "Dissipation Floor" = "#4daf4a"),
                    name = "Lifecycle Landmarks") +
  
  # 2. Break the COLOR scale memory for the vertical lines
  new_scale_color() +
  geom_vline(data = hurricane_landmarks, 
             aes(xintercept = Estimate, color = Feature), 
             linetype = "dashed", size = 1, show.legend = FALSE) +
  scale_color_manual(values = c("Peak Intensity (RMW)" = "#e41a1c", 
                                "Max Weakening Rate" = "#377eb8", 
                                "Dissipation Floor" = "#4daf4a"))+ggtitle(" CAT5 Hurricane Critical point predictions from a cubic LMM(slope + intercept) 95%CIs") + labs(caption="sample mean in blue, confidence intervals for cubic critical points derived by delta method with CR2 vcovs.Cor(slope~intercept=0.3)",y="Windspeed(kt)",x="6 Hour observation Interval")+theme_light()
mixed_predictions <- ggpredict(remodelintslope3full, terms = c("Period [all]", "name"), type = "random")

plot(mixed_predictions) + theme_minimal()+ggtitle("Best Linear Unbiased Predictions(FULL RE) for CAT 5 storm windspeeds+sample mean")+ 
  stat_summary(
    aes(x = Period, y = value),
    fun = mean, 
    geom = "line", 
    color = "darkblue", 
    linewidth = 1.1,
    data=panel,inherit.aes = FALSE)+geom_hline(yintercept =0)+labs(y="Windspeed(kt)",x="6 Hour observation Period",caption ="Corr(Intercept~SlopePeriod)=0.86,Corr(Intercept~SlopePeriod^2)=-0.91,Corr(Intercept~SlopePeriod^3)=-0.86" )+theme_light()
modelsummary(remodelintslope3full,vcov = list(clubSandwich::vcovCR(remodelintslope3full,cluster = panel$name,type = "CR2")))
# 1. Setup the framework for the 36 new storms
# Let's assume you have a vector of the time periods you want to simulate across.
# Replace 'unique(panel$Period)' with whatever your sequence of periods is.
periods <- unique(panel$Period) 

new_data <- expand.grid(
  Period = periods,
  name = paste0("S_", 1:300) # Creating 36 brand new storm identifiers
)

# 2. Run the simulation
# allow.new.levels = TRUE permits simulating for 'names' not in the original data.
# re.form = NA tells lme4 to sample new random effects for the new groups.
simulated_values <- simulate(
  remodelintslope3full, 
  nsim = 1, 
  seed = 123,               # Seed for reproducibility
  newdata = new_data, 
  allow.new.levels = TRUE, 
  re.form = NA
)

# 3. Add the simulated windspeeds to the data frame
new_data$sim_windspeed <- simulated_values[[1]]

# 4. Enforce the 0 lower bound (truncate negative windspeeds)
new_data$sim_windspeed <- pmax(0, new_data$sim_windspeed)
ggplot(new_data,aes(x=Period,y=sim_windspeed,colour = factor(name)))+geom_line() +geom_smooth(method = "loess", aes(group = 1), colour = "darkblue", size = 1.2, se = FALSE)+theme(legend.position = "none")+ggtitle("300 simulated storm paths from an orthogonal  cubic random coefficients LMM trained on 6 historic CAT 5s")
summary(new_data)
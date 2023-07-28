#R RStudio 2022.07.2 Build 576
#R 4.1.2
#Tongue end/Bourne Eau Predator-prey interactions data set

#Required libraries####

library(tidyverse)
library(dunn.test)
library(car)
library(ggeffects)
library(DHARMa)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(forcats)
library(cowplot)

#1 Load data and prepare ####

##tongue_1.csv - Primary data set: sample period, month, PRE, PRE duration, attack in  PRE, number of attacks, foraging rate, predator size, predator species
## shoal size, shoal density, shoal area, switching rate
#perform data profiling on this data and plot interesting relationships
tongue_1=read_csv("./data/tongue_1.csv")
##tongue_2 - Secondary data set: attack duration, shoal behavior response, shoal behavior duration, shoal change (density, area)
tongue_2=read_csv("./data/tongue_2.csv")

#1.1 Set up labels and factors####

# Create a lookup table for variable labels
labels_table <- list(
  month = c('October', 'November', 'December'),
  sample = c('Dawn','Day','Dusk','Night'),
  action_pred = c('No attack', 'Attack'),
  sp_pred = c('Phalacrocorax carbo', 'Esox lucius'),
  shoal_beh = c('No response', 'Flee', 'Flee \n (into weedscreen)', 'Flee \n (away from weedscreen)', 
                'Avoid', 'Avoid \n (into weedscreen)', 'Avoid \n (away from weedscreen)'),
  attack_in_pre = c('No predator attack \n during PRE', 'Predator attack \n during PRE'),
  response = c('Flee', 'Avoid')
)

#1.1.1 convert_to_factors function####

#Function takes two arguments - data (DF to be converted) and labels_table (variable labels)
#Loop over each var in labels_table, check if current var is present in DF, if present extract levels
#Retrieve labels, convert var to factors, end the if statement, end the loop, return the DF

convert_to_factors <- function(data, labels_table) {
  for (var in names(labels_table)) {
    if (var %in% names(data)) {
      levels <- unique(data[[var]])
      labels <- labels_table[[var]]
      data[[var]] <- factor(data[[var]], levels = levels, labels = labels)
    }
  }
  return(data)
}

#Convert factor variables of tongue_1 and tongue_2
tongue_1 <- convert_to_factors(tongue_1, labels_table)
tongue_2 <- convert_to_factors(tongue_2, labels_table)

#2 Data exploration and summary statistics ####
#2.1 Normality testing####
shapiro.test(tongue_1$pre_dur)
qqnorm(tongue_1$pre_dur, pch = 1, frame = FALSE)
qqline(tongue_1$pre_dur, col = "steelblue", lwd = 2)

shapiro.test(tongue_1$n_attack)
qqnorm(tongue_1$n_attack, pch = 1, frame = FALSE)
qqline(tongue_1$n_attack, col = "steelblue", lwd = 2)

shapiro.test(tongue_1$foraging_rate)
qqnorm(tongue_1$foraging_rate, pch = 1, frame = FALSE)
qqline(tongue_1$foraging_rate, col = "steelblue", lwd = 2)

shapiro.test(tongue_1$switch_rate)
qqnorm(tongue_1$switch_rate, pch = 1, frame = FALSE)
qqline(tongue_1$switch_rate, col = "steelblue", lwd = 2)

shapiro.test(tongue_1$shoal_size)
qqnorm(tongue_1$shoal_size, pch = 1, frame = FALSE)
qqline(tongue_1$shoal_size, col = "steelblue", lwd = 2)

shapiro.test(tongue_1$shoal_density)
qqnorm(tongue_1$shoal_density, pch = 1, frame = FALSE)
qqline(tongue_1$shoal_density, col = "steelblue", lwd = 2)

#Data is not normally distrusted
#data is also multinomial (month, sample period) with bimodal elements (attacks, pred_sp) so is very unlikely to achieve normality


#2.2 PRE metadata####
#(sample periods, months, predator species, attack during)

#PRE count by species, % of PRE
tongue_1  %>% group_by(sp_pred) %>%  summarise(n = n())%>% mutate(rel.freq = paste0(round(100 * n/sum(n), 0), "%"))
#PRE count by month, species, % of PRE
tongue_1  %>%  group_by(month,sp_pred) %>% summarise(n = n())%>%mutate(rel.freq = paste0(round(100 * n/sum(n), 0), "%"))
#PRE count by sample, species, % of PRE
tongue_1  %>%group_by(sample,sp_pred) %>%  summarise(n = n())%>%mutate(rel.freq = paste0(round(100 * n/sum(n), 0), "%"))
#PRE duration by species
tongue_1  %>% group_by(sp_pred) %>% summarise(med = median(pre_dur),min = min(pre_dur),max = max(pre_dur),IQR = IQR(pre_dur))
#PRE duration by month, species
tongue_1  %>% group_by(month,sp_pred) %>% summarise(med = median(pre_dur), min = min(pre_dur),max = max(pre_dur), IQR = IQR(pre_dur))
#PRE duration by sample, species
tongue_1  %>% group_by(sample,sp_pred) %>%  summarise(med = median(pre_dur),min = min(pre_dur), max = max(pre_dur),IQR = IQR(pre_dur))
#PRE duration by attack status, pike
tongue_1  %>% filter(sp_pred=="Esox lucius")%>%group_by(action_pred) %>% summarise( med = median(pre_dur),min = min(pre_dur),max = max(pre_dur),IQR = IQR(pre_dur))

#PRE duration statistical

#between species
wilcox.test(tongue_1$pre_dur~tongue_1$sp_pred)
#W = 847, p-value = 0.0001112

#sample period
#corm
tongue_1 %>% filter(sp_pred=="Phalacrocorax carbo")%>% kruskal.test(data=.,pre_dur~sample)
#Kruskal-Wallis chi-squared = 0.59326, df = 2, p-value = 0.7433

#pike
tongue_1 %>% filter(sp_pred=="Esox lucius")%>%kruskal.test(data=.,pre_dur~sample)
#Kruskal-Wallis chi-squared = 14.619, df = 3, p-value = 0.002173

#post-hoc testing
tongue_1 %>% filter(sp_pred=="Esox lucius")%>%dunn.test(x=.$pre_dur, g=.$sample)

#between attack/no attack pike
tongue_1 %>% filter(sp_pred=="Esox lucius")%>%wilcox.test(data=.,pre_dur~attack_in_pre)
#W = 1946, p-value = 0.04029

#2.3 Predator metadata####
#(size, attacks, foraging rate)

#Predator size
tongue_1  %>% group_by(sp_pred) %>% summarise(med = median(size_pred),min = min(size_pred),max = max(size_pred),IQR = IQR(size_pred))
#Number of attacks predator species
tongue_1  %>% group_by(sp_pred) %>% summarise(sum=sum(n_attack),n=n())%>%mutate(rel.freq = paste0(round(100 * n/sum(n), 0), "%"))
#Number of PRE with attack, pike
tongue_1  %>% filter(sp_pred=="Esox lucius")%>% group_by(attack_in_pre) %>% summarise(sum=sum(n_attack), n=n())%>%mutate(rel.freq = paste0(round(100 * n/sum(n), 0), "%"))
#Predator foraging rate, cormorant
tongue_1  %>% filter(sp_pred=="Phalacrocorax carbo")%>%summarise(med = median(foraging_rate),min = min(foraging_rate),max = max(foraging_rate),IQR = IQR(foraging_rate))
#Predator foraging rate, cormorant, month
tongue_1  %>% filter(sp_pred=="Phalacrocorax carbo")%>%group_by(month)%>%summarise(med = median(foraging_rate),min = min(foraging_rate),max = max(foraging_rate), IQR = IQR(foraging_rate))
#Predator foraging rate, cormorant, sample
tongue_1  %>% filter(sp_pred=="Phalacrocorax carbo")%>%group_by(sample)%>%summarise( med = median(foraging_rate), min = min(foraging_rate), max = max(foraging_rate),IQR = IQR(foraging_rate))
#Predator foraging rate, pike
tongue_1  %>% filter(sp_pred=="Esox lucius")%>% filter(attack_in_pre=="Predator attack \n during PRE")%>% summarise(med = median(foraging_rate),min = min(foraging_rate),max = max(foraging_rate),IQR = IQR(foraging_rate))
##Time until attack pike PRE
tongue_1  %>% filter(sp_pred=="Esox lucius")%>% filter(attack_in_pre=="Predator attack \n during PRE")%>%  summarise( med = median(tte_pred), min = min(tte_pred), max = max(tte_pred),IQR = IQR(tte_pred))
###Attack duration
tongue_2  %>%  filter(action_pred=="Attack")%>%group_by(sp_pred)%>%summarise(med = median(action_dur),min = min(action_dur), max = max(action_dur),IQR = IQR(action_dur))

#statistical test
tongue_2  %>%  filter(action_pred=="Attack")%>%wilcox.test(data=.,action_dur~sp_pred)
#W = 1599.5, p-value = 0.007516

#Foraging rate statistical
#Correlation testing
#Add ID ROW
tongue_1$ID<-1:nrow(tongue_1)

#Corm, foraging_rate, row ID
#provide an alternative way of testing correlation
tongue_1%>%filter(sp_pred=="Phalacrocorax carbo")%>% filter(action_pred=="Attack")%>%
  cor.test(data=.,.$ID, .$foraging_rate,  method="spearman", exact=FALSE)
ggscatter(tongue_1%>%filter(sp_pred=="Phalacrocorax carbo")%>% filter(action_pred=="Attack") ,x = "ID", y = "foraging_rate", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ID", ylab = "Foraging rate")
#p-value = 0.7678  rho 0.06085365 

#Pike, foraging_rate, row ID
tongue_1%>%filter(sp_pred=="Esox lucius")%>% filter(action_pred=="Attack")%>%
  cor.test(data=.,.$ID, .$foraging_rate,  method="spearman", exact=FALSE)
ggscatter(tongue_1%>%filter(sp_pred=="Esox lucius")%>% filter(action_pred=="Attack") ,x = "ID", y = "foraging_rate", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ID", ylab = "Foraging rate")
#p-value = 0.115  rho -0.3543564

#Foraging_rate between species
wilcox.test(tongue_1$foraging_rate~tongue_1$sp_pred)
#W = 2916, p-value = 4.249e-13

#Pike
#sample period
tongue_1  %>% filter(sp_pred=="Esox lucius")%>% 
  kruskal.test(data=., foraging_rate~sample)
#Kruskal-Wallis chi-squared = 13.666, df = 3, p-value = 0.003397

#post-hoc testing
tongue_1  %>% filter(sp_pred=="Esox lucius")%>% 
  dunn.test(x=.$foraging_rate, g=.$sample)

#Corm
#sample period
tongue_1  %>% filter(sp_pred=="Phalacrocorax carbo")%>% 
  kruskal.test(data=., foraging_rate~sample)
#Kruskal-Wallis chi-squared = 0.52065, df = 2, p-value = 0.7708

#post-hoc testing
tongue_1  %>% filter(sp_pred=="Phalacrocorax carbo")%>% 
  dunn.test(x=.$foraging_rate, g=.$sample)

#2.4 Prey shoal metadata####

#Prey shoal size
tongue_1  %>% filter(!is.na(shoal_size))%>%summarise(sum=sum(shoal_size))
#Prey shoal size, month, sample
tongue_1  %>% filter(!is.na(shoal_size))%>%group_by(month,sample)%>%summarise(sum=sum(shoal_size))
#Prey shoal density, month, sample
tongue_1  %>% filter(!is.na(shoal_density))%>%group_by(month, sample)%>%summarise( med = median(shoal_density), min = min(shoal_density),max = max(shoal_density), IQR = IQR(shoal_density))
#Prey shoal density, species, attack behaviour
tongue_1  %>% filter(!is.na(shoal_density))%>%group_by(sp_pred, action_pred)%>%summarise(med = median(shoal_density),min = min(shoal_density),max = max(shoal_density),IQR = IQR(shoal_density))

#Prey shoal statistical

#shoal size months
tongue_1  %>% filter(!is.na(shoal_size))%>%
  kruskal.test(data=., shoal_size~month)
#Kruskal-Wallis chi-squared = 4.4675, df = 2, p-value = 0.1071

#shoal size sample period
tongue_1  %>% filter(!is.na(shoal_size))%>%
  kruskal.test(data=., shoal_size~sample)
#Kruskal-Wallis chi-squared = 42.692, df = 3, p-value = 2.861e-09

#shoal density between species
tongue_1%>%
  wilcox.test(data=., shoal_density~sp_pred)
#W = 1672.5, p-value = 0.01204

#shoal density pike attack, no attack
tongue_1%>%filter(sp_pred=="Esox lucius")%>%
  wilcox.test(data=., shoal_density~attack_in_pre)
#W = 548.5, p-value = 0.0001161

#Prey shoal behavioural responses
#Areal response of prey, predator species

#cormorant
tongue_2  %>%filter(!is.na(change_area))%>%filter(sp_pred=="Phalacrocorax carbo")%>%group_by(action_pred) %>% summarise( med = median(change_area),min = min(change_area),max = max(change_area),IQR = IQR(change_area))
#pike
tongue_2  %>%filter(!is.na(change_area))%>%filter(sp_pred=="Esox lucius")%>% group_by(action_pred) %>% summarise(med = median(change_area),min = min(change_area),max = max(change_area),IQR = IQR(change_area))

#Density response of prey, predator species

#cormorant
tongue_2  %>%filter(!is.na(change_density))%>%filter(sp_pred=="Phalacrocorax carbo")%>%group_by(action_pred) %>% summarise( med = median(change_density), min = min(change_density), max = max(change_density),IQR = IQR(change_density))
#pike
tongue_2  %>%filter(!is.na(change_density))%>%filter(sp_pred=="Esox lucius")%>%group_by(action_pred) %>% summarise(med = median(change_density), min = min(change_density),max = max(change_density),IQR = IQR(change_density))

#Prey response statistical

#pike areal response, attack, no attack
tongue_2 %>%filter(sp_pred=="Esox lucius")%>%
  wilcox.test(data=., change_area~action_pred)
#W = 1403.5, p-value = 3.789e-05

#pike density response, attack, no attack
tongue_2 %>%filter(sp_pred=="Esox lucius")%>%
  wilcox.test(data=., change_density~action_pred)
#W = 792, p-value = 0.301

#areal response species
tongue_2 %>%filter(action_pred=="Attack")%>%
  wilcox.test(data=., change_area~sp_pred, exact=FALSE)
#W = 510.5, p-value = 0.1614

#density response species
#areal response species
tongue_2 %>%filter(action_pred=="Attack")%>%
  wilcox.test(data=., change_density~sp_pred, exact=FALSE)
#W = 290, p-value = 0.0439

#Prey behavioural response categories 

#shoal behaviors
tongue_2 %>% filter(!is.na(shoal_beh))%>%
  group_by(shoal_beh)%>%
  summarise(
    n = n())%>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 0), "%"))

#shoal behaviors, main response categories
tongue_2 %>% filter(!is.na(response))%>%
  group_by(response)%>%
  summarise(
    n = n())%>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 0), "%"))

#shoal behaviors, cormorants
tongue_2 %>% filter(!is.na(shoal_beh))%>%filter(sp_pred=="Phalacrocorax carbo")%>%
  group_by(shoal_beh)%>%
  summarise(
    n = n())%>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 0), "%"))

#shoal behaviors, pike
tongue_2 %>% filter(sp_pred=="Esox lucius")%>%
  group_by(action_pred,shoal_beh)%>%
  summarise(
    n = n())%>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 0), "%"))

#Behaviour duration, sub categories
tongue_2  %>%filter(!is.na(beh_response_dur))%>%
  group_by(shoal_beh) %>% 
  summarise(
    med = median(beh_response_dur),
    min = min(beh_response_dur),
    max = max(beh_response_dur),
    IQR = IQR(beh_response_dur))

#Behaviour duration, main categories
tongue_2  %>%filter(!is.na(beh_response_dur))%>%group_by(response) %>% summarise( med = median(beh_response_dur), min = min(beh_response_dur),max = max(beh_response_dur), IQR = IQR(beh_response_dur))

#behaviour duration across categories
tongue_2  %>% filter(!is.na(beh_response_dur))%>%
  kruskal.test(data=., beh_response_dur~shoal_beh)
#Kruskal-Wallis chi-squared = 116.75, df = 6, p-value < 2.2e-16

#post-hoc testing
tongue_2  %>% filter(!is.na(beh_response_dur))%>%
  dunn.test(x=.$beh_response_dur, g=.$shoal_beh)

#Proportion of flee responses in cormorant attacks
prop.test(x=c(12, 26, 8), n=c(46, 46, 46),  conf.level = 0.95,)

#Proportion of flee responses in pike attacks
prop.test(x=c(10, 17, 22), n=c(49, 49, 49))

#Proportion of avoid responses in pike no attacks
prop.test(x=c(19, 75, 8, 17), n=c(119, 119, 119, 119))


#2.4.1 Weed screen switch rate####

#Primary consideration when handling data: Ensure weed screen switching is dependent on predator behavior.
#I.E., cannot simply compare attacks to no attack based on initial prey response, as no attack could include attacks which happened after first response.
# Create new data set that:
# Excludes all NA values from switch_rate (important to double check that genuine 0's are correctly recorded)
# Includes switch_rate from the initial attack behavior only to ensure comparable between predator species
# Pike switch rate for 'no attack' MUST include switch_rate from only events with NO ATTACKS throughout PRE

#Pike
switch_1 <-tongue_1  %>%filter(!is.na(switch_rate))%>%filter(sp_pred=="Esox lucius")%>%filter(action_pred=="Attack")
switch_2 <-tongue_1  %>%filter(!is.na(switch_rate))%>%filter(sp_pred=="Esox lucius")%>%filter(attack_in_pre=="No predator attack \n during PRE")
#Corm
switch_3 <-tongue_1  %>%filter(!is.na(switch_rate))%>%filter(sp_pred=="Phalacrocorax carbo")

#Convert column group from factor to numeric for row binding
#as.character first required to convert factor to character, and then to numeric
switch_1$sp_pred = as.numeric(as.character(switch_1$sp_pred))
#rename group to 1
switch_1["sp_pred"] <- 2
switch_3$sp_pred = as.numeric(as.character(switch_3$sp_pred))
#rename group to 1
switch_3["sp_pred"] <- 1

##modified for figure creation, code will need adjusting accordingly
switch_df <- rbind(switch_1, switch_2, switch_3)

#Descriptive statistics
switch_df  %>%group_by(sp_pred, action_pred)%>%summarise(med = median(switch_rate),IQR = IQR(switch_rate), n = n())

#Switch rate pike, corm attack

switch_df %>%filter(action_pred=="Attack")%>%
  wilcox.test(data=.,switch_rate~sp_pred,exact = FALSE)
#W = 277, p-value = 0.7077

#Switch rate pike, attack, no attack
switch_df %>%filter(sp_pred=="Esox lucius")%>%
  wilcox.test(data=.,switch_rate~action_pred,exact = FALSE)
#W = 290.5, p-value = 0.004362

###Correlation testing###

#Corm, PRE duration, switch rate
switch_df %>%filter(sp_pred=="Phalacrocorax carbo")%>%
  cor.test(data=.,.$pre_dur, .$switch_rate,  method="spearman", exact=FALSE)
#p-value = 0.001746 rs 0.5837079 

ggscatter(switch_df %>%filter(sp_pred=="Phalacrocorax carbo"), x = "pre_dur", y = "switch_rate", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PRE Duration", ylab = "Switrch rate")

#Pike, PRE duration, switch rate (only during attacks to keep comparability with corm)
switch_df %>%filter(sp_pred=="Esox lucius")%>%filter(action_pred=="Attack")%>%
  cor.test(data=.,.$pre_dur, .$switch_rate,  method="spearman", exact=FALSE)
#p-value = 0.03617 rs 0.4830374 

ggscatter(switch_df %>%filter(sp_pred=="Esox lucius")%>%filter(action_pred=="Attack"), x = "pre_dur", y = "switch_rate", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PRE Duration", ylab = "Switrch rate")

###PRE duration of 500 a clear outlier. Remove from switch_df data frame and rerun correlation testing
#remove outlier (PRE duration 500)
switch_df <-switch_df[-c(19),]

#repeat correlation testing for foraging rate and switch rate

#Corm, foraging_rate, switch rate
switch_df %>%filter(sp_pred=="Phalacrocorax carbo")%>%
  cor.test(data=.,.$foraging_rate, .$switch_rate,  method="spearman", exact=FALSE)
#p-value = 2.57e-06 rs 0.7804289 

ggscatter(switch_df %>%filter(sp_pred=="Phalacrocorax carbo"), x = "foraging_rate", y = "switch_rate", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "foraging_rate", ylab = "Switrch rate")

#Pike, foraging_rate, switch rate (only during attacks to keep comparability with corm)
switch_df %>%filter(sp_pred=="Esox lucius")%>%filter(action_pred=="Attack")%>%
  cor.test(data=.,.$foraging_rate, .$switch_rate,  method="spearman", exact=FALSE)
#p-value = 0.286 rs 0.2581295 

ggscatter(switch_df %>%filter(sp_pred=="Esox lucius")%>%filter(action_pred=="Attack"), x = "foraging_rate", y = "switch_rate", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "foraging_rate", ylab = "Switrch rate")

#3 GLM ####
#Prepare data set for modelling - only interested in PREs which only contain attacks, 
#thus limiting cormorants to 26 PREs, pike to 19 PREs
# Response variable is a continuous numeric value

model_df <- rbind(switch_1, switch_3)

#Set up factors
model_df$month <- factor(model_df$month, labels = c("October", "November", "December"))
model_df$sample <- factor(model_df$sample, labels = c("Dawn","Day","Dusk"))
model_df$sp_pred <- factor(model_df$sp_pred, labels=c("P carbo", "E lucius"))

model_df  %>%
  group_by(sp_pred, action_pred)%>%
  summarise(
    med = median(switch_rate),
    IQR = IQR(switch_rate),
    n = n())

#check for missing data
colSums(is.na(model_df))
#not in any of the variables used

#Check distribution of data
fit1 <- lm(model_df$switch_rate ~ model_df$foraging_rate)
hist(residuals(fit1))

qqnorm(residuals(fit1))
qqline(residuals(fit1))
shapiro.test(residuals(fit1))

#Response variable significantly differs from a normal distribution
#From here, standard LM with Gaussian distribution not possible
#Use Generalised Linear Model and consider gamma distribution, or a log linked Gaussian

#Use Brown & Forsythe test to check for Homogeneity of Variance
leveneTest(model_df$switch_rate,
            model_df$foraging_rate,
            location = c("median"),
            trim.alpha = 0.25)

#F value Pr(>F) 0.9507 0.4251, does not deviate from homogeneity

#check balance of data across categorical variables
table(model_df$switch_rate) # many more observations of 1.6 and 3.3
table(model_df$foraging_rate) # many more observations of 1.6 and 3.3
table(model_df$sp_pred) # 26 corm, 20 pike
table(model_df$sample) # 22 dawn, 9 day, 15 dusk
table(model_df$month) # 19 October, 8 November, 19 December

#Check % of zeros in response variable
sum(model_df$switch_rate == 0,
    na.rm = TRUE) * 100 / nrow(model_df)
#10% of switch_rates are 0, unlikely zero-inflation

##check for multi-collinearity in continuous variables
cor.test(model_df$foraging_rate, model_df$shoal_density, method="spearman", exact=FALSE)
plot(model_df$foraging_rate, model_df$pre_dur)

#Alternatively check Variance Inflation Factor co-lineartiy (VIF >3)
#VIF requires package car
vif(glm(switch_rate ~ sample + sp_pred + month + foraging_rate + shoal_density,
        family = poisson,
        data = model_df))

#No evidence for co-linearity
###Explore plots to determine relationships of interest
plot(model_df$switch_rate~ model_df$foraging_rate) #positive linear
plot(model_df$switch_rate~ model_df$shoal_density) #positive linear
plot(model_df$switch_rate~ model_df$sample) #dusk possibly highest, maybe interacts with sp_pred
plot(model_df$switch_rate~ model_df$sp_pred)
plot(model_df$switch_rate~ model_df$pre_dur) #strong positive linear relationship, 1 outlier present
plot(model_df$switch_rate~ model_df$month) #November lowest, but probably reduced cormorant events 

##remove outlier of PRE dur 500s
model_df <-model_df[-c(19),]

#Of these, a model containing foraging rate and duration would be interesting
#Now consider interactions 

ggplot(model_df, aes(x=sp_pred, y=switch_rate)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE) +facet_grid(~sample)

#strength of interaction between switch_rate and dusk is stronger in cormorant events
#Not enough data points to correctly test this interaction 
#The data exploration showed
#continuous non-negative response variable with 5 values at 0
# No NAs
# 1 outlier in predictor variable pre_dur, removed
# Non-normally distributed response variable, but homogeneous 
# 10% zeros in response
# No collinearity issues
# Weak interactions, insufficient data to model
# Probable independence of response variable - each PRE was a unique event

#Remove all unnecessary variables from data frame
model_df<-subset(model_df, select=-c(date,month,sample, time,pre_id, pre,attack_in_pre,n_pred, 
                                     size_pred,action_pred,tte_pred,n_attack,shoal_size,shoal_density,
                                     shoal_area,change_size,change_density,change_area,
                                     shoal_beh,tte_shoal_end))

#Given that thee are zeros in the response variable, these will need to be offset to allow the model to run
# non-positive y values not allowed in Gaussian and Gamma models
# add 0.0000000001 to response
model_df[c(2,5,11,13,22), "switch_rate"] <- 0.0000000001

model1_foraging <- glm(switch_rate ~ foraging_rate*sp_pred,family=Gamma(),
                    data = model_df)
summary (model1_foraging)

plot(ggpredict(model1_foraging, terms = c("foraging_rate","sp_pred[P carbo]")))
plot(ggpredict(model1_foraging, terms = c("foraging_rate","sp_pred[E lucius]")))

#Model outputs suggest a very poor fit with a gamma distribution
#Predicts foraging rate decreases switch_rate, does not make sense in context of observed data
#Poor estimate of species interaction with very high p values
#reject this model and consider alternative distribution 

model1_duration <- glm(switch_rate ~ pre_dur*sp_pred,family=Gamma(),
                       data = model_df)
summary (model1_duration)

#Similar results to previous model, poor prediction and alternative distribution required

#3.1 GLM species interaction####

#consider if using Gaussian with link = log will improve model fit

model2_foraging<- glm(switch_rate ~ foraging_rate*sp_pred,family=gaussian(link="log"),
                    data = model_df)
summary (model2_foraging)

#plot predicted values at 0:8, stepped by 0.1 = ~80 predicted observations
plot(ggpredict(model2_foraging, terms = c("foraging_rate[0:8, by=0.1]","sp_pred[P carbo]")))
plot(ggpredict(model2_foraging, terms = c("foraging_rate[0:8, by=0.1]","sp_pred[E lucius]")))

#Model estimates make more sense in the context of observed data. P
#positive relationship between foraging rate and weed screen switch rate
#Effect reduced in pike events compared to cormorant events

model2_duration<- glm(switch_rate ~ pre_dur*sp_pred,family=gaussian(link="log"),
                      data = model_df)
summary (model2_duration)

#Improved over gamma distribution. Model predictions much more sensible
#plot predicted values at 0:8, stepped by 0.1 = ~80 predicted observations
plot(ggpredict(model2_duration, terms = c("pre_dur[0:150, by=0.1]","sp_pred[P carbo]")))
plot(ggpredict(model2_duration, terms = c("pre_dur[0:150, by=0.1]","sp_pred[E lucius]")))

#Check model fit with DHARMA
#model2_foraging

fittedModel1 <- model2_foraging
simuout1 <- simulateResiduals(fittedModel = fittedModel1)
plot(simuout1)

#Plotted vs residuals ok
testDispersion(simuout1)
#Dispersion close to 1
#model2_duration

fittedModel2 <- model2_duration
simuout2 <- simulateResiduals(fittedModel = fittedModel2)
plot(simuout2)

#Plotted vs residuals deviance detected, but accepted
testDispersion(simuout2)
#Dispersion close to 1

#Adding a random effect of month or sample period would not improve the models given the variance determined earlier

#plot predictions
plot(ggpredict(model2_duration, terms = c("pre_dur[0:150, by=1]","sp_pred[P carbo]")))
plot(ggpredict(model2_duration, terms = c("pre_dur[0:150, by=1]","sp_pred[E lucius]")))
plot(ggpredict(model2_foraging, terms = c("foraging_rate[0:8, by=0.1]","sp_pred[P carbo]")))
plot(ggpredict(model2_foraging, terms = c("foraging_rate[0:8, by=0.1]","sp_pred[E lucius]")))

#3.2 GLM no species####

#Revisiting GLM analysis RE: T. Reid Nelson comments
##Remove species interaction
#this is model presenetd in paper

model2.1_foraging<- glm(switch_rate ~ foraging_rate+sp_pred,family=gaussian(link="log"),
                        data = model_df)

summary (model2.1_foraging)
table1 <- summary (model2.1_foraging)$coefficients
write.csv(table1, file = "model2.1_foraging.csv")

#AIC marginally reduced by 2 from removing species interaction

plot(ggpredict(model2.1_foraging, terms = c("foraging_rate[0:6, by=0.1]","sp_pred")))
plot(ggpredict(model2.1_foraging, terms = c("foraging_rate[0:8, by=0.1]","sp_pred[E lucius]")))

model2.1_duration<- glm(switch_rate ~ pre_dur+sp_pred,family=gaussian(link="log"),
                        data = model_df)

summary (model2.1_duration)
table2 <- summary (model2.1_duration)$coefficients
write.csv(table2, file = "model2.1_duration.csv")

#AIC marginally increased by 3 from removing species interaction

plot(ggpredict(model2.1_duration, terms = c("pre_dur[0:100, by=1]","sp_pred")))
plot(ggpredict(model2.1_duration, terms = c("pre_dur[0:150, by=1]","sp_pred[E lucius]")))

#Check model fit with DHARMA
#model2.1_foraging

fittedModel1 <- model2.1_foraging
simuout1 <- simulateResiduals(fittedModel = fittedModel1)
plot(simuout1)

#Plotted vs residuals ok
testDispersion(simuout1)
#Dispersion close to 1
#model2.1_duration

fittedModel2 <- model2.1_duration
simuout2 <- simulateResiduals(fittedModel = fittedModel2)
plot(simuout2)
testDispersion(simuout2)
#Dispersion close to 1

#4 Produce figures ####

#create colour vector

colours <- c("steelblue1", "grey70")

#4.1 Create theme function####

theme_JN <- function(base_size=10){ 
  theme_grey() %+replace%
    theme(
      axis.text = element_text(colour="black"),
      axis.title = element_text(colour="black"),
      axis.ticks = element_line(colour="black"),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(colour = "black",fill = NA),
      panel.spacing.x = unit(12, "pt")
    ) 
}

#4.2 Figure 1 (Figure 2 in paper): N PRE, facet by predator species, month and fill by sample period ####

prebar <- ggplot(tongue_1, aes(x=factor(attack_in_pre), y=pre, fill=factor(sample, labels=c("Dawn","Day","Dusk", "Night"))))+
  geom_bar(color="black",stat='identity',position=position_stack(reverse=TRUE), size=1)+
  geom_bar(stat='identity',position=position_stack(reverse=TRUE))+
  labs(x = "Predator attack behaviour", y = expression ("Predation Related Event (PRE) count"), fill="Sample period")+
  scale_y_continuous(breaks = seq(0, 50, 10),limits=c(0,55), expand=c(0,0))+
  scale_x_discrete(labels=c("No attack","Attack"),expand=c(0.6,0))+
  scale_fill_manual(values=c("#d1495b", "#00798c", "#edae49", "#2e4057"))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.4, "lines"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(angle = 90,hjust=0.5,colour="black", size = 10),
        strip.text.x = element_text(face = "italic"),
        legend.position = c(.85, 0.85),
        legend.title = element_blank(),
        legend.key.height = unit(0.4,"cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text(size=8),
        legend.background = element_rect(colour = 'black', fill = 'white')) +
  facet_grid(factor(month, labels = c("October", "November", "December"))~fct_rev(factor(sp_pred, labels = c("Phalacrocorax carbo (n=27)", "Esox lucius (n=120)"))),scale="free") +coord_flip()
prebar

ggsave(filename="pre_main_bar.svg", plot=prebar,device = "svg",units="cm",width=14,height=12, dpi=600)



#4.4 Figure 2 (Figure 3 in paper): Prey areal response, Prey density response, weed screen switch rate ####

#Require 3 separate plots, arranged into grid using plot_grid (cowplot)

#GRID 1: Prey areal response
####Box plot change in area no attack pike, change in area attack both species

#Create new data frame
area_1 <-tongue_2  %>%filter(!is.na(change_area))%>%filter(sp_pred=="Esox lucius")
#Check class. Stored as factor earlier
class(area_1$action_pred)
#Convert back to numeric for this data frame
area_1$action_pred<- as.numeric(area_1$action_pred)
#Pike No attack = 1, Attack = 2
area_2 <-tongue_2  %>%filter(!is.na(change_area))%>%filter(sp_pred=="Phalacrocorax carbo")
#Convert action_pred to numeric, label 3
area_2["action_pred"] <- 3
class(area_2$action_pred)
area_df <- rbind(area_1, area_2)
##Pike No attack = 1, Attack = 2, Corm attack = 3

###Create summary statistics table
stat.test1 <- area_df %>%
  wilcox_test(change_area ~ action_pred) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()
stat.test1

#adding xy values
stat.test1 <- stat.test1 %>% add_xy_position(x = "action_pred")
#removing unwanted rows
stat.test1 <-stat.test1[-c(2),]
#adjusting y position manually 
stat.test1["y.position"] <- c(4.5,4.5)

#Create x labels

x_titles <- c(
  expression(atop(italic("E. lucius"), "no attack")),
  expression(atop(italic("E. lucius"), "attack")),
  expression(atop(italic("P. carbo"), "attack")))

#GRID 1: Prey areal response ggplot

changearea <- ggplot(area_df, aes(x = factor(action_pred), y = change_area))+ 
  stat_boxplot(geom = "errorbar", width=0.05, position = position_dodge(0.9))+
  geom_boxplot(outlier.shape=NA, coef = 0, width=0.5, position = position_dodge(0.9)) +
  scale_y_continuous(breaks = seq(-5, 4, 1),limits=c(-5.5, 5)) +
  scale_x_discrete(breaks = c(1, 2, 3), labels = x_titles)+
  labs(y = expression (atop("Prey areal response", ~ (Delta~"area m" ^2))))+
  scale_color_grey() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,1,0,1), units="mm"),
        axis.text.y=element_text(colour="black", size = 10)) +
  scale_fill_brewer(palette = "Greys")+
  stat_pvalue_manual(stat.test1)
changearea

#GRID 2: Prey density response

####Box plot change in density no attack pike, change in density attack both species

###Create summary statistics table
stat.test2 <- area_df %>%
  wilcox_test(change_density ~ action_pred) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()
stat.test2

#adding xy values
stat.test2 <- stat.test2 %>% add_xy_position(x = "action_pred")
#removing unwanted rows
stat.test2 <-stat.test2[-c(2),]
#adjusting y position manually
stat.test2["y.position"] <- c(31,31)

changedensity <- ggplot(area_df, aes(x = factor(action_pred), y = change_density))+ 
  stat_boxplot(geom = "errorbar", width=0.05, position = position_dodge(0.9))+
  geom_boxplot(outlier.shape=NA, coef = 0, width=0.5, position = position_dodge(0.9)) +
  scale_y_continuous(breaks = seq(-25, 25, 5),limits=c(-28, 33)) +
  scale_x_discrete(breaks = c(1, 2, 3), labels = x_titles)+
  labs(y = expression (atop("Prey density response", ~ (Delta~"density m" ^2))))+
  scale_color_grey() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,1,0,1), units="mm"),
        axis.text.y=element_text(colour="black", size = 10)) +
  scale_fill_brewer(palette = "Greys")+
  stat_pvalue_manual(stat.test2)
changedensity

#GRID 3: Prey weed screen switch rate

####Box plot switch rate no attack pike, switch rate attack both species

#convert attack factor to numerical for figure 
switch_df$action_pred<- as.numeric(switch_df$action_pred)
#Pike No attack = 1, Attack = 2
#Convert corm attacks labelled as 2, to 3
switch_df[c(71:96), "action_pred"] <- 3
##Pike No attack = 1, Attack = 2, Corm attack = 3

###Create summary statistics table
stat.test3 <- switch_df %>%
  wilcox_test(switch_rate ~ action_pred, paired=FALSE) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()
stat.test3

#adding xy values
stat.test3 <- stat.test3 %>% add_xy_position(x = "action_pred")
#removing unwanted rows
stat.test3 <-stat.test3[-c(2),]
#adjusting y position manually
stat.test3["y.position"] <- c(11,11)

switch_rate_plot <- ggplot(switch_df, aes(x = factor(action_pred), y = switch_rate))+ 
  stat_boxplot(geom = "errorbar", width=0.05, position = position_dodge(0.9))+
  geom_boxplot(outlier.shape=NA, coef = 0, width=0.5, position = position_dodge(0.9)) +
  scale_y_continuous(breaks = seq(0, 10, 2),limits=c(0, 12)) +
  scale_x_discrete(breaks = c(1, 2, 3), labels = x_titles)+
  labs(x = "Predator species and attack behaviour", y = expression (atop("Switch rate",~("switches?minute" ^-1))))+
  scale_color_grey() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,1,1,1), units="mm"),#top, right, bottom, left
        axis.text.y=element_text(colour="black", size = 10),
        axis.text.x=element_text(colour="black", size = 10)) +
  scale_fill_brewer(palette = "Greys")+
  stat_pvalue_manual(stat.test3)
switch_rate_plot

#Arrange GRID 1, GRID 2, GRID 3 with cowplot, add plot labels

pre_response_bind <-plot_grid(changearea, changedensity, switch_rate_plot,
                              ncol = 1, nrow = 3, rel_heights = c(3.25,3,3.75),align = "v") + 
  draw_plot_label(label = c("a)", "b)", "c)"), 
                  size = 10,
                  x = c(0.18, 0.18, 0.18), 
                  y = c(0.98, 0.65, 0.35)) 
pre_response_bind

ggsave(filename="pre_response_bind.svg", plot=pre_response_bind, device = "svg",units="cm", width=10,height=14, dpi=600)

#4.5 Figure 3 (Figure 4 in paper) GLM####

#Get predicted model fit line using ggpredict and store as dataframe

modforc_dur <-ggpredict(model2.1_duration, terms = c("pre_dur[0:150, by=1]","sp_pred[P carbo]"))
modforp_dur <-ggpredict(model2.1_duration, terms = c("pre_dur[0:150, by=1]","sp_pred[E lucius]"))
modforc <-ggpredict(model2.1_foraging, terms = c("foraging_rate[0:8, by=0.1]","sp_pred[P carbo]"))
modforp <-ggpredict(model2.1_foraging, terms = c("foraging_rate[0:8, by=0.1]","sp_pred[E lucius]"))

#Convert column group from factor to numeric for row binding
#as.character first required to convert factor to character, and then to numeric
modforc_dur$group = as.numeric(as.character(modforc_dur$group))
#rename group to 1
modforc_dur["group"] <- 1
#repeat for pike model
modforp_dur$group = as.numeric(as.character(modforp_dur$group))
modforp_dur["group"] <- 2
#bind data frames
predator_PRE_dur_mod <- bind_rows(modforc_dur, modforp_dur)

#Follow same process for attack rate

modforc$group = as.numeric(as.character(modforc$group))
modforc["group"] <- 1
modforp$group = as.numeric(as.character(modforp$group))
modforp["group"] <- 2
predator_attack_mod <- bind_rows(modforc, modforp)

PREduration_mod <-ggplot(predator_PRE_dur_mod, aes(x=x, y=predicted, color=as.factor(group),fill=as.factor(group)))+
  geom_jitter(model_df, height = 0, width = 0, shape = 21, inherit.aes = FALSE,mapping=aes(x=pre_dur, y=switch_rate, fill=as.factor(sp_pred)))+
  geom_line(lwd=1)+
  geom_ribbon(aes(x=x,ymin=conf.low, ymax=conf.high),alpha=0.3, colour="black",linetype=0)+
  scale_fill_manual(values=colours)+
  scale_colour_manual(values=colours)+
  scale_y_continuous(breaks = seq(0, 14, 2) ,limits=c(0,17), expand=c(0.02,0), oob = scales::squish) +
  scale_x_continuous(breaks = seq(0,150,30), limits=c(0,150), expand=c(0,0),  oob = scales::squish) +
  labs(x = expression ("PRE Duration s" ^-1))+
  coord_cartesian(ylim=c(0,14),xlim=c(0,150), clip="off") +
  theme_JN()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(5.5,10,5.5,5.5, "pt"))
PREduration_mod 

attack_mod <-ggplot(predator_attack_mod, aes(x=x, y=predicted, color=as.factor(group),fill=as.factor(group)))+
  geom_jitter(model_df, shape = 21, height = 0, width = 0, inherit.aes = FALSE,mapping=aes(x=foraging_rate, y=switch_rate,fill=as.factor(sp_pred)))+
  geom_line(lwd=1)+
  geom_ribbon(aes(x=x,ymin=conf.low, ymax=conf.high),alpha=0.3,linetype=0)+
  scale_fill_manual(values=colours)+
  scale_colour_manual(values=colours)+
  scale_y_continuous(breaks = seq(0, 14, 2) ,limits=c(0,22), expand=c(0.02,0), oob = scales::squish) +
  scale_x_continuous(breaks = seq(0,8,1), limits=c(0,10), expand=c(0,0), oob = scales::squish) +
  coord_cartesian(ylim=c(0,14),xlim=c(0,8),clip = "off") +
  labs(x = expression ("Attack rate" ~ ("attacks·minute" ^-1)), y = expression("Switch rate" ~ ("switches·minute" ^-1)))+
  theme_JN()+
  theme(legend.position = "none",
        plot.margin = margin(5.5,10,5.5,5.5, "pt"))
attack_mod

combined_mod <-plot_grid(attack_mod, PREduration_mod,
                            ncol = 2, nrow = 1, rel_widths = c(5.3,4.7),align = "h") + 
  draw_plot_label(label = c("a)", "b)", "CI = 95%"), 
                  size = 10,
                  x = c(0.10, 0.55, 0.8), 
                  y = c(0.97, 0.97, 0.97)) 
combined_mod
ggsave(filename="combined_mod.svg", plot=combined_mod, device = "svg",units="cm", width=14,height=7)


#4.3 Figure 4 (Figure 5 in paper): N behavior responses, facet by predator species and fill by attack behavior####

shoalbehbar <- ggplot(tongue_2%>%filter(!is.na(shoal_beh)), aes(x=shoal_beh, fill=action_pred))+
  geom_bar(color="black",stat='count',position="stack")+
  labs(x = "Prey shoal behaviour response", y = "Behaviour count")+
  scale_y_continuous(breaks = seq(0, 80, 20),limits=c(0,85), expand=c(0,0))+
  scale_fill_manual(values=c("forestgreen", "red4"))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.4, "lines"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.position = c(.85, 0.85),
        legend.title = element_blank(),
        legend.key.height = unit(0.4,"cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text(size=8),
        legend.background = element_rect(colour = 'black', fill = 'white'),
        strip.text.x = element_text(face = "italic"))+
  facet_grid(~fct_rev(factor(sp_pred)))+coord_flip()
shoalbehbar

ggsave(filename="shoal_beh_bar.svg", plot=shoalbehbar,device = "svg",units="cm",width=14,height=8, dpi=600)


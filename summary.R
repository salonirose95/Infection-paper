library(readxl)
library(plyr)
library(xlsx)
library(fmsb)
library(ggplot2)
library(tidyr)
library(grid)
library(reshape2)
library(dplyr)
library(readxl)
library(stringr)
library(matrixStats)
library(gplots)
library(ggsignif)
library(ggpubr)
library(gridExtra)
library(ggsci)
library(extrafont)
library(FSA)
library(survival)
library(survminer)
library(betareg)
library(multcomp)
library(gamlss)
library(glmmADMB)
library(DHARMa)



df <-read_excel("Data/Compiled data/Figure 3(Competitve mating assay)-data.xlsx", sheet = "Serratia-CMA-Infected+PBSfemale")

sum(str_count(na.omit(df$`first copulated`), pattern = "C"))
sum(str_count(na.omit(df$`first copulated`), pattern = "E"))
setwd("/Users/ssr978/Documents/Revision Infection")
df <- read_excel("Data/Compiled data/ECC15-Cs-male-diff time points.xlsx",sheet="total")

dunnTest(CoD ~ Treatment, df,method = "bonferroni")



df <- read_excel("Data/Compiled data/Figure S3-Dahomey-females.xlsx", 
                 sheet = "S.aureus")

df <- read_excel("~/Documents/Revision Infection/Data/Compiled data/ECC15-Cs-male-diff time points.xlsx", 
                 sheet = "total")
df <- read_excel("Data/Compiled data/Serratia-Cs-male-diff time points.xlsx", 
                 sheet = "total")


summarydf <-df %>% group_by(Treatment) %>%  summarise(Median = median(CoL, na.rm = TRUE), IQR=IQR(CoL, na.rm = TRUE),Sample_size=length(CI)) %>% as.data.frame()

kruskal.test(CI ~ Treatment,df)
a1<-dunnTest(CI ~ Treatment, df,method = "bonferroni")
a2<-a1$res

df <-DF

summarydf <-df %>% group_by(Treatment) %>%  summarise(Median = median(`Walking speed`, na.rm = TRUE), IQR=IQR(`Walking speed`, na.rm = TRUE),Sample_size=length(`Walking speed`)) %>% as.data.frame()
summarydf <-df %>% group_by(Treatment) %>%  summarise(Median = median(`Distance travelled`, na.rm = TRUE), IQR=IQR(`Distance travelled`, na.rm = TRUE),Sample_size=length(`Distance travelled`)) %>% as.data.frame()
kruskal.test(`Walking speed` ~ Treatment,df)
a1<-dunnTest(`Walking speed` ~ Treatment, df,method = "bonferroni")
a2<-a1$res


a <- df[df$Treatment=="Uninfected",]
b <- df[df$Treatment== "PBS",]
c <- df[df$Treatment== "S.aureus",]
pairwise.fisher.test(c(sum(is.na(a$CoL)), sum(is.na(b$CoL)), sum(is.na(c$CoL))), c(nrow(a), nrow(b), nrow(c)),p.adjust.method="bonferroni")

a1 <- (1- sum(is.na(a$CoL))/nrow(a)) *100
b1 <- (1- sum(is.na(b$CoL))/nrow(b))*100 
c1 <- (1- sum(is.na(c$CoL))/nrow(c))*100 

df <- read_excel("Data/Compiled data/Figure S2-Dahomey-males.xlsx", 
                 sheet = "L.monocytogenes")
df1 <- subset(df, select = -c(2,4,5))
model1 <- gamlss(`CI-old` ~ factor(Treatment)+ Time, data =df1,family=BEINF)
a1<-summary(model1)

model1$residuals
emm1 <- emmeans(model1, pairwise ~ Treatment,adjust = "bonf")
a2<- as.data.frame(emm1$contrasts)
df <- read_excel("~/Documents/Revision Infection/Data/Compiled data/Figure S1-CS-Survival-curve.xlsx", 
                 sheet = "Microccocus-female")
survdiff(Surv(Time, Event) ~ Treatment, data=df)
a1 <-(pairwise_survdiff(Surv(Time, Event) ~ Treatment, df, p.adjust.method = "bonferroni"))

a2<- a1$p.value

df <- read_excel("Data/Compiled data/serratia-locomotion.xlsx", 
                 sheet = "total")

model1 <- glm(df$`Distance travelled` ~ factor(Treatment)+ Time, data =df)
a1<-summary(model1)


cox <- coxph(Surv(Time, Event) ~ Treatment, data=df)    
summary(cox)
df<-read_excel("~/Documents/Revision Infection/Data/Compiled data/Figure 2(CS-female infected)-data.xlsx", 
               sheet = "L.monocytogenes")
df <- read_excel("Data/Compiled data/Figure 4(Immune system activation)-data.xlsx", 
                 sheet = "fatbody-toll-females")

df <- read_excel("Data/Compiled data/Figure S2-Dahomey-males.xlsx", 
                 sheet = "L.monocytogenes")

df <- read_excel("~/Documents/Revision Infection/Data/Compiled data/Figure S3-Dahomey-females.xlsx", 
                 sheet = "L.mono")I


g <- df[df$Treatment=="Uninfected",]
h <- df[df$Treatment== "PBS",]
i <- df[df$Treatment== "CO2",]

g1 <- (1- sum(is.na(g$CoL))/nrow(g)) *100
h1 <- (1- sum(is.na(h$CoL))/nrow(h))*100 
i1 <- (1- sum(is.na(i$CoL))/nrow(i))*100 



df <- read_excel("Data/Compiled data/Figure S1-CS-Survival-curve.xlsx", 
                 sheet = "Microccocus-male")
ss <-df %>% group_by(Treatment) %>%  summarise(Sample_size=length(Treatment)) %>% as.data.frame()





library(emmeans)
df<-read_excel("~/Documents/Revision Infection/Data/Compiled data/Figure 1(CS-male infected)-data.xlsx", 
               sheet = "L.monocytogenes")




model2<- glm(Mating_Sucess ~ factor(Treatment),df, family="quasibinomial")
summary(model2)
emmeans(model2, pairwise ~ Treatment)





model1 <- gamlss(CI_frac ~ factor(Treatment)+ OD+Time+OD*Time, data =df1,family=BEINF)
summary(model1)
contr <- matrix(0, nrow = 3, ncol = length(coef(model1)))
colnames(contr) <- names(coef(model1))
rownames(contr) <- c("med - low", "high - low", "high - med")
contr[, 2:3] <- rbind(c(1, 0), c(0, 1), c(-1, 1))
contr[, 1:5]

pw<-glht(model2, linfct = contr)
summary(pw)

model2 <- lm(CI_frac ~ factor(Treatment), data =df)
plot(fitted(model1), residuals(model1))
coef(model1)
 

a1 <-glmmadmb(df1$`CI-old` ~ Treatment, data =df1)
summary(a1)


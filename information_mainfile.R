################################################################################
################################################################################
## Data Preparation and Analysis for 
## Impact on info channel on vaccination decision for CoV-2
## (part of Study SOEP-RKI-2)
##
## Sep. 2022
## SZinn
## joint work with Susanne Jordan & Sarah Jane Böttger (RKI)
##
################################################################################
################################################################################

rm(list=ls())
library(haven)
library(naniar)
library(mice)
library(miceadds)
library(survey)
library(sandwich)
library(lmtest)
library(stargazer)
library(ggplot2)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Read Data
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

setwd("C:\\Users\\Freddie\\Documents\\RKI-Studie\\Welle_2\\DATA")
DAT <- read_dta("CoMobu2_isced_LOCF.dta")
setwd("C:\\Users\\Freddie\\Documents\\RKI-Studie\\Welle_2\\DATA\\RKISOEP2_Gewichte und Dokumentation")
G <- read_dta("RKISOEP2_phrf_21_v2.dta")
table(DAT$pid %in% G$pid)
DAT <- merge(DAT, G[,-1], by="pid", all.x=TRUE)

# Ab 18J
table(DAT$agegrp17C, exclude=NULL) # lowest group <18y
DAT <- DAT[DAT$agegrp17C %in% c(1:4),] # N=10448 with age higher than 17y
table(DAT$agegrp17C, exclude=NULL) 

# Transform IDs to avoid numerical issues
DAT <- DAT[order(DAT$pid, DAT$hid),]
DAT$id <- 1:nrow(DAT)
hidds <- as.numeric(unique(DAT$hid))
hidm <- cbind(sort(hidds), 1:length(hidds))
colnames(hidm) <- c("hid","idh")
dim(DAT)
DAT <- merge(DAT, hidm, by="hid", all.x=TRUE)
DAT <- DAT[order(DAT$id, DAT$idh),]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Data Preparation / Variables
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Focal Variables
# ------------------------------------------------------------------------------
# Mind einmal geimpft
table(DAT$pcovimpf1, DAT$pcovimpf_n_kurz, exclude=NULL)
table(DAT[is.na(DAT$pcovimpf1),"pcovimpf_n_kurz"], exclude=NULL)
DAT[is.na(DAT$pcovimpf1), "pcovimpf1"] <- DAT[is.na(DAT$pcovimpf1),"pcovimpf_n_kurz"]
table(DAT$pcovimpf1, exclude=NULL)
DAT$vacc1 <- ifelse(is.na(DAT$pcovimpf1), NA, ifelse(DAT$pcovimpf1 %in% 1, 1, 0))

# Variable *Grundimmunisierung* zum Zeitpunkt der Erhebung
# 2x geimpft || 1x geimpft und genesen
table(DAT$IPcovB_g, exclude=NULL)
DAT$GrundImm <- ifelse(is.na(DAT$IPcovB_g), NA, ifelse(DAT$IPcovB_g %in% 1, 1, 0))
table(DAT$GrundImm, exclude=NULL) # miss 4.4%

# Impfentscheidung v1
table(DAT$prki2iabsi, exclude=NULL)
table(DAT$GrundImm, DAT$prki2iabsi, exclude=NULL)
DAT1 <- DAT[DAT$GrundImm %in% c(1,NA),]
DAT1$ImpEntPos <- ifelse(is.na(DAT1$GrundImm), NA, 1)
DAT2 <- DAT[DAT$GrundImm %in% 0,]
DAT2$ImpEntPos <- ifelse(DAT2$prki2iabsi %in% c(1,2), 1, ifelse(DAT2$prki2iabsi %in% c(3,4,5), 0, NA))
DAT <- rbind.data.frame(DAT1, DAT2)
table(DAT$ImpEntPos, exclude=NULL)
#    0    1   <NA> 
#   499 9233  716 
table(DAT$GrundImm,DAT$ImpEntPos, exclude=NULL) # miss 6.9%

# Impfentscheidung v2
table(DAT$pcovimpf1, exclude=NULL)
table(DAT$pcovimpf1, DAT$prki2iabsi, exclude=NULL)
DAT1 <- DAT[DAT$pcovimpf1 %in% c(1,NA),] 
DAT1$ImpEntPos_v2 <- ifelse(is.na(DAT1$pcovimpf1), NA, 1)
DAT2 <- DAT[DAT$pcovimpf1 %in% 2,]
DAT2$ImpEntPos_v2 <- ifelse(DAT2$prki2iabsi %in% c(1,2), 1, ifelse(DAT2$prki2iabsi %in% c(3,4,5) | DAT2$pcovimpf1 %in% 2, 0, NA))
DAT <- rbind.data.frame(DAT1, DAT2)
table(DAT$ImpEntPos_v2, exclude=NULL)
#   0    1    <NA> 
#  1737 7939  772 
table(DAT$pcovimpf1,DAT$ImpEntPos_v2, exclude=NULL) # miss 7.4%

# Wichtigkeit Informationsquelle 
table(DAT$prki2iueb01, exclude=NULL) # Familie / Freunde
DAT$infoFam <- ifelse(DAT$prki2iueb01 %in% -4, NA, DAT$prki2iueb01)
table(DAT$prki2iueb02, exclude=NULL) # Ärzte
DAT$infoDoc <- ifelse(DAT$prki2iueb02 %in% -4, NA, DAT$prki2iueb02)
table(DAT$prki2iueb03, exclude=NULL) # Apotheke
DAT$infoApo <- ifelse(DAT$prki2iueb03 %in% -4, NA, DAT$prki2iueb03)
table(DAT$prki2iueb04, exclude=NULL) # TV/Radio
DAT$infoOef <- ifelse(DAT$prki2iueb04 %in% -4, NA, DAT$prki2iueb04)
table(DAT$prki2iueb05, exclude=NULL) # Zeitung / Online news
DAT$infoNews <- ifelse(DAT$prki2iueb05 %in% -4, NA, DAT$prki2iueb05)
table(DAT$prki2iueb06, exclude=NULL) # soc. Media
DAT$infoSozM <- ifelse(DAT$prki2iueb06 %in% -4, NA, DAT$prki2iueb06)
table(DAT$prki2iueb08, exclude=NULL) # Adm
DAT$infoAdm <- ifelse(DAT$prki2iueb08 %in% -4, NA, DAT$prki2iueb08)
table(DAT$prki2iueb09, exclude=NULL) # Krankenkassen
DAT$infoKK <- ifelse(DAT$prki2iueb09 %in% -4, NA, DAT$prki2iueb09)
table(DAT$prki2iueb10, exclude=NULL) # Gesundheitsportale
DAT$infoGP <- ifelse(DAT$prki2iueb10 %in% -4, NA, DAT$prki2iueb10)

# ------------------------------------------------------------------------------
# Controls
# ------------------------------------------------------------------------------
# Education
table(DAT$pgisced11_LOCF, DAT$pgisced11, exclude=NULL)
DAT$edu <- ifelse(is.na(DAT$pgisced11_LOCF), NA, 
             ifelse(DAT$pgisced11_LOCF %in% c(6,7,8), "high", 
                    ifelse(DAT$pgisced11_LOCF %in% c(3,4,5), "med", "low")))
table(DAT$edu, exclude=NULL)

# Gender
table(DAT$psex, exclude=NULL)
DAT$sexFem <- ifelse(DAT$psex %in% 2, 1, 0)

# Subj. Health
DAT$sjHealth <- ifelse(is.na(DAT$ple0008) | DAT$ple0008 %in% c(-4), NA, DAT$ple0008) 
table(DAT$sjHealth, exclude=NULL)

# Subj. Inf
table(DAT$prki2infor2, exclude=NULL)
DAT$sjInf <- ifelse(is.na(DAT$prki2infor2) | DAT$prki2infor2 %in% c(-4), NA, DAT$prki2infor2) 
table(DAT$sjInf, exclude=NULL)

# Worries
table(DAT$prki2krsor, exclude=NULL)
DAT$sjWorr <- ifelse(is.na(DAT$prki2krsor) | DAT$prki2krsor %in% c(-4), NA, DAT$prki2krsor) 
table(DAT$sjWorr, exclude=NULL)

# Time Dummy
table(DAT$datj_n2, exclude=NULL)

# How many people in HHs
XX <- aggregate(DAT$id, by=list(DAT$idh),FUN=length); table(XX$x) 

# ------------------------------------------------------------------------------
# Select Focals and Controls
# ------------------------------------------------------------------------------
D <- DAT[, c("id", "idh", "phrf_full", 
             "vacc1", "GrundImm", "ImpEntPos", "ImpEntPos_v2",
             "infoFam", "infoDoc", "infoApo", "infoOef", "infoNews", "infoSozM", "infoAdm", "infoKK", "infoGP",
             "edu", "sexFem", "agegrp17C", 
             "sjHealth", "sjInf", "sjWorr",  
             "datj_n2"
             )]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Check Missing Pattern
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Study missingness pattern
missP <- md.pattern(D, plot=F)
round(missP[nrow(missP),]/nrow(D)*100,2) # 7.5% miss on mind. 1x vaccinated
table(complete.cases(D))/nrow(D) # 77% complete cases, loose 23% due to NA

# Plot missingness pattern
gg_miss_upset(D)

# Test for MCAR
D.test <- D[, !(colnames(D)%in% "edu")]
D.test$eduLow <- ifelse(is.na(D$edu), NA, ifelse(D$edu %in% "low", 1,0))
D.test$eduMed <- ifelse(is.na(D$edu), NA, ifelse(D$edu %in% "med", 1,0))
mcar_test(D.test[,-c(1:3)]) # H_0 is MCAR: rejected

set.seed(123, sample.kind = "Rejection") 
D.imp <- as.matrix(D.test)
imp <- mice(D.imp, m=1, maxit=1)
meth <- imp$method
meth[meth!=""] <- "cart"
pred <- imp$predictorMatrix
pred[, c("id", "idh")] <- 0
imp <- mice(D.imp, predictorMatrix = pred, m=30, seed=2388)

# ------------------------------------------------------------------------------
# A. Models: one vaccination only 
# - logit model, weighted, cluster robust std err, imputed 
# ------------------------------------------------------------------------------
# for each channel, weigthed and unweighted, logit
setwd("C:\\Users\\Freddie\\Documents\\RKI-Studie\\Welle_2\\Syntax")
source("models.R")
D_cc <- complete(imp, action=1)

# 1. Info Fam
eq <- as.formula(vacc1 ~ infoFam + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                  + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res1_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res1 <- fitModel(imp, eq, D_cc, we=FALSE)

# 2. Info Doc
eq <- as.formula(vacc1 ~ infoDoc + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res2_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res2 <- fitModel(imp, eq, D_cc, we=FALSE)

# 3. Info Apo
eq <- as.formula(vacc1 ~ infoApo + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res3_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res3 <- fitModel(imp, eq, D_cc, we=FALSE)

# 4. Info infoOef
eq <- as.formula(vacc1 ~ infoOef + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res4_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res4 <- fitModel(imp, eq, D_cc, we=FALSE)

# 5. Info infoNews
eq <- as.formula(vacc1 ~ infoNews + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res5_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res5 <- fitModel(imp, eq, D_cc, we=FALSE)

# 6. Info infoSozM
eq <- as.formula(vacc1 ~ infoSozM + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res6_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res6 <- fitModel(imp, eq, D_cc, we=FALSE)

# 7. Info infoAdm
eq <- as.formula(vacc1 ~ infoAdm + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res7_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res7 <- fitModel(imp, eq, D_cc, we=FALSE)

# 8. Info infoKK
eq <- as.formula(vacc1 ~ infoKK + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res8_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res8 <- fitModel(imp, eq, D_cc, we=FALSE)

# 9. Info infoGP
eq <- as.formula(vacc1 ~ infoGP + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res9_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res9 <- fitModel(imp, eq, D_cc, we=FALSE)

# Plot it
vars <- c("info channel", "edu low vs. high", "edu med vs. high", "female", 
          "35-49 vs. 18-34", "50-64 vs. 18-34", "65+ vs. 18-34",
          "subj. health", "subj. inform", "worried Corona", "2022 vs. 2021")
resA <- rbind.data.frame(res1[-1,], res1_ww[-1,])
resA$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res1)-1))
resA$Category <- "Family/Friends"
resA$Variables <- rep(vars, 2)
resB <- rbind.data.frame(res2[-1,], res2_ww[-1,])
resB$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res2)-1))
resB$Category <- "MedDoc"
resB$Variables <- rep(vars, 2)
resC <- rbind.data.frame(res3[-1,], res3_ww[-1,])
resC$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res3)-1))
resC$Category <- "Pharm"
resC$Variables <- rep(vars, 2)
resD <- rbind.data.frame(res4[-1,], res4_ww[-1,])
resD$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res4)-1))
resD$Category <- "PublicMat"
resD$Variables <- rep(vars, 2)
resE <- rbind.data.frame(res5[-1,], res5_ww[-1,])
resE$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res5)-1))
resE$Category <- "NewsP"
resE$Variables <- rep(vars, 2)
resF <- rbind.data.frame(res6[-1,], res6_ww[-1,])
resF$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res6)-1))
resF$Category <- "SocMed"
resF$Variables <- rep(vars, 2)
resG <- rbind.data.frame(res7[-1,], res7_ww[-1,])
resG$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res7)-1))
resG$Category <- "AdmMat"
resG$Variables <- rep(vars, 2)
resH <- rbind.data.frame(res8[-1,], res8_ww[-1,])
resH$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res8)-1))
resH$Category <- "HealthIns"
resH$Variables <- rep(vars, 2)
resI <- rbind.data.frame(res9[-1,], res9_ww[-1,])
resI$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res9)-1))
resI$Category <- "HealthWeb"
resI$Variables <- rep(vars, 2)

res <- rbind.data.frame(resA,resB,resC,resD,resE,resF,resG,resH,resI)
res$Estimate <- as.numeric(res$Estimate)
res$CI_low <- as.numeric(res$CI_low)
res$CI_up <- as.numeric(res$CI_up)
res$Category <- as.factor(res$Category)
res$Variables <- factor(res$Variables, levels=rev(vars))

zp <- ggplot(res,aes(colour=Type)) 
zp <- zp + geom_pointrange(aes(x = Variables, y = exp(Estimate), ymin = exp(CI_low), ymax = exp(CI_up)),
                             lwd = 0.1, position = position_dodge(width = 1/2)) 
zp <- zp + facet_grid(cols=vars(Category)) + scale_color_manual(values=c("darkred", "darkblue"))
zp <- zp + geom_hline(yintercept=1, color = "grey10")
zp <- zp + xlab("") + ylab("Odd-Ratio") 
zp <- zp + ggtitle("At Least Once Vaccinated") 
zp <- zp + coord_flip()
zp

resR <- res[res$Variables %in% "info channel",]
zp <- ggplot(resR,aes(colour=Type)) 
zp <- zp + geom_pointrange(aes(x = Category, y = exp(Estimate), ymin = exp(CI_low), ymax = exp(CI_up)),
                           lwd = 0.7, position = position_dodge(width = 1/2)) 
zp <- zp + geom_hline(yintercept=1, color = "grey10")
zp <- zp + xlab("") + ylab("Odd-Ratio") 
zp <- zp + ggtitle("At Least Once Vaccinated (under controls)") 
zp <- zp + coord_flip()
zp

# ------------------------------------------------------------------------------
# B. Models: basic immunisation
# - logit model, weighted, cluster robust std err, imputed 
# ------------------------------------------------------------------------------
# for each channel, weigthed and unweighted, logit

# 1. Info Fam
eq <- as.formula(GrundImm ~ infoFam + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res1_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res1 <- fitModel(imp, eq, D_cc, we=FALSE)

# 2. Info Doc
eq <- as.formula(GrundImm ~ infoDoc + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res2_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res2 <- fitModel(imp, eq, D_cc, we=FALSE)

# 3. Info Apo
eq <- as.formula(GrundImm ~ infoApo + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res3_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res3 <- fitModel(imp, eq, D_cc, we=FALSE)

# 4. Info infoOef
eq <- as.formula(GrundImm ~ infoOef + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res4_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res4 <- fitModel(imp, eq, D_cc, we=FALSE)

# 5. Info infoNews
eq <- as.formula(GrundImm ~ infoNews + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res5_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res5 <- fitModel(imp, eq, D_cc, we=FALSE)

# 6. Info infoSozM
eq <- as.formula(GrundImm ~ infoSozM + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res6_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res6 <- fitModel(imp, eq, D_cc, we=FALSE)

# 7. Info infoAdm
eq <- as.formula(GrundImm ~ infoAdm + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res7_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res7 <- fitModel(imp, eq, D_cc, we=FALSE)

# 8. Info infoKK
eq <- as.formula(GrundImm ~ infoKK + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res8_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res8 <- fitModel(imp, eq, D_cc, we=FALSE)

# 9. Info infoGP
eq <- as.formula(GrundImm ~ infoGP + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res9_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res9 <- fitModel(imp, eq, D_cc, we=FALSE)

# Plot it
vars <- c("info channel", "edu low vs. high", "edu med vs. high", "female", 
          "35-49 vs. 18-34", "50-64 vs. 18-34", "65+ vs. 18-34",
          "subj. health", "subj. inform", "worried Corona", "2022 vs. 2021")
resA <- rbind.data.frame(res1[-1,], res1_ww[-1,])
resA$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res1)-1))
resA$Category <- "Family/Friends"
resA$Variables <- rep(vars, 2)
resB <- rbind.data.frame(res2[-1,], res2_ww[-1,])
resB$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res2)-1))
resB$Category <- "MedDoc"
resB$Variables <- rep(vars, 2)
resC <- rbind.data.frame(res3[-1,], res3_ww[-1,])
resC$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res3)-1))
resC$Category <- "Pharm"
resC$Variables <- rep(vars, 2)
resD <- rbind.data.frame(res4[-1,], res4_ww[-1,])
resD$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res4)-1))
resD$Category <- "PublicMat"
resD$Variables <- rep(vars, 2)
resE <- rbind.data.frame(res5[-1,], res5_ww[-1,])
resE$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res5)-1))
resE$Category <- "NewsP"
resE$Variables <- rep(vars, 2)
resF <- rbind.data.frame(res6[-1,], res6_ww[-1,])
resF$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res6)-1))
resF$Category <- "SocMed"
resF$Variables <- rep(vars, 2)
resG <- rbind.data.frame(res7[-1,], res7_ww[-1,])
resG$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res7)-1))
resG$Category <- "AdmMat"
resG$Variables <- rep(vars, 2)
resH <- rbind.data.frame(res8[-1,], res8_ww[-1,])
resH$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res8)-1))
resH$Category <- "HealthIns"
resH$Variables <- rep(vars, 2)
resI <- rbind.data.frame(res9[-1,], res9_ww[-1,])
resI$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res9)-1))
resI$Category <- "HealthWeb"
resI$Variables <- rep(vars, 2)

res <- rbind.data.frame(resA,resB,resC,resD,resE,resF,resG,resH,resI)
res$Estimate <- as.numeric(res$Estimate)
res$CI_low <- as.numeric(res$CI_low)
res$CI_up <- as.numeric(res$CI_up)
res$Category <- as.factor(res$Category)
res$Variables <- factor(res$Variables, levels=rev(vars))

zp <- ggplot(res,aes(colour=Type)) 
zp <- zp + geom_pointrange(aes(x = Variables, y = exp(Estimate), ymin = exp(CI_low), ymax = exp(CI_up)),
                           lwd = 0.1, position = position_dodge(width = 1/2)) 
zp <- zp + facet_grid(cols=vars(Category)) + scale_color_manual(values=c("darkred", "darkblue"))
zp <- zp + geom_hline(yintercept=1, color = "grey10")
zp <- zp + xlab("") + ylab("Odd-Ratio") 
zp <- zp + ggtitle("Basic Immunisation") 
zp <- zp + coord_flip()
zp

resR <- res[res$Variables %in% "info channel",]
zp <- ggplot(resR,aes(colour=Type)) 
zp <- zp + geom_pointrange(aes(x = Category, y = exp(Estimate), ymin = exp(CI_low), ymax = exp(CI_up)),
                           lwd = 0.7, position = position_dodge(width = 1/2)) 
zp <- zp + geom_hline(yintercept=1, color = "grey10")
zp <- zp + xlab("") + ylab("Odd-Ratio") 
zp <- zp + ggtitle("Basic Immunisation (under controls)") 
zp <- zp + coord_flip()
zp

# ------------------------------------------------------------------------------
# C. Models: basic immunisation & willingness for vaccination
# - logit model, weighted, cluster robust std err, imputed 
# ------------------------------------------------------------------------------
# for each channel, weigthed and unweighted, logit

# 1. Info Fam
eq <- as.formula(ImpEntPos ~ infoFam + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res1_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res1 <- fitModel(imp, eq, D_cc, we=FALSE)

# 2. Info Doc
eq <- as.formula(ImpEntPos ~ infoDoc + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res2_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res2 <- fitModel(imp, eq, D_cc, we=FALSE)

# 3. Info Apo
eq <- as.formula(ImpEntPos ~ infoApo + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res3_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res3 <- fitModel(imp, eq, D_cc, we=FALSE)

# 4. Info infoOef
eq <- as.formula(ImpEntPos ~ infoOef + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res4_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res4 <- fitModel(imp, eq, D_cc, we=FALSE)

# 5. Info infoNews
eq <- as.formula(ImpEntPos ~ infoNews + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res5_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res5 <- fitModel(imp, eq, D_cc, we=FALSE)

# 6. Info infoSozM
eq <- as.formula(ImpEntPos ~ infoSozM + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res6_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res6 <- fitModel(imp, eq, D_cc, we=FALSE)

# 7. Info infoAdm
eq <- as.formula(ImpEntPos ~ infoAdm + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res7_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res7 <- fitModel(imp, eq, D_cc, we=FALSE)

# 8. Info infoKK
eq <- as.formula(ImpEntPos~ infoKK + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res8_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res8 <- fitModel(imp, eq, D_cc, we=FALSE)

# 9. Info infoGP
eq <- as.formula(ImpEntPos ~ infoGP + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res9_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res9 <- fitModel(imp, eq, D_cc, we=FALSE)

# Plot it
vars <- c("info channel", "edu low vs. high", "edu med vs. high", "female", 
          "35-49 vs. 18-34", "50-64 vs. 18-34", "65+ vs. 18-34",
          "subj. health", "subj. inform", "worried Corona", "2022 vs. 2021")
resA <- rbind.data.frame(res1[-1,], res1_ww[-1,])
resA$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res1)-1))
resA$Category <- "Family/Friends"
resA$Variables <- rep(vars, 2)
resB <- rbind.data.frame(res2[-1,], res2_ww[-1,])
resB$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res2)-1))
resB$Category <- "MedDoc"
resB$Variables <- rep(vars, 2)
resC <- rbind.data.frame(res3[-1,], res3_ww[-1,])
resC$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res3)-1))
resC$Category <- "Pharm"
resC$Variables <- rep(vars, 2)
resD <- rbind.data.frame(res4[-1,], res4_ww[-1,])
resD$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res4)-1))
resD$Category <- "PublicMat"
resD$Variables <- rep(vars, 2)
resE <- rbind.data.frame(res5[-1,], res5_ww[-1,])
resE$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res5)-1))
resE$Category <- "NewsP"
resE$Variables <- rep(vars, 2)
resF <- rbind.data.frame(res6[-1,], res6_ww[-1,])
resF$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res6)-1))
resF$Category <- "SocMed"
resF$Variables <- rep(vars, 2)
resG <- rbind.data.frame(res7[-1,], res7_ww[-1,])
resG$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res7)-1))
resG$Category <- "AdmMat"
resG$Variables <- rep(vars, 2)
resH <- rbind.data.frame(res8[-1,], res8_ww[-1,])
resH$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res8)-1))
resH$Category <- "HealthIns"
resH$Variables <- rep(vars, 2)
resI <- rbind.data.frame(res9[-1,], res9_ww[-1,])
resI$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res9)-1))
resI$Category <- "HealthWeb"
resI$Variables <- rep(vars, 2)

res <- rbind.data.frame(resA,resB,resC,resD,resE,resF,resG,resH,resI)
res$Estimate <- as.numeric(res$Estimate)
res$CI_low <- as.numeric(res$CI_low)
res$CI_up <- as.numeric(res$CI_up)
res$Category <- as.factor(res$Category)
res$Variables <- factor(res$Variables, levels=rev(vars))

zp <- ggplot(res,aes(colour=Type)) 
zp <- zp + geom_pointrange(aes(x = Variables, y = exp(Estimate), ymin = exp(CI_low), ymax = exp(CI_up)),
                           lwd = 0.1, position = position_dodge(width = 1/2)) 
zp <- zp + facet_grid(cols=vars(Category)) + scale_color_manual(values=c("darkred", "darkblue"))
zp <- zp + geom_hline(yintercept=1, color = "grey10")
zp <- zp + xlab("") + ylab("Odd-Ratio") 
zp <- zp + ggtitle("Positive Vaccination Decision") 
zp <- zp + coord_flip()
zp

resR <- res[res$Variables %in% "info channel",]
zp <- ggplot(resR,aes(colour=Type)) 
zp <- zp + geom_pointrange(aes(x = Category, y = exp(Estimate), ymin = exp(CI_low), ymax = exp(CI_up)),
                           lwd = 0.7, position = position_dodge(width = 1/2)) 
zp <- zp + geom_hline(yintercept=1, color = "grey10")
zp <- zp + xlab("") + ylab("Odd-Ratio") 
zp <- zp + ggtitle("Positive Vaccination Decision (under controls)") 
zp <- zp + coord_flip()
zp

# ------------------------------------------------------------------------------
# D. Models: vaccinated at least once & willingness for vaccination
# - logit model, weighted, cluster robust std err, imputed 
# ------------------------------------------------------------------------------
# for each channel, weigthed and unweighted, logit

# 1. Info Fam
eq <- as.formula(ImpEntPos_v2 ~ infoFam + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res1_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res1 <- fitModel(imp, eq, D_cc, we=FALSE)

# 2. Info Doc
eq <- as.formula(ImpEntPos_v2 ~ infoDoc + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res2_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res2 <- fitModel(imp, eq, D_cc, we=FALSE)

# 3. Info Apo
eq <- as.formula(ImpEntPos_v2 ~ infoApo + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res3_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res3 <- fitModel(imp, eq, D_cc, we=FALSE)

# 4. Info infoOef
eq <- as.formula(ImpEntPos_v2 ~ infoOef + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res4_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res4 <- fitModel(imp, eq, D_cc, we=FALSE)

# 5. Info infoNews
eq <- as.formula(ImpEntPos_v2 ~ infoNews + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res5_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res5 <- fitModel(imp, eq, D_cc, we=FALSE)

# 6. Info infoSozM
eq <- as.formula(ImpEntPos_v2 ~ infoSozM + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res6_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res6 <- fitModel(imp, eq, D_cc, we=FALSE)

# 7. Info infoAdm
eq <- as.formula(ImpEntPos_v2 ~ infoAdm + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res7_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res7 <- fitModel(imp, eq, D_cc, we=FALSE)

# 8. Info infoKK
eq <- as.formula(ImpEntPos_v2 ~ infoKK + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res8_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res8 <- fitModel(imp, eq, D_cc, we=FALSE)

# 9. Info infoGP
eq <- as.formula(ImpEntPos_v2 ~ infoGP + eduLow + eduMed + sexFem + as.factor(agegrp17C) + 
                   + sjHealth + sjInf + sjWorr + as.factor(datj_n2))
res9_ww <- fitModel(imp, eq, D_cc, we=TRUE)
res9 <- fitModel(imp, eq, D_cc, we=FALSE)

# Plot it
vars <- c("info channel", "edu low vs. high", "edu med vs. high", "female", 
          "35-49 vs. 18-34", "50-64 vs. 18-34", "65+ vs. 18-34",
          "subj. health", "subj. inform", "worried Corona", "2022 vs. 2021")
resA <- rbind.data.frame(res1[-1,], res1_ww[-1,])
resA$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res1)-1))
resA$Category <- "Family/Friends"
resA$Variables <- rep(vars, 2)
resB <- rbind.data.frame(res2[-1,], res2_ww[-1,])
resB$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res2)-1))
resB$Category <- "MedDoc"
resB$Variables <- rep(vars, 2)
resC <- rbind.data.frame(res3[-1,], res3_ww[-1,])
resC$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res3)-1))
resC$Category <- "Pharm"
resC$Variables <- rep(vars, 2)
resD <- rbind.data.frame(res4[-1,], res4_ww[-1,])
resD$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res4)-1))
resD$Category <- "PublicMat"
resD$Variables <- rep(vars, 2)
resE <- rbind.data.frame(res5[-1,], res5_ww[-1,])
resE$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res5)-1))
resE$Category <- "NewsP"
resE$Variables <- rep(vars, 2)
resF <- rbind.data.frame(res6[-1,], res6_ww[-1,])
resF$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res6)-1))
resF$Category <- "SocMed"
resF$Variables <- rep(vars, 2)
resG <- rbind.data.frame(res7[-1,], res7_ww[-1,])
resG$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res7)-1))
resG$Category <- "AdmMat"
resG$Variables <- rep(vars, 2)
resH <- rbind.data.frame(res8[-1,], res8_ww[-1,])
resH$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res8)-1))
resH$Category <- "HealthIns"
resH$Variables <- rep(vars, 2)
resI <- rbind.data.frame(res9[-1,], res9_ww[-1,])
resI$Type <- as.factor(rep(c("unweighted","weighted"), each=nrow(res9)-1))
resI$Category <- "HealthWeb"
resI$Variables <- rep(vars, 2)

res <- rbind.data.frame(resA,resB,resC,resD,resE,resF,resG,resH,resI)
res$Estimate <- as.numeric(res$Estimate)
res$CI_low <- as.numeric(res$CI_low)
res$CI_up <- as.numeric(res$CI_up)
res$Category <- as.factor(res$Category)
res$Variables <- factor(res$Variables, levels=rev(vars))

zp <- ggplot(res,aes(colour=Type)) 
zp <- zp + geom_pointrange(aes(x = Variables, y = exp(Estimate), ymin = exp(CI_low), ymax = exp(CI_up)),
                           lwd = 0.1, position = position_dodge(width = 1/2)) 
zp <- zp + facet_grid(cols=vars(Category)) + scale_color_manual(values=c("darkred", "darkblue"))
zp <- zp + geom_hline(yintercept=1, color = "grey10")
zp <- zp + xlab("") + ylab("Odd-Ratio") 
zp <- zp + ggtitle("Positive Vaccination Decision") 
zp <- zp + coord_flip()
zp

resR <- res[res$Variables %in% "info channel",]
zp <- ggplot(resR,aes(colour=Type)) 
zp <- zp + geom_pointrange(aes(x = Category, y = exp(Estimate), ymin = exp(CI_low), ymax = exp(CI_up)),
                           lwd = 0.7, position = position_dodge(width = 1/2)) 
zp <- zp + geom_hline(yintercept=1, color = "grey10")
zp <- zp + xlab("") + ylab("Odd-Ratio") 
zp <- zp + ggtitle("Positive Vaccination Decision (under controls)") 
zp <- zp + coord_flip()
zp

setwd("C:\\Users\\Freddie\\Documents\\RKI-Studie\\Paper\\Informationen")
write.csv2(res, "regression_vaccPosDec_allChannels_completeModels_weightedUnweighted.txt", row.names = FALSE)
write.csv2(resR, "regression_vaccPosDec_allChannels_weightedUnweighted.txt", row.names = FALSE)



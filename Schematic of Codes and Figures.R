library(ggplot2)
library(grid)
library(BWQS)
library(bkmr)
library(gtsummary)

d <- as.data.frame(read.csv(".....csv")) # Read the data
d_set1 <- d[d$imp==1,]  # choose one of the imputed dataset

# create liver injury outcome

d_set1_comp <- d_set1[!(d_set1$helixid %in% unique(c(d_set1$helixid[is.na(d_set1$alt)], d_set1$helixid[is.na(d_set1$ast)], d_set1$helixid[is.na(d_set1$ggt)]))),] 
d_set1_comp$liver_injury <- rep(NA_character_, dim(d_set1_comp)[1])
d_set1_comp$liver_injury <- ifelse((d_set1_comp$alt > quantile(d_set1_comp$alt, 0.9)) | (d_set1_comp$ast > quantile(d_set1_comp$ast, 0.9)) | (d_set1_comp$ggt > quantile(d_set1_comp$ggt, 0.9)), "1", '0')
d_set1_comp$liver_injury <- factor(d_set1_comp$liver_injury, levels = c("0","1"))

# Convert child age from days to years

d_set1_comp$hs_child_age_days_none <- d_set1_comp$hs_child_age_days_none/365

# Function to convert exposures to quartiles 

convert.to.quantile <- function(x){q = 4; as.integer(cut(x, quantile(x, probs=0:q/q), include.lowest=TRUE))}

# Summary data

summary_data <- as.data.frame(d_set1_comp[,c('h_mbmi_none','h_age_none','h_parity_none', 'h_edumc_none',"h_cohort", 'h_pecl_none',"liver_injury","alt","ast","ggt","hs_child_age_days_none",
                                             "e3_sex_none","hs_c_bmi_none")])
tbl_summary(summary_data, by = liver_injury, statistic = list(all_continuous() ~ "{mean} ({sd}) [{median} ({p25}, {p75})]"), 
            digits = list(all_continuous() ~ c(2, 2), all_categorical() ~ c(0,2)), percent = "column", missing_text = "Missing") %>% 
  add_n() %>% add_p(test = list(all_categorical() ~ "fisher.test",all_continuous() ~ "wilcox.test", all_dichotomous() ~ "lme4"), pvalue_fun = function(x) style_pvalue(x, digits = 3)) %>% bold_labels()  %>% 
  italicize_labels() %>%
  modify_header(label = "**Variable**") %>%  add_overall() %>% add_stat_label() %>%
  bold_p(t = 0.05)

################################################################################

# Code for Figure 1

prenatal_exposures <- d_set1_comp[,c('hs_dde_madj_log2','hs_ddt_madj_log2','hs_hcb_madj_log2','hs_pcb118_madj_log2',	'hs_pcb138_madj_log2','hs_pcb153_madj_log2','hs_pcb170_madj_log2','hs_pcb180_madj_log2','hs_pbde153_madj_log2',	'hs_pbde47_madj_log2',
                                     "hs_bpa_madj_log2", 'hs_oxbe_madj_log2', 'hs_trcs_madj_log2',
                                     'hs_bupa_madj_log2','hs_etpa_madj_log2','hs_mepa_madj_log2','hs_prpa_madj_log2',
                                     'hs_mehp_madj_log2','hs_meohp_madj_log2', 'hs_mehhp_madj_log2',
                                     'hs_mecpp_madj_log2', 'hs_mbzp_madj_log2', 'hs_ohminp_madj_log2','hs_oxominp_madj_log2',
                                     'hs_mep_madj_log2','hs_mibp_madj_log2','hs_mnbp_madj_log2',
                                     'hs_dep_madj_log2','hs_detp_madj_log2','hs_dmp_madj_log2','hs_dmtp_madj_log2',
                                     "hs_pfhxs_m_log2","hs_pfna_m_log2","hs_pfoa_m_log2","hs_pfos_m_log2","hs_pfunda_m_log2",
                                     "hs_as_m_log2","hs_cd_m_log2","hs_co_m_log2","hs_cs_m_log2","hs_cu_m_log2",
                                     "hs_hg_m_log2","hs_mn_m_log2","hs_mo_m_log2","hs_pb_m_log2")]

prenatal_exposures_bar <- data.frame(Chemicals = colnames(prenatal_exposures), Groups = c(rep("OC Pesticides",3), rep("PCBs",5), rep("PBDEs",2),rep("Phenols",3),rep("Parabens",4),rep("HMWPs", 7),rep("LMWPs", 3),rep("OP Pesticides",4), rep("PFAS",5), rep("Metals",9)), Mean_log_conc = apply(prenatal_exposures, 2, function(x) {mean(x)}), SE = apply(prenatal_exposures, 2, function(x) {sd(x)/sqrt(length(x))}))

prenatal_exposures_bar$Groups <- factor(prenatal_exposures_bar$Groups, levels = c("OC Pesticides","PCBs","PBDEs",
                                                                                  "Phenols","Parabens","HMWPs","LMWPs","OP Pesticides","PFAS","Metals"))

prenatal_exposures_bar$Chemicals <- factor(prenatal_exposures_bar$Chemicals, levels = colnames(prenatal_exposures))

xnames <- rep(NA_character_,length(prenatal_exposures_bar$Chemicals))
for(i in 1:length(xnames)){
  xnames[i] <- toupper(strsplit(as.character(prenatal_exposures_bar$Chemicals),"_")[[i]][2])
}
xnames[13] <- "TCS"
xnames[22] <- "MBzP"
xnames[23] <- "OHMiNP"
xnames[24] <- "OXOMiNP"
xnames[25] <- "MiBP"
xnames[26] <- "MnBP"
xnames[32] <- "PFHxS"
xnames[36] <- "PFUnDA"
xnames[37] <- "As"
xnames[38] <- "Cd"
xnames[39] <- "Co"
xnames[40] <- "Cs"
xnames[41] <- "Cu"
xnames[42] <- "Hg"
xnames[43] <- "Mn"
xnames[44] <- "Mo"
xnames[45] <- "Pb"
xnames[27] <- "MnBP"
xnames[26] <- "MiBP"
xnames[25] <- "MEP"


p <- (ggplot(prenatal_exposures_bar, aes(x=Chemicals, y=Mean_log_conc, fill=Groups)) + 
        geom_bar(stat="identity", color = "black", width = 0.7, position=position_dodge()) + 
        geom_errorbar(aes(ymin=Mean_log_conc-2*SE, ymax=Mean_log_conc+2*SE), width=.2, position=position_dodge(.9)) 
      + labs(x ="",y = "Mean (log2 Concentrations)", tag  = "(A)") +  scale_fill_brewer(palette="Paired"))  +coord_flip() + theme_bw()

pm = p + theme(axis.text.x= element_text(face="bold", angle = 90, size = 10),
               axis.text.y = element_text(size = 9,face = "bold"),
               plot.tag = element_text(size = 14,face = "bold"),
               axis.title=element_text(size=12,face="bold"),plot.margin = unit(c(0,0,0,0), "lines")) + scale_x_discrete(breaks=prenatal_exposures_bar$Chemicals,
                                                                                                                        labels=xnames) + guides(fill = guide_legend(title = ""))

SoilSciGuylabs <- c(rep("OC Pesticides",3), rep("PCBs",5), rep("PBDEs",2),rep("Phenols",3),rep("Parabens",4),rep("HMWPs", 7),rep("LMWPs", 3),rep("OP Pesticides",4), rep("PFAS",5), rep("Metals",9))
prechemplot <- ggplot(melt(cor(prenatal_exposures)), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.7, width=0.7) +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  theme_bw() +
  coord_equal() +
  labs(x="",y="",fill="Correlation", tag = "(B)") +
  theme(axis.text.x=element_text(size=9, angle=50, vjust=1, hjust=1, 
                                 margin=margin(-3,0,0,0), face = "bold"),
        axis.ticks = element_blank(),
        plot.tag = element_text(size = 14,face = "bold"),
        axis.text.y=element_text(size=9, margin=margin(0,-2,0,0), face = "bold"),
        panel.grid.major=element_blank(), plot.margin = margin(0,-5,0,-5),
        legend.text = element_text(angle = 90, hjust = 1, vjust = 1)) +
  ggplot2::scale_y_discrete(labels = xnames) + ggplot2::scale_x_discrete(labels = xnames)

P <- ggarrange(pm, prechemplot, ncol = 2, nrow = 1, common.legend = F, legend = "bottom")
p <- P + theme_bw()


################################################################################

# Bayesian Weighted Quantile Sum regression

## Schematic of the main code

bwqs(liver_injury ~ h_mbmi_none + h_age_none  + h_parity_none + h_edumc_none + h_cohort + hs_child_age_days_none +e3_sex_none, 
     mix_name = c(Set of exposures), data = prenatal_exposures, q = 4,  chains = 4, 
     c_int = c(0.025,0.975), family = "binomial") ## For continuous outcome, the family needs to be changed

## Schematic of the Figure

d1 <- as.data.frame(read.csv("......../bwqs_pre_OC_Pesticide.csv"))
d2 <- as.data.frame(read.csv("......../bwqs_PCBs.csv"))
d3 <- as.data.frame(read.csv("......../bwqs_PBDs.csv"))
d4 <- as.data.frame(read.csv("......../bwqs_Phenols.csv"))
d5 <- as.data.frame(read.csv("......../bwqs_Parabens.csv"))
d6 <- as.data.frame(read.csv("......../bwqs_High_weight_Phthalates.csv"))
d7 <- as.data.frame(read.csv("......../bwqs_Low_weight_Phthalates.csv"))
d8 <- as.data.frame(read.csv("......../bwqs_OP_pesticides.csv"))
d9 <- as.data.frame(read.csv("......../bwqs_PFAS.csv"))
d10 <- as.data.frame(read.csv("......../bwqs_metals.csv"))

r <- rbind(d1[d1$X == "beta1",c("mean","X2.5.","X97.5.")],
           d2[d2$X == "beta1",c("mean","X2.5.","X97.5.")],
           d3[d3$X == "beta1",c("mean","X2.5.","X97.5.")],
           d4[d4$X == "beta1",c("mean","X2.5.","X97.5.")],
           d5[d5$X == "beta1",c("mean","X2.5.","X97.5.")],
           d6[d6$X == "beta1",c("mean","X2.5.","X97.5.")],
           d7[d7$X == "beta1",c("mean","X2.5.","X97.5.")],
           d8[d8$X == "beta1",c("mean","X2.5.","X97.5.")],
           d9[d9$X == "beta1",c("mean","X2.5.","X97.5.")],
           d10[d10$X == "beta1",c("mean","X2.5.","X97.5.")]
)

r <- apply(r, 2, function(x){exp(x)})
r <- round(as.data.frame(r),3)
r$Group = c("OC Pesticides","PCBs","PBDEs","Phenols","Parabens","HMWPs","LMWPs","OP Pesticides","PFAS",
            "Metals")
r <- r[,c(4,1,2,3)]
colnames(r) = c("Group","Odds Ratio","95% CI:Lower","95% CI:Upper")
datatable(r, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

label <- c("OC Pesticides","PCBs","PBDEs","Phenols","Parabens","HMWPs","LMWPs","OP Pesticides","PFAS",
           "Metals")
mean  <-  r$`Odds Ratio`
lower <- r$`95% CI:Lower`
upper <- r$`95% CI:Upper`

df <- data.frame(label, mean, lower, upper)
df$label <- factor(df$label, levels=(df$label))
dff <- df
dff$or <- paste0(round(dff$mean,2)," (",round(dff$lower,2),", ",round(dff$upper,2),")")
dff <- dff[,c("label","or")]
dff$num <- rep(1, 10)


tbl <- (ggplot(dff, aes(x=num, y=label))  + geom_text(label=dff$or, size = 5,fontface = "bold" ) + xlab("") + ylab("") + theme_classic2() 
        + guides(x = "none") + theme(axis.ticks.x = element_blank(),
                                     axis.text.x = element_blank(),axis.text.y = element_text(size=15,face = "bold") )) 
fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(size = 1.5) + 
  geom_hline(yintercept=1, lty=2) +  
  scale_y_continuous(trans='log2') + guides(y ="none") +
  coord_flip() +  
  xlab("") + ylab("")  +  theme(plot.title=element_text(size=14,face="bold"),
                                axis.text.x=element_text(size=12,face="bold"),
                                axis.text.y = element_text(size=10,face = "bold"),
                                axis.title=element_text(size=12,face="bold"),
                                plot.tag = element_text(size = 14,face = "bold"))
fp <- fp + scale_x_discrete(expand = expansion(mult = c(0, 0.06),
                                               add = c(0.3, 0)))
fp <- ggpubr::ggarrange(tbl,fp,  widths = c(1,1),
                        ncol = 2,nrow = 1,common.legend = T, legend = "none")

fp <- annotate_figure(fp, top = text_grob("(A) Liver Injury Risk: Forest Plots", 
                                          color = "black", face = "bold", size = 14), fig.lab.pos = c("top.left"))

prenatal_exposures_bar <- data.frame(Chemicals = c(d1$X[10:nrow(d1)],d2$X[14:nrow(d2)],d3$X[10:nrow(d3)],d4$X[10:nrow(d4)],
                                                   d5$X[10:nrow(d5)],d6$X[10:nrow(d6)],d7$X[10:nrow(d7)],d8$X[10:nrow(d8)],
                                                   d9$X[10:nrow(d9)],d10$X[10:nrow(d10)]), 
                                     
                                     Groups = rep(r$Group, c(length(d1$X[10:nrow(d1)]),length(d2$X[14:nrow(d2)]),length(d3$X[10:nrow(d3)]),
                                                             length(d4$X[10:nrow(d4)]),length(d5$X[10:nrow(d5)]),length(d6$X[10:nrow(d6)]),
                                                             length(d7$X[10:nrow(d7)]),length(d8$X[10:nrow(d8)]),length(d9$X[10:nrow(d9)]),
                                                             length(d10$X[10:nrow(d10)]))),
                                     
                                     Mean_log_conc = c(d1[d1$X %in% d1$X[10:nrow(d1)], c("mean")],
                                                       d2[d2$X %in% d2$X[14:nrow(d2)], c("mean")],
                                                       d3[d3$X %in% d3$X[10:nrow(d3)], c("mean")],
                                                       d4[d4$X %in% d4$X[10:nrow(d4)], c("mean")],
                                                       d5[d5$X %in% d5$X[10:nrow(d5)], c("mean")],
                                                       d6[d6$X %in% d6$X[10:nrow(d6)], c("mean")],
                                                       d7[d7$X %in% d7$X[10:nrow(d7)], c("mean")],
                                                       d8[d8$X %in% d8$X[10:nrow(d8)], c("mean")],
                                                       d9[d9$X %in% d9$X[10:nrow(d9)], c("mean")],
                                                       d10[d10$X %in% d10$X[10:nrow(d10)], c("mean")]),
                                     
                                     lower = c(d1[d1$X %in% d1$X[10:nrow(d1)], c("X2.5.")],
                                               d2[d2$X %in% d2$X[14:nrow(d2)], c("X2.5.")],
                                               d3[d3$X %in% d3$X[10:nrow(d3)], c("X2.5.")],
                                               d4[d4$X %in% d4$X[10:nrow(d4)], c("X2.5.")],
                                               d5[d5$X %in% d5$X[10:nrow(d5)], c("X2.5.")],
                                               d6[d6$X %in% d6$X[10:nrow(d6)], c("X2.5.")],
                                               d7[d7$X %in% d7$X[10:nrow(d7)], c("X2.5.")],
                                               d8[d8$X %in% d8$X[10:nrow(d8)], c("X2.5.")],
                                               d9[d9$X %in% d9$X[10:nrow(d9)], c("X2.5.")],
                                               d10[d10$X %in% d10$X[10:nrow(d10)], c("X2.5.")]),
                                     
                                     upper = c(d1[d1$X %in% d1$X[10:nrow(d1)], c("X97.5.")],
                                               d2[d2$X %in% d2$X[14:nrow(d2)], c("X97.5.")],
                                               d3[d3$X %in% d3$X[10:nrow(d3)], c("X97.5.")],
                                               d4[d4$X %in% d4$X[10:nrow(d4)], c("X97.5.")],
                                               d5[d5$X %in% d5$X[10:nrow(d5)], c("X97.5.")],
                                               d6[d6$X %in% d6$X[10:nrow(d6)], c("X97.5.")],
                                               d7[d7$X %in% d7$X[10:nrow(d7)], c("X97.5.")],
                                               d8[d8$X %in% d8$X[10:nrow(d8)], c("X97.5.")],
                                               d9[d9$X %in% d9$X[10:nrow(d9)], c("X97.5.")],
                                               d10[d10$X %in% d10$X[10:nrow(d10)], c("X97.5.")]))


prenatal_exposures_bar$Chemicals <- factor(prenatal_exposures_bar$Chemicals, levels = prenatal_exposures_bar$Chemicals)
prenatal_exposures_bar$Groups <- factor(prenatal_exposures_bar$Groups, levels =  c("OC Pesticides","PCBs","PBDEs","Phenols","Parabens","HMWPs","LMWPs",
                                                                                   "OP Pesticides","PFAS","Metals"))
prenatal_exposures_bar$CI <-paste0( round(prenatal_exposures_bar$Mean_log_conc,3)," (", round(prenatal_exposures_bar$lower,3),", ",round(prenatal_exposures_bar$upper,3),")")

p <- (ggplot(prenatal_exposures_bar, aes(x=Chemicals, y=Mean_log_conc, fill = Groups)) + 
        geom_bar(stat="identity", color = "black",width = 0.7, position=position_dodge()) + 
        geom_errorbar(aes(ymin=Mean_log_conc, ymax=upper), width=.5,
                      position=position_dodge(.9))  +
        labs(title = "(B) Liver Injury Risk: Estimated Weights", 
             x = "", 
             y = "Estimated Weights Relative to each mixture\nassociations (95% CrIs)") 
      +  scale_fill_brewer(palette="Paired") + theme_bw())

pm <- p + theme(plot.title=element_text(size=14,face="bold"),
                axis.text.x=element_text(face="bold", angle = 90, vjust = 0.5, size = 10),
                plot.tag = element_text(size = 14,face = "bold"),
                axis.title=element_text(size=12,face="bold"),
                strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                legend.position="bottom") + scale_x_discrete(breaks=prenatal_exposures_bar$Chemicals,labels=xnames)

pm = (pm  + geom_segment(aes(x = 0.5, y = 1/3, xend = 3.5, yend = 1/3), color = "black",lwd = 1) 
      +  geom_segment(aes(x = 3.5, y = 1/5, xend = 8.5, yend = 1/5), color = "black", lwd = 1)
      +  geom_segment(aes(x = 8.5, y = 1/2, xend = 10.5, yend = 1/2), color = "black", lwd = 1)
      +  geom_segment(aes(x = 10.5, y = 1/3, xend = 13.5, yend = 1/3), color = "black", lwd = 1)
      +  geom_segment(aes(x = 13.5, y = 1/4, xend = 17.5, yend = 1/4), color = "black", lwd = 1)
      +  geom_segment(aes(x = 17.5, y = 1/7, xend = 24.5, yend = 1/7), color = "black", lwd = 1)
      +  geom_segment(aes(x = 24.5, y = 1/3, xend = 27.5, yend = 1/3), color = "black", lwd = 1)
      +  geom_segment(aes(x = 27.5, y = 1/4, xend = 31.5, yend = 1/4), color = "black", lwd = 1)
      +  geom_segment(aes(x = 31.5, y = 1/5, xend = 36.5, yend = 1/5), color = "black", lwd = 1)
      +  geom_segment(aes(x = 36.5, y = 1/9, xend = 45.5, yend = 1/9), color = "black", lwd = 1))

liver_bwqs <- ggpubr::ggarrange(fp,pm,
                                ncol = 2,nrow = 1,common.legend = T, legend = "none")


################################################################################

# Bayesian Kernel Machine Regression

## Plot of predictor-response function for seperate chemical

Z <- prenatal_exposures[,c(Group of chemicals)]
y <- outcome
X <- set of covariates
fitkm_chemicals <- kmbayes(y = y, Z = as.data.frame((Z)), X = as.data.frame((X)), iter = 2000, family = "binomial",verbose = T, varsel = TRUE) ## For continuous outcome, the family needs to be changed

pred.resp.univar <- PredictorResponseUnivar(fit = fitkm_chemicals)

ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") +  facet_wrap(~ variable) + ylab("h(z)")

## Estimated posterior inclusion probabilities

ExtractPIPs(fitkm_chemicals)

### Overall Risk summaries

risks.overall_chemicals <- OverallRiskSummaries(fit = fitkm_chemicals, y = y, Z = Z, X = X, 
                                                    qs = seq(0.05, 0.95, by = 0.05), q.fixed = 0.5, method = "approx")

ggplot(risks.overall_chemicals, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_pointrange()


## Figure

d1 <- as.data.frame(read.csv("....../BKMR/risks.overall_OC_Pesticides.csv"))
d2 <- as.data.frame(read.csv("....../BKMR/risks.overall_PCBs.csv"))
d3 <-as.data.frame(read.csv("....../BKMR/risks.overall_PBDs.csv"))
d4 <- as.data.frame(read.csv("....../BKMR/risks.overall_Phenols.csv"))
d5 <- as.data.frame(read.csv("....../BKMR/risks.overall_Parabens.csv"))
d6 <- as.data.frame(read.csv("....../BKMR/risks.overall_High_weight_Phthalates.csv"))
d7 <- as.data.frame(read.csv("....../BKMR/risks.overall_Low_weight_Phthalates.csv"))
d8 <- as.data.frame(read.csv("....../BKMR/risks.overall_OP_pesticides.csv"))
d9 <- as.data.frame(read.csv("....../BKMR/risks.overall_PFAS.csv"))
d10 <- as.data.frame(read.csv("....../BKMR/risks.overall_Metals.csv"))

r <- data.frame(quantile = rep(c("25th","50th","75th"),10), 
                
                Estimate = c(d1$est[d1$quantile == 0.25],d1$est[d1$quantile == 0.50],d1$est[d1$quantile == 0.75],
                             d2$est[d2$quantile == 0.25],d2$est[d2$quantile == 0.50],d2$est[d2$quantile == 0.75],
                             d3$est[d3$quantile == 0.25],d3$est[d3$quantile == 0.50],d3$est[d3$quantile == 0.75],
                             d4$est[d4$quantile == 0.25],d4$est[d4$quantile == 0.50],d4$est[d4$quantile == 0.75],
                             d5$est[d5$quantile == 0.25],d5$est[d5$quantile == 0.50],d5$est[d5$quantile == 0.75],
                             d6$est[d6$quantile == 0.25],d6$est[d6$quantile == 0.50],d6$est[d6$quantile == 0.75],
                             d7$est[d7$quantile == 0.25],d7$est[d7$quantile == 0.50],d7$est[d7$quantile == 0.75],
                             d8$est[d8$quantile == 0.25],d8$est[d8$quantile == 0.50],d8$est[d8$quantile == 0.75],
                             d9$est[d9$quantile == 0.25],d9$est[d9$quantile == 0.50],d9$est[d9$quantile == 0.75],
                             d10$est[d10$quantile == 0.25],d10$est[d10$quantile == 0.50],d10$est[d10$quantile == 0.75]), 
                
                
                SD = c(d1$sd[d1$quantile == 0.25],d1$sd[d1$quantile == 0.50],d1$sd[d1$quantile == 0.75],
                       d2$sd[d2$quantile == 0.25],d2$sd[d2$quantile == 0.50],d2$sd[d2$quantile == 0.75],
                       d3$sd[d3$quantile == 0.25],d3$sd[d3$quantile == 0.50],d3$sd[d3$quantile == 0.75],
                       d4$sd[d4$quantile == 0.25],d4$sd[d4$quantile == 0.50],d4$sd[d4$quantile == 0.75],
                       d5$sd[d5$quantile == 0.25],d5$sd[d5$quantile == 0.50],d5$sd[d5$quantile == 0.75],
                       d6$sd[d6$quantile == 0.25],d6$sd[d6$quantile == 0.50],d6$sd[d6$quantile == 0.75],
                       d7$sd[d7$quantile == 0.25],d7$sd[d7$quantile == 0.50],d7$sd[d7$quantile == 0.75],
                       d8$sd[d8$quantile == 0.25],d8$sd[d8$quantile == 0.50],d8$sd[d8$quantile == 0.75],
                       d9$sd[d9$quantile == 0.25],d9$sd[d9$quantile == 0.50],d9$sd[d9$quantile == 0.75],
                       d10$sd[d10$quantile == 0.25],d10$sd[d10$quantile == 0.50],d10$sd[d10$quantile == 0.75]),
                
                Group = rep(c("OC Pesticides","PCBs","PBDEs","Phenols","Parabens","HMWPs","LMWPs","OP Pesticides","PFAS","Metals"), each = 3))

r$lower <- r$Estimate - 1.96*r$SD
r$upper <- r$Estimate + 1.96*r$SD
r$Group <- factor(r$Group, levels = c("OC Pesticides","PCBs","PBDEs","Phenols","Parabens","HMWPs","LMWPs","OP Pesticides","PFAS","Metals"))
r$Estimate <- round(r$Estimate,3)

fp <- ggplot(data=r, aes(x=Group, y=Estimate, ymin=lower, ymax=upper, fill = quantile, color = quantile)) +
  labs(title = "(A) Liver Injury Risk: Forest Plot") +
  geom_linerange(size=1.5,position=position_dodge(width = 0.75)) +
  geom_hline(yintercept=0, size = 2, color = "black") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,9.5), lty =1, lwd = 0.7) +
  geom_point(size=4, shape=21, colour="black",position=position_dodge(width = 0.75)) +
  xlab("") + ylab("Beta for Group Association\n (95% Confidence Interval)") +
  theme_bw() + scale_x_discrete(breaks=c("OC Pesticides","PCBs","PBDEs","Phenols","Parabens",
                                         "HMWPs","LMWPs","OP Pesticides","PFAS","Metals"), 
                                labels = c("OC Pesticides","PCBs","PBDEs","Phenols","Parabens",
                                           "HMWPs","LMWPs","OP Pesticides","PFAS","Metals")) + 
  theme(plot.title=element_text(size=14,face="bold"), 
        axis.text.x=element_text(size = 12, face="bold", angle = 90, vjust = 0), 
        axis.text.y = element_text(size = 12, face = "bold"), axis.title=element_text(size=12,face="bold"), 
        legend.position = "bottom",
        plot.tag = element_text(size = 14,face = "bold"))  + coord_flip() 


d1 <-  as.data.frame(read.csv("....../BKMR/pip_OC_Pesticides.csv"))
d2 <-  as.data.frame(read.csv("....../BKMR/pip_PCBs.csv"))
d3 <-  as.data.frame(read.csv("....../BKMR/pip_PBDs.csv"))
d4 <-  as.data.frame(read.csv("....../BKMR/pip_Phenols.csv"))
d5 <-  as.data.frame(read.csv("....../BKMR/pip_Parabens.csv"))
d6 <-  as.data.frame(read.csv("....../BKMR/pip_High_weight_Phthalates.csv"))
d7 <-  as.data.frame(read.csv("....../BKMR/pip_Low_weight_Phthalates.csv"))
d8 <-  as.data.frame(read.csv("....../BKMR/pip_OP_Pesticides.csv"))
d9 <-  as.data.frame(read.csv("....../BKMR/pip_PFAS.csv"))
d10 <-  as.data.frame(read.csv("....../BKMR/pip_Metals.csv"))

prenatal_exposures_bar <- data.frame(Chemicals = c(d1$variable,d2$variable,d3$variable,d4$variable,
                                                   d5$variable,d6$variable,d7$variable,d8$variable,d9$variable,d10$variable), 
                                     Groups = c(rep("OC Pesticides",3), rep("PCBs",5), rep("PBDEs",2),rep("Phenols",3),rep("Parabens",4),
                                                rep("HMWPs", 7),rep("LMWPs", 3),rep("OP Pesticides",4), rep("PFAS",5), rep("Metals",9)),
                                     Mean_log_conc = c(d1$PIP,d2$PIP,d3$PIP,d4$PIP,d5$PIP,
                                                       d6$PIP,d7$PIP,d8$PIP,d9$PIP,d10$PIP))


prenatal_exposures_bar$Chemicals <- factor(prenatal_exposures_bar$Chemicals, levels = prenatal_exposures_bar$Chemicals)
prenatal_exposures_bar$Groups <- factor(prenatal_exposures_bar$Groups, levels = c("OC Pesticides","PCBs","PBDEs",
                                                                                  "Phenols","Parabens","HMWPs","LMWPs","OP Pesticides","PFAS","Metals"))


prenatal_exposures_bar$Mean_log_conc <- c(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "OC Pesticides"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "OC Pesticides"]),
                                          
                                          prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "PCBs"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "PCBs"]),
                                          
                                          prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "PBDEs"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "PBDEs"]),
                                          
                                          prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "Phenols"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "Phenols"]),
                                          
                                          prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "Parabens"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "Parabens"]),
                                          
                                          prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "HMWPs"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "HMWPs"]),
                                          
                                          prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "LMWPs"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "LMWPs"]),
                                          
                                          prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "OP Pesticides"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "OP Pesticides"]),
                                          
                                          prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "PFAS"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "PFAS"]),
                                          
                                          prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "Metals"]/sum(prenatal_exposures_bar$Mean_log_conc[prenatal_exposures_bar$Groups == "Metals"]))


prenatal_exposures_bar$Mean_log_conc <- round(prenatal_exposures_bar$Mean_log_conc,3)


p <- (ggplot(prenatal_exposures_bar, aes(x=Chemicals, y=Mean_log_conc, fill=Groups)) + 
        geom_bar(stat="identity", color = "black", width = 0.7, position=position_dodge()) + 
        labs(title = "(B) Liver Injury Risk:\nScaled Posterior Inclusion Probability", 
             x = "Prenatal Chemical Exposures", y = "") 
      +  scale_fill_brewer(palette="Paired") + theme_bw())

pm = (p  + geom_segment(aes(x = 0.5, y = 1/3, xend = 3.5, yend = 1/3), color = "black",lwd = 1) 
      +  geom_segment(aes(x = 3.5, y = 1/5, xend = 8.5, yend = 1/5), color = "black", lwd = 1)
      +  geom_segment(aes(x = 8.5, y = 1/2, xend = 10.5, yend = 1/2), color = "black", lwd = 1)
      +  geom_segment(aes(x = 10.5, y = 1/3, xend = 13.5, yend = 1/3), color = "black", lwd = 1)
      +  geom_segment(aes(x = 13.5, y = 1/4, xend = 17.5, yend = 1/4), color = "black", lwd = 1)
      +  geom_segment(aes(x = 17.5, y = 1/7, xend = 24.5, yend = 1/7), color = "black", lwd = 1)
      +  geom_segment(aes(x = 24.5, y = 1/3, xend = 27.5, yend = 1/3), color = "black", lwd = 1)
      +  geom_segment(aes(x = 27.5, y = 1/4, xend = 31.5, yend = 1/4), color = "black", lwd = 1)
      +  geom_segment(aes(x = 31.5, y = 1/5, xend = 36.5, yend = 1/5), color = "black", lwd = 1)
      +  geom_segment(aes(x = 36.5, y = 1/9, xend = 45.5, yend = 1/9), color = "black", lwd = 1))


pm <- pm + theme(plot.title=element_text(size=14,face="bold"),
                 axis.text.x=element_text(face="bold", angle = 90, vjust = 0.5, size = 10),
                 axis.title=element_text(size=12,face="bold"),
                 plot.tag = element_text(size = 14,face = "bold"),
                 strip.text.y = element_text(hjust=0,vjust = 0,angle=180,face="bold"),
                 legend.position="bottom") + scale_x_discrete(breaks=prenatal_exposures_bar$Chemicals,labels=xnames)


# Generalized Linear Mixed Effect Regression (Main exposures are turned to quartiles)

glmer(liver_injury ~  LOG TRANSFORMED EXPOSURE + h_mbmi_none + h_age_none + h_parity_none + h_edumc_none + hs_child_age_days_none + h_cohort
      +e3_sex_none +  (1|h_cohort), data = prenatal_exposures, family = "binomial") ## For continuous outcome, the family needs to be changed


## Figure

linear_table <- read.csv(".../linear_table_ck18.csv", header = T)
linear_table$Chemicals <- xnames
linear_table$Chemicals <- factor(linear_table$Chemicals, levels = linear_table$Chemicals)

label <- linear_table$Chemicals
mean  <-  linear_table$Value
lower <- linear_table$lower_ci
upper <- linear_table$upper_ci

df <- data.frame(label, mean, lower, upper)
df$label <- factor(df$label, levels=(df$label))
df$Group <- c(rep("OC Pesticides",3), rep("PCBs",5), rep("PBDEs",2),rep("Phenols",3),rep("Parabens",4),rep("HMWPs", 7),rep("LMWPs", 3),rep("OP Pesticides",4), rep("PFAS",5), rep("Metals",9))
df$Group <- factor(df$Group, levels = c("OC Pesticides","PCBs",
                                        "PBDEs", "Phenols","Parabens","HMWPs","LMWPs","OP Pesticides","PFAS","Metals"))

df$mean <- round(df$mean,3)
df$CI <- paste0("(",round(df$lower,3),", ",round(df$upper,3),")")

fp_ck18 <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper, color = Group)) +
  geom_pointrange(size = 1) + labs(title = "CK-18 (IU/L):\nLinear Mixed Effect Models") +
  geom_hline(yintercept=0, lty=2) +  scale_color_brewer(palette="Paired") + 
  coord_flip() +  
  xlab("Prenatal Chemical Exposures") + ylab("Beta Coefficient\n(95% Confidence Interval)") +
  theme_bw()  +  theme(plot.title=element_text(size=12,face="bold"),
                       axis.text.x=element_text(face="bold", size = 11),
                       axis.text.y = element_text(face = "bold"),
                       axis.title=element_text(size=12,face="bold"))

fp_ck18 <- fp_ck18 + geom_vline(xintercept = c(3.5,8.5,10.5,13.5, 17.5, 24.5, 27.5, 31.5, 36.5),colour = "black",lwd = 1)
fp_ck18 <- fp_ck18 + theme(legend.title = element_blank(),
                           legend.spacing.y = unit(0, "mm"), 
                           panel.border = element_rect(colour = "black", fill=NA),
                           aspect.ratio = 1, 
                           legend.background = element_blank(),
                           legend.box.background = element_rect(colour = "black"))

linear_table <- read.csv(".../linear_table_liver_injury.csv", header = T)
linear_table$chemicals <- xnames
linear_table$chemicals <- factor(linear_table$chemicals, levels = linear_table$chemicals)

label <- linear_table$chemicals
mean  <-  exp(linear_table$Estimate)
lower <- exp(linear_table$lower_ci)
upper <- exp(linear_table$upper_ci)

df <- data.frame(label, mean, lower, upper)
df$label <- factor(df$label, levels=(df$label))
df$Group <- c(rep("OC Pesticides",3), rep("PCBs",5), rep("PBDEs",2),rep("Phenols",3),rep("Parabens",4),rep("HMWPs", 7),rep("LMWPs", 3),rep("OP Pesticides",4), rep("PFAS",5), rep("Metals",9))
df$Group <- factor(df$Group, levels = c("OC Pesticides","PCBs",
                                        "PBDEs", "Phenols","Parabens","HMWPs","LMWPs","OP Pesticides","PFAS","Metals"))


fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper, color = Group)) +
  geom_pointrange(size = 1) + 
  scale_color_brewer(palette="Paired") +
  labs(title = "Liver Injury Risk:\nGeneralized Linear Mixed Effect Models") +
  geom_hline(yintercept=1, lty=2)  + 
  coord_flip() +  
  xlab("Prenatal Chemical Exposures") + ylab("Odds Ratio\n(95% Confidence Interval)") +
  theme_bw()  +  theme(plot.title=element_text(size=12,face="bold"),
                       axis.text.x=element_text(face="bold", size = 11),
                       axis.text.y = element_text(face = "bold"),
                       axis.title=element_text(size=12,face="bold")) 

fp <- fp + geom_vline(xintercept = c(3.5,8.5,10.5,13.5, 17.5, 24.5, 27.5, 31.5, 36.5),colour = "black",lwd = 1)
fp <- fp + theme(legend.title = element_blank(),
                 legend.spacing.y = unit(0, "mm"), 
                 panel.border = element_rect(colour = "black", fill=NA),
                 aspect.ratio = 1, 
                 legend.background = element_blank(),
                 legend.box.background = element_rect(colour = "black"))
fp_mix <- ggpubr::ggarrange(fp,fp_ck18,
                            ncol = 2,nrow = 1,common.legend = T, legend = "bottom")


library(tidycensus)
library(tigris)
library(dplyr)
library(tidyr)
library(spdep)
library(RSpectra)
library(Matrix)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(dbarts)
source('MCMC.R')
acs_20_vars = load_variables(
  year = 2021, 
  "acs5",
  cache = TRUE
)

keep_vars <- c(
  "B01001_001", # Total population
  "B01001_002", # Male population
  "B01001_026", # Female population
  "B17001_002", # Below poverty level
  "B19013_001", # Median household income
  "B19058_001", # Food stamp/SNAP count
  "B25077_001", # Median home value
  "B15003_017", # High school graduates
  "B15003_022", # Bachelor's degree
  "B15003_023", # Master's degree
  "B15003_024", # Professional degree
  "B15003_025", # Doctorate degree
  "B23025_005", # Unemployment count
  "B23025_003", # Labor force total
  "B02001_002", # White alone
  "B02001_003", # Black or African American alone
  "B02001_005", # Asian alone
  "B03001_003", # Hispanic or Latino origin
  "B19001_002", # Households with income <10k
  "B19001_017"  # Households with income >200k
)

acs_dat <- get_acs(
  geography="tract",
  state="CA",
  year=2021,
  survey="acs5",
  variables=keep_vars
)

acs_dat <- acs_dat %>% pivot_wider(names_from="variable", values_from=c("estimate", "moe")) %>%
  mutate(SNAPRate=estimate_B19058_001/estimate_B01001_001 ,
         PovRate=estimate_B17001_002/estimate_B01001_001, 
         Male=estimate_B01001_002/estimate_B01001_001, 
         SNAP=estimate_B19058_001/estimate_B01001_001, 
         HomeV=estimate_B25077_001, 
         HS=estimate_B15003_017/estimate_B01001_001, 
         BS=estimate_B15003_022/estimate_B01001_001,
         UNEMP=estimate_B23025_005/estimate_B01001_001,
         MedInc=estimate_B19013_001, 
         MedIncSE=(moe_B19013_001/1.645)) %>% drop_na()
acs_shape <- tracts(cb=F, state="CA", year=2021)
acs_shape <- acs_shape %>% filter(GEOID %in% acs_dat$GEOID)
acs_dat <- acs_dat %>% filter(GEOID %in% acs_shape$GEOID)

racevars <- c(White = "P2_005N", 
              Black = "P2_006N", 
              Asian = "P2_008N", 
              Hispanic = "P2_002N")

demographics <- get_decennial(
  geography = "tract",
  state="CA",
  variables = racevars,
  summary_var = "P2_001N",
  year = 2020
) 

demographics <- demographics %>% mutate(pct=value/summary_value) %>%
  dplyr::select(-c(value, summary_value))%>%
  pivot_wider(names_from="variable", values_from="pct")

acs_dat <- acs_dat %>% left_join(demographics, by="GEOID")

acs_dat <- acs_dat %>% filter(!is.na(MedIncSE)) %>% filter(GEOID %in% acs_shape$GEOID)
acs_dat <- acs_dat[order(match(acs_dat$GEOID,acs_shape$GEOID)),]
acs_shape <- acs_shape %>% filter(GEOID %in% acs_dat$GEOID)

### EDA ###

p1 <- ggplot(acs_dat, aes(x=PovRate, y=log(MedInc)))+
  geom_point(alpha=0.2)+
  theme_classic()+
  xlab("Poverty Rate")+
  ylab("Log Median Income")+
  geom_smooth(se=F)

p2 <- ggplot(acs_dat, aes(x=SNAPRate, y=log(MedInc)))+
  geom_point(alpha=0.2)+
  theme_classic()+
  xlab("SNAP Rate")+
  ylab("Log Median Income")+
  geom_smooth(se=F)

p3 <- ggplot(acs_dat, aes(x=HomeV, y=log(MedInc)))+
  geom_point(alpha=0.2)+
  theme_classic()+
  xlab("Median Home Value")+
  ylab("Log Median Income")+
  geom_smooth(se=F)

p4 <- ggplot(acs_dat, aes(x=BS, y=log(MedInc)))+
  geom_point(alpha=0.2)+
  theme_classic()+
  xlab("Pct. BS")+
  ylab("Log Median Income")+
  geom_smooth(se=F)

plot_grid(p1, p2, p3, p4, nrow=2)
ggsave("eda.png", dpi=600)


ggplot(acs_dat, aes(x=PovRate, y=Asian, color=log(MedInc)))+
  geom_point(alpha=0.2)+
  scale_color_viridis_c()

### Model Fit ###

set.seed(1)
iter <- 2000
burn <- 500



X <- cbind(acs_dat$PovRate, 
           acs_dat$SNAPRate, 
           acs_dat$Male, 
           acs_dat$HomeV, 
           acs_dat$HS, 
           acs_dat$BS, 
           acs_dat$UNEMP, 
           acs_dat$White, 
           acs_dat$Black, 
           acs_dat$Hispanic, 
           acs_dat$Asian)
colnames(X) <- c("PovRate", "SNAPRate", "Male", "HomeV", "HS", "BS",
                 "UNEMP", "White", "Black", "Hispanic", "Asian")
Yobs <- acs_dat$MedInc
system.time(modFH <- FH_Fit(log(Yobs), X, (acs_dat$MedIncSE)^2/Yobs^2, iter=iter, burn=burn))
system.time(modFHBART <- BARTFH(y=log(Yobs), x=cbind(1,scale(X)), D=(acs_dat$MedIncSE)^2/Yobs^2, n.iter=iter, n.burn=burn))

BART_theta <- t(modFHBART$f_x+modFHBART$u)

predFH <- rowMeans(exp(modFH$Preds))
predBARTFH <- rowMeans(exp(BART_theta))

seFH <- apply(exp(modFH$Preds), 1, sd)
seBARTFH <- apply(exp(BART_theta), 1, sd)



##### Plots

PlotDF <- acs_shape 
PlotDF$FH <- predFH; PlotDF$`BART-FH` <- predBARTFH; PlotDF$Direct <- acs_dat$MedInc
PlotDF <- PlotDF %>% pivot_longer(14:16, names_to="Model", values_to="Median Income")
PlotDF$Model <- factor(PlotDF$Model, levels=c("Direct", "FH", "BART-FH"))

ggplot(PlotDF)+
  geom_sf(linewidth=0, aes(fill=`Median Income`))+
  facet_wrap(~Model, nrow=1)+
  theme_map()+
  scale_fill_viridis_c()
ggsave("CA_ests.jpg", dpi=600)


PlotDF <- acs_shape 
PlotDF$FH <- seFH; PlotDF$`BART-FH` <- seBARTFH; PlotDF$Direct <- acs_dat$MedIncSE
PlotDF <- PlotDF %>% pivot_longer(14:16, names_to="Model", values_to="Standard Error")
PlotDF$Model <- factor(PlotDF$Model, levels=c("Direct", "FH", "BART-FH"))

ggplot(PlotDF)+
  geom_sf(linewidth=0, aes(fill=`Standard Error`))+
  facet_wrap(~Model, nrow=1)+
  theme_map()+
  scale_fill_viridis_c(trans="sqrt")
ggsave("CA_se.jpg", dpi=600)


## scatterplots


PlotDF <- acs_shape 
PlotDF$FH <- predFH; PlotDF$`BART-FH` <- predBARTFH; PlotDF$Direct <- acs_dat$MedInc
PlotDF <- PlotDF %>% pivot_longer(14:15, names_to="Model", values_to="Median Income Estimate")
PlotDF$Model <- factor(PlotDF$Model, levels=c("Direct", "FH", "BART-FH"))


ggplot(PlotDF)+
  geom_point(alpha=0.3, aes(x=`Median Income Estimate`, y=Direct))+
  facet_wrap(~Model)+
  geom_abline(slope=1, color="red")+
  xlab("Model-based Estimate")+
  ylab("Direct Estimate")+
  theme_classic()
ggsave("CA_est_Scatter.jpg", dpi=600)





PlotDF <- acs_shape 
PlotDF$FH <- seFH; PlotDF$`BART-FH` <- seBARTFH; PlotDF$Direct <- acs_dat$MedIncSE
PlotDF <- PlotDF %>% pivot_longer(14:15, names_to="Model", values_to="SE")
PlotDF$Model <- factor(PlotDF$Model, levels=c("Direct", "FH", "BART-FH"))

ggplot(PlotDF)+
  geom_point(alpha=0.3, aes(x=SE, y=Direct))+
  facet_wrap(~Model)+
  geom_abline(slope=1, color="red")+
  xlab("Model-based Standard Error")+
  ylab("Direct Estimate Standard Error")+
  theme_classic()
ggsave("CA_SE_Scatter.jpg", dpi=600)





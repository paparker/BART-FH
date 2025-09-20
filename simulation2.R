library(tidycensus)
library(tigris)
library(dplyr)
library(tidyr)
library(spdep)
library(RSpectra)
library(Matrix)
library(dbarts)
library(spdep)
library(xtable)
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


####################
set.seed(1)
nrep <- 200
iter <- 2000
burn <- 500



lowFH <- highFH <- lowBART <- highBART <- lowNN100 <- highNN100 <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)
lowDir <- highDir <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)
low50FH <- high50FH <- low50BART <- high50BART <- low50NN100 <- high50NN100 <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)
low50Dir <- high50Dir <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)
predFH <- predDir <- predBART <- predNN100 <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)
xbFH <-  xbBART <- xbNN100 <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)

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

for(i in 1:nrep){
  Yobs <- rnorm(nrow(X), mean=acs_dat$MedInc, sd=acs_dat$MedIncSE)
  Yobs[Yobs < 0] <- 100
  ## fit models
  modFH <- FH_Fit(log(Yobs), X, (acs_dat$MedIncSE)^2/Yobs^2, iter=iter, burn=burn)
  modBART <- BARTFH(y=log(Yobs), x=cbind(1,scale(X)), D=(acs_dat$MedIncSE)^2/Yobs^2, n.iter=iter, n.burn=burn)
  modNN100 <- FH_RNN_Fit(Y=log(Yobs), X=cbind(1,scale(X)), S2=(acs_dat$MedIncSE)^2/Yobs^2, nh=100, iter=iter, burn=burn)
  
  lowFH[,i] <- apply(exp(modFH$Preds), 1, quantile, probs=0.025)
  highFH[,i] <- apply(exp(modFH$Preds), 1, quantile, probs=0.975)
  low50FH[,i] <- apply(exp(modFH$Preds), 1, quantile, probs=0.25)
  high50FH[,i] <- apply(exp(modFH$Preds), 1, quantile, probs=0.75)
  predFH[,i] <- rowMeans(exp(modFH$Preds))
  xbFH[,i] <- rowMeans(exp(modFH$XB))
  
  
  BART_theta <- t(modBART$f_x+modBART$u)
  lowBART[,i] <- apply(exp(BART_theta), 1, quantile, probs=0.025)
  highBART[,i] <- apply(exp(BART_theta), 1, quantile, probs=0.975)
  low50BART[,i] <- apply(exp(BART_theta), 1, quantile, probs=0.25)
  high50BART[,i] <- apply(exp(BART_theta), 1, quantile, probs=0.75)
  predBART[,i] <- rowMeans(exp(BART_theta))
  xbBART[,i] <- rowMeans(exp(t(modBART$f_x)))
  
  
  lowNN100[,i] <- apply(exp(modNN100$Preds), 1, quantile, probs=0.025)
  highNN100[,i] <- apply(exp(modNN100$Preds), 1, quantile, probs=0.975)
  low50NN100[,i] <- apply(exp(modNN100$Preds), 1, quantile, probs=0.25)
  high50NN100[,i] <- apply(exp(modNN100$Preds), 1, quantile, probs=0.75)
  predNN100[,i] <- rowMeans(exp(modNN100$Preds))
  xbNN100[,i] <- rowMeans(exp(modNN100$XB))
  
  
  lowDir[,i] <- Yobs - 1.96*acs_dat$MedIncSE
  highDir[,i] <- Yobs + 1.96*acs_dat$MedIncSE
  low50Dir[,i] <- Yobs - 0.67*acs_dat$MedIncSE
  high50Dir[,i] <- Yobs + 0.67*acs_dat$MedIncSE
  predDir[,i] <- Yobs
  
  print(i)
}



######### Results


## Coverage rate

cr50 <- c(mean(acs_dat$MedInc < high50Dir & acs_dat$MedInc > low50Dir),
          mean(acs_dat$MedInc < high50FH & acs_dat$MedInc > low50FH),
          mean(acs_dat$MedInc < high50NN100 & acs_dat$MedInc > low50NN100),
          mean(acs_dat$MedInc < high50BART & acs_dat$MedInc > low50BART))


cr95 <- c(mean(acs_dat$MedInc < highDir & acs_dat$MedInc > lowDir),
          mean(acs_dat$MedInc < highFH & acs_dat$MedInc > lowFH),
          mean(acs_dat$MedInc < highNN100 & acs_dat$MedInc > lowNN100),
          mean(acs_dat$MedInc < highBART & acs_dat$MedInc > lowBART))



## Interval score (Gneiting and Raftery 2007)
IS50 <- c(mean((high50Dir -low50Dir) + 2/.05*(low50Dir - acs_dat$MedInc)*(acs_dat$MedInc<low50Dir) + 
                 2/.05*(acs_dat$MedInc-high50Dir)*(acs_dat$MedInc>high50Dir)),
          mean((high50FH -low50FH) + 2/.05*(low50FH - acs_dat$MedInc)*(acs_dat$MedInc<low50FH) + 
                 2/.05*(acs_dat$MedInc-high50FH)*(acs_dat$MedInc>high50FH)),
          mean((high50NN100 -low50NN100) + 2/.05*(low50NN100 - acs_dat$MedInc)*(acs_dat$MedInc<low50NN100) + 
                 2/.05*(acs_dat$MedInc-high50NN100)*(acs_dat$MedInc>high50NN100)),
          mean((high50BART -low50BART) + 2/.05*(low50BART - acs_dat$MedInc)*(acs_dat$MedInc<low50BART) + 
                 2/.05*(acs_dat$MedInc-high50BART)*(acs_dat$MedInc>high50BART)))


IS95 <- c(mean((highDir -lowDir) + 2/.05*(lowDir - acs_dat$MedInc)*(acs_dat$MedInc<lowDir) + 
                 2/.05*(acs_dat$MedInc-highDir)*(acs_dat$MedInc>highDir)),
          mean((highFH -lowFH) + 2/.05*(lowFH - acs_dat$MedInc)*(acs_dat$MedInc<lowFH) + 
                 2/.05*(acs_dat$MedInc-highFH)*(acs_dat$MedInc>highFH)),
          mean((highNN100 -lowNN100) + 2/.05*(lowNN100 - acs_dat$MedInc)*(acs_dat$MedInc<lowNN100) + 
                 2/.05*(acs_dat$MedInc-highNN100)*(acs_dat$MedInc>highNN100)),
          mean((highBART -lowBART) + 2/.05*(lowBART - acs_dat$MedInc)*(acs_dat$MedInc<lowBART) + 
                 2/.05*(acs_dat$MedInc-highBART)*(acs_dat$MedInc>highBART)))

## MSE
mse <- c((mean((acs_dat$MedInc - predDir)^2)),
         (mean((acs_dat$MedInc - predFH)^2)),
         (mean((acs_dat$MedInc - predNN100)^2)),
         (mean((acs_dat$MedInc - predBART)^2)))

mse <- mse/mse[1]

## Rel. Abs Bias
ab <- c(mean(abs(acs_dat$MedInc - rowMeans(predDir))/acs_dat$MedInc),
        mean(abs(acs_dat$MedInc - rowMeans(predFH))/acs_dat$MedInc),
        mean(abs(acs_dat$MedInc - rowMeans(predNN100))/acs_dat$MedInc),
        mean(abs(acs_dat$MedInc - rowMeans(predBART))/acs_dat$MedInc))




######### Results (Quantiles)


## Coverage rate
q_cr50 <- rbind(quantile(rowMeans(acs_dat$MedInc < high50Dir & acs_dat$MedInc > low50Dir), probs=c(0.25, 0.75)),
                quantile(rowMeans(acs_dat$MedInc < high50FH & acs_dat$MedInc > low50FH), probs=c(0.25, 0.75)),
                quantile(rowMeans(acs_dat$MedInc < high50NN100 & acs_dat$MedInc > low50NN100), probs=c(0.25, 0.75)),
                quantile(rowMeans(acs_dat$MedInc < high50BART & acs_dat$MedInc > low50BART), probs=c(0.25, 0.75)))

q_cr95 <- rbind(quantile(rowMeans(acs_dat$MedInc < highDir & acs_dat$MedInc > lowDir), probs=c(0.25, 0.75)),
                quantile(rowMeans(acs_dat$MedInc < highFH & acs_dat$MedInc > lowFH), probs=c(0.25, 0.75)),
                quantile(rowMeans(acs_dat$MedInc < highNN100 & acs_dat$MedInc > lowNN100), probs=c(0.25, 0.75)),
                quantile(rowMeans(acs_dat$MedInc < highBART & acs_dat$MedInc > lowBART), probs=c(0.25, 0.75)))

## Interval score (Gneiting and Raftery 2007)
q_IS50 <- rbind(quantile(rowMeans((highDir -lowDir) + 2/.05*(lowDir - acs_dat$MedInc)*(acs_dat$MedInc<lowDir) + 
                                    2/.05*(acs_dat$MedInc-highDir)*(acs_dat$MedInc>highDir)), probs=c(0.25, 0.75)),
                quantile(rowMeans((highFH -lowFH) + 2/.05*(lowFH - acs_dat$MedInc)*(acs_dat$MedInc<lowFH) + 
                                    2/.05*(acs_dat$MedInc-highFH)*(acs_dat$MedInc>highFH)), probs=c(0.25, 0.75)),
                quantile(rowMeans((highNN100 -lowNN100) + 2/.05*(lowNN100 - acs_dat$MedInc)*(acs_dat$MedInc<lowNN100) + 
                                    2/.05*(acs_dat$MedInc-highNN100)*(acs_dat$MedInc>highNN100)), probs=c(0.25, 0.75)),
                quantile(rowMeans((highBART -lowBART) + 2/.05*(lowBART - acs_dat$MedInc)*(acs_dat$MedInc<lowBART) + 
                                    2/.05*(acs_dat$MedInc-highBART)*(acs_dat$MedInc>highBART)), probs=c(0.25, 0.75)))


q_IS95 <- rbind(quantile(rowMeans((high50Dir -low50Dir) + 2/.05*(low50Dir - acs_dat$MedInc)*(acs_dat$MedInc<low50Dir) + 
                                    2/.05*(acs_dat$MedInc-high50Dir)*(acs_dat$MedInc>high50Dir)), probs=c(0.25, 0.75)),
                quantile(rowMeans((high50FH -low50FH) + 2/.05*(low50FH - acs_dat$MedInc)*(acs_dat$MedInc<low50FH) + 
                                    2/.05*(acs_dat$MedInc-high50FH)*(acs_dat$MedInc>high50FH)), probs=c(0.25, 0.75)),
                quantile(rowMeans((high50NN100 -low50NN100) + 2/.05*(low50NN100 - acs_dat$MedInc)*(acs_dat$MedInc<low50NN100) + 
                                    2/.05*(acs_dat$MedInc-high50NN100)*(acs_dat$MedInc>high50NN100)), probs=c(0.25, 0.75)),
                quantile(rowMeans((high50BART -low50BART) + 2/.05*(low50BART - acs_dat$MedInc)*(acs_dat$MedInc<low50BART) + 
                                    2/.05*(acs_dat$MedInc-high50BART)*(acs_dat$MedInc>high50BART)), probs=c(0.25, 0.75)))

## MSE
q_mse <- rbind(quantile(rowMeans((acs_dat$MedInc - predDir)^2)/rowMeans((acs_dat$MedInc - predDir)^2), probs=c(0.25, 0.75)),
               quantile(rowMeans((acs_dat$MedInc - predFH)^2)/rowMeans((acs_dat$MedInc - predDir)^2), probs=c(0.25, 0.75)),
               quantile(rowMeans((acs_dat$MedInc - predNN100)^2)/rowMeans((acs_dat$MedInc - predDir)^2), probs=c(0.25, 0.75)),
               quantile(rowMeans((acs_dat$MedInc - predBART)^2)/rowMeans((acs_dat$MedInc - predDir)^2), probs=c(0.25, 0.75)))



## Abs Bias
q_ab <- rbind(quantile((abs(acs_dat$MedInc - rowMeans(predDir))/acs_dat$MedInc), probs=c(0.25, 0.75)),
              quantile((abs(acs_dat$MedInc - rowMeans(predFH))/acs_dat$MedInc), probs=c(0.25, 0.75)),
              quantile((abs(acs_dat$MedInc - rowMeans(predNN100))/acs_dat$MedInc), probs=c(0.25, 0.75)),
              quantile((abs(acs_dat$MedInc - rowMeans(predBART))/acs_dat$MedInc), probs=c(0.25, 0.75)))


## Table
fmt_cell <- function(mean_vec, qmat, digits_mean = 3, digits_q = 3, scale = 1, suffix = "") {
  stopifnot(length(mean_vec) == nrow(qmat))
  out <- character(length(mean_vec))
  for (i in seq_along(mean_vec)) {
    m  <- formatC(scale * mean_vec[i], digits = digits_mean, format = "f")
    q1 <- formatC(scale * qmat[i, 1],  digits = digits_q,   format = "f")
    q3 <- formatC(scale * qmat[i, 2],  digits = digits_q,   format = "f")
    out[i] <- paste0(m, " (", q1, ", ", q3, ")", suffix)
  }
  out
}

estimators <- c("Direct","FH","NN100","BART-FH")

tab <- data.frame(
  Estimator   = estimators,
  MSE         = fmt_cell(mse,     q_mse, digits_mean = 3, digits_q = 3),
  Rel_Abs_Bias= fmt_cell(ab,      q_ab,  digits_mean = 3, digits_q = 3),
  Cov_Rate_50    = fmt_cell(100*cr50,  100*q_cr50, digits_mean = 1, digits_q = 1, suffix = " \\%"),
  Int_Score_50   = fmt_cell(IS50/10000, q_IS50/10000, digits_mean = 3, digits_q = 3),
  Cov_Rate_95    = fmt_cell(100*cr95,  100*q_cr95, digits_mean = 1, digits_q = 1, suffix = " \\%"),
  Int_Score_95   = fmt_cell(IS95/10000, q_IS95/10000, digits_mean = 3, digits_q = 3)
  , check.names = FALSE
)
saveRDS(tab, "sim2tab.rds")


print(
  xtable(tab, caption = "Point metrics with IQR in parentheses by estimator.",
         label = "tab:metrics_iqr"),
  include.rownames = FALSE,
  sanitize.text.function = identity  # keep parentheses and \% as-is
)




### Plot

plotDF <- data.frame(FH=rowMeans((acs_dat$MedInc - predFH)^2),
                     BART=rowMeans((acs_dat$MedInc - predBART)^2),
                     SE_group=case_when(
                       acs_dat$MedIncSE <= quantile(acs_dat$MedIncSE, 0.33) ~ "Low Standard Error",
                       acs_dat$MedIncSE > quantile(acs_dat$MedIncSE, 0.33) & acs_dat$MedIncSE <= quantile(acs_dat$MedIncSE, 0.67) ~ "Medium Standard Error",
                       acs_dat$MedIncSE > quantile(acs_dat$MedIncSE, 0.67) ~ "High Standard Error"
                     ))
plotDF$SE_group <- factor(plotDF$SE_group,
                   levels = c("Low Standard Error", "Medium Standard Error", "High Standard Error"),
                   labels = c("Low Standard Error", "Medium Standard Error", "High Standard Error"))
plotDF <- plotDF %>% pivot_longer(1:2, names_to ="Model", values_to = "MSE")
ggplot(plotDF, aes(x=Model, y=log(MSE), fill=Model))+
  facet_wrap(~SE_group)+
  geom_boxplot()+
  ylab("Log MSE")
ggsave("mse_boxplot.jpg", dpi=600)

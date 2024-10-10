library(stats)
library(ggplot2)
library(glmnet)
library('variancePartition')
library(tidyverse)
library(lme4)
library(mgcv)
library(VGAM)
library(nlme)

setwd("C:/Users/m139105/Documents/digipath/inputs")
df <- read.csv('expression.csv', row.names='gene')
epi_vals <- read.csv('variance_input_df.csv', row.names='sid')

# Assuming 'response' is your outcome variable and 'cont_var1', 'cont_var2' are continuous predictors
customModel <- lme(df ~ epi_vals$combined_inflam_perc + epi_vals$combined_lymph_perc,data=epi_vals,
                   correlation = corExp(form = ~ epi_vals$combined_inflam_perc + epi_vals$combined_lymph_perc, nugget = TRUE))

#model <- vgam(df ~ epi_vals$combined_inflam_perc + epi_vals$combined_lymph_perc,family = vectorff())


epi_vals$combined_total_lymph_perc <- epi_vals$combined_inflam_perc + epi_vals$combined_lymph_perc

form <- ~ epi_manual + lymph_manual
varPart <- fitExtractVarPartModel( df, form, epi_vals )
plotVarPart(varPart)
varPart.df <- as.data.frame(varPart)
vp <- varPart.df%>%arrange(-epi_manual)
#arrange(vpvarPart,epi_manual)
#plotVarPart( vp )
#plotPercentBars( vp[1:10,] )

deep_form <- ~ epi_pan_perc + inflam_perc
deep_varPart <- fitExtractVarPartModel( df, deep_form, epi_vals )
plotVarPart(deep_varPart)
deep_varPart.df <- as.data.frame(deep_varPart)
deep_vp <- deep_varPart.df%>%arrange(-epi_pan_perc)

deep_mon_form <- ~ epi_mon_perc + lymph_perc
deep_mon_varPart <- fitExtractVarPartModel( df, deep_mon_form, epi_vals )
plotVarPart(deep_mon_varPart)
deep_mon_varPart.df <- as.data.frame(deep_mon_varPart)
deep_mon_vp <- deep_mon_varPart.df%>%arrange(-epi_mon_perc)

deep_combined_form <- ~ combined_epi_perc  + combined_total_lymph_perc
deep_combined_varPart <- fitExtractVarPartModel( df, deep_combined_form, epi_vals )
plotVarPart(deep_combined_varPart)
deep_combined_varPart.df <- as.data.frame(deep_combined_varPart)
deep_combined_vp <- deep_combined_varPart.df%>%arrange(-combined_epi_perc)

deconv_combined_form <- ~ scaden_epi  + scaden_lymph
deconv_combined_varPart <- fitExtractVarPartModel( df, deconv_combined_form, epi_vals )
plotVarPart(deconv_combined_varPart)
deconv_combined_varPart.df <- as.data.frame(deconv_combined_varPart)
deconv_combined_vp <- deconv_combined_varPart.df%>%arrange(-scaden_epi)

data_frame_merge_prep1 <- merge(vp, deep_vp,
                          by = 'row.names')
colnames(data_frame_merge_prep1)[1] <- "RowName"
row.names(data_frame_merge_prep1) <- data_frame_merge_prep1$RowName

data_frame_merge_prep2 <- merge(data_frame_merge_prep1, deep_mon_vp,
                          by = 'row.names')
colnames(data_frame_merge_prep2)[1] <- "RowName"
row.names(data_frame_merge_prep2)<- data_frame_merge_prep2$RowName

names(data_frame_merge_prep2)[names(data_frame_merge_prep2) == "Residuals.x"] <- "manual_residuals"
names(data_frame_merge_prep2)[names(data_frame_merge_prep2) == "Residuals.y"] <- "pan_residuals"
names(data_frame_merge_prep2)[names(data_frame_merge_prep2) == "Residuals"] <- "mon_residuals"

data_frame_merge_prep3 <- merge(data_frame_merge_prep2, deep_combined_vp,
                                by = 'row.names')
colnames(data_frame_merge_prep3)[1] <- "RowName"
row.names(data_frame_merge_prep3)<- data_frame_merge_prep3$RowName
names(data_frame_merge_prep3)[names(data_frame_merge_prep3) == "Residuals"] <- "combined_residuals"


data_frame_merge_prep4 <- merge(data_frame_merge_prep3, deconv_combined_vp,
                                by = 'row.names')
colnames(data_frame_merge_prep3)[1] <- "RowName"
row.names(data_frame_merge_prep4)<- data_frame_merge_prep4$RowName
names(data_frame_merge_prep4)[names(data_frame_merge_prep4) == "Residuals"] <- "scaden_residuals"

data_frame_merge <- data_frame_merge_prep4[, -c(2, 3, 4)]

write.csv(data_frame_merge, file="variance_explained_df_scaden_included.csv")

p <- ggplot(data_frame_merge, aes(x=epi_manual, y=epi_pan_perc)) + geom_point() +
  geom_smooth(method=lm,
              color="darkred", fill="blue") + labs(x = "Var explained by Manual Epi proportion",
                                                   y= " Var explained by DL Pan Epi proportion") +
  coord_fixed() + lims(x = c(0,1),y = c(0,1)) + geom_abline(intercept = 0, slope = 1, linetype = "dotted")

ggsave("deep_vs_manual_epi_pan.png",p)

p <- ggplot(data_frame_merge, aes(x=epi_manual, y=epi_mon_perc)) + geom_point() +
  geom_smooth(method=lm,
              color="darkred", fill="blue") + labs(x = "Var explained by Manual Epi proportion",
                                                   y= " Var explained by DL Mon Epi proportion") +
  coord_fixed() + lims(x = c(0,1),y = c(0,1)) + geom_abline(intercept = 0, slope = 1, linetype = "dotted")

ggsave("deep_vs_manual_epi_mon.png",p)

p <- ggplot(data_frame_merge, aes(x=til, y=inflam_perc)) + geom_point() +
  geom_smooth(method=lm,
              color="darkred", fill="blue") + labs(x = "Var explained by Manual TILs",
                                                   y= " Var explained by DL Inflammatory proportion") +
  coord_fixed() + lims(x = c(0,1),y = c(0,1)) + geom_abline(intercept = 0, slope = 1, linetype = "dotted")
ggsave("til_vs_inflam.png",p)

p <- ggplot(data_frame_merge, aes(x=til, y=lymph_perc)) + geom_point() +
  geom_smooth(method=lm,
              color="darkred", fill="blue") + labs(x = "Var explained by Manual TILs",
                                                   y= " Var explained by DL Lymph proportion") +
  coord_fixed() + lims(x = c(0,1),y = c(0,1)) + geom_abline(intercept = 0, slope = 1, linetype = "dotted")
ggsave("til_vs_lymph.png",p)

p <- ggplot(data_frame_merge, aes(x=Residuals.x, y=Residuals.y)) + geom_point() +
  geom_smooth(method=lm,
              color="darkred", fill="blue") + labs(x = "Manual Residual",
                                                   y= "Deep Residual") +
  coord_fixed() + lims(x = c(0,1),y = c(0,1)) + geom_abline(intercept = 0, slope = 1, linetype = "dotted")
ggsave("deep_vs_manual_residual.png",p)

merged_df_manual <- data.frame(row.names = data_frame_merge$Row.names)
merged_df_manual$epi_variance <- c(data_frame_merge$epi_manual)
merged_df_manual$inflam_til_variance <- c(data_frame_merge$til)
merged_df_manual$residual_variance <- c(data_frame_merge$Residuals.x)
merged_df_manual$type <- "manual"

merged_df_deep <- data.frame(row.names = data_frame_merge$Row.names)
merged_df_deep$epi_variance <- c(data_frame_merge$epi_pan_perc)
merged_df_deep$inflam_til_variance <- c(data_frame_merge$inflam_perc)
merged_df_deep$residual_variance <- c(data_frame_merge$Residuals.y)
merged_df_deep$type <- "deep"

merged_df <- rbind(merged_df_manual, merged_df_deep)

p<-ggplot(merged_df, aes(x=residual_variance, fill=type, color=type)) +
  geom_histogram(position="identity", alpha=0.5)
ggsave("deep_vs_manual_residual_hist.png",p)

p<-ggplot(merged_df, aes(x=epi_variance, fill=type, color=type)) +
  geom_histogram(position="identity", alpha=0.5)
ggsave("deep_vs_manual_epi_hist.png",p)

p<-ggplot(merged_df, aes(x=inflam_til_variance, fill=type, color=type)) +
  geom_histogram(position="identity", alpha=0.5)
ggsave("deep_vs_manual_inflam_hist.png",p)

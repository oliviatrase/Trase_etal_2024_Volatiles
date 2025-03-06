## Clear Environment
rm(list=ls())

#### Import Libraries ####
library('lme4')
library('ggplot2')
library('dplyr')
library('tidyr')
library('plyr')
library('DHARMa')
library('emmeans')
library('sjPlot')
library('ggeffects')
library('ggpubr')
library('gridExtra')
library('patchwork')

#### Formatting ####
windowsFonts(Times=windowsFont("Times New Roman"))
options(dplyr.summarise.inform = FALSE)

#### Import Data ####
df <- read.csv("https://raw.githubusercontent.com/oliviatrase/Trase_etal_2025_Volatiles/refs/heads/main/20230628_EPNchoices.csv") 

#### Clean Data ####
df <- subset(df,Variety=="MastersChoice")

## Create Orientation variable combining which arm was in front and which side had WCR added
df$Orientation <- ifelse((df$Arm == 'front') & (df$WCR == 'yes'), 'orientation1',
                         ifelse((df$Arm == 'front') & (df$WCR == 'no'), 'orientation2',
                                ifelse((df$Arm == 'back') & (df$WCR == 'yes'), 'orientation2',
                                       ifelse((df$Arm == 'back') & (df$WCR == 'no'), 'orientation1',NA))))

## Melt data frame such that each observation is a pair of plants rather than a single plant
df_agg <- df %>% 
  dplyr::select(-c(Rep,Arm,Sample)) %>%
  pivot_wider(
    names_from = WCR, 
    values_from = c(EPN,Mass),
    values_fn = sum
  )

#### Calculations ####
## Calculate Chemotaxis Index
df_agg$ChemIndex<-(df_agg$EPN_yes-df_agg$EPN_no)/(df_agg$EPN_yes+df_agg$EPN_no)

## Calculate Percentage of EPN toward WCR (PercYes) and away from WCR (PercNo)
df_agg$PercYes<-(df_agg$EPN_yes)/(df_agg$EPN_yes+df_agg$EPN_no)
df_agg$PercNo<-(df_agg$EPN_no)/(df_agg$EPN_yes+df_agg$EPN_no)

## Calculate plant biomass per pair (total and difference between the two)
df_agg$Mass_total = df_agg$Mass_yes + df_agg$Mass_no
df_agg$Mass_diff = df_agg$Mass_yes - df_agg$Mass_no

## Calculate total EPN emerged per plant pair
df_agg$EPN_total = df_agg$EPN_no + df_agg$EPN_yes

#### Data Cleaning Part II ####
## Remove pairs for which less than 10 EPN emerged
df_agg <- df_agg[(df_agg$EPN_total>10),] # removes 4 pairs

## Remove pairs for which plant biomass data was not taken
df_agg <- df_agg %>% drop_na(Mass_total) # removes 8 pairs

## Check number of observations in each year/cover crop treatment
df_agg %>% 
  dplyr::group_by(SoilYear, SoilType) %>% 
  dplyr::summarise(n = n())

## Average # and standard deviation of EPN emerged
mean(df_agg$EPN_total, na.rm =T) # 160.19
sd(df_agg$EPN_total)/sqrt(length((df_agg$EPN_total))) #10.55
# 160 +/- 11

#### Statistics ####
## Run linear mixed effects model to determine whether cover crop, year, or plant biomass predicts total EPN emergence
## Using date of experiment and orientation on the table as random effects
EPNtotal_lme <- nlme::lme(EPN_total ~ SoilType + SoilYear + Mass_diff, data = df_agg, random = ~1|Orientation/Block)
car::qqPlot(resid(EPNtotal_lme)) ## looks okay, one outlier
summ_EPNtotal_lme <- summary(EPNtotal_lme)
print(paste('DF = ',summ_EPNtotal_lme$tTable[,3][1])) # residual degrees of freedom
car::Anova(EPNtotal_lme, type = "II") # no significance

## Run linear mixed effects model to determine whether cover crop, year, or plant biomass predicts chemotaxis index
## using date of experiment and orientation as random effects
chem_model <- nlme::lme(ChemIndex ~ SoilType  + SoilYear +Mass_diff, random = ~ 1|Orientation/Block, data=df_agg)
car::qqPlot(resid(chem_model)) # looks good
summary(chem_model)
car::Anova(chem_model, type = "II") # no significance

#### Plot Total EPN Emergence ####
ggplot(df_agg, aes( x=SoilType, y=EPN_total, fill = Variety))+
  scale_fill_manual(values = c("#878787")) + 
  geom_boxplot()+
  theme_classic()+
  theme(legend.position = "none",
        text=element_text(family="Times"))+
  facet_wrap(~SoilYear)

#### Plot Chemotaxis Indices ####
ggplot(df_agg,aes(x = ChemIndex, y = SoilType, fill = Variety)) +
  geom_violin(aes(),show.legend = FALSE,size=1)+
  scale_fill_manual(values = c("#878787")) + 
  stat_summary(fun=mean, colour="black", geom="point",
  shape=18, size=3, show.legend=FALSE) +
  xlab("Chemotaxis Index")+
  ylab(NULL)+
  xlim(-1.25,1.25)+
  geom_vline(xintercept = 0, linetype="dashed", color = "black",size=1)+
  theme_classic()+theme(axis.text.y=element_text(size=16,color="black"),
                        axis.text.x=element_text(size=10,color="black"),
                        axis.line.y = element_blank(),
                        axis.title=element_text(size=14),
                        strip.text.x = element_text(size = 16),
                        text=element_text(family="Times"))+
  facet_wrap(~SoilYear)

#### Predicted EPN count ####

## Loop through each cover crop treatment and year
prediction_plot_list <- list()
for (i in 1:4){
  cc <- unique(df$SoilType)[i]
  prediction_plot_list1 <- list()
  for (j in 1:2){
    year <- unique(df$SoilYear)[j]
    df_subset <- df[(df$SoilType==cc)&(df$SoilYear==year),]
    
    ## get unique plant pair
    df_subset$Pair <- paste(df_subset$Block,df_subset$Order)
    
    ## run generalized linear mixed model with poisson distribution to predict EPN count
    EPNcount_model <- glmer(EPN ~ WCR+Mass + (1|Arm) + (1|Block) + (1|Pair), family = poisson, data = df_subset)
    # print(car::qqPlot(resid(EPNcount_model))) # looks good
    
    ## get model estimates of fixed effects, including standard error and pvalues
    EPNcount_model_summary <- summary(EPNcount_model)
    
    WCR_estimate <- round(EPNcount_model_summary$coefficients[2,1],2)
    WCR_SE <- round(EPNcount_model_summary$coefficients[2,2],2)
    # WCR_pvalue <- round(EPNcount_model_summary$coefficients[2,4],4)
    WCR_pvalue1 <- EPNcount_model_summary$coefficients[2,4]
    if (WCR_pvalue1 < 0.0001){
      WCR_pvalue <- "<0.0001"
    } else if (WCR_pvalue1 > 0.05) {
      WCR_pvalue <- round(WCR_pvalue1,2)
    } else {
      WCR_pvalue <- format(round(WCR_pvalue1,4), scientific = FALSE)
    }
    
    Mass_estimate <- round(EPNcount_model_summary$coefficients[3,1],2)
    Mass_SE <- round(EPNcount_model_summary$coefficients[3,2],2)
    # Mass_pvalue <- round(EPNcount_model_summary$coefficients[3,4],4)
    Mass_pvalue1 <- EPNcount_model_summary$coefficients[3,4]
    if (Mass_pvalue1 < 0.0001){
      Mass_pvalue <- "<0.0001"
    } else if (Mass_pvalue1 > 0.05) {
      Mass_pvalue <- round(Mass_pvalue1,2)
    } else {
      Mass_pvalue <- round(Mass_pvalue1,4)
    }
    
    ## Put model terms into formula for visualization
    WCR_formula <- bquote(paste(beta[WCR],"=",.(WCR_estimate),"(SE = ",.(WCR_SE),", ",italic(p) ," = ",.(WCR_pvalue),")"))
    Mass_formula <- bquote(paste(beta[Biomass],"=",.(Mass_estimate),"(SE = ",.(Mass_SE),", ",italic(p) ," = ",.(Mass_pvalue),")"))
    
    ## Get predicted values for visualization
    EPNcount_prediction <- ggpredict(EPNcount_model, terms = c("Mass","WCR"),ci.lvl = 0.95)
    
    ## Plot
    predict_plot<-ggplot(EPNcount_prediction, aes(x, predicted)) + 
      geom_line(aes(linetype=group, color=group),size=1) +
      geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2) + 
      scale_fill_manual(values = c("#006769","#d16002")) + 
      scale_colour_manual(values = c("#006769","#d16002")) +
      xlab("                                                                   Predicted maize biomass (g)")+
      ylab("                                    Predicted EPN count")+
      ylim(0,500)+
      xlim(0.1,0.8)+
      theme_classic()+
      theme(text=element_text(family="Times"),
            legend.position = "none")
    
    ## Add formulas to plot
    predict_plot <- predict_plot + 
      annotate("text", x = c(0.1), y = c(490), label = WCR_formula, family="Times", size = 3,hjust = 0)+
      annotate("text", x = c(0.1), y = c(440), label = Mass_formula, family="Times", size = 3,hjust = 0)
    
    ## Append plots to list
    prediction_plot_list1[[j]] <- predict_plot
  }
  prediction_plot_list[[i]] <- prediction_plot_list1
}

## Flatten the nested list into a single list
flat_plot_list <- unlist(prediction_plot_list, recursive = FALSE)

## Display Plots in a Grid
nrow <- 4
ncol <- 2

for (i in seq_along(flat_plot_list)) {
  row_idx <- ceiling(i / ncol)  # Determine row number
  col_idx <- (i - 1) %% ncol + 1  # Determine column number
  
  ## Remove y-axis text only for the second column (keep it for column 1)
  if (col_idx == 2) {
    flat_plot_list[[i]] <- flat_plot_list[[i]] + theme(axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())
  }
  ## Remove x-axis text for all rows except the last one (row 4)
  if (row_idx != nrow) {
    flat_plot_list[[i]] <- flat_plot_list[[i]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  ## Remove y-axis labels for all but third plot
  if (row_idx != 3) {
    flat_plot_list[[i]] <- flat_plot_list[[i]] +theme(axis.title.y = element_blank())
  }
}

final_plot <- wrap_plots(flat_plot_list, ncol = ncol, nrow = nrow)

## Plot predictions and add more plot labels
y_labels <- c("Fallow","","Radish","" ,"Triticale","" ,"Pea","")
x_labels <- c("2022","2023")
final_plot + 
  plot_annotation(
    tag_levels = list(y_labels),  # Left-side labels
    tag_prefix = "", tag_suffix = "",
    title = paste(x_labels, collapse = "                                  "))&  # Spaced out manually for alignment
    theme(plot.tag.position = "left",
      plot.tag = element_text(size = 14, face = "bold", hjust = 0, family = "Times"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.35, family = "Times"),
      plot.margin = unit(c(0.25, 0.25, 0.5, 0.25),"cm")
    )

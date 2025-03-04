## clear environment
rm(list=ls())

#### Import libraries ####
library(rcartocolor)
library(stringr)
library(reshape)
library(dplyr)
library(ggplot2)
library(ggpattern)
library(MASS)
library(ggfortify)
library(vegan)
library(lmerTest)
library(lme4)
library(emmeans)
library(multcomp)
library(caret)
library(randomForest)
library(RColorBrewer )
library(glmmTMB)
library(DHARMa)
library(extrafont)
library(gridExtra)
library(ecoCopula)
library(sjPlot)
library(xtable)
library(rempsyc)
library(ggsci)
library(pals)
# library(statmod)

#### Formatting ####
windowsFonts(Times=windowsFont("Times New Roman"))
options(dplyr.summarise.inform = FALSE)

#### Functions ####
`%ni%` <- Negate(`%in%`)

#### Import data ####
Master <- read.csv("https://raw.githubusercontent.com/oliviatrase/Trase_etal_2024_Volatiles/refs/heads/main/Master.csv")
Volatiles <- read.csv("https://raw.githubusercontent.com/oliviatrase/Trase_etal_2024_Volatiles/refs/heads/main/20240816_Volatiles2.csv")

#### Getting 3-carene data ####
voc_check_df <-  read.csv("https://raw.githubusercontent.com/oliviatrase/Trase_etal_2024_Volatiles/refs/heads/main/20240530_Unknowns_NateLib.csv")
# get values that are likely 3-carene
voc_check_df <- subset(voc_check_df, (Component.RT>8.88) & (Component.RT < 8.91))
# remove values that are obviously not 3-carene
voc_check_df$Compound.Name <- as.character(voc_check_df$Compound.Name)
voc_check_df <- voc_check_df %>% 
  dplyr::filter(!grepl('Octane', Compound.Name))
voc_check_df <- voc_check_df %>% 
  dplyr::filter(!grepl('Undecane', Compound.Name))
voc_check_df <- voc_check_df %>% 
  dplyr::filter(!grepl('Standard', File.Name))
# Rename compounds to 3-carene (checked spectra on mass hunter)
voc_check_df$Compound_new <- "delta-3-Carene"
# fix file.names for merger
voc_check_df$File.Name <- gsub("\\.D","",voc_check_df$File.Name)
voc_check_df$File.Name <- paste0("X",voc_check_df$File.Name)
voc_check_df$File.Name <- gsub("(.+?_.+?)_.*" ,"\\1", voc_check_df$File.Name)

voc_check_df_final <- voc_check_df %>% 
  group_by(Compound_new, File.Name) %>%
  summarise(max_area=max(Component.Area)) %>%
  reshape2::dcast(Compound_new ~ File.Name, value.var="max_area")

voc_check_df_final$CAS <- "13466-78-9"
list_of_cols <- setdiff(colnames(Volatiles), colnames(voc_check_df_final))
voc_check_df_final[list_of_cols] <- NA
# combine with Volatiles data
Volatiles_check <- rbind(Volatiles, voc_check_df_final)


#### Clean data ####
# Get rid of data not needed
Master <- subset(Master,(SoilYear!=2023)&(SoilType!="ConvFallow"))

# change row names
rownames(Master)<-Master$Sample
Master$Sample <- NULL
rownames(Volatiles_check)<-Volatiles_check$Compound_new

# change NAs to 0s
Volatiles_check[is.na(Volatiles_check)] <- 0

#### Create Matrix ####
Volatile_matrix <- as.data.frame(t(Volatiles_check[4:468]))
Volatile_matrix$Sample.Name <- rownames(Volatile_matrix)
Volatile_matrix$Sample <- word(Volatile_matrix$Sample.Name, 2, sep = "_")
rownames(Volatile_matrix)<- Volatile_matrix$Sample
Volatile_matrix$Sample.Name<-NULL
Volatile_matrix$Sample <- NULL
Volatile_matrix[is.na(Volatile_matrix)] <- 0

#### Normalize using Acetic Acid nonyl ester ####

# 5 microliters of 40 ng/ul acetic acid nonyl ester and octanol
# 200 ng acetic acid nonyl ester in 5 ul
# divide peak area of acetic acid noyl ester by 200, get peak area per ng
# divide all other peaks by that number to get ng

Volatile_matrix$Multiplier <- Volatile_matrix$`Acetic acid, nonyl ester`/200
Volatile_matrix_norm <- Volatile_matrix/Volatile_matrix$Multiplier
hist(Volatile_matrix$`Acetic acid, nonyl ester`) 
  # Approximately normal ranging from 500,000 to 2,000,000 peak area

#### Remove contaminants/standards ####
Volatile_matrix_norm <- as.data.frame(Volatile_matrix_norm)
Volatile_matrix_norm$`Acetic acid, nonyl ester`<-NULL
Volatile_matrix_norm$Multiplier <- NULL
length(colnames(Volatile_matrix_norm)) 
  # 103 compounds

# remove contamination
Volatile_matrix_norm <- Volatile_matrix_norm %>% 
  dplyr::select(-contains(c("silox","oluene","ylene","phthalate","omosalate",
                            "Noise","maybe","uinone","naphth","tyrene","diisobutyrate")))

Volatile_matrix_norm$`3-Carene` <- NULL
length(colnames(Volatile_matrix_norm)) 
# removed 18 compounds

#### Analyze Total Concentration of VOCs ####
Volatile_matrix_norm$TotalVOC <- rowSums(Volatile_matrix_norm)

# Check normality
hist(Volatile_matrix_norm$TotalVOC) # not normally distributed
hist(log(Volatile_matrix_norm$TotalVOC)) # more normal

# merge totals with metadata
Volatiles_total <- subset(Volatile_matrix_norm, select=c("TotalVOC"))
Treatments_df <- subset(Master,select=c("SoilType","WCR","Variety","Block"))
Volatiles_total2 <- merge(Treatments_df, Volatiles_total, by=0)
rownames(Volatiles_total2)<-Volatiles_total2$Row.names
Volatiles_total2$Row.names <-NULL
Volatiles_total2<-na.omit(Volatiles_total2)

## Statistics ##
# Linear Model using log transformation, Block as fixed effect
total_linear_model <- lm((log(TotalVOC))~WCR+SoilType+Variety+Block, data=Volatiles_total2)
AIC(total_linear_model) # 1072
summary(total_linear_model) 
  # WCR and Block significant
hist(resid(total_linear_model)) # looks good
plot(resid(total_linear_model)) # looks good
car::qqPlot(resid(total_linear_model)) # good enough
  # linear model works, will just add Block as a random variable in a mixed model

# Linear Mixed Model using log transformation
total_lme_model <- nlme::lme(log(TotalVOC)~WCR+Variety+SoilType,random=~1|Block, data=Volatiles_total2)
AIC(total_lme_model) # 1006
summary(total_lme_model) # WCR significant
car::qqPlot(resid(total_lme_model)) # looks good
  # WCR significant
hist(resid(total_lme_model)) # looks good
plot(resid(total_lme_model)) # looks good

AOV_total <- anova(total_lme_model) # WCR significant. CC marginally significant (0.11)
AOV_total
AOV_total$`p-value`[2] # WCR pval 1.86e-08
  # WCR super significant

# For loop looking at differences within each variety 
var_list <- list()
VAR_DF_list <- list()
var_df_list <- list()
for (j in 1:2){
  var <- unique(Volatiles_total2$Variety)[j]
  print("nlme::lme for log transformed total VOC by variety")
  print(var)
  var_df <- subset(Volatiles_total2, Variety==var)
  
  ## run statistical model
  model_var <- nlme::lme(log(TotalVOC)~WCR+SoilType,random=~1|Block, data=var_df)
  
  ## check assumptions
  # qqnorm(model_var,abline = c(0, 1))
  # print(summary(model_var))
  
  ## print sumamry
  AOV_var <- anova(model_var)
  print(AOV_var)

  ## run pairwise post-hoc tests
  var_emmeans <- emmeans::emmeans(model_var, list(pairwise ~ c(WCR,SoilType)), adjust = "tukey")
  var_letters<-multcomp::cld(object = var_emmeans$`emmeans of WCR, SoilType`, Letters = letters,decreasing=T)
  var_letters<-as.data.frame(var_letters)
  var_letters <- subset(var_letters,select=c("WCR","SoilType",".group"))

  ## if there are no differences at all, replace tukey letters with blank
  if (length(unique(var_letters$.group))==1){
    var_letters$.group <- ""
  }
  
  ## creat dataframe with post-hoc results
  colnames(var_letters)<- c("WCR","SoilType","var_letters")
  var_letters$Variety <- var
  
  ## append resulting dataframe to a list
  var_df_list[[j]] <-var_letters
 
  ## Within each variety, look for differences between cover crop treatments within 
  ## WCR-infested treatments and within control treatments
  WCR_CC_df_list <-list()
  for (k in 1:2){
    wcr <- unique(var_df$WCR)[k]
    wcr_df <- subset(var_df, WCR==wcr)
    
    ## run model using log transformed total VOC
    model_soiltype <- nlme::lme(log(TotalVOC)~SoilType,random=~1|Block, data=wcr_df)
    summ_model_soiltype <- as.data.frame(summary(model_soiltype)$tTable)
    
    ## print results
    print("nlme::lme for log transformed total VOC by variety and WCR status")
    print(paste("DF:",summ_model_soiltype$DF[1]))
    print(var)
    print(paste(wcr,"WCR"))
    print(anova(model_soiltype))
  
    ## run post-hoc tests
    soiltype_emmeans <- emmeans(model_soiltype, list(pairwise ~ c(SoilType)), adjust = "tukey")
    soiltype_letters<-multcomp::cld(object = soiltype_emmeans$emmeans, Letters = letters,decreasing=T)
    soiltype_letters<-as.data.frame(soiltype_letters)
    soiltype_letters <- subset(soiltype_letters,select=c("SoilType",".group"))
    
    ## if there are no difference at all, replace tukey letters with blank
    if (length(unique(soiltype_letters$.group))==1){
      soiltype_letters$.group <- ""
    }
    
    ## append post-hoc results to a dataframe
    colnames(soiltype_letters)<- c("SoilType","CCletters")
    soiltype_letters$WCR <- wcr
    
    ## append dataframe to a list
    WCR_CC_df_list[[k]] <-soiltype_letters
    
  }
  
  ## combine list of post-hoc results into one larger dataframe
  WCR_CC_df <- do.call(rbind,WCR_CC_df_list)
  WCR_CC_df$Variety <- var
  VAR_DF_list[[j]] <- WCR_CC_df
  
  # Differences between WCR vs no WCR within each cover crop
  soiltype <- list()
  pair_pval <- list()
  for (i in 1:5){
    pair <- unique(var_df$SoilType)[i]
    pair_df <- subset(var_df,SoilType==pair)
    
    ## run model
    model_pair <- nlme::lme(log(TotalVOC)~WCR,random=~1|Block, data=pair_df)
    summ_model_pair <- as.data.frame(summary(model_pair)$tTable)
    pair_aov <- anova(model_pair)
    
    ## print results
    print(var)
    print(pair)
    print(paste("DF:",summ_model_pair$DF[1]))
    print(pair_aov)
    
    ## append results to a list
    pair_pval[[i]]<-pair_aov$`p-value`[2]
    soiltype[[i]]<-pair
    
  }
  
  ## create dataframe from results
  pair_pval_df <- data.frame("SoilType"=unlist(soiltype),"Pvalue" = unlist(pair_pval))
  pair_pval_df$Variety <- var
  
  ## append results to a list
  var_list[[j]]<- pair_pval_df
  
}

## create larger dataframe from list of results
var_df_letters <- do.call(rbind,var_df_list) 
VAR_DF <- do.call(rbind,VAR_DF_list)

# Add significance stars
sig_df <- do.call(rbind, var_list)
sig_df$WCR <- "no"

sig_df$sig <- NA
sig_df <- sig_df %>% mutate(sig = case_when(
  sig_df$Pvalue <0.001 ~ '***',
  sig_df$Pvalue <0.01 ~ '**',
  sig_df$Pvalue <0.05 ~ '*',
  TRUE ~ sig))

# For loop looking at differences between variety and cover crop within each wcr treatment
CC_DF_list <- list()
for (i in 1:2){
  cc <- unique(Volatiles_total2$WCR)[i]
  cc_df <- subset(Volatiles_total2,WCR==cc)
  
  ## run model
  model_cc <- nlme::lme(log(TotalVOC)~SoilType*Variety,random=~1|Block, data=cc_df)
  summ_model_cc <- as.data.frame(summary(model_cc)$tTable)
  AOV_cc <- car::Anova(model_cc,type="II")
  
  ## print results
  print("nlme::lme for log transformed total VOC by WCR treatment")
  print(paste(cc,"WCR"))
  print(paste("DF:",summ_model_cc$DF[1]))
  print(AOV_cc)
  
  ## run post hoc tests
  cc_emmeans <- emmeans(model_cc, list(pairwise ~ c(SoilType,Variety)), adjust = "tukey")
  cc_letters<-multcomp::cld(object = cc_emmeans$emmeans, Letters = letters,decreasing=T)
  cc_letters<-as.data.frame(cc_letters)
  cc_letters <- subset(cc_letters,select=c("SoilType","Variety",".group"))
  
  ## if there are no differences, replace tukey letters with blank
  if (length(unique(cc_letters$.group))==1){
    cc_letters$.group <- ""
  }
  
  ## create dataframe from results and append to list
  colnames(cc_letters)<- c("SoilType","Variety","soilvar_letters")
  cc_letters$WCR <- cc
  CC_DF_list[[i]] <-cc_letters
  
}

## create larger dataframe from list of results
CC_DF<-do.call(rbind,CC_DF_list)


#### Plot Total VOC Concentration ####

# Aggregate by mean and standard error
Volatiles_avg <- Volatiles_total2 %>% dplyr::group_by(Variety,SoilType,WCR) %>% dplyr::summarise(Total_mean = mean(TotalVOC),
                                                       Total_sd = sd(TotalVOC),
                                                       n = n(),
                                                       Total_se = Total_sd/sqrt(n))

# Add in statistical data
Volatiles_avg<-merge(Volatiles_avg, sig_df, by=c("Variety","SoilType","WCR"), all.x=T)
Volatiles_avg<-merge(Volatiles_avg, VAR_DF, by=c("Variety","SoilType","WCR"), all.x=T)
Volatiles_avg<-merge(Volatiles_avg, CC_DF, by=c("Variety","SoilType","WCR"), all.x=T)
Volatiles_avg<-merge(Volatiles_avg,var_df_letters,by=c("Variety","SoilType","WCR"), all.x=T)

# Variety vs Variety
TotalVOC_plot <- ggplot(Volatiles_avg,aes(x=SoilType, y=Total_mean,fill=WCR)) +
  geom_bar(stat="identity",position='dodge',width=0.9,color="black") +
  geom_errorbar(aes(ymin=(Total_mean-Total_se), ymax=(Total_mean+Total_se)), 
                width=0.4,
                position = position_dodge(0.9), 
                colour="black",
                size=1) +
  scale_fill_manual(values=c('white','darkgray'))+ 
  scale_color_manual(values=c('white','darkgray'))+ 
  scale_x_discrete(labels=c("CCFallow" = "Fallow", "CCPea" = "Pea","CCRadish"="Radish","CCTriticale"="Triticale","CCMix"="Mixture"))+ 
  ggtitle("") + 
  labs(x=" ",y="ng VOC (mean ± se)",family="Times")+
  ylim(0,1.2*(max(Volatiles_avg$Total_mean)+max(Volatiles_avg$Total_se)))+
  theme_classic()+ 
  theme(axis.text.x= element_text(size= 16,angle=45,vjust=1,hjust=1),
        axis.text.y= element_text(size= 12),axis.title.y= element_text(size= 16),
        strip.text.x = element_text(size = 20),text=element_text(family="Times"))

## add significant statistical results to plot
TotalVOC_plot + 
  geom_text(aes(x = SoilType, label = sig, y = (Total_mean+Total_se)),
            position=position_dodge(.9), size=8,vjust=-0.25,hjust=-0.005,family="Times")+
  facet_wrap(~Variety)


#### Individual Volatiles - data clean up ####
## remove outliers

### Remove compounds that appear in less than three samples in one block ###
## get dataframe where rows are each sample and columns are volatiles
VOC_check <- merge(Treatments_df,Volatile_matrix_norm[,1:(length(colnames(Volatile_matrix_norm))-1)],by=0)
rownames(VOC_check)<- VOC_check$Row.names
VOC_check$Row.names<-NULL
VOC_check$TotalVOC<-NULL

## create list of compounds
VOCs <- colnames(VOC_check[,5:length(colnames(VOC_check))])

## loop through each volatile compound and get the number of occurrences of that compound within each block
compounds_to_remove <- list()
for (i in 1:length(VOCs)){#
  compound <- VOCs[i]
  df <- subset(VOC_check,select = c("Block","Variety","SoilType","WCR",compound))
  df <- data.frame(df)
  colnames(df) <- c("Block","Variety","SoilType","WCR","compound")
  
  df_agg <- df %>% dplyr::group_by(Block) %>% dplyr::summarise(count= sum(compound > 0))
  counts <- df_agg$count
  
  # remove any compounds for which there are less than three occurrences in any block
  if (any(counts<3)){
    compounds_to_remove[[i]]<-compound
    next
  }
}

length(VOCs) # 85 compounds before
length(unlist(compounds_to_remove)) # 23 compounds that don't appear more than twice in one block

## remove unevenly distributed/rare compounds
Volatiles_norm_subset1 <- subset(VOC_check,select = names(VOC_check) %ni% unlist(compounds_to_remove))

### remove volatiles that occur in less than 3 samples ###
cols_to_remove <- names(Volatiles_norm_subset1)[sapply(Volatiles_norm_subset1, function(col) sum(col > 0) < 3)]
Volatiles_norm_subset2 <- Volatiles_norm_subset1[, !(names(Volatiles_norm_subset1) %in% cols_to_remove)]

### remove samples that contain fewer than 4 volatiles ###
rows_to_remove <- which(rowSums(Volatiles_norm_subset2 > 0) < 3)
Volatile_matrix_norm2 <- Volatiles_norm_subset2[ !(names(Volatiles_norm_subset2) %in% rows_to_remove),]



## optional write csv of new volatile matrix
# write.csv(Volatile_matrix_norm2,"Volatiles_normalized2.csv")

#### Volatile Composition ####

Volatile_matrix_norm2$Block <- as.factor(Volatile_matrix_norm2$Block)
numVOC <- length(colnames(Volatile_matrix_norm2))

## account for blocking factor using setBlocks, using 9999 permutations
perm<-how(nperm=9999)
setBlocks(perm)<-with(Volatile_matrix_norm2[5:numVOC], Volatile_matrix_norm2$Block)

## PERMANOVA
maov_total <- adonis2(as.matrix(Volatile_matrix_norm2[5:numVOC])~WCR+SoilType+Variety, data=Volatile_matrix_norm2, permutations = perm) 
maov_total

## calculate dissimilarity indices for all samples
dist_total <- vegdist(Volatile_matrix_norm2[5:numVOC])

## calculate dispersion between WCR and control treatments
betadisper_total_wcr <- betadisper(dist_total,Volatile_matrix_norm2$WCR)
pmod_wcr <- permutest(betadisper_total_wcr, permutations = perm, pairwise = F)
pmod_wcr # F(1,385) = 0.63, p = 0.39

## calculate dispersion between cover crop treatments
betadisper_total_cc <- betadisper(dist_total,Volatile_matrix_norm2$SoilType)
pmod_cc <- permutest(betadisper_total_cc, permutations = perm, pairwise = F)
pmod_cc # F(4,382) = 1.61, p = 0.07

## calculate dispersion between varieties
betadisper_total_var <- betadisper(dist_total,Volatile_matrix_norm2$Variety)
pmod_var <- permutest(betadisper_total_var, permutations = perm, pairwise = F)
pmod_var # F(1,385) = 0.8, p = 0.31

## Calculate PERMANOVA and dispersion By Cover Crop ##
## FALLOW
fallow_all <- subset(Volatile_matrix_norm2,SoilType=='CCFallow')
perm<-how(nperm=9999)
setBlocks(perm)<-with(fallow_all, Block)
fallow_matrix <- fallow_all[5:length(colnames(fallow_all))]
fallow_matrix<-as.matrix(fallow_matrix)
maov_fallow <- adonis2(fallow_matrix~WCR+Variety, data=fallow_all, permutations = perm) 
maov_fallow # WCR: pseudo-F(1,75) = 2.18, R2 = 2.8%, p = 0.004; 
            # Variety: pseudo-F(1,75) = 0.21, R2 = 0.27%, p = 0.99

dist_fallow <- vegdist(fallow_matrix)

betadisper_fallow_wcr <- betadisper(dist_fallow,fallow_all$WCR)
pmod_f_wcr <- permutest(betadisper_fallow_wcr, permutations = perm, pairwise = F)
pmod_f_wcr # not significant

betadisper_fallow_var <- betadisper(dist_fallow,fallow_all$Variety)
pmod_f_var <- permutest(betadisper_fallow_var, permutations = perm, pairwise = F)
pmod_f_var # not significant

## PEA
pea_all <- subset(Volatile_matrix_norm2,SoilType=='CCPea')
perm<-how(nperm=9999)
setBlocks(perm)<-with(pea_all, Block)
pea_matrix <- pea_all[5:length(colnames(pea_all))]
pea_matrix<-as.matrix(pea_matrix)
maov_pea <- adonis2(pea_matrix~WCR+Variety, data=pea_all, permutations = perm) 
maov_pea  # WCR: pseudo-F(1,75) = 2.23, R2 = 2.9%, p = 0.01; 
          # Variety: pseudo-F(1,75) = 1.11, R2 = 1.4%, p = 0.15

dist_pea <- vegdist(pea_matrix)

betadisper_pea_wcr <- betadisper(dist_pea,pea_all$WCR)
pmod_p_wcr <- permutest(betadisper_pea_wcr, permutations = perm, pairwise = F)
pmod_p_wcr # not significant

betadisper_pea_var <- betadisper(dist_pea,pea_all$Variety)
pmod_p_var <- permutest(betadisper_pea_var, permutations = perm, pairwise = F)
pmod_p_var # Variety dispersion significant: pval = 0.039

## Triticale
triticale_all <- subset(Volatile_matrix_norm2,SoilType=='CCTriticale')
perm<-how(nperm=9999)
setBlocks(perm)<-with(triticale_all, Block)
triticale_matrix <- triticale_all[5:length(colnames(triticale_all))]
triticale_matrix<-as.matrix(triticale_matrix)
maov_triticale <- adonis2(triticale_matrix~WCR+Variety, data=triticale_all, permutations = perm) 
maov_triticale  # WCR: pseudo-F(1,74) = 3.05, R2 = 3.9%, p = 0.0004; 
                # Variety: pseudo-F(1,74) = 0.72, R2 = 0.9%, p = 0.33

dist_triticale <- vegdist(triticale_matrix)

betadisper_triticale_wcr <- betadisper(dist_triticale,triticale_all$WCR)
pmod_t_wcr <- permutest(betadisper_triticale_wcr, permutations = perm, pairwise = F)
pmod_t_wcr # not significant

betadisper_triticale_var <- betadisper(dist_triticale,triticale_all$Variety)
pmod_t_var <- permutest(betadisper_triticale_var, permutations = perm, pairwise = F)
pmod_t_var # not significant

## Radish
radish_all <- subset(Volatile_matrix_norm2,SoilType=='CCRadish')
perm<-how(nperm=9999)
setBlocks(perm)<-with(radish_all, Block)
radish_matrix <- radish_all[5:length(colnames(radish_all))]
radish_matrix<-as.matrix(radish_matrix)
maov_radish <- adonis2(radish_matrix~WCR+Variety, data=radish_all, permutations = perm) 
maov_radish   # WCR: pseudo-F(1,75) = 2.14, R2 = 2.8%, p = 0.002; 
              # Variety: pseudo-F(1,75) = 0.7, R2 = 0.9%, p = 0.33

dist_radish <- vegdist(radish_matrix)

betadisper_radish_wcr <- betadisper(dist_radish,radish_all$WCR)
pmod_r_wcr <- permutest(betadisper_radish_wcr, permutations = perm, pairwise = F)
pmod_r_wcr # not significant

betadisper_radish_var <- betadisper(dist_radish,radish_all$Variety)
pmod_r_var <- permutest(betadisper_radish_var, permutations = perm, pairwise = F)
pmod_r_var # not significant

## Mixture
mix_all <- subset(Volatile_matrix_norm2,SoilType=='CCMix')
perm<-how(nperm=9999)
setBlocks(perm)<-with(mix_all, Block)
mix_matrix <- mix_all[5:length(colnames(mix_all))]
mix_matrix<-as.matrix(mix_matrix)
maov_mix <- adonis2(mix_matrix~WCR+Variety, data=mix_all, permutations = perm) 
maov_mix  # WCR: pseudo-F(1,73) = 2.32, R2 = 3.1%, p = 0.006; 
          # Variety: pseudo-F(1,73) = 0.48, R2 = 0.6%, p = 0.66

dist_mix <- vegdist(mix_matrix)

betadisper_mix_wcr <- betadisper(dist_mix,mix_all$WCR)
pmod_m_wcr <- permutest(betadisper_mix_wcr, permutations = perm, pairwise = F)
pmod_m_wcr # not significant

betadisper_mix_var <- betadisper(dist_mix,mix_all$Variety)
pmod_m_var <- permutest(betadisper_mix_var, permutations = perm, pairwise = F)
pmod_m_var # not significant


## by Variety ##
# MASTERS CHOICE
MC_all <- subset(Volatile_matrix_norm2,Variety=='MastersChoice')
perm<-how(nperm=9999)
setBlocks(perm)<-with(MC_all, Block)
MC_matrix <- MC_all[5:length(colnames(MC_all))]
MC_matrix<-as.matrix(MC_matrix)
maov_MC <- adonis2(MC_matrix~WCR+SoilType, data=MC_all, permutations = perm) 
maov_MC   # WCR: pseudo-F(1,188) = 4.63, R2 = 2.4%, p = 0.0001; 
          # SoilType: pseudo-F(4,188) = 0.89, R2 = 1.8%, p = 0.11

dist_MC <- vegdist(MC_matrix)

betadisper_MC_wcr <- betadisper(dist_MC,MC_all$WCR)
pmod_MC_wcr <- permutest(betadisper_MC_wcr, permutations = perm, pairwise = F)
pmod_MC_wcr # not significant

betadisper_MC_var <- betadisper(dist_MC,MC_all$SoilType)
pmod_MC_var <- permutest(betadisper_MC_var, permutations = perm, pairwise = F)
pmod_MC_var # CC pval = 0.011

# BLUE RIVER ORGANIC
BRO_all <- subset(Volatile_matrix_norm2,Variety=='BlueRiverOrganic')
perm<-how(nperm=9999)
setBlocks(perm)<-with(BRO_all, Block)
BRO_matrix <- BRO_all[5:length(colnames(BRO_all))]
BRO_matrix<-as.matrix(BRO_matrix)
maov_BRO <- adonis2(BRO_matrix~WCR+SoilType, data=BRO_all, permutations = perm) 
maov_BRO  # WCR: pseudo-F(1,187) = 5.67, R2 = 2.9%, p = 0.0001; 
          # SoilType: pseudo-F(4,187) = 0.58, R2 = 1.2%, p = 0.7

dist_BRO <- vegdist(BRO_matrix)

betadisper_BRO_wcr <- betadisper(dist_BRO,BRO_all$WCR)
pmod_BRO_wcr <- permutest(betadisper_BRO_wcr, permutations = perm, pairwise = F)
pmod_BRO_wcr # not significant

betadisper_BRO_var <- betadisper(dist_BRO,BRO_all$SoilType)
pmod_BRO_var <- permutest(betadisper_BRO_var, permutations = perm, pairwise = F)
pmod_BRO_var # not significant

#### Chemical Class Analysis ####

## Import Class Data and clean #
Classes <- read.csv("https://raw.githubusercontent.com/oliviatrase/Trase_etal_2024_Volatiles/refs/heads/main/20240520_Classes.csv")

## create data frame where each compound is a row and each sample is a column
## merge with compound class data
Volatiles_Class2 <- as.data.frame(t(Volatile_matrix_norm2[,5:(length(colnames(Volatile_matrix_norm2)))]))
Volatiles_Class2$Compound <- rownames(Volatiles_Class2)
# intersect(Classes$Compound, Volatiles_Class2$Compound)
Volatiles_Class <- merge(Classes,Volatiles_Class2,by="Compound",all.y=T)
Volatiles_Class$Compound <- NULL
Volatiles_Class$Compound_new <- NULL

## Group compounds by chemical class and get the average 
Volatiles_Classes <- Volatiles_Class %>% dplyr::group_by(Class) %>% dplyr::summarise(across(!matches("Class"),\(x) mean(x, na.rm = TRUE)))
Volatiles_Classes<- as.data.frame(Volatiles_Classes)
## remove NA values
Volatiles_Classes<-Volatiles_Classes[!is.na(Volatiles_Classes$Class),]
## clean up columns
rownames(Volatiles_Classes) <- Volatiles_Classes$Class
Volatiles_Classes$Class <- NULL
## transpose so chemical classes are columns and samples are rows
Volatiles_Classes <- as.data.frame(t(Volatiles_Classes))

## Sort classes based on concentration
class_order_df <- data.frame("Class" = colnames(Volatiles_Classes), "sums" = colMeans(Volatiles_Classes))
class_order_df[order(class_order_df$sums),]$Class
Order_of_classes <- c("Unknown","MCFA","Isothiocyanate","Phenol","Alkene","Ketone","Benzothiazole","Aldehyde",
                      "Salicylate","Sesquiterpene/Sesquiterpenoid","Quinone","Furan",
                      "Alkane","Ester","Alcohol","Monoterpene/Monoterpenoid")

## create new dataframe containing chemical class info and metadata
Volatiles_norm_IDs <- Volatile_matrix_norm2[,1:4]
Volatiles_Classes_merged <- merge(Volatiles_norm_IDs,Volatiles_Classes,by=0)
rownames(Volatiles_Classes_merged) <- Volatiles_Classes_merged$Row.names
Volatiles_Classes_merged$Row.names <- NULL
Volatiles_Classes_merged$Block <- as.factor(Volatiles_Classes_merged$Block)

#### Classes: Stacked bar plots ####

## Get average concentration for each class 
Volatile_Classes_total <-Volatiles_Classes_merged %>% dplyr::group_by(Variety,SoilType,WCR) %>% dplyr::summarise(across(is.numeric, ~ as.numeric(mean(.))))
Volatile_Classes_total <- as.data.frame(Volatile_Classes_total)

## reclassify compounds for visualization
Volatile_Classes_melt <- reshape::melt(Volatile_Classes_total, id.var = c("Variety","SoilType","WCR"),variable.name = "Class")
Volatile_Classes_melt$Class<-NA
Volatile_Classes_melt <- Volatile_Classes_melt %>% mutate(Class = case_when(
  Volatile_Classes_melt$variable == "MCFA" ~ "Other",
  Volatile_Classes_melt$variable == "Isothiocyanate" ~ "Other",
  Volatile_Classes_melt$variable == "Phenol" ~ "Other",
  Volatile_Classes_melt$variable == "Sesquiterpene" ~ "Sesquiterpene/Sesquiterpenoid",
  Volatile_Classes_melt$variable == "Sesquiterpenoid" ~ "Sesquiterpene/Sesquiterpenoid",
  Volatile_Classes_melt$variable == "Monoterpene" ~ "Monoterpene/Monoterpenoid",
  Volatile_Classes_melt$variable == "Monoterpenoid" ~ "Monoterpene/Monoterpenoid",
  TRUE ~ variable))
Volatile_Classes_melt$Class <- factor(Volatile_Classes_melt$Class,
                                      levels = c("Unknown","Other","Alkene","Ketone","Benzothiazole","Aldehyde",
                                                 "Salicylate","Sesquiterpene/Sesquiterpenoid","Quinone","Furan",
                                                 "Alkane","Ester","Alcohol","Monoterpene/Monoterpenoid"))

## create Trt variable (combination of cover crop and WCR treatment) for visualization
Volatile_Classes_melt$Trt <- paste(Volatile_Classes_melt$SoilType,Volatile_Classes_melt$WCR)

## get appropriate colors for visualization
stacked_bar_colors <- c("#222222","#F2F3F4","#1E88E5","#D81B60", # black, gray, sapphire blue, deep pink
                        "#33A02C","#E69F00","#6A3D9A","#F0E442", # emerald, golden orange, purple, yellow
                        "#56B4E9","#B22222","#008080","#CC79A7", # cyan, red, teal, pink
                        "#0067A5") # blue

## Stacked bar plot of average chemical class concentrations
ggplot(Volatile_Classes_melt, aes(fill=Class,y=value,x=Trt, pattern=WCR))+
  # scale_alpha_discrete(range = c(0.6, 1))+
  # geom_bar(position="stack",stat="identity")+
  geom_bar_pattern(position = "fill",stat="identity",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  # scale_fill_carto_d(name = "Class",type = "diverging", palette = "Vivid", direction = -1) +
  # scale_fill_manual(values = category20_d3("default")(16)[1:16]) +
  scale_fill_manual(values=stacked_bar_colors)+
  scale_pattern_manual(values = c(no = "none", yes = "stripe")) +
  facet_wrap(~Variety)+
  theme_classic()+
  theme(axis.text.x= element_text(size= 12,angle=45,vjust=1,hjust=1),
        axis.text.y= element_text(size= 12),
        axis.title.y= element_text(size= 12),
        strip.text.x = element_text(size = 20),
        text = element_text(family = "Times"))
 
#### Class Bar Plots ####
## working on changing the stats in the manuscript
Classes_for_model <- colnames(Volatiles_Classes_merged[5:20])

## loop through each chemical class, run statistical tests, and plot
Class_plots <- list()
Class_plot_names <-list()
WCR_pval_list2 <- list()
class_names_list <-list()
for (i in 1:length(Classes_for_model)){
  ## get the comppound name
  chem_class <- Classes_for_model[i]
  class_names_list[[i]] <- chem_class
  
  df <- subset(Volatiles_Classes_merged, select = c("Block","Variety","SoilType","WCR",chem_class))
  df <- data.frame(df)
  colnames(df) <- c("Block","Variety","SoilType","WCR","chem_class")
  df$Block <- as.factor(df$Block)

  ## run a negative binomial glm
  model1 <- glmmTMB((chem_class)~Variety+SoilType+WCR+Block,
                    # ziformula = ~1,
                    data=df,
                    family=nbinom2(),
                    control = glmmTMBControl(rank_check = "adjust"))
  # print(DHARMa::simulateResiduals(model1,plot=T))
  
  AOV_model1 <-(car::Anova(model1, type = "II"))
  model1_df <- as.data.frame(summary(model1)$AICtab)$`summary(model1)$AICtab`[5]
  # print(chem_class)
  # print(model1_df)
  # print(AOV_model1)

  ## assign p-values to values
  Var_pval <- AOV_model1$`Pr(>Chisq)`[1]
  CC_pval <- AOV_model1$`Pr(>Chisq)`[2]
  WCR_pval <- AOV_model1$`Pr(>Chisq)`[3]
  
  ## put WCR pvalues in a list
  WCR_pval_list2[[i]] <- WCR_pval
  
  ## get summary data for each chem class
  df_agg <- df %>% dplyr::group_by(Variety,SoilType,WCR) %>% dplyr::summarise(mean = mean(chem_class),
                                                                              sd = sd(chem_class),
                                                                              n = n(),
                                                                              count = sum(chem_class >0),
                                                                              se = sd/sqrt(n))
  ## loop through WCR treatment and run stats 
  for (n in 1:2){
    wcr_2 <- unique(df$WCR)[n]
    df_wcr <- subset(df, (WCR==wcr_2))
    
    ## run negative binomia model
    model_wcr <- glmmTMB((chem_class)~SoilType+Variety+Block,
                            # ziformula = ~1,
                            data=df_wcr,
                            family=nbinom2(),
                            control = glmmTMBControl(rank_check = "adjust"))
    # plot(DHARMa::simulateResiduals(model_varwcr))
    # print(wcr_2)
    # print(summary(model_wcr)$AIC[5])
    # print(car::Anova(model_wcr,type = "II"))
    ## run post hoc tests
    emmeans_wcr <- emmeans(model_wcr, ~c(SoilType,Variety),adjust="none")
    letters_wcr <- multcomp::cld(object = emmeans_wcr, Letters = letters,decreasing=TRUE)
    df_wcr_letters <- subset(letters_wcr,select=c("SoilType","Variety",".group"))
    colnames(df_wcr_letters)<-c("SoilType","Variety","wcr_letters") 
    df_wcr_letters$WCR <- wcr_2
    
    # print(df_wcr_letters)
  }
  ## loop through each cover crop treatment and run stats (negative binomial)
  diffs_df_list2<-list()
  for (k in 1:5){
    cc <- unique(df$SoilType)[k]
    cc_df <- subset(df, SoilType == cc)
    
    ## negative binomial 
    model_cc <- glmmTMB((chem_class)~Variety+WCR+Block,
                        # ziformula = ~1,
                        data=cc_df,
                        family=nbinom2(),
                        control = glmmTMBControl(rank_check = "adjust"))
    
    ## post hoc test compound
    emmeans_var <- emmeans(model_cc, ~WCR|Variety,adjust="none")
    pairs_emm <- subset(as.data.frame(pairs(emmeans_var)),select=c("Variety","p.value"))
    pairs_emm$WCR <- "no"
    
    emmeans_var1 <- emmeans(model_cc, ~c(WCR,Variety),adjust="none")
    letters_var <- multcomp::cld(object = emmeans_var1, Letters = letters,decreasing=TRUE)
    var_diffs_df <- subset(as.data.frame(letters_var),select = c("Variety","WCR",".group"))
    var_diffs_df$.group<-str_trim(var_diffs_df$.group)
    
    ## if there is no difference, replace tukey letters with blank
    if (length(unique(var_diffs_df$.group))==1){
      var_diffs_df$.group <- ""
    }
    
    ## put results in data frame with significance stars
    var_diffs_df$SoilType <- cc
    var_diffs_df<-merge(var_diffs_df,pairs_emm,by=c("Variety","WCR"),all.x=T)
    var_diffs_df$sig<-NA
    var_diffs_df <- var_diffs_df %>% mutate(sig = case_when(
      var_diffs_df$p.value <0.001 ~ '***',
      var_diffs_df$p.value <0.01 ~ '**',
      var_diffs_df$p.value <0.05 ~ '*',
      # var_diffs_df$p.value > 0.05 ~ 'n.s.',
      TRUE ~ sig))
    
    ## append results data frame to list
    diffs_df_list2[[k]]<- var_diffs_df
  }
  
  var_diffs_df2 <- do.call(rbind, diffs_df_list2)
  df_agg2<-merge(df_agg, var_diffs_df2,by=c("SoilType","Variety","WCR"))
  
  ## loop through each variety and WCR treatment to get differences between cover crop treatments
  df_varwcr_list <-list()
  for (l in 1:2){
    variety_1 <- unique(df$Variety)[l]
    df_varwcr_list1 <-list()
    for (m in 1:2){
      wcr_1 <- unique(df$WCR)[m]
      df_varwcr <- subset(df, ((Variety==variety_1)&(WCR==wcr_1)))
      
      ## run negative binomia model
      model_varwcr <- glmmTMB((chem_class)~SoilType+Block,
                              # ziformula = ~1,
                              data=df_varwcr,
                              family=nbinom2(),
                              control = glmmTMBControl(rank_check = "adjust"))
      # plot(DHARMa::simulateResiduals(model_varwcr))
      
      ## run post hoc tests
      emmeans_varwcr <- emmeans(model_varwcr, ~SoilType,adjust="none")
      letters_varwcr <- multcomp::cld(object = emmeans_varwcr, Letters = letters,decreasing=TRUE)
      df_varwcr_letters <- subset(letters_varwcr,select=c("SoilType",".group"))
      colnames(df_varwcr_letters)<-c("SoilType","varwcr_letters") 
      
      ## if there are no differences, replace tukey letters with blanks
      if (length(unique(df_varwcr_letters$varwcr_letters))==1){
        df_varwcr_letters$varwcr_letters <- ""
      }
      
      df_varwcr_letters$WCR <- wcr_1
      
      ## append results dataframe to list
      df_varwcr_list1[[m]] <-(df_varwcr_letters)
    }
    
    ## combind lists of results dataframes
    df_varwcr_letters2 <- do.call(rbind,df_varwcr_list1)
    df_varwcr_letters2$Variety <- variety_1
    df_varwcr_list[[l]]<-df_varwcr_letters2
  }
  
  ## combind lists of results dataframes
  df_varwcr_letters3<-do.call(rbind,(df_varwcr_list))
  
  ## merge results dataframe with summary dataframe
  df_agg3<-merge(df_agg2, df_varwcr_letters3,by=c("SoilType","Variety","WCR"))
  
  # print(subset(df_agg3,WCR == "yes",select = c("SoilType","Variety","WCR","varwcr_letters")))
  
  ### plot results ###
  
  ## set y limit for each plot
  ylimit <- 2*(max(df_agg3$mean))
  
  ## set post hoc letter height for each plot
  letter_height1 <- 1.75*(max(df_agg3$mean))
  
  # windowsFonts(Times=windowsFont("Times New Roman"))
  Var_plot <- ggplot(df_agg3,aes(x=SoilType, y=mean,fill=Variety,alpha=WCR)) +
    scale_alpha_discrete(range = c(0.4, 1))+
    geom_bar(stat="identity",position='dodge',width=0.9,color="black",size=1) +
    geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)),
                  width=0.4,
                  position = position_dodge(0.9),
                  colour="black",
                  size=1) +
    
    ## plot post hoc letters within each cover crop treatment
    geom_text(aes(label=varwcr_letters,x=SoilType, y=letter_height1, color=Variety,alpha=WCR),
              position=position_dodge(.9),
              vjust = -1.5, size=6,family="Times") +
    
    scale_fill_manual(values=c('royalblue','coral3'))+
    scale_color_manual(values=c('royalblue','coral3'))+
    scale_x_discrete(labels=c("CCFallow" = "Fallow", "CCPea" = "Pea","CCRadish"="Radish","CCTriticale"="Triticale","CCMix"="Mixture"))+
    
    ggtitle(chem_class) +
    labs(x=" ",y="ng Chemical Class (mean ± se)")+
    
    theme_classic()+
    theme(axis.text.x= element_text(size= 14),
          axis.text.y= element_text(size= 12),
          axis.title.y= element_text(size= 12),
          strip.text.x = element_text(size = 20),
          text = element_text(family = "Times"))
  
  Var_plot <- Var_plot+
    ## plot differences between cover crop treatments
    geom_text(aes(label=.group,x=SoilType, y=mean+se),
              color='black',
              position=position_dodge(.9),
              vjust = -1.5, size=4,family="Times") +
    ## plot differences between WCR treatments 
    geom_text(aes(label=sig,x=SoilType, y=mean+se,alpha="yes"),
              color='black',
              position=position_dodge(.9),
              vjust = -2.5, size=5,family="Times")
  
  ## save each bar graph as a pdf
  # ggsave(filename = paste('D:/OLIVIA PHD/Volatiles/Figures/20250225_class-figures/20250225_', chem_class, '.pdf', sep = ''),
  #        plot = Var_plot,
  #        height = 7, width = 7)
  Class_plots[[i]]<-Var_plot
  Class_plot_names[[i]]<-chem_class
}

# Class_plots

#### Classes Heatmap ####

## Remove chemical classes categorized as 'other' or 'unknown' for visualization
Volatiles_Classes_merged2 <- Volatiles_Classes_merged
Volatiles_Classes_merged2$MCFA <- NULL
Volatiles_Classes_merged2$Isothiocyanate <- NULL
Volatiles_Classes_merged2$Phenol <- NULL
Volatiles_Classes_merged2$Unknown <- NULL

## loop through each chemical class, subset by variety and cover crop, and run negative binomial to assess differences in WCR vs control
VOCs_subset <- colnames(Volatiles_Classes_merged2)[5:length(colnames(Volatiles_Classes_merged2))]
class_names_list <-list()
class_wcr_pval_df_list <- list()
for (i in 1:length(VOCs_subset)){
  ## subset data by chemical class
  chem_class <- VOCs_subset[i]
  class_names_list[[i]] <- chem_class
  df <- subset(Volatiles_Classes_merged2, select = c("Block","Variety","SoilType","WCR",chem_class))
  df <- data.frame(df)
  colnames(df) <- c("Block","Variety","SoilType","WCR","chem_class")

  cc_wcr_pval_df_list <- list()
  ## loop through varieties
  for (j in 1:2){
    var <- unique(df$Variety)[j]
    df_var <- subset(df, Variety == var)
    
    diff_list <-list()
    CC_list <-list()
    WCR_pval_list<-list()
    
    ## loop through cover crop treatments
    for (k in 1:5){
      cc <- unique(df_var$SoilType)[k]
      df_cc <- subset(df_var,SoilType==cc)
      df_cc$chem_class<-as.numeric(df_cc$chem_class)
      
      ## get difference between WCR and control
      df_cc_agg <- df_cc %>% dplyr::group_by(WCR) %>% dplyr::summarise(mean = mean(chem_class,na.rm=T))
      diff <- df_cc_agg$mean[1] - df_cc_agg$mean[2] # positive, no wcr greater than yes wcr
      diff_list[[k]]<-diff
      
      ## run negative binomial model
      model_wcr <- glmmTMB((chem_class)~WCR+Block,
                           # ziformula = ~1,
                           data=df_cc,
                           family=nbinom2(),
                           control = glmmTMBControl(rank_check = "adjust"))
      # print(DHARMa::simulateResiduals(model_wcr,plot=T))
      
      ## append WCR pvalue to a list
      WCR_pval <- summary(model_wcr)$coefficients$cond[,4][2]
      
      CC_list[[k]]<- cc
      WCR_pval_list[[k]] <- WCR_pval
    }
    wcr_pval_df <- data.frame("SoilType"=unlist(CC_list),"WCR_pval"=unlist(WCR_pval_list),"Diff"=unlist(diff_list))
    wcr_pval_df$Variety <- var
    
    cc_wcr_pval_df_list[[j]] <-wcr_pval_df
  }
  cc_wcr_pval_df<- do.call(rbind,cc_wcr_pval_df_list)
  cc_wcr_pval_df$Class <- chem_class
  
  class_wcr_pval_df_list[[i]] <- cc_wcr_pval_df
}

## combine all chemical class results into one data frame
class_wcr_pval_df<-do.call(rbind, class_wcr_pval_df_list)

## assign an alpha value (shade darker or lighter) depending on the p value for the visualization
class_wcr_pval_df$sig <- NA
class_wcr_pval_df <- class_wcr_pval_df %>% mutate(sig = case_when(
  class_wcr_pval_df$WCR_pval <0.001 ~ 1,
  class_wcr_pval_df$WCR_pval <0.01 ~ 0.75,
  class_wcr_pval_df$WCR_pval <0.05 ~ 0.5,
  class_wcr_pval_df$WCR_pval > 0.05 ~ 0.1,
  TRUE ~ sig))

## assign a color depending on the difference (positive vs negative) between the WCR vs the control
class_wcr_pval_df$color <- NA
class_wcr_pval_df <- class_wcr_pval_df %>% mutate(color = case_when(
  class_wcr_pval_df$Diff > 0 ~ 'WCR less',
  class_wcr_pval_df$Diff < 0 ~ 'WCR greater',
  TRUE ~ color))

## Benjamini & Hochberg adjustment (false discovery rate) 
class_wcr_pval_df$padj <- p.adjust(class_wcr_pval_df$WCR_pval, method = "BH", n = length(class_wcr_pval_df$WCR_pval))

## assign significance stars to the adjusted p values
class_wcr_pval_df$sig_adj <- NA
class_wcr_pval_df <- class_wcr_pval_df %>% mutate(sig_adj = case_when(
  class_wcr_pval_df$padj <0.001 ~ "***",
  class_wcr_pval_df$padj <0.01 ~ "**",
  class_wcr_pval_df$padj <0.05 ~ "*",
  class_wcr_pval_df$padj > 0.05 ~ NA,
  TRUE ~ sig_adj))

class_wcr_pval_df$Class <- factor(class_wcr_pval_df$Class,
                                      levels = c("Sesquiterpene","Sesquiterpenoid","Quinone",
                                                 "Monoterpene","Monoterpenoid","Ketone","Furan",
                                                 "Ester","Benzothiazole","Alkene","Alkane","Aldehyde","Alcohol"))

ggplot(class_wcr_pval_df, aes(x = SoilType, y = Class, fill = color, alpha = sig)) +
  geom_tile(height = 1) +
  geom_text(aes(label=sig_adj))+
  labs(title = "Heatmap")+
  scale_fill_manual(values=c('coral3','royalblue'))+
  scale_color_manual(values=c('coral3','royalblue'))+
  scale_x_discrete(labels=c("CCFallow" = "Fallow", "CCPea" = "Pea","CCRadish"="Radish","CCTriticale"="Triticale","CCMix"="Mixture"))+
  theme_classic() +
  facet_grid(cols=vars(Variety),drop=T,scales='free',space="free_y")+
  theme(axis.text.x= element_text(size= 10,angle=45,vjust=1,hjust=1),
        text = element_text(family="Times"))

#### Individual Compounds ####

## get list of volatiles remaining
VOCs_subset <- colnames(Volatile_matrix_norm2[,5:length(colnames(Volatile_matrix_norm2))])
total_col <- length(colnames(Volatile_matrix_norm2))
VOC_num <-length(VOCs_subset) # 62 volatile remain

## loop through each compound, run statistical tests, and plot
Compound_plots <- list()
Compound_plot_names <-list()
WCR_pval_list <- list()
compound_names_list <-list()
for (i in 1:VOC_num){
  ## get the compound name
  compound <- colnames(Volatile_matrix_norm2)[-c(1:4)][i]
  compound_names_list[[i]] <- compound
  
  ## subset data frame
  df <- subset(Volatile_matrix_norm2, select = c("Block","Variety","SoilType","WCR",compound))
  df <- data.frame(df)
  colnames(df) <- c("Block","Variety","SoilType","WCR","compound")
  df$Block <- as.factor(df$Block)

  ## run a negative binomial glm
  model1 <- glmmTMB((compound)~Variety+SoilType+WCR+Block,
                    # ziformula = ~1,
                    data=df,
                    family=nbinom2(),
                    control = glmmTMBControl(rank_check = "adjust"))
  # print(DHARMa::simulateResiduals(model1,plot=T))
  AOV_model1 <-(car::Anova(model1, type = "II"))
  # print(AOV_model1)

  ## assign p-values to values
  Var_pval <- AOV_model1$`Pr(>Chisq)`[1]
  CC_pval <- AOV_model1$`Pr(>Chisq)`[2]
  WCR_pval <- AOV_model1$`Pr(>Chisq)`[3]
  
  ## put WCR pvalues in a list
  WCR_pval_list[[i]] <- WCR_pval
  
  ## get summary data for each compound
  df_agg <- df %>% dplyr::group_by(Variety,SoilType,WCR) %>% dplyr::summarise(mean = mean(compound),
                                                            sd = sd(compound),
                                                            n = n(),
                                                            count = sum(compound >0),
                                                            se = sd/sqrt(n))
  ## loop through each cover crop treatment and run stats (negative binomial)
  diffs_df_list2<-list()
  for (k in 1:5){
    cc <- unique(df$SoilType)[k]
    cc_df <- subset(df, SoilType == cc)
    
    ## negative binomial 
    model_cc <- glmmTMB((compound)~Variety+WCR+Block,
                        # ziformula = ~1,
                        data=cc_df,
                        family=nbinom2(),
                        control = glmmTMBControl(rank_check = "adjust"))
    
    ## post hoc test
    emmeans_var <- emmeans(model_cc, ~WCR|Variety,adjust="none")
    pairs_emm <- subset(as.data.frame(pairs(emmeans_var)),select=c("Variety","p.value"))
    pairs_emm$WCR <- "no"
    
    emmeans_var1 <- emmeans(model_cc, ~c(WCR,Variety),adjust="none")
    letters_var <- multcomp::cld(object = emmeans_var1, Letters = letters,decreasing=TRUE)
    var_diffs_df <- subset(as.data.frame(letters_var),select = c("Variety","WCR",".group"))
    var_diffs_df$.group<-str_trim(var_diffs_df$.group)
    
    ## if there is no difference, replace tukey letters with blank
    if (length(unique(var_diffs_df$.group))==1){
      var_diffs_df$.group <- ""
    }
    
    ## put results in data frame with significance stars
    var_diffs_df$SoilType <- cc
    var_diffs_df<-merge(var_diffs_df,pairs_emm,by=c("Variety","WCR"),all.x=T)
    var_diffs_df$sig<-NA
    var_diffs_df <- var_diffs_df %>% mutate(sig = case_when(
      var_diffs_df$p.value <0.001 ~ '***',
      var_diffs_df$p.value <0.01 ~ '**',
      var_diffs_df$p.value <0.05 ~ '*',
      # var_diffs_df$p.value > 0.05 ~ 'n.s.',
      TRUE ~ sig))
    
    ## append results data frame to list
    diffs_df_list2[[k]]<- var_diffs_df
  }
  
  var_diffs_df2 <- do.call(rbind, diffs_df_list2)
  df_agg2<-merge(df_agg, var_diffs_df2,by=c("SoilType","Variety","WCR"))
  
  ## loop through each variety and WCR treatment to get differences between cover crop treatments
  df_varwcr_list <-list()
  for (l in 1:2){
    variety_1 <- unique(df$Variety)[l]
    df_varwcr_list1 <-list()
    for (m in 1:2){
      wcr_1 <- unique(df$WCR)[m]
      df_varwcr <- subset(df, ((Variety==variety_1)&(WCR==wcr_1)))
      
      ## run negative binomia model
      model_varwcr <- glmmTMB((compound)~SoilType+Block,
                              # ziformula = ~1,
                              data=df_varwcr,
                              family=nbinom2(),
                              control = glmmTMBControl(rank_check = "adjust"))
      # plot(DHARMa::simulateResiduals(model_varwcr))
      
      ## run post hoc tests
      emmeans_varwcr <- emmeans(model_varwcr, ~SoilType,adjust="none")
      letters_varwcr <- multcomp::cld(object = emmeans_varwcr, Letters = letters,decreasing=TRUE)
      df_varwcr_letters <- subset(letters_varwcr,select=c("SoilType",".group"))
      colnames(df_varwcr_letters)<-c("SoilType","varwcr_letters") 
      
      ## if there are no differences, replace tukey letters with blanks
      if (length(unique(df_varwcr_letters$varwcr_letters))==1){
        df_varwcr_letters$varwcr_letters <- ""
      }
      
      df_varwcr_letters$WCR <- wcr_1
      
      ## append results dataframe to list
      df_varwcr_list1[[m]] <-(df_varwcr_letters)
    }
    
    ## combind lists of results dataframes
    df_varwcr_letters2 <- do.call(rbind,df_varwcr_list1)
    df_varwcr_letters2$Variety <- variety_1
    df_varwcr_list[[l]]<-df_varwcr_letters2
  }
  
  ## combind lists of results dataframes
  df_varwcr_letters3<-do.call(rbind,(df_varwcr_list))
  
  ## merge results dataframe with summary dataframe
  df_agg3<-merge(df_agg2, df_varwcr_letters3,by=c("SoilType","Variety","WCR"))
  
  ### plot results ###
  
  ## set y limit for each plot
  ylimit <- 2*(max(df_agg3$mean))
  
  ## set post hoc letter height for each plot
  letter_height1 <- 1.75*(max(df_agg3$mean))

  # windowsFonts(Times=windowsFont("Times New Roman"))
  Var_plot <- ggplot(df_agg3,aes(x=SoilType, y=mean,fill=Variety,alpha=WCR)) +
    scale_alpha_discrete(range = c(0.4, 1))+
    geom_bar(stat="identity",position='dodge',width=0.9,color="black",size=1) +
    geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)),
                  width=0.4,
                  position = position_dodge(0.9),
                  colour="black",
                  size=1) +
    
    ## plot post hoc letters within each cover crop treatment
    geom_text(aes(label=varwcr_letters,x=SoilType, y=letter_height1, color=Variety,alpha=WCR),
              position=position_dodge(.9),
              vjust = -1.5, size=6,family="Times") +
  
    scale_fill_manual(values=c('royalblue','coral3'))+
    scale_color_manual(values=c('royalblue','coral3'))+
    scale_x_discrete(labels=c("CCFallow" = "Fallow", "CCPea" = "Pea","CCRadish"="Radish","CCTriticale"="Triticale","CCMix"="Mixture"))+
    
    ggtitle(compound) +
    labs(x=" ",y="ng Compound (mean ± se)")+
    
    theme_classic()+
    theme(axis.text.x= element_text(size= 14),
          axis.text.y= element_text(size= 12),
          axis.title.y= element_text(size= 12),
          strip.text.x = element_text(size = 20),
          text = element_text(family = "Times"))
  
  Var_plot <- Var_plot+
    ## plot differences between cover crop treatments
    geom_text(aes(label=.group,x=SoilType, y=mean+se),
              color='black',
              position=position_dodge(.9),
              vjust = -1.5, size=4,family="Times") +
    ## plot differences between WCR treatments 
    geom_text(aes(label=sig,x=SoilType, y=mean+se,alpha="yes"),
              color='black',
              position=position_dodge(.9),
              vjust = -2.5, size=5,family="Times")
  
  ## save each bar graph as a pdf
  # ggsave(filename = paste('D:/OLIVIA PHD/Volatiles/Figures/20250204_compound-figures/20250204_', compound, '.pdf', sep = ''),
  #        plot = Var_plot,
  #        height = 7, width = 7)
  # print(Var_plot)
  Compound_plots[[i]]<-Var_plot
  Compound_plot_names[[i]]<-compound
}
## name compound plots
names(Compound_plots) <- Compound_plot_names
## View any plot of interest
# Compound_plots[2]

## Combine WCR results into dataframe
WCR_pval_df <- data.frame("Compound"=unlist(compound_names_list),"WCR_pval"=unlist(WCR_pval_list))

## pull out compounds with significant differences between WCR treatments
WCR_sig_compounds <- subset(WCR_pval_df,WCR_pval<0.05)['Compound']
# 51 compounds where WCR is significant


#### Individual Heatmap ####
VOCs_subset <- colnames(Volatile_matrix_norm2[,5:length(colnames(Volatile_matrix_norm2))])

## loop through compounds
compound_wcr_pval_df_list <- list()
compound_names_list <- list()
for (i in 1:length(VOCs_subset)){
  compound <- VOCs_subset[i]
  compound_names_list[[i]] <- compound
  df <- subset(Volatile_matrix_norm2, select = c("Block","Variety","SoilType","WCR",compound))
  df <- data.frame(df)
  colnames(df) <- c("Block","Variety","SoilType","WCR","compound")

  ## loop through each variety and cover crop, and run a model looking at differences between WCR and no WCR treatment
  cc_wcr_pval_df_list <- list()
  for (j in 1:2){
    ## loop through variety
    var <- unique(df$Variety)[j]
    df_var <- subset(df, Variety == var)
    
    diff_list <-list()
    CC_list <-list()
    WCR_pval_list<-list()
    for (k in 1:5){
      ## loop through cover crop
      cc <- unique(df_var$SoilType)[k]
      df_cc <- subset(df_var,SoilType==cc)
      
      ## get average compound concentration for WCR vs control
      df_cc_agg <- df_cc %>% dplyr::group_by(WCR) %>% dplyr::summarise(mean = mean(compound))
      
      ## calculate whether control or WCR treatment has greater compound concentration
      diff <- df_cc_agg$mean[1] - df_cc_agg$mean[2] # positive, no wcr greater than yes wcr
      diff_list[[k]]<-diff
      
      # model_wcr <- nlme::lme(sqrt(compound)~WCR,random=~1|Block,data=df_cc)
      # print(car::qqPlot(resid(model_wcr))) # this model wasn't as good
      
      ## run negative binomial model to compare WCR and control for each treatment
      model_wcr <- glmmTMB((compound)~WCR+Block,
                          # ziformula = ~1,
                          data=df_cc,
                          family=nbinom2(),
                          control = glmmTMBControl(rank_check = "adjust"))
      # print(DHARMa::simulateResiduals(model_wcr,plot=T))
      
      ## extract model p value
      sum_model_wcr <- summary(model_wcr)
      ttable_model <- as.data.frame(sum_model_wcr$coefficients$cond)
      WCR_pval <- ttable_model$`Pr(>|z|)`[2]
      
      ## append data to lists
      CC_list[[k]]<- cc
      WCR_pval_list[[k]] <- WCR_pval
    }
    ## concat lists into data frame
    wcr_pval_df <- data.frame("SoilType"=unlist(CC_list),"WCR_pval"=unlist(WCR_pval_list),"Diff"=unlist(diff_list))
    wcr_pval_df$Variety <- var
    cc_wcr_pval_df_list[[j]] <-wcr_pval_df
  }
  ## concat results into dataframe
  cc_wcr_pval_df<- do.call(rbind,cc_wcr_pval_df_list)
  cc_wcr_pval_df$Compound <- compound
  compound_wcr_pval_df_list[[i]] <- cc_wcr_pval_df
}

compound_wcr_pval_df<-do.call(rbind, compound_wcr_pval_df_list)

## assign alpha (color) values depending on significance
compound_wcr_pval_df$sig <- NA
compound_wcr_pval_df <- compound_wcr_pval_df %>% mutate(sig = case_when(
  compound_wcr_pval_df$WCR_pval <0.001 ~ 1,
  compound_wcr_pval_df$WCR_pval <0.01 ~ 0.75,
  compound_wcr_pval_df$WCR_pval <0.05 ~ 0.5,
  compound_wcr_pval_df$WCR_pval > 0.05 ~ 0.1,
  TRUE ~ sig))

## assign color (warm vs cool) based on whether compound is found in higher concentrations in control vs WCR
compound_wcr_pval_df$color <- NA
compound_wcr_pval_df <- compound_wcr_pval_df %>% mutate(color = case_when(
  compound_wcr_pval_df$Diff > 0 ~ 'WCR less',
  compound_wcr_pval_df$Diff < 0 ~ 'WCR greater',
  TRUE ~ color))

## Benjamini & Hochberg adjustment (false discovery rate)
compound_wcr_pval_df$padj <- p.adjust(compound_wcr_pval_df$WCR_pval, method = "BH", n = length(compound_wcr_pval_df$WCR_pval))
compound_wcr_pval_df$sig_adj <- NA
compound_wcr_pval_df <- compound_wcr_pval_df %>% mutate(sig_adj = case_when(
  compound_wcr_pval_df$padj <0.001 ~ "***",
  compound_wcr_pval_df$padj <0.01 ~ "**",
  compound_wcr_pval_df$padj <0.05 ~ "*",
  compound_wcr_pval_df$padj > 0.05 ~ NA,
  TRUE ~ sig_adj))

## Get classifications 
## Import Class Data
Classes <- read.csv("https://raw.githubusercontent.com/oliviatrase/Trase_etal_2024_Volatiles/refs/heads/main/20240520_Classes.csv")
compound_wcr_pval_df1 <- merge(Classes,compound_wcr_pval_df,by="Compound",all.y=T)

## Create classification for 'Other'
compound_wcr_pval_df1$Class2<-NA
compound_wcr_pval_df1 <- compound_wcr_pval_df1 %>% mutate(Class2 = case_when(
  compound_wcr_pval_df1$Class == "MCFA" ~ "Other",
  compound_wcr_pval_df1$Class == "Isothiocyanate" ~ "Other",
  compound_wcr_pval_df1$Class == "Phenol" ~ "Other",
  TRUE ~ Class))

## refactor compound classes for visualization
compound_wcr_pval_df1$Class2 <- factor(compound_wcr_pval_df1$Class2,
                                      levels=c("Sesquiterpene","Sesquiterpenoid","Salicylate","Quinone",
                                               "Monoterpene","Monoterpenoid","Ketone","Furan",
                                               "Ester","Benzothiazole","Alkene","Alkane","Aldehyde","Alcohol",
                                               "Unknown","Other"))
## remove unknowns
classes_to_remove <- c("Unknown")
compound_wcr_pval_df1 <- compound_wcr_pval_df1[compound_wcr_pval_df1$Class2 %ni% classes_to_remove,]

## create heatmap from results
compound_wcr_pval_df1$Compound<-as.factor(compound_wcr_pval_df1$Compound)
ggplot(compound_wcr_pval_df1, aes(x = SoilType, y = Compound_new, fill = color, alpha = sig)) +
  geom_tile(height = 1) +
  geom_text(aes(label=sig_adj))+
  scale_fill_manual(values=c('coral3','royalblue'))+
  scale_color_manual(values=c('coral3','royalblue'))+
  scale_x_discrete(labels=c("CCFallow" = "Fallow", 
                            "CCPea" = "Pea",
                            "CCRadish"="Radish",
                            "CCTriticale"="Triticale",
                            "CCMix"="Mixture"))+
  theme_classic() +
  facet_grid(vars(Class2),vars(Variety),scales='free',space="free_y")+
  theme(strip.text.y.right = element_text(angle = 0,family="Times"),strip.background = element_blank(),
        axis.text.x= element_text(size= 10,angle=45,vjust=1,hjust=1,family="Times"),
        axis.text.y= element_text(size=10,family="Times"),
        text=element_text(family="Times"))

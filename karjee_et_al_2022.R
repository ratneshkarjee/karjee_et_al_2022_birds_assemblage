# Measuring bird diversity 
# load required packages
library(BiodiversityR)
library(tidyverse)
library(viridis)
library(xlsx)



#Clean Environment
rm(list = ls())


# Choose file directory
setwd("D:/bird_Review/ProofReading")
bird <- read.xlsx("bird_assem.xlsx", sheetIndex = 2) # sheet 2 contains data
nrow(bird) # Check rows 
colnames(bird) # check columns

# Find uniqueness in rows
unique(bird$Season) # For Habitat only
unique(bird$Sub.Plot) # Sub plot


#Make a new data frame for each habitat.
RA <- filter(bird, Habitat =='RA') %>% select(ap:gi)
WB <- filter(bird, Habitat =='WB') %>% select(ap:gi)
CW <- filter(bird, Habitat =='CW') %>% select(ap:gi) 
SF <- filter(bird, Habitat =='SF') %>% select(ap:gi)

                          ############################
#========================== MEASURING BIODIVERSITY ============================#
                          ###########################

# Calculate indices-based richness 
comu.sp = specpool(bird[, 4:120])
sitwise = specpool(bird[, 4:120], pool = bird$Habitat)



#Overview of estimated richness
comu.sp 
sitwise


#Make a data frame of the estimated richness
df <- rbind.data.frame(comu.sp, sitwise)
head(df)

# Export the table
write.xlsx(x= df, file =  "estimatedRichness.xlsx")



# Plot richness estimator
par(mar=c(5,5,5,5), mfrow=c(1,5))
rich.es = poolaccum(bird[, 4:120])
jpeg(file = "estimator_S.jpeg", width = 15,
     height = 10, units = "cm", res = 800)
plot(rich.es, xlab= "size", ylab="richness")
dev.off()



                      ######################################
#=====================SPECIES AREA CURVE FOR SITES/Habitats====================#
                      ######################################

#Calculate Species Area Curve
# We simply generate rarefaction for each habitat
sacRA <- specaccum(RA, method = 'rarefaction')
sacWB <- specaccum(WB, method = 'rarefaction')
sacCW <- specaccum(CW, method = 'rarefaction')
sacSF <- specaccum(SF, method = 'rarefaction')


# Export SAC Plot
jpeg(file = "Sitewise_SAC2.jpeg", width = 17, 
     height = 15, units = "cm", res = 800)
par(mar=c(5,5,5,5), mfrow=c(2,2))

plot(sacRA, main = "Residential Area", 
     font.main=4,xlab="Individuals", ylab="Rarefraction",
     xvar = "individuals",  add = FALSE, col = "#CC0000", 
     ci.type = "polygon", ci.lty=0, lwd=1,ci.col = "red")

plot(sacWB, main = "Water Bodies", 
     font.main=4,xlab="Individuals", ylab="Rarefraction",
     xvar = "individuals",  add = FALSE, col = "navy", 
     ci.type = "polygon", ci.lty=0, lwd=1,ci.col = "royal blue")

plot(sacCW, main = "Cropland & Wasteland", 
     font.main=4,xlab="Individuals", ylab="Rarefraction",
     xvar = "individuals",  add = FALSE, col = "red", 
     ci.type = "polygon", ci.lty=0, lwd=1,ci.col = "deep pink")

plot(sacSF, main = "Sal Forest", 
     font.main=4,xlab="Individuals", ylab="Rarefraction",
     xvar = "individuals",  add = FALSE, col = "lime green", 
     ci.type = "polygon", ci.lty=0, lwd=1,ci.col = "green")

dev.off()





# Data rearrangement....
df1 <- bird %>% select(-c(1,3)) # Remove unwanted column



# Generate habitat-wise abundance table
df2 <- df1 %>% gather(bird.com, count, 
                      ap:gi) %>% group_by(Habitat) %>% 
  summarise(sum_of_counts = sum(count)) %>% ungroup() %>% 
  mutate(total_count = sum(sum_of_counts))
head(df2)



# Prepare species composition table
# Make separate table for each site
cRA <- colSums(RA)
cWB <- colSums(WB)
cCW <- colSums(CW)
cSF <- colSums(SF)

# Prepare data for beta diversity 
# Unite all the column into a single data frame
df3 <- as.data.frame(rbind(RA=cRA, WB=cWB, CW=cCW, SF=cSF))
betaDiv <- t(df3)
betaDiv <- as.data.frame(betaDiv)
write.xlsx(x= betaDiv, file = "betaDIV.xlsx") # Export as xlsx
betaDiv <-read.xlsx("betaDIV.xlsx", sheetIndex = 1)
#================================E-------N-------D=============================



                         ###########################
#======================== CALCULATE BETA DIVERSITY ============================#
                        ###########################
# Check the data set
head(betaDiv)
bDiv <- betaDiv %>% select(RA:SF)
tbDiv = t(bDiv)
dist(tbDiv, method = "euclidean")
vegdist(tbDiv, method = "bray", diag = TRUE, upper = FALSE)
vegdist(tbDiv, method = "gower")
vegdist(tbDiv, method = "morisita")
designdist(tbDiv, method = "1-(A+B-2*J)/(A+B-J)", 
           terms = "binary", name = "Jaccard")

#================================== E*N*D* ====================================#



                             #########################
#============================  Rank abundance curve  ==========================#
                            #########################

# Prepare rank abundance matrix [L25-L29]
AbuRA = as.data.frame(rankabundance(RA))
AbuWB = as.data.frame(rankabundance(WB))
AbuCW = as.data.frame(rankabundance(CW))
AbuSF = as.data.frame(rankabundance(SF))


# Export rank-abundance plot
# Overlapping Rank Abundance Plot 
jpeg(file = "rankabu_overlap.jpeg", width = 20, 
     height = 15, units = "cm", res = 300)
par(mfrow = c(1, 1))

rankabunplot(AbuWB,scale='abundance', 
             main = NA, 
             font.main=4, col = "royal blue", lwd= 2, pch=1, 
             add=F,specnames=c(1),ylim = c(0,461),xlim = c(0,117)) 
rankabunplot(AbuRA,scale='abundance', 
             main = NA, 
             font.main=4, col = "red", lwd= 2, pch=1, 
             add=T,specnames=c(1), ylim = c(0,461),xlim = c(0,117))
rankabunplot(AbuCW,scale='abundance', 
             main = NA, 
             font.main=4, col = "deep pink", lwd= 2, pch=1, 
             add=T,specnames=c(1),ylim = c(0,461),xlim = c(0,117))
rankabunplot(AbuSF,scale='abundance', 
             main = NA, 
             font.main=4, col = "green", lwd= 2, pch=1, 
             add=T,specnames=c(1),ylim = c(0,461),xlim = c(0,117))

legend("topright", 
       legend = c("Water Bodies", "Residential Area",
                  "Cropland & Wasteland", "Sal Forest" ),
       col = c("royal blue", "red", "deep pink ",  "green" ),
       text.col = "black", horiz = F, bty = "n", 
       pt.cex = .5, cex = .5, lwd=3)

dev.off()



# Site-wise Rank Abundance Plot
jpeg(file = "rankabu_sitewise.jpeg", width = 15, 
     height = 15, units = "cm", res = 300)
par(mfrow=c(2,2))
rankabunplot(AbuRA,scale='abundance', 
             main = "Residential Area", 
             font.main=4, col = "red", lwd= 2, pch=1, 
             addit=FALSE,specnames=c(1:5), add=T, ylim = c(0, 191)) 
rankabunplot(AbuWB,scale='abundance', 
             main = "Water Bodies", 
             font.main=4, col = "royal blue", lwd= 2, pch=1, 
             addit=FALSE,specnames=c(1:5), add=T, ylim = c(0,461)) 
rankabunplot(AbuCW,scale='abundance', 
             main = "Cropland & Wasteland", 
             font.main=4, col = "deep pink", lwd= 2, pch=1, 
             addit=FALSE,specnames=c(1:5), add=T, ylim = c(0,147))
rankabunplot(AbuSF,scale='abundance', 
             main = "Sal Forest", 
             font.main=4, col = "green", lwd= 2, pch=1, 
             addit=FALSE,specnames=c(1:5), add=T, ylim = c(0, 121))
dev.off()
#===============================#E---N---D#====================================#





                              ########################
#============================= BIODIVERSITY ANALYSIS ==========================# 
                             ########################



# Residential Area
hRA = diversityresult(RA,y=NULL, index="Shannon",
                      method="pooled", digits=4)
hRA

# Water Bodies
hWB = diversityresult(WB,y=NULL, index="Shannon",
                      method="pooled", digits=4)
hWB


# For Cropland & Wasteland
hCW = diversityresult(CW,y=NULL, index="Shannon", 
                      method="pooled", digits=4)
hCW


# For Forest
hSF = diversityresult(SF,y=NULL, index="Shannon", 
                      method="pooled", digits=4)
hSF


#Calculate Inverse Simpson
# For Residential Area
cdRA = diversityresult(RA, y=NULL, index="Simpson", 
                       method="pooled", digits=3)
cdRA

# For Waterbodies
cdWB = diversityresult(WB, y=NULL, index="Simpson", 
                       method="pooled", digits=3)
cdWB

# For Cropland & Wasteland
cdCW = diversityresult(CW, y=NULL, index="Simpson", 
                       method="pooled", digits=3)
cdCW

# For Forest
cdSF = diversityresult(SF, y=NULL, index="Simpson", 
                       method="pooled", digits=3)
cdSF


#=============================Inverse Simpson=============================#

# For Residential Area
ivRA = diversityresult(RA, y=NULL, index="inverseSimpson", 
                       method="pooled", digits=3)
ivRA


# For Waterbodies
ivWB = diversityresult(WB, y=NULL, index="inverseSimpson", 
                       method="pooled", digits=3)
ivWB


# For Cropland & Wasteland
ivCW = diversityresult(CW, y=NULL, index="inverseSimpson", 
                       method="pooled", digits=3)
ivCW


# For Forest
ivFR = diversityresult(SF, y=NULL, index="inverseSimpson", 
                       method="pooled", digits=3)

ivFR



#==================================JEvenness==============================#


# For Residential Area
jRA = diversityresult(RA, y=NULL, index="Jevenness",
                      method="pooled", digits=3)
jRA



# For Waterbodies
jWB = diversityresult(WB, y=NULL, index="Jevenness",
                      method="pooled", digits=3)
jWB


# For Cropland & Wasteland
jCW = diversityresult(CW, y=NULL, index="Jevenness",
                      method="pooled", digits=3)
jCW


# For Forest
jFR = diversityresult(SF, y=NULL, index="Jevenness",
                      method="pooled", digits=3)
jFR


#======================EEvenness==================#

# For Residential Area
eRA = diversityresult(RA, y=NULL, index="Eevenness",
                      method="pooled", digits=3)
eRA


# For Waterbodies
eWB = diversityresult(WB, y=NULL, index="Eevenness",
                      method="pooled", digits=3)
eWB

# For Cropland & Wasteland
eCW = diversityresult(CW, y=NULL, index="Eevenness",
                      method="pooled", digits=3)
eCW

# For Forest
eFR = diversityresult(SF, y=NULL, index="Eevenness",
                      method="pooled", digits=3)
eFR



                            ###################
#===========================   N - M - D - S  ===============================#
                            ###################


library(ggplot2) # used for plotting graph
library(viridis) # used for color 


# Import data. [Please be sure about your working directory]
head(bird) # to check the structure of data frame


# Now split the data into two set
# First for env variable and second for community
bird.env <- bird[,1:4]
bird.com <- bird[,5:121]



# Now we perform NMDS
nmds <- metaMDS(bird.com,"bray")
nmds # To check the stress value 


# You can print your graph 
#To save your graph follow this
jpeg(file = "stressplot_4.06.21.jpg", width = 20, 
     height = 20, units = "cm", res = 800)
stressplot(nmds)
dev.off()


# Now, extract NMDS value
# Create data frame with extracted NMDS values
# Next step, "sites' and 'species" are default argument
# Use the scores function from vegan/Biodiversity to extract NMDS value
# Extract value for species
sp.score <- as.data.frame(scores(nmds, "species"))  
sp.score$species <- rownames(sp.score)#Assign values against species
head(sp.score)  


# Extract value for Site/Habitat/Group
df.score <- as.data.frame(scores(nmds, "sites"))
df.score <- cbind(as.data.frame(df.score), Habitat = bird.env$Habitat)
head(df.score) # To check extracted values against site




#Calculate group centroid ..... 
df.cent <- aggregate(cbind(NMDS1, NMDS2) ~ Habitat, data = df.score, FUN = mean)
sp.cent <- aggregate(cbind(NMDS1, NMDS2) ~ species, data=sp.score, FUN = mean)


#Plot NMDS 
ggplot()+ geom_jitter(data = df.score, aes(NMDS1, NMDS2,color=Habitat))+ 
  stat_density_2d(data=df.score, 
                  aes(NMDS1, NMDS2, color =Habitat, fill=Habitat),
                  geom = "polygon", alpha=.3) +
  geom_text(data = sp.score, aes(NMDS1, NMDS2, label=species), size=3) +
  scale_fill_viridis(discrete=TRUE) + 
  scale_color_viridis(discrete=TRUE) + theme_bw()+ coord_equal()

#Add xy coordination line in NMDS plot
ggplot()+ geom_jitter(data = df.score, aes(NMDS1, NMDS2,color=Habitat))+ 
  stat_density_2d(data=df.score, 
                  aes(NMDS1, NMDS2, color =Habitat, fill=Habitat),
                  geom = "polygon", alpha=.3) +
  geom_text(data = sp.score, aes(NMDS1, NMDS2, label=species), size=3) +
  geom_vline(xintercept = c(0), color = "black", linetype = 2) +
  geom_hline(yintercept = c(0), color = "black", linetype = 2) +
  scale_fill_viridis(discrete=TRUE) + 
  scale_color_viridis(discrete=TRUE) + theme_bw()+ coord_equal() +
  theme(panel.background = element_blank())

# Remove overlapping text
ggplot()+ geom_jitter(data = df.score, aes(NMDS1, NMDS2,color=Habitat))+ 
  stat_density_2d(data=df.score, 
                  aes(NMDS1, NMDS2, color =Habitat, fill=Habitat),
                  geom = "polygon", alpha=.3) +
  geom_vline(xintercept = c(0), color = "black", linetype = 2) +
  geom_hline(yintercept = c(0), color = "black", linetype = 2) +
  scale_fill_viridis(discrete=TRUE) + 
  scale_color_viridis(discrete=TRUE) + theme_bw()+ coord_equal() +
  theme(panel.background = element_blank())+
  geom_text_repel(data = sp.score, 
                  aes(NMDS1, NMDS2,label=species),size = 3)



# Print photo
jpeg(file = "NMDS_final123.jpeg", width = 20, 
     height = 15, units = "cm", res = 800)
par(mfrow = c(1, 1))
dev.off()


# Anova for the bird 
fit <- adonis2(bird.com ~ Habitat, data=bird.env, 
               permutations=999, method="bray")
fit


#Check assumption of homogeneity of multivariate dispersion
distances_data <- vegdist(bird.com)
aov <-anova(betadisper(distances_data, bird.env$Habitat))
aov

summary(aov)



#============================ SIGNIFICANT TEST ==============================#
#Required following package
library(ggpubr) 

# Create a new data frame
Habitat <- bird.env[-c(1,3)] # Habitat column
head(Habitat)


# Analysis bird richness for each sampling effort
Richness <-apply(bird.[,-1]>0,1,sum)
Richness <- as.data.frame(Richness)
head(Richness)


# Analysis bird abundance for each sampling effort
Abundance <-apply(bird.com[,-1],1,sum)
Abundance <-as.data.frame(Abundance)
head(Abundance)



# Create new data frame correspondence to each site
anov.hab <- cbind(Habitat, Richness, Abundance)
head(anov.hab)


# Test Sig for independence var vs dependence var

#for richness
ric.aov <- aov(Richness ~ Habitat, data = anov.hab)
summary(ric.aov)# Summary of the analysis
TukeyHSD(ric.aov) 

# for abundance
abu.aov <- aov(Abundance ~ Habitat, data = anov.hab)
summary(abu.aov)# Summary of the analysis
TukeyHSD(abu.aov) 




# Plot graphics
# for richness
jpeg(file = "anova.jpeg", width = 20, 
     height = 15, units = "cm", res = 800)
par(mfrow = c(1, 1))

rs <- ggplot(anov.hab, aes(x=Habitat, y=Richness, fill=Habitat, color=Habitat)) +
  geom_violin(trim=FALSE, alpha=.25) +
  geom_boxplot(width = .3, alpha=.4, fatten=5) +
  geom_jitter(shape=16, position=position_jitter(0.4),
              alpha=.4, size=5) + 
  stat_compare_means(method = "anova",label.y = 40)+      
  scale_fill_viridis(discrete=TRUE) + scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_text(face = 1, size=10), 
        strip.text = element_text(face=1, size=9),
        axis.text=element_text(face=1),
        axis.title = element_text(face=1),
        plot.title = element_text(face = 1, hjust = 0.5,size=13))


# For abundance
ab <-ggplot(anov.hab, aes(x=Habitat, y=Abundance, fill=Habitat, color=Habitat)) +
  geom_violin(trim=FALSE, alpha=.25) +
  geom_boxplot(width = .3, alpha=.4, fatten=5) +
  geom_jitter(shape=16, position=position_jitter(0.4),
              alpha=.4, size=5) + 
  stat_compare_means(method = "anova",label.y = 180)+      
  scale_fill_viridis(discrete=TRUE) + scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_text(face = 1, size=10), 
        strip.text = element_text(face=1, size=9),
        axis.text=element_text(face=1),
        axis.title = element_text(face=1),
        plot.title = element_text(face = 1, hjust = 0.5,size=13))

ggarrange(rs, ab,ncol = 2, nrow = 1, 
          common.legend = TRUE, legend="right")

dev.off()



                    ########################################
#================= Analysis species sharing common habitat ==================# 
                  #########################################

# We already have bird composition table
# Make it a binary data frame
betaD <- read.xlsx("betaDIV.xlsx", sheetIndex = 1)
df4 <- betaDiv %>% mutate_if(is.numeric, ~1 * (. > 0))
head(df4)
df5 <- rownames_to_column(df4)
colnames(df5)[1] <- "Species"
head(df5)
write.xlsx(x= df4, file = "binary_bird.xlsx")


# It requires only binary data (1,0)
# we already have data (df5)

df6 <- df5 %>% mutate(RACW = ifelse(RA + CW == 2,1,0),
                      RASF = ifelse(RA + SF == 2,1,0), 
                      RAWB = ifelse(RA + WB == 2,1,0),
                      CWSF = ifelse(CW + SF == 2,1,0),
                      CWWB = ifelse(CW + WB == 2,1,0),
                      SFWB = ifelse(SF + WB == 2,1,0),
                      RACWSF = ifelse(RA + CW + SF == 3,1,0),
                      RACWWB = ifelse(RA + CW + WB == 3,1,0),
                      RASFWB = ifelse(RA + SF + WB == 3,1,0),
                      CWSFWB = ifelse(CW + SF + WB == 3,1,0),
                      RACWSFWB = ifelse(RA + CW + SF + WB == 4,1,0)) %>%
  summarise(RACW = sum(RACW), RASF = sum(RASF), RAWB = sum(RAWB),
            CWSF = sum(CWSF), CWWB = sum(CWWB), SFWB = sum(SFWB),
            RACWSF = sum(RACWSF), RACWWB = sum(RACWWB),
            RASFWB = sum(RASFWB), CWSFWB = sum(CWSFWB),
            RACWSFWB = sum(RACWSFWB))


#============================= Venn Diagram ==================================#
library(VennDiagram)
grid.newpage()
jpeg(file = "VennFinal.jpeg", width = 17, 
     height = 15, units = "cm", res = 300)
par(mfrow = c(1, 1))
VD4<- draw.quad.venn(area1 = 56, area2 = 64, area3 = 54,
                     area4 = 37, n12 = 39, n13 = 26, n14 = 23,
                     n23 = 19, n24 = 19, n34 = 13, n123 = 16,
                     n124 = 16, n134 = 13, n234 = 9, n1234 = 9,
                     category = c("RA", "CW", "SF", "WB"),
                     fill = c( "orange", "red", "green", "blue"), 
                     lty = "dashed", lwd = 1.75, cex = 1, cat.cex = 1.20, 
                     cat.col = c("orange", "red", "green", "blue"))
dev.off()


#============================= HEAT MAP =============================#
library(superheat)
fd <- read.xlsx("fd_guild.xlsx", sheetIndex = 1)
fdc <- fd %>% column_to_rownames("Guild")

# plot
superheat(fdc,pretty.order.rows = TRUE,
          pretty.order.cols = TRUE, 
          scale = TRUE, row.dendrogram = TRUE,
          grid.hline.col = "white",
          grid.vline.col = "white",
          row.title = "Guild", row.title.size = 5,
          column.title = "Habitat", column.title.size = 5)

jpeg(file = "heatCluster.jpeg", width = 17, 
     height = 15, units = "cm", res = 800)
par(mar=c(5,5,5,5), mfrow=c(1,1))
dev.off()


#---------------------------- T H A N K S---------------------------#
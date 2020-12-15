#############################################################
####             Fitness ~ Genetic Diversity             ####
####                                                     ####
####                 Thibaut Capblancq                   ####
####                     3/10/2020                       ####
#############################################################

library(ggplot2)
library(reshape2)
library(rasterVis)
library(sf)
library(geoR)
library(fields)
library(maps)
library(sp)
library(raster)
library(ggpubr)
library(rgdal)
library(rgeos)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(gtable)
library(grid)
library(sjPlot)


##################
####   DATA   ####

TAB <- read.table("./FitnessTraits_GeneticParameters_RedSpruce.txt", header = T)

##################


##############################
####   Traits variation   ####

### Trait ~ family ###

# Order by latitude mean in the pop
mean <- lapply(split(TAB$Latitude, f = TAB$Population), function(x) mean(x, na.rm=T))
number <- lapply(split(TAB$Latitude, f = TAB$Population), length)
new_values <- lapply(1:length(mean), function(x) rep(mean[[x]], number[[x]]))
names(new_values) <- names(mean)
TAB$Lat_mean <- unlist(new_values)
TAB <- TAB[order(TAB$Lat_mean),]

TAB_height <- melt(TAB[,c(2,5,11)])
TAB_weight <- melt(TAB[,c(2,5,8)])
TAB_germ <- melt(TAB[,c(2,5,9)])
TAB_surv <- melt(TAB[,c(2,5,10)])
TAB_fitness <- melt(TAB[,c(2,5,12)])

table <-rbind(TAB_weight, TAB_germ, TAB_surv, TAB_height, TAB_fitness)
table$population <- as.character(table$Population)
table$population <- factor(table$population, levels = unique(table$population))
table$dummy <- "Family"
  
p1 <- ggplot(data = table, aes(x=population, y=value)) +
  geom_boxplot(fill="darkgrey", outlier.alpha = 0.1) +
  facet_grid(rows = vars(variable), cols = vars(dummy), scales="free", space = "free_x") +
  xlab("") +
  ylab("") +
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 7),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Traits ~ Region ###

TAB_boxplot <- rbind(data.frame(Region=TAB[,5], variable=TAB[,8],trait=rep("SeedWeight", nrow(TAB))), data.frame(Region=TAB[,5], variable=TAB[,9],trait=rep("Germination", nrow(TAB))), data.frame(Region=TAB[,5], variable=TAB[,10],trait=rep("Survival", nrow(TAB))), data.frame(Region=TAB[,5], variable=TAB[,11],trait=rep("Height", nrow(TAB))), data.frame(Region=TAB[,5], variable=TAB[,12],trait=rep("Fitness", nrow(TAB))))
TAB_boxplot$dummy <- "Region" 

p2 <- ggplot(TAB_boxplot, aes(Region, variable)) +
  facet_grid(rows = vars(trait), cols = vars(dummy), scales="free_y") +
  geom_violin(fill="lightgrey", colour = NA) +
  geom_boxplot(width = 0.1, fill="black") +
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 7), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf("./FitnessTraits.pdf", width = 12, height = 7)
ggarrange(p1, p2, ncol = 2, widths = c(10,4), align = "hv")
dev.off()

#############################


##########################################
####   Genetic parameters variation   ####

### DATA ###

TAB_2 <- data.frame(Longitude = -TAB$Longitude, Latitude = TAB$Latitude, Population = TAB$Population, genetic_diversity = TAB$Genetic_Diversity, genetic_load = TAB$Genetic_Load, homozygosity = TAB$Population_Homozygosity, Famhomozygosity = TAB$Family_Homozygosity, Fitness = TAB$Fitness, PC1 = TAB$PC1, PC2 = TAB$PC2)
TAB_pop <- aggregate(TAB_2, by = list(TAB_2$Population), function(x) mean(x, na.rm = T))[,-c(1,4)]

### Kriging - Figure ###

# Red spruce range
range <- readOGR("./picerube.shp")

# Create grid of prediction points: 
sp1<-seq(min(TAB_2[,1])-2,max(TAB_2[,1])+2,length=1000) 
sp2<-seq(min(TAB_2[,2])-2,max(TAB_2[,2])+2,length=1000) 
sp<-expand.grid(sp1,sp2)

# Perform ordinary Kriging (value of cov.pars and nugget are copied from mle output): 
pred_GD<-krige.conv(data=TAB_pop$genetic_diversity,coords=TAB_pop[,1:2],locations=sp,krige=krige.control(cov.model="gaussian", cov.pars=c(20,5),nugget=1)) 
pred_GL<-krige.conv(data=TAB_pop$genetic_load,coords=TAB_pop[,1:2],locations=sp,krige=krige.control(cov.model="gaussian", cov.pars=c(20,5),nugget=1)) 
pred_Hom<-krige.conv(data=TAB_pop$homozygosity,coords=TAB_pop[,1:2],locations=sp,krige=krige.control(cov.model="gaussian", cov.pars=c(20,5),nugget=1)) 
pred_Fam<-krige.conv(data=TAB_pop$Famhomozygosity,coords=TAB_pop[,1:2],locations=sp,krige=krige.control(cov.model="gaussian", cov.pars=c(20,5),nugget=1)) 

# Rasterize the prediction
GT.s <- points2grid(SpatialPoints(as.matrix(sp)))
reorder <- as.vector(matrix(1:nrow(sp), nc=slot(GT.s, "cells.dim")[2])[,slot(GT.s, "cells.dim")[2]:1])
SGDF.s_GD <- SpatialGridDataFrame(grid=GT.s, data=as.data.frame(pred_GD[1:2])[reorder,])
r_GD<-raster(SGDF.s_GD)
SGDF.s_GL <- SpatialGridDataFrame(grid=GT.s, data=as.data.frame(pred_GL[1:2])[reorder,])
r_GL<-raster(SGDF.s_GL)
SGDF.s_Hom <- SpatialGridDataFrame(grid=GT.s, data=as.data.frame(pred_Hom[1:2])[reorder,])
r_Hom<-raster(SGDF.s_Hom)
SGDF.s_Fam <- SpatialGridDataFrame(grid=GT.s, data=as.data.frame(pred_Fam[1:2])[reorder,])
r_Fam<-raster(SGDF.s_Fam)

# Crop red spruce range only
x_GD<-mask(r_GD, range)
x_GL<-mask(r_GL, range)
x_Hom<-mask(r_Hom, range)
x_Fam<-mask(r_Fam, range)

# Converting the predicted rasters into data.frame 
map.p_GD <- as.data.frame(rasterToPoints(x_GD))
map.p_GD$variable <- rep("Genetic Diversity", nrow(map.p_GD))
map.p_GL <- as.data.frame(rasterToPoints(x_GL))
map.p_GL$variable <- rep("Genetic Load", nrow(map.p_GL))
map.p_Hom <- as.data.frame(rasterToPoints(x_Hom))
map.p_Hom$variable <- rep("Population Homozygosity", nrow(map.p_Hom))
map.p_Fam <- as.data.frame(rasterToPoints(x_Fam))
map.p_Fam$variable <- rep("Family Homozygosity", nrow(map.p_Fam))

# Administrative boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))
lakes <- cbind(lakes, st_coordinates(st_centroid(lakes)))
provinces <- c("Québec", "Ontario", "New Brunswick", "Saskatchewan")
canada <- getData("GADM",country="CAN",level=1)
ca.provinces <- canada[canada$NAME_1 %in% provinces,]

# Plotting
p1 <- ggplot() + 
  geom_sf(data = world, fill=gray(.95), size=0.3) +
  geom_raster(data = map.p_GD, aes(x=x, y=y, fill=predict)) +
  geom_point(data = TAB_pop, aes(x=Longitude, y=Latitude), size = 0.1) +
  geom_sf(data = lakes, fill="#A6CAE0", size=0.1) +
  geom_sf(data = states, fill=NA, size=0.1, colour = gray(.4)) +
  geom_path(data=ca.provinces, aes(x=long,y=lat,group=group), size = 0.1, colour = gray(.4)) +
  coord_sf(xlim = c(-85, -65), ylim = c(33, 48), expand = FALSE) +
  facet_wrap(~variable) +
  scale_fill_gradient(low = "yellow", high="red") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.85, 0.2), legend.title = element_blank(), legend.key.size = unit(0.35, "cm"), legend.box.background = element_rect(colour = "black"), legend.spacing.x = unit(0.2, 'cm'))

p2 <- ggplot() + 
  geom_sf(data = world, fill=gray(.95), size=0.3) +
  geom_raster(data = map.p_GL, aes(x=x, y=y, fill=predict)) +
  geom_point(data = TAB_pop, aes(x=Longitude, y=Latitude), size = 0.1) +
  geom_sf(data = lakes, fill="#A6CAE0", size=0.1) +
  geom_sf(data = states, fill=NA, size=0.1, colour = gray(.4)) +
  geom_path(data=ca.provinces, aes(x=long,y=lat,group=group), size = 0.1, colour = gray(.4)) +
  coord_sf(xlim = c(-85, -65), ylim = c(33, 48), expand = FALSE) +
  facet_wrap(~variable) +
  scale_fill_gradient(low = "yellow", high="red") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.85, 0.2), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.title = element_blank(), legend.key.size = unit(0.35, "cm"), legend.box.background = element_rect(colour = "black"), legend.spacing.x = unit(0.2, 'cm'))

p3 <- ggplot() + 
  geom_sf(data = world, fill=gray(.95), size=0.3, colour = "black") +
  geom_raster(data = map.p_Hom, aes(x=x, y=y, fill=predict)) +
  geom_point(data = TAB_pop, aes(x=Longitude, y=Latitude), size = 0.1) +
  geom_sf(data = lakes, fill="#A6CAE0", size=0.1) +
  geom_sf(data = states, fill=NA, size=0.1, colour = gray(.4)) +
  geom_path(data=ca.provinces, aes(x=long,y=lat,group=group), size = 0.1, colour = gray(.4)) +
  coord_sf(xlim = c(-85, -65), ylim = c(33, 48), expand = FALSE) +
  facet_wrap(~variable) +
  scale_fill_gradient(low = "yellow", high="red") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.85, 0.2), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.title = element_blank(), legend.key.size = unit(0.35, "cm"), legend.box.background = element_rect(colour = "black"), legend.spacing.x = unit(0.2, 'cm'))

p4 <- ggplot() + 
  geom_sf(data = world, fill=gray(.95), size=0.3) +
  geom_raster(data = map.p_Fam, aes(x=x, y=y, fill=predict)) +
  geom_point(data = TAB_pop, aes(x=Longitude, y=Latitude), size = 0.1) +
  geom_sf(data = lakes, fill="#A6CAE0", size=0.1) +
  geom_sf(data = states, fill=NA, size=0.1, colour = gray(.4)) +
  geom_path(data=ca.provinces, aes(x=long,y=lat,group=group), size = 0.1, colour = gray(.4)) +
  coord_sf(xlim = c(-85, -65), ylim = c(33, 48), expand = FALSE) +
  facet_wrap(~variable) +
  scale_fill_gradient(low = "yellow", high="red") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.85, 0.2), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank(), legend.key.size = unit(0.35, "cm"), legend.box.background = element_rect(colour = "black"), legend.spacing.x = unit(0.2, 'cm'))

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)

g = rbind(cbind(g4, g3, size = "last"), cbind(g1, g2, size = "last"), size = "last")

pdf("./Krigging_maps_gaussian_v2.pdf", width = 10)
grid.newpage()
grid.draw(g)
dev.off()

##########################################


######################################
####   Correlation among traits   ####

### Models ###

# Seed weight
summary(lm(TAB$SeedWeight~TAB$Germination))
summary(lm(TAB$SeedWeight~TAB$Survival))
summary(lm(TAB$SeedWeight~TAB$Height))
summary(lm(TAB$SeedWeight~TAB$Fitness))

# Germination
summary(lm(TAB$Germination~TAB$Survival))
summary(lm(TAB$Germination~TAB$Height))
summary(lm(TAB$Germination~TAB$Fitness))

# Survival
summary(lm(TAB$Survival~TAB$Height))
summary(lm(TAB$Survival~TAB$Fitness))

# Height
summary(lm(TAB$Height~TAB$Fitness))

### Figure ###

# Seed weigth
TAB_Seedweight_Height <- data.frame(variable=rep("Height", nrow(TAB)), y=TAB[,"SeedWeight"], x=TAB[,"Height"])
TAB_Seedweight_Sur <- data.frame(variable=rep("Survival", nrow(TAB)), y=TAB[,"SeedWeight"], x=TAB[,"Survival"])
TAB_Seedweight_Ger <- data.frame(variable=rep("Germination", nrow(TAB)), y=TAB[,"SeedWeight"], x=TAB[,"Germination"])
TAB_Seedweight_OF <- data.frame(variable=rep("Overall Fitness", nrow(TAB)), y=TAB[,"SeedWeight"], x=TAB[,"Fitness"])
TAB_Seedweight <- rbind(TAB_Seedweight_Sur, TAB_Seedweight_Ger, TAB_Seedweight_Height, TAB_Seedweight_OF)
TAB_Seedweight$variable <- factor(as.character(TAB_Seedweight$variable), levels = c("Survival", "Germination", "Height", "Overall Fitness"))
TAB_Seedweight$dummy <- "Seed Weight"

p1 <- ggplot(TAB_Seedweight, aes(x=x, y=y)) +
  geom_point(shape=19, size = 0.3) + 
  geom_smooth(method="lm") + 
  facet_grid(rows = vars(dummy), cols = vars(variable), scales="free_x") +
  theme_bw(base_size = 12) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank())

# Germination
TAB_Germination_SW <- data.frame(variable=rep("Seed Weight", nrow(TAB)), y=TAB[,"Germination"], x=TAB[,"SeedWeight"])
TAB_Germination_Sur <- data.frame(variable=rep("Survival", nrow(TAB)), y=TAB[,"Germination"], x=TAB[,"Survival"])
TAB_Germination_Height <- data.frame(variable=rep("Height", nrow(TAB)), y=TAB[,"Germination"], x=TAB[,"Height"])
TAB_Germination_OF <- data.frame(variable=rep("Overall Fitness", nrow(TAB)), y=TAB[,"Germination"], x=TAB[,"Fitness"])
TAB_Germination <- rbind(TAB_Germination_SW, TAB_Germination_Sur, TAB_Germination_Height, TAB_Germination_OF)
TAB_Germination$variable <- factor(as.character(TAB_Germination$variable), levels = c("Seed Weight", "Survival", "Height", "Overall Fitness"))
TAB_Germination$dummy <- "Germination"

p2 <- ggplot(TAB_Germination, aes(x=x, y=y)) +
  geom_point(shape=19, size = 0.3) + 
  geom_smooth(method="lm") + 
  facet_grid(rows = vars(dummy), cols = vars(variable), scales="free_x") +
  theme_bw(base_size = 12) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank())

# Survival
TAB_Survival_SW <- data.frame(variable=rep("Seed Weight", nrow(TAB)), y=TAB[,"Survival"], x=TAB[,"SeedWeight"])
TAB_Survival_Height <- data.frame(variable=rep("Height", nrow(TAB)), y=TAB[,"Survival"], x=TAB[,"Height"])
TAB_Survival_Ger <- data.frame(variable=rep("Germination", nrow(TAB)), y=TAB[,"Survival"], x=TAB[,"Germination"])
TAB_Survival_OF <- data.frame(variable=rep("Overall Fitness", nrow(TAB)), y=TAB[,"Survival"], x=TAB[,"Fitness"])
TAB_Survival <- rbind(TAB_Survival_SW, TAB_Survival_Ger, TAB_Survival_Height, TAB_Survival_OF)
TAB_Survival$variable <- factor(as.character(TAB_Survival$variable), levels = c("Seed Weight", "Germination", "Height", "Overall Fitness"))
TAB_Survival$dummy <- "Survival"

p3 <- ggplot(TAB_Survival, aes(x=x, y=y)) +
  geom_point(shape=19, size = 0.3) + 
  geom_smooth(method="lm") + 
  facet_grid(rows = vars(dummy), cols = vars(variable), scales="free_x") +
  theme_bw(base_size = 12) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank())

# Heigth
TAB_Height_SW <- data.frame(variable=rep("Seed Weight", nrow(TAB)), y=TAB[,"Height"], x=TAB[,"SeedWeight"])
TAB_Height_Sur <- data.frame(variable=rep("Survival", nrow(TAB)), y=TAB[,"Height"], x=TAB[,"Survival"])
TAB_Height_Ger <- data.frame(variable=rep("Germination", nrow(TAB)), y=TAB[,"Height"], x=TAB[,"Germination"])
TAB_Height_OF <- data.frame(variable=rep("Overall Fitness", nrow(TAB)), y=TAB[,"Height"], x=TAB[,"Fitness"])
TAB_Height <- rbind(TAB_Height_SW, TAB_Height_Sur, TAB_Height_Ger, TAB_Height_OF)
TAB_Height$variable <- factor(as.character(TAB_Height$variable), levels = c("Seed Weight", "Survival", "Germination", "Overall Fitness"))
TAB_Height$dummy <- "Height"

p4 <- ggplot(TAB_Height, aes(x=x, y=y)) +
  geom_point(shape=19, size = 0.3) + 
  geom_smooth(method="lm") + 
  facet_grid(rows = vars(dummy), cols = vars(variable), scales="free_x") +
  theme_bw(base_size = 12) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# Fitness
TAB_Fitness_SW <- data.frame(variable=rep("Seed Weight", nrow(TAB)), y=TAB[,"Fitness"], x=TAB[,"SeedWeight"])
TAB_Fitness_Sur <- data.frame(variable=rep("Survival", nrow(TAB)), y=TAB[,"Fitness"], x=TAB[,"Survival"])
TAB_Fitness_Ger <- data.frame(variable=rep("Germination", nrow(TAB)), y=TAB[,"Fitness"], x=TAB[,"Germination"])
TAB_Fitness_Height <- data.frame(variable=rep("Height", nrow(TAB)), y=TAB[,"Fitness"], x=TAB[,"Height"])
TAB_Fitness <- rbind(TAB_Fitness_SW, TAB_Fitness_Sur, TAB_Fitness_Ger, TAB_Fitness_Height)
TAB_Fitness$variable <- factor(as.character(TAB_Fitness$variable), levels = c("Seed Weight", "Survival", "Germination", "Height"))
TAB_Fitness$dummy <- "Overall Fitness"

p5 <- ggplot(TAB_Fitness, aes(x=x, y=y)) +
  geom_point(shape=19, size = 0.3) +
  geom_smooth(method="lm") +
  facet_grid(rows = vars(dummy), cols = vars(variable), scales="free_x") +
  theme_bw(base_size = 12) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# All

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)
g5 <- ggplotGrob(p5)

g = rbind(g1, g2, g3, g4, g5, size = "first")

pdf("./Linear_regression.pdf", width = 7, height = 7)
grid.newpage()
grid.draw(g)
dev.off()

######################################


#####################################################
####   Traits ~ Genetic parameters associations  ####

### Models ###

modelGermination <- lm(Germination ~ Population_Homozygosity + Family_Homozygosity + Genetic_Diversity + Genetic_Load + SeedWeight + Region + Region*Population_Homozygosity + Region*Family_Homozygosity + Region*Genetic_Diversity + Region*Genetic_Load, data = TAB)
summary(modelGermination)
modelSurvival <- lm(Survival ~ Population_Homozygosity + Family_Homozygosity + Genetic_Diversity + Genetic_Load + SeedWeight + Region + Region*Population_Homozygosity + Region*Family_Homozygosity + Region*Genetic_Diversity + Region*Genetic_Load, data=TAB)
summary(modelSurvival)
modelHeight <- lm(Height ~ Population_Homozygosity + Family_Homozygosity + Genetic_Diversity + Genetic_Load + SeedWeight + Region + Region*Population_Homozygosity + Region*Family_Homozygosity + Region*Genetic_Diversity + Region*Genetic_Load, data=TAB)
summary(modelHeight)
modelFitness <- lm(Fitness ~ Population_Homozygosity + Family_Homozygosity + Genetic_Diversity + Genetic_Load + SeedWeight + Region + Region*Population_Homozygosity + Region*Family_Homozygosity + Region*Genetic_Diversity + Region*Genetic_Load, data=TAB)
summary(modelFitness)

### Figure - Regression coefficients for Fitness ###

p1 <- plot_model(modelGermination, show.values = TRUE, value.offset = .5, group.terms = c(1,1,1,1,2,1,2,1,1,1,1,1,1,1,1),colors = c("black", "orange"), value.size = 2.3, show.intercept = F, title = "", vline.color = "lightgrey") + 
  facet_wrap(~"Germination") + 
  geom_text(aes(x=2.3, y=60), label = "Overall model:", color = "grey40", hjust = 0, size = 2.7, family = "Times") +
  geom_text(aes(x=1.8, y=60), label = "pval < 2.2e-16", color = "grey40", hjust = 0, size = 2.7, family = "Times") +
  geom_text(aes(x=1.2, y=60), label = paste("R^2 == ", "0.279"), color = "grey40", parse = T, hjust = 0, size = 2.7, family = "Times") +
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.title.x=element_blank(), panel.grid = element_blank())
p2 <- plot_model(modelSurvival, show.values = TRUE, value.offset = .5, group.terms = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),colors = c("black", "orange"), value.size = 2.3, show.intercept = F, title = "", vline.color = "lightgrey") + 
  facet_wrap(~"Survival") + 
  geom_text(aes(x=2, y=20), label = "Overall model:", color = "grey40", hjust = 0, size = 2.7, family = "Times") +
  geom_text(aes(x=1.4, y=20), label = "pval = 0.57", color = "grey40", hjust = 0, size = 2.7, family = "Times") +
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
p3 <- plot_model(modelHeight, show.values = TRUE, value.offset = .5, group.terms = c(1,1,1,2,2,2,1,1,2,1,1,1,2,2,1),colors = c("black", "orange"), value.size = 2.3, show.intercept = F, title = "", vline.color = "lightgrey") + 
  facet_wrap(~"Height") + 
  geom_text(aes(x=2.3, y=750), label = "Overall model:", color = "grey40", hjust = 0, size = 2.7, family = "Times") +
  geom_text(aes(x=1.8, y=750), label = "pval = 1e-13", color = "grey40", hjust = 0, size = 2.7, family = "Times") +
  geom_text(aes(x=1.2, y=750), label = paste("R^2 == ", "0.233"), color = "grey40", parse = T, hjust = 0, size = 2.7, family = "Times") +
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
p4 <- plot_model(modelFitness, show.values = TRUE, value.offset = .5, group.terms = c(1,1,1,1,2,1,1,1,1,1,1,2,1,2,1), colors = c("black", "orange"), value.size = 2.3, show.intercept = F, title = "", vline.color = "lightgrey") + 
  facet_wrap(~"Overall Fitness") + 
  geom_text(aes(x=2.3, y=2400), label = "Overall model:", color = "grey40", hjust = 0, size = 2.7, family = "Times") +
  geom_text(aes(x=1.8, y=2400), label = "pval = 2.8e-15", color = "grey40", hjust = 0, size = 2.7, family = "Times") +
  geom_text(aes(x=1.2, y=2400), label = paste("R^2 == ", "0.253"), color = "grey40", parse = T, hjust = 0, size = 2.7, family = "Times") +
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

pdf("./Coefficients_Traits.pdf", height = 3.5, width = 12)
ggarrange(p1, p2, p3, p4, ncol = 4, widths = c(1.76,1,1,1))
dev.off()

### Figure - Marginal effects for Fitness ###

theme_set(theme_bw(base_size = 10))
p1 <- plot_model(modelFitness, type = "pred", terms = c("Population_Homozygosity", "Region"), title = "", colors = c("black", "black", "black")) + 
  facet_wrap(~"Population Homozygosity") + 
  ylim(-20,40) + 
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.title.x=element_blank(), panel.grid = element_blank())
p2 <- plot_model(modelFitness, type = "pred", terms = c("Family_Homozygosity", "Region"), title = "", axis.title = "", colors = c("black", "black", "black")) + 
  facet_wrap(~"Family Homozygosity") + 
  ylim(-20,40) +  
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.title.x=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
p3 <- plot_model(modelFitness, type = "pred", terms = c("Genetic_Diversity", "Region"), title = "", axis.title = "", colors = c("black", "black", "black")) + 
  facet_wrap(~"Genetic Diversity") + 
  ylim(-20,40) + 
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.title.x=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
p4 <- plot_model(modelFitness, type = "pred", terms = c("Genetic_Load", "Region"), title = "", axis.title = "", colors = c("black", "black", "black")) + 
  facet_wrap(~"Genetic Load") + 
  ylim(-20,40) + 
  theme_bw(base_size = 10, base_family = "Times") +
  theme(axis.title.x=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

pdf("./Marginal_effect_Fitness.pdf", height = 3.5, width = 12)
ggarrange(p1,p2,p3,p4, ncol = 4, common.legend = T, legend = "right", widths = c(1.1,1,1,1))
dev.off()

#####################################################

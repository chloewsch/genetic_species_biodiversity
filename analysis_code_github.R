#### Code accompanying: Schmidt, Dray, & Garroway 
# Genetic and species-level biodiversity patterns are linked by demography and ecological opportunity
# https://onlinelibrary.wiley.com/doi/abs/10.1111/evo.14407
# Dryad data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.g79cnp5qv

## Load packages:
library(tidyverse)

# dbMEMs and variation partitioning
library(vegan)
library(adespatial)

# Structural equation models (SEM):
library(lme4)
library(lmerTest)
library(piecewiseSEM)

# SAR (simultaneous autoregressive) models
library(sf)
library(spdep)
library(spatialreg)

# Plots
library(rnaturalearth)
library(viridis)
library(extrafont)

# Data:
sgdata <- read.csv("Schmidt_genetic_environments_dryad.csv", head = T)

#### dbMEM ####
## Extract coordinates
xy <- dplyr::select(sgdata, lon, lat)

# check for linear gradient
anova(lm(sgdata$poplvlsd ~ ., data=xy)) # species: yes
anova(lm(sgdata$gene_diversity ~ ., data=xy)) # genetic: yes

# Detrend data
sp.det <- resid(lm(sgdata$poplvlsd ~ ., data=xy))
gd.det <- resid(lm(sgdata$gene_diversity ~ ., data=xy))

# dbMEM
mammal.dbmem <- dbmem(xy, silent=FALSE)

## pull the spatial weights matrix out of the dbmem object
mlistw <- attributes(mammal.dbmem)$listw

# convert to dataframe
mammal.dbmem <- as.data.frame(mammal.dbmem)
length(attributes(mammal.dbmem)$values) # # of MEMs detected

# Global significance test
memlm <- lm(gd.det ~., mammal.dbmem)
summary(memlm) # genetic: yes

memlms <- lm(sp.det ~., mammal.dbmem)
summary(memlms) # species: yes

## MEM selection
## gene diversity
gr2da <- RsquareAdj(memlm)$adj.r.squared
# forward selection
gmemfwd <- forward.sel(gd.det, as.matrix(mammal.dbmem), 
                       adjR2thresh = gr2da)
# sort & extract selected MEMs
gmems <- sort(gmemfwd[,2])
gmem.red <- mammal.dbmem[,gmems]

## species diversity
sr2da <- RsquareAdj(memlms)$adj.r.squared
# forward selection
smemfwd <- forward.sel(sp.det, as.matrix(mammal.dbmem), 
                       adjR2thresh = sr2da)
# sort & extract selected MEMs
smems <- sort(smemfwd[,2])
smem.red <- mammal.dbmem[,smems]

## Shared MEMs
shared <- gmems[gmems %in% smems]
shared.df <- mammal.dbmem[,shared]

## Moran's I for each MEM to assign cut-off for broad scale Moran's I > 0.25
moranig <- moran.randtest(as.matrix(gmem.red), mlistw)
morani_gd <- moranig$obs[moranig$obs > 0.25] # these are Moran's I values

moranis <- moran.randtest(as.matrix(smem.red), mlistw)
morani_sp <- moranis$obs[moranis$obs > 0.25]

gmem.broad <- mammal.dbmem[,gmems[c(1:length(morani_gd))]]
smem.broad <- mammal.dbmem[,smems[c(1:length(morani_sp))]]


## Predicted values for maps
fit <- lm(scale(sgdata$gene_diversity)~., data = gmem.broad) 
fittedGD <- predict(fit)

fitsp <- lm(scale(sgdata$poplvlsd)~., data = smem.broad)
fittedSP <- predict(fitsp)

#### Variation partitioning ####

# Variation explained by MEMs (total spatial variation)
sp_mem <- lm(scale(sgdata$poplvlsd) ~., data = smem.red) # all MEMs selected for species
gd_mem <- lm(scale(sgdata$gene_diversity) ~., data= gmem.red) # all MEMs selected for genetic

# Variation explained by shared MEMs (shared spatial variation)
shareGD <- lm(scale(sgdata$gene_diversity) ~., data= shared.df)
shareSR <- lm(scale(sgdata$poplvlsd) ~., data= shared.df)

# GD
(gATS <- summary(gd_mem)$adj.r.squared) # total spatial variation
(gANS <- 1- gATS) # nonspatial variation
(gASS <- summary(shareGD)$adj.r.squared) # shared spatial variation
gASS/gATS # % shared spatial of total spatial variation
(gANSS <- gATS - gASS) # non-shared spatial variation


# SD: 
(ATS <- summary(sp_mem)$adj.r.squared) # total spatial variation
(ANS <- 1- ATS) # nonspatial variation
(ASS <- summary(shareSR)$adj.r.squared) # shared spatial variation
ASS/ATS # % shared spatial of total spatial variation
(ANSS <- ATS - ASS) # non-shared spatial variation

#### SEM ####
sgdata.sem <- sgdata
sgdata.sem$species <- as.factor(sgdata.sem$species)

# scale and center all variables:
sgdata.sem$genetic_diversity <- scale(sgdata.sem$gene_diversity)
sgdata.sem$species_richness <- scale(sgdata.sem$poplvlsd)
sgdata.sem$energy <- scale(sgdata.sem$PET_mean)
sgdata.sem$human_popden <- scale(log(sgdata.sem$popden+1))
sgdata.sem$body_size <- scale(log(sgdata.sem$adult_mass_g))
sgdata.sem$heterogeneity_5SI <- scale(sgdata.sem$lc_5kSI)
sgdata.sem$heterogeneity_20SI <- scale(sgdata.sem$lc_20kSI)
sgdata.sem$heterogeneity_50SI <- scale(sgdata.sem$lc_50kSI)
sgdata.sem$heterogeneity_100SI <- scale(sgdata.sem$lc_100kSI)

#### SEM 5k km2 ####
M0.5k <- psem(lm(body_size ~ energy, data = sgdata.sem),
              lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_5SI + (1|species), data = sgdata.sem),
              lm(species_richness ~ body_size + human_popden + heterogeneity_5SI, data = sgdata.sem),
              body_size %~~% human_popden, ## site level can't predict range level
              body_size %~~% heterogeneity_5SI
)

summary(M0.5k) # check model diagnostics & d-sep tests for additional suggested links

## Add energy -> SR
M1.5k <- psem(lm(body_size ~ energy, data = sgdata.sem),
              lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_5SI + (1|species), data = sgdata.sem),
              lm(species_richness ~ body_size + human_popden + heterogeneity_5SI + energy, data = sgdata.sem),
              body_size %~~% human_popden, ## site level can't predict range level
              body_size %~~% heterogeneity_5SI
)

summary(M1.5k)

# Add GD -> SR
M2.5k <- psem(lm(body_size ~ energy, data = sgdata.sem),
              lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_5SI + (1|species), data = sgdata.sem),
              lm(species_richness ~ body_size + human_popden + heterogeneity_5SI + energy + genetic_diversity, data = sgdata.sem),
              body_size %~~% human_popden, ## site level can't predict range level
              body_size %~~% heterogeneity_5SI
)

summary(M2.5k) # model acceptable

#### SEM 20k km2 ####
M0.20k <- psem(lm(body_size ~ energy, data = sgdata.sem),
               lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_20SI + (1|species), data = sgdata.sem),
               lm(species_richness ~ body_size + human_popden + heterogeneity_20SI, data = sgdata.sem),
               body_size %~~% human_popden, ## site level can't predict range level
               body_size %~~% heterogeneity_20SI
)

summary(M0.20k)

## Add energy -> SR
M1.20k <- psem(lm(body_size ~ energy, data = sgdata.sem),
               lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_20SI + (1|species), data = sgdata.sem),
               lm(species_richness ~ body_size + human_popden + heterogeneity_20SI + energy, data = sgdata.sem),
               body_size %~~% human_popden, ## site level can't predict range level
               body_size %~~% heterogeneity_20SI
)

summary(M1.20k)

## Add GD -> SR
M2.20k <- psem(lm(body_size ~ energy, data = sgdata.sem),
               lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_20SI + (1|species), data = sgdata.sem),
               lm(species_richness ~ body_size + human_popden + heterogeneity_20SI + energy + genetic_diversity, data = sgdata.sem),
               body_size %~~% human_popden, ## site level can't predict range level
               body_size %~~% heterogeneity_20SI
)

summary(M2.20k) # model acceptable

#### SEM 50k km2 ####
M0.50k <- psem(lm(body_size ~ energy, data = sgdata.sem),
               lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_50SI + (1|species), data = sgdata.sem),
               lm(species_richness ~ body_size + human_popden + heterogeneity_50SI, data = sgdata.sem),
               body_size %~~% human_popden, ## site level can't predict range level
               body_size %~~% heterogeneity_50SI
)

summary(M0.50k)

## Add energy -> SR
M1.50k <- psem(lm(body_size ~ energy, data = sgdata.sem),
               lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_50SI + (1|species), data = sgdata.sem),
               lm(species_richness ~ body_size + human_popden + heterogeneity_50SI + energy, data = sgdata.sem),
               body_size %~~% human_popden, ## site level can't predict range level
               body_size %~~% heterogeneity_50SI
)

summary(M1.50k)

## Add GD -> SR
M2.50k <- psem(lm(body_size ~ energy, data = sgdata.sem),
               lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_50SI + (1|species), data = sgdata.sem),
               lm(species_richness ~ body_size + human_popden + heterogeneity_50SI + energy + genetic_diversity, data = sgdata.sem),
               body_size %~~% human_popden, ## site level can't predict range level
               body_size %~~% heterogeneity_50SI
)

summary(M2.50k) # model acceptable

#### SEM 100k km2 ####
M0.100k <- psem(lm(body_size ~ energy, data = sgdata.sem),
                lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_100SI + (1|species), data = sgdata.sem),
                lm(species_richness ~ body_size + human_popden + heterogeneity_100SI, data = sgdata.sem),
                body_size %~~% human_popden, ## site level can't predict range level
                body_size %~~% heterogeneity_100SI
)

summary(M0.100k)

## Add energy -> SR
M1.100k <- psem(lm(body_size ~ energy, data = sgdata.sem),
                lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_100SI + (1|species), data = sgdata.sem),
                lm(species_richness ~ body_size + human_popden + heterogeneity_100SI + energy, data = sgdata.sem),
                body_size %~~% human_popden, ## site level can't predict range level
                body_size %~~% heterogeneity_100SI
)

summary(M1.100k)

## Add GD -> SR
M2.100k <- psem(lm(body_size ~ energy, data = sgdata.sem),
                lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_100SI + (1|species), data = sgdata.sem),
                lm(species_richness ~ body_size + human_popden + heterogeneity_100SI + energy + genetic_diversity, data = sgdata.sem),
                body_size %~~% human_popden, ## site level can't predict range level
                body_size %~~% heterogeneity_100SI
)

summary(M2.100k) # model acceptable

#### Testing residual spatial autocorrelation ####
# Chosen model: 5k km2
xy <- dplyr::select(sgdata.sem, lon, lat) ## same as in MEM section
xysf <- st_as_sf(xy, coords = c("lon", "lat"), crs=4326)

# Generate neighborhood matrix
nb1 <- dnearneigh(xysf, 0, 582) # minimum distance (km) to not get an empty neighbor set

BSres <- residuals(lm(body_size ~ energy, data = sgdata.sem))
moran.test(BSres, nb2listw(nb1)) # yes, I = 0.4

GDres <- residuals(lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_5SI + (1|species), data = sgdata.sem))
moran.test(GDres, nb2listw(nb1)) # yes, I = 0.01 (low)

SRres <- residuals(lm(species_richness ~ body_size + human_popden + heterogeneity_5SI + energy, data = sgdata.sem))
moran.test(SRres, nb2listw(nb1)) # yes, I = 0.43

##### Autoregressive (SAR) sem ####
# 5k
M2.5ks <- psem(spatialreg::lagsarlm(body_size ~ energy, listw = nb2listw(nb1), data = sgdata.sem),
               lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_5SI + (1|species), data = sgdata.sem),
               spatialreg::lagsarlm(species_richness ~ body_size + human_popden + heterogeneity_5SI + energy + genetic_diversity, 
                                    listw = nb2listw(nb1), data = sgdata.sem),
               body_size %~~% human_popden, ## site level can't predict range level
               body_size %~~% heterogeneity_5SI
)
summary(M2.5ks)

# 20k
M2.20ks <- psem(spatialreg::lagsarlm(body_size ~ energy, listw = nb2listw(nb1), data = sgdata.sem),
                lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_20SI + (1|species), data = sgdata.sem),
                spatialreg::lagsarlm(species_richness ~ body_size + human_popden + heterogeneity_20SI + energy + genetic_diversity, 
                                     listw = nb2listw(nb1), data = sgdata.sem),
                body_size %~~% human_popden, ## site level can't predict range level
                body_size %~~% heterogeneity_20SI
)
summary(M2.20ks)

# 50k
M2.50ks <- psem(spatialreg::lagsarlm(body_size ~ energy, listw = nb2listw(nb1), data = sgdata.sem),
                lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_50SI + (1|species), data = sgdata.sem),
                spatialreg::lagsarlm(species_richness ~ body_size + human_popden + heterogeneity_50SI + energy + genetic_diversity, 
                                     listw = nb2listw(nb1), data = sgdata.sem),
                body_size %~~% human_popden, ## site level can't predict range level
                body_size %~~% heterogeneity_50SI
)
summary(M2.50ks)

# 100k
M2.100ks <- psem(spatialreg::lagsarlm(body_size ~ energy, listw = nb2listw(nb1), data = sgdata.sem),
                 lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_100SI + (1|species), data = sgdata.sem),
                 spatialreg::lagsarlm(species_richness ~ body_size + human_popden + heterogeneity_100SI + energy + genetic_diversity, 
                                      listw = nb2listw(nb1), data = sgdata.sem),
                 body_size %~~% human_popden, ## site level can't predict range level
                 body_size %~~% heterogeneity_100SI
)
summary(M2.100ks)

#### Variation explained (Nagelkerke R2 for SAR models) ####
bsmod <- lagsarlm(body_size ~ energy, listw = nb2listw(nb1), data = sgdata.sem)
summary(bsmod, Nagelkerke = TRUE)

srmod <- lagsarlm(species_richness ~ body_size + human_popden + heterogeneity_5SI + energy, 
                  listw = nb2listw(nb1), data = sgdata.sem)
summary(srmod, Nagelkerke = TRUE)

#### SEM 5k km2, AET ####
# Repeat structural equation models using actual evapotranspiration as a measure of
# energy availability instead of PET
sgdata.sem$water <- scale(sgdata.sem$AET_mean)

M3.5ksA <- psem(lagsarlm(body_size ~ water, listw = nb2listw(nb1), data = sgdata.sem),
                lmer(genetic_diversity ~ body_size + human_popden + heterogeneity_5SI + (1|species), data = sgdata.sem),
                lagsarlm(species_richness ~ body_size + human_popden + heterogeneity_5SI + water + genetic_diversity, 
                         listw = nb2listw(nb1), data = sgdata.sem),
                body_size %~~% human_popden, ## site level can't predict range level
                body_size %~~% heterogeneity_5SI
)
summary(M3.5ksA)

bsmodw <- lagsarlm(body_size ~ water, listw = nb2listw(nb1), data = sgdata.sem)
summary(bsmodw, Nagelkerke = TRUE)

srmodw <- lagsarlm(species_richness ~ body_size + human_popden + heterogeneity_5SI + water, 
                   listw = nb2listw(nb1), data = sgdata.sem)
summary(srmodw, Nagelkerke = TRUE)

#### Effect of heterogeneity on population divergence (FST) ####
sgdata.sem.fst <- sgdata.sem[!is.na(sgdata.sem$global_fst),] # remove rows where FST = NA

## Extract coordinates 
xy.fst <- dplyr::select(sgdata.sem.fst, lon, lat)

# Detrend data
anova(lm(sgdata.sem.fst$global_fst ~ ., data=xy.fst)) # FST related to lat + lon
fst.det <- resid(lm(sgdata.sem.fst$global_fst ~ ., data=xy.fst))

# dbMEM (to control for isolation by distance)
mammal.dbmem.fst <- as.data.frame(dbmem(xy.fst, silent=FALSE))

# Global significance test
mem.fst <- lm(fst.det ~., mammal.dbmem.fst)
summary(mem.fst) # yes
fstr2da <- RsquareAdj(mem.fst)$adj.r.squared

# forward selection
fstmemfwd <- forward.sel(fst.det, as.matrix(mammal.dbmem.fst),
                         adjR2thresh = fstr2da)
# sort & extract selected MEMs
fstmems <- sort(fstmemfwd[,2])
fstmem.red <- mammal.dbmem.fst[,fstmems]

summary(lm(sgdata.sem.fst$lc_5kSI ~., data=fstmem.red)) ## MEMs 2 and 42 are related to heterogeneity

# regression:
hetmod <- lmer(scale(sgdata.sem.fst$global_fst) ~ scale(sgdata.sem.fst$lc_5kSI) + 
                 MEM1 + MEM14 + MEM27 + MEM31 + MEM47 + MEM48 + MEM68 + MEM87 + MEM92 + MEM124 +
                 MEM126 + MEM143 + MEM170 + MEM186 + MEM194 + (1|sgdata.sem.fst$species), data=fstmem.red)
plot(residuals(hetmod))
hist(residuals(hetmod))
summary(hetmod)

#### Plots ####
world <- ne_countries(scale="medium")
canadausa <- world[world$sovereignt=="Canada" | world$sovereignt=="United States of America" & world$type !="Dependency",]
rm(world)

##### Fig. 2: Predicted maps, variation partitioning, environment maps #####
# Predicted genetic diversity
gdpoint <- ggplot() + geom_polygon(data = canadausa,
                                   aes(x = long, y = lat, group = group),
                                   fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  coord_map("conic", lat0=40, orientation=c(90, 0, -90), xlim=c(-140,-50)) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  geom_point(aes(x=sgdata$lon, y=sgdata$lat, fill=fittedGD), shape=21, size=0.2) + ## 2 layers of points because I had to hack the legend to have just "min" and "max". You can't see this layer. This legend is "fill" (shown) as opposed to "color" (hidden)
  scale_fill_viridis(option = "inferno", breaks=c(0,1),labels=c("low","high"),
                     limits=c(0,1)) +
  labs(fill="") +
  geom_point(aes(x=sgdata$lon, y=sgdata$lat, color=fittedGD), show.legend = FALSE, size=2) +
  scale_color_viridis(option = "inferno") +
  guides(fill = guide_colourbar(ticks = FALSE)) +
  labs(title = "Genetic diversity") +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"))

# Predicted species richness
sppoint <- ggplot() + geom_polygon(data = canadausa,
                                   aes(x = long, y = lat, group = group),
                                   fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  coord_map("conic", lat0=40, orientation=c(90, 0, -90), xlim=c(-140,-50)) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  geom_point(aes(x=sgdata$lon, y=sgdata$lat, fill=fittedSP), shape = 21, size=0.2) +
  scale_fill_viridis(option = "inferno", breaks=c(0,1),labels=c("low","high"),
                     limits=c(0,1)) +
  guides(fill = guide_colourbar(ticks = FALSE)) +
  labs(fill="") +
  geom_point(aes(x=sgdata$lon, y=sgdata$lat, color=fittedSP), show.legend = FALSE, size=2) +
  scale_color_viridis(option = "inferno") + 
  labs(title = "Species richnes") +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"))

# Variation partitioning #
variable <- c(rep("gene_diversity", 3), rep("species_diversity", 3))
variation <- factor(rep(c("shared_spatial", "non_shared_spatial", "var_non_spatial"), 2), 
                    levels = c("shared_spatial", "non_shared_spatial", "var_non_spatial"))
value <- c(gASS, gANSS, gANS, ASS, ANSS, ANS)
sdgdR2tab <- data.frame(variable, variation, value)

pal <- c("#782144", "#d35f8d", "#e9afc6")

q <- ggplot(sdgdR2tab, aes(fill = variation, y=value, x=variable)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = pal,
                    name= "variation", labels = c("shared spatial", "non-shared spatial", "non-spatial")) +
  scale_x_discrete(labels=c("genetic diversity", "species richness")) +
  labs(title = "", x= "", y= "proportion of variation") +
  theme(text=element_text(size=18,  family="Lato Black"))


## elevation 
elevras <- raster("elevation/w001001.adf")
crs(elevras) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
elevras <- aggregate(elevras, fact=2)
elev_crop <- crop(elevras, canadausa)
elevCU <- mask(elev_crop, canadausa)
elevnao <- elevCU[!is.na(values(elevCU))]
elevnaocoords <- coordinates(elevCU)[!is.na(values(elevCU)),]
## convert to df for ggplot
elev_spdf <- cbind.data.frame(elevnaocoords, elevnao)

elevmap <- ggplot()+
  geom_polygon(data = canadausa,
               aes(x = long, y = lat, group = group),
               fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  coord_map("conic", lat0=40, orientation=c(90, 0, -90), xlim=c(-140,-50)) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +
  geom_tile(data=elev_spdf, aes(x=x, y=y, fill=elevnao), show.legend = FALSE) +
  scale_fill_viridis(option = "inferno") + 
  labs(title = "Elevation") +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"))

## PET
PET <- raster("et0_yr/et0_yr.tif")
PETras <- aggregate(PET, fact=10)
PET_crop <- crop(PETras, canadausa)
PETCU <- mask(PET_crop, canadausa)
PETnao <- PETCU[!is.na(values(PETCU))]
PETnaocoords <- coordinates(PETCU)[!is.na(values(PETCU)),]
## convert to df for ggplot
PET_df <- cbind.data.frame(PETnaocoords, PETnao)

PETmap <- ggplot()+
  geom_polygon(data = canadausa,
               aes(x = long, y = lat, group = group),
               fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  coord_map("conic", lat0=40, orientation=c(90, 0, -90), xlim=c(-140,-50)) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +
  geom_tile(data=PET_df, aes(x=x, y=y, fill=PETnao), show.legend = FALSE) +
  scale_fill_viridis(option = "inferno") + ## inferno for mammals, "D" default for amphibs
  labs(title = "Potential evapotranspiration") +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"))

## Human population density
popdenras1 <- raster("PopDen1.asc")
popden <- aggregate(popdenras1, fact=2)
popdennao <- values(popden)[!is.na(values(popden))]
popdennaocoords <- coordinates(popden)[!is.na(values(popden)),]
## convert to df for ggplot
popden_df <- cbind.data.frame(popdennaocoords, popdennao)
names(popden_df) <- c("lon", "lat", "popden")
popden_df$logpopden <- log(popden_df$popden+1)

popmap <- ggplot()+
  geom_polygon(data = canadausa,
               aes(x = long, y = lat, group = group),
               fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  coord_map("conic", lat0=40, orientation=c(90, 0, -90), xlim=c(-140,-50)) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +
  geom_tile(data=popden_df, aes(x=lon, y=lat, fill=logpopden), show.legend = FALSE) +
  scale_fill_viridis(option = "inferno") +
  labs(title = "Human population density") +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"))

## Heterogeneity map: Land cover 250m
landcov250 <- raster("landcover2010/NA_LandCover_2010_25haMMU.tif")
lc_agg <- aggregate(landcov250,fact=30,modal)
canadausa.sf <- st_as_sf(canadausa)
canusa <- st_transform(canadausa.sf, crs = crs(lc_agg))
lc_crop <- crop(lc_agg, canusa)
lcCU <- mask(lc_crop, canusa)
lc_wgs <- projectRaster(lcCU, crs = crs(canadausa))

# make dataframe for ggplot
lc_df <- rasterToPoints(lc_wgs)
lc_df <- data.frame(lc_df)
names(lc_df) <- c("lon", "lat", "lctype")
lc_df$lctype <- as.factor(lc_df$lctype)

lcmap <- ggplot()+
  geom_polygon(data = canadausa,
               aes(x = long, y = lat, group = group),
               fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  coord_map("conic", lat0=40, orientation=c(90, 0, -90), xlim=c(-140,-50)) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +
  geom_tile(data=lc_df, aes(x=lon, y=lat, fill = lctype), show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(title = "Habitat heterogeneity") +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"))

##### Raw data maps #####
gdraw <- ggplot() + geom_polygon(data = canadausa,
                                 aes(x = long, y = lat, group = group),
                                 fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  coord_map("conic", lat0=40, orientation=c(90, 0, -90), xlim=c(-140,-50)) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  labs(color="") +
  geom_point(aes(x=sgdata$lon, y=sgdata$lat, color=sgdata$gene_diversit), size=2) +
  scale_color_viridis(option = "inferno") +
  guides(color = guide_colourbar(ticks = FALSE)) +
  labs(title = "Genetic diversity (raw data)") +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"))

sppraw <- ggplot() + geom_polygon(data = canadausa,
                                  aes(x = long, y = lat, group = group),
                                  fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  coord_map("conic", lat0=40, orientation=c(90, 0, -90), xlim=c(-140,-50)) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  guides(color = guide_colourbar(ticks = FALSE)) +
  labs(color="") +
  geom_point(aes(x=sgdata$lon, y=sgdata$lat, color=sgdata$poplvlsd), size=2) +
  scale_color_viridis(option = "inferno") + 
  labs(title = "Species richness (raw data)") +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"))
rm(list=ls())  #Erase workspace
library(sf)
library(sp)
library(rstudioapi)
library(dplyr)
library(ggplot2)

setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

####Prepping and loading the surveillance data and shapefile####
it_map1 <- st_read("./Italy Shapefile/ITA_adm2.shp") #this map is used for masking the projections of the model
it_map <- map_data("italy") #this map is used for visualisation at the end 
sampled_map <- it_map1[c(1 : 31, 36 : 40, 57 : 63, 72 : 88, 101, 102),]

basilicata_S_locations <- read.csv("./BASILICATA_EC_S.csv")
basilicata_G_locations <- read.csv("./BASILICATA_EC_G.csv")
basilicata_C_locations <- read.csv("./New data/CATTLE_BUFFALO_FARMS_COORDINATES.csv")
basilicata_C_locations <- basilicata_C_locations[-c(177, 2208, 707, 1336, 2290, 2289, 1237 , 2216, 2232, 2221), ] #removing anomalous locations
salerno_locations <- read.csv("./New data/SHEEP_FARMS_COORDINATES_SALERNO_2015.csv")
salerno_locations[, 5] <- "Sheep" 
colnames(salerno_locations) <- c("Id", "Ec", "Longitude", "Latitude", "Species")

basilicata_C_locations[is.na(basilicata_C_locations)] <- 0
basilicata_C_locations$Latitude <- as.numeric(basilicata_C_locations$Latitude)
basilicata_C_locations$Longitude <- as.numeric(basilicata_C_locations$Longitude)
basilicata_C_locations = subset(basilicata_C_locations, select = -c(slaughtering.date))
basilicata_G_locations = subset(basilicata_G_locations, select = -c(Municipality))
basilicata_S_locations = subset(basilicata_S_locations, select = -c(Municipality))
salerno_locations = subset(salerno_locations, select = -c(Municipality))
colnames(basilicata_C_locations) <- c("Id", "Latitude", "Longitude", "Species", "Ec")
colnames(basilicata_G_locations) <- c("Id", "Latitude", "Longitude", "Species", "Ec")
colnames(basilicata_S_locations) <- c("Id", "Latitude", "Longitude", "Species", "Ec")


#Combining all the data into one data frame
basilicata_all_locations <- rbind(basilicata_S_locations, basilicata_G_locations, basilicata_C_locations, salerno_locations)
basilicata_all_locations$Latitude <- as.numeric(basilicata_all_locations$Latitude)
basilicata_all_locations$Longitude <- as.numeric(basilicata_all_locations$Longitude)

#Aggregating the data with the same location
d <- group_by(basilicata_all_locations, Latitude, Longitude) %>% 
  summarize(
    total = n(),
    positive = sum(Ec),
    prev = if (positive >= 1){1}
    else{0},
    no_species = toString(unique(Species))
            )

write.csv(basilicata_all_locations, "./all samples.csv")
write.csv(d, "./Unique_locations.csv")

####The spatial model####
library(leaflet)
library(raster)

#Creating a raster of the area where prevalence will be projected onto
r <- raster(ncol=240, nrow=240, xmn=6.6, xmx=18.6, ymn=35.4, ymx=47.4)
r[] <- 1:length(r)
r <- crop(r, extent(sampled_map))
r <- mask(r, sampled_map)


#Plotting the locations of samples using leaflet
pal <- colorFactor(c("darkgreen", "red"), domain = c(0, 1))

leaflet(d) %>%  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~Longitude, lat = ~Latitude, color = ~pal(prev)) %>%
  addLegend("bottomright", pal = pal, values = ~prev, title = "Farms tested positive for Ec") %>%
  addScaleBar(position = c("bottomleft"))

d$prev <- as.character(d$prev)

#Plotting the locations of samples using ggplot2
sample_locations <- ggplot(mapping=aes(x=long,y=lat))+
  geom_point(data = d , aes(x = Longitude, y = Latitude, col = prev), size = 4, 
             shape = 20,)+
  scale_color_discrete(breaks = c("0", "1"))+
  scale_color_manual(values = c("0" = "darkgreen", "1" = "red"), labels = c("Negative farm", "Positive farm"), name = "Farm location")+
  geom_polygon(data=it_map,mapping=aes(group=group),
               fill=NA,color='black')+
  coord_map()+
  theme_bw()+
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Farm Locations")+
  theme(legend.position = c(0.8,0.8), legend.key.size = unit(2.4, 'cm'), legend.key.height = unit(1, 'cm'),
        legend.title = element_text(size=28), legend.text = element_text(size=26),
        plot.title = element_text(hjust = 0.5, face = "bold", size = (34)), axis.title=element_text(size=30))

png("All data sample locations.png", width = 1200, height = 1200,)
sample_locations
dev.off()

d$prev <- as.numeric(d$prev)

#Building triangulation mesh over the sampled regions
library(INLA)
coo <- cbind(d$Longitude, d$Latitude)
mesh <- inla.mesh.2d(loc = coo, max.edge = c(0.5, 5), cutoff = 0.01)

mesh$n

plot(mesh)
points(coo, col = "red")

#Building the SPDE model on the mesh
spde <- inla.spde2.matern(mesh = mesh, alpha = 0.01)

#Index set
indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

#Projection matrix
A <- inla.spde.make.A(mesh = mesh, loc = coo)

#Prediction locations - area where the model will predict data 
dp <- rasterToPoints(r)
dim(dp)

coop <- dp[, c("x", "y")]

#Projector matrix
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

#stack for estimation stk.e
stk.e <- inla.stack(tag = "est",
                    data = list(y = d$positive, numtrials = d$total),
                    A = list(1, A),
                    effects = list(data.frame(b0 = rep(1, times = nrow(d))), s = indexs))

#stack for prediction stk.p
stk.p <- inla.stack(tag = "pred",
                    data = list(y = NA, numtrials = NA),
                    A = list(1, Ap),
                    effects = list(data.frame(b0 = rep(1, times = nrow(dp))), s = indexs))

#stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

#Formula for model
formula <- y ~ 0 + b0 + f(s, model = spde)

#Running the model
res <- inla(formula, family = "binomial", Ntrials = numtrials,
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))

summary(res)
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]

r_prev_mean <- rasterize(x = coop, y = r, field = prev_mean, fun = mean)

save.image("INLAData_all_new.RData")

#### Model output visualisation ####
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
load("INLAData_all_new.RData")


library(dplyr)
library(leaflet)
library(raster)
library(INLA)
library(ggplot2)
library(viridis)

#Plotting the mean prevalence 
r_prev_mean <- as.data.frame(r_prev_mean, xy = TRUE)
colnames(r_prev_mean) <- c("long", "lat", "prev")

m <- ggplot(it_map, aes(x=long,y=lat))+
  geom_polygon(aes(group = group), fill = "grey70", colour = "white")+
  geom_tile(data=r_prev_mean,mapping=aes(fill=prev))+
  ggtitle("A")+
  scale_fill_viridis(option = "magma",na.value = "transparent", labels=scales::percent)+
  #scale_fill_gradient2(low="blue",mid= "yellow" ,high="red",limits = c(0.14,0.79), na.value = "transparent")+
  geom_polygon(aes(group = group), fill = NA, colour = "white")+
  coord_map()+
  theme_bw()+
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Mean probability of infection of farms")+
  theme(legend.position = c(0.83,0.83), legend.key.size = unit(2.8, 'cm'), legend.key.height = unit(2, 'cm'),
        legend.title = element_text(size=42), legend.text = element_text(size=38),
        plot.title = element_text(hjust = 0.5, face = "bold", size = (40)), axis.title=element_text(size=42))

png("All data prev mean.png", width = 2000, height = 2000,)
m
dev.off()


pal <- colorNumeric("viridis", c(0.1, 0.8), na.color = "transparent")

leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_mean, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_mean), title = "Relative Risk") %>%
  addScaleBar(position = c("bottomleft"))

#Plotting the lower limit of the mean (95% CI)

r_prev_ll <- rasterize(x = coop, y = r, field = prev_ll, fun = mean)

r_prev_ll <- as.data.frame(r_prev_ll, xy = TRUE)
colnames(r_prev_ll) <- c("long", "lat", "prev")

mll <- ggplot(it_map, aes(x=long,y=lat))+
  geom_polygon(aes(group = group), fill = "grey70", colour = "white")+
  geom_tile(data=r_prev_ll,mapping=aes(fill=prev))+
  ggtitle("C")+
  scale_fill_viridis(option = "magma",na.value = "transparent", labels=scales::percent, limit = c(0.05, 0.93))+
  #scale_fill_gradient2(low="blue",mid= "yellow" ,high="red",limits = c(0.14,0.79), na.value = "transparent")+
  geom_polygon(aes(group = group), fill = NA, colour = "white")+
  coord_map()+
  theme_bw()+
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Lower limit probability of \ninfection of farms")+
  theme(legend.position = c(0.83,0.83), legend.key.size = unit(2.8, 'cm'), legend.key.height = unit(2, 'cm'),
        legend.title = element_text(size=42), legend.text = element_text(size=38),
        plot.title = element_text(hjust = 0.5, face = "bold", size = (40)), axis.title=element_text(size=42))

png("All data prev mean ll.png", width = 2000, height = 2000,)
mll
dev.off()

leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_ll, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_ll), title = "LL") %>%
  addScaleBar(position = c("bottomleft"))

#Plotting the upper limit of the mean (95% CI)

r_prev_ul <- rasterize(x = coop, y = r, field = prev_ul, fun = mean)

r_prev_ul <- as.data.frame(r_prev_ul, xy = TRUE)
colnames(r_prev_ul) <- c("long", "lat", "prev")

mul <- ggplot(it_map, aes(x=long,y=lat))+
  geom_polygon(aes(group = group), fill = "grey70", colour = "white")+
  geom_tile(data=r_prev_ul,mapping=aes(fill=prev))+
  ggtitle("B")+
  scale_fill_viridis(option = "magma",na.value = "transparent", labels=scales::percent, limit = c(0.05, 0.93))+
  #scale_fill_gradient2(low="blue",mid= "yellow" ,high="red",limits = c(0.14,0.79), na.value = "transparent")+
  geom_polygon(aes(group = group), fill = NA, colour = "white")+
  coord_map()+
  theme_bw()+
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Upper limit probability of \ninfection of farms")+
  theme(legend.position = c(0.83,0.83), legend.key.size = unit(2.8, 'cm'), legend.key.height = unit(2, 'cm'),
        legend.title = element_text(size=42), legend.text = element_text(size=38),
        plot.title = element_text(hjust = 0.5, face = "bold", size = (40)), axis.title=element_text(size=42))

png("All data prev mean ul.png", width = 2000, height = 2000,)
mul
dev.off()

leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_ul, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_ul), title = "UL") %>%
  addScaleBar(position = c("bottomleft"))

#Plotting the mean and sd of the spatial effect

rang <- apply(mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(mesh,
                            xlim = rang[, 1], ylim = rang[, 2],
                            dims = c(300, 300)
)

mean_s <- inla.mesh.project(proj, res$summary.random$s$mean)
sd_s <- inla.mesh.project(proj, res$summary.random$s$sd)

df <- expand.grid(x = proj$x, y = proj$y)
df$mean_s <- as.vector(mean_s)
df$sd_s <- as.vector(sd_s)

library(cowplot)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

plot_grid(gmean, gsd)

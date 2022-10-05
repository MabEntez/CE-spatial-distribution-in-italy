rm(list=ls())  #Erase workspace
library(sf)
library(sp)
library(rstudioapi)
library(dplyr)
library(ggplot2)

setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

####Prepping and loading the surveillance data and shapefile####
it_map1 <- st_read("./Italy Shapefile/ITA_adm2.shp")  #this map is used for masking the projections of the model
it_map <- map_data("italy") #this map is used for visualisation at the end 
bas_map <- it_map1[c(11,12,22),]

basilicata_all_locations <- read.csv("./BASILICATA_EC_S.csv")
basilicata_all_locations <- basilicata_all_locations[, c("Longitude",  "Latitude", "Ec")]
salerno_locations <- read.csv("./New data/SHEEP_FARMS_COORDINATES_SALERNO_2015.csv")
salerno_locations <- rename(salerno_locations, Ec = Positive, Latitude = POINT_Y, Longitude = POINT_X)
salerno_locations <- salerno_locations[, c("Longitude",  "Latitude", "Ec")]

all_sheep_locations <- basilicata_all_locations
all_sheep_locations <- rbind(basilicata_all_locations, salerno_locations)

#Aggregating the data with the same location
d <- group_by(all_sheep_locations, Latitude, Longitude) %>% 
  summarize(
    total = n(),
            positive = sum(Ec),
            prev = if (positive >= 1){1}
            else{0}
            )

####The spatial model####
library(leaflet)
library(raster)

#Creating a raster of the area where prevalence will be projected onto
r <- raster(ncol=240, nrow=120, xmn=13.7, xmx=16.9, ymn=39.8, ymx=41.6)
r[] <- 1:length(r)
r <- crop(r, extent(bas_map))
r <- mask(r, bas_map)

#Plotting the locations of samples using leaflet
pal <- colorFactor(c("darkgreen", "red"), domain = c(0, 1))
leaflet(d) %>%  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~Longitude, lat = ~Latitude, color = ~pal(prev)) %>%
  addLegend("bottomright", pal = pal, values = ~prev, title = "Farms tested positive for Ec") %>%
  addScaleBar(position = c("bottomleft"))

#Building triangulation mesh over the sampled regions
library(INLA)
coo <- cbind(d$Longitude, d$Latitude)
mesh <- inla.mesh.2d(loc = coo, max.edge = c(0.1, 5), cutoff = 0.01)

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

save.image("INLAData_S_bascamp.RData")

#### Model output visualisation ####
load("INLAData_S_bascamp.RData")

library(rstudioapi)
library(dplyr)
library(ggplot2)
library(leaflet)
library(viridis)

#Plotting the mean prevalence 
r_prev_mean <- as.data.frame(r_prev_mean, xy = TRUE)
colnames(r_prev_mean) <- c("long", "lat", "prev")

m.s <- ggplot(it_map, aes(x=long,y=lat))+
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
  theme(legend.position = c(0.83,0.87), legend.key.size = unit(2, 'cm'), legend.key.height = unit(1.8, 'cm'),
        legend.title = element_text(size=34), legend.text = element_text(size=32),
        plot.title = element_text(hjust = 0.5, face = "bold", size = (40)), axis.title=element_text(size=36))

png("Sheep prev mean.png", width = 2000, height = 2000,)
m.s
dev.off()

pal <- colorNumeric("viridis", c(0.69, 0.77), na.color = "transparent")

leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_mean, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_mean), title = "Relative Risk") %>%
  addScaleBar(position = c("bottomleft"))

#Plotting the lower limit of the mean (95% CI)
r_prev_ll <- rasterize(x = coop, y = r, field = prev_ll, fun = mean)

leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_ll, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_ll), title = "LL") %>%
  addScaleBar(position = c("bottomleft"))

#Plotting the upper limit of the mean (95% CI)
r_prev_ul <- rasterize(x = coop, y = r, field = prev_ul, fun = mean)

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

library(viridis)
library(cowplot)
library(ggplot2)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

plot_grid(gmean, gsd)

rm(list=ls())  #Erase workspace
library(sf)
library(sp)
library(rstudioapi)
library(dplyr)
library(ggplot2)

setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

####Prepping and loading the surveillance data and shapefile####
#Loading shapefiles
it_map1 <- st_read("./Italy Shapefile/ITA_adm2.shp")  #this map is used for masking the projections of the model
it_map <- map_data("italy") #this map is used for visualisation at the end 
sampled_map <- it_map1[c(1 : 22 ,36 : 40,57 : 63),]

#Aggregating the data with the same location
d <- group_by(basilicata_B_locations, Latitude, Longitude) %>% 
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
r <- raster(ncol=240, nrow=240, xmn=6.6, xmx=18.6, ymn=35.4, ymx=47.4)
r[] <- 1:length(r)
r <- crop(r, extent(sampled_map))
r <- mask(r, sampled_map)
pal <- colorFactor(c("darkgreen", "red"), domain = c(0, 1))
leaflet(d) %>%  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~Longitude, lat = ~Latitude, color = ~pal(prev)) %>%
  addLegend("bottomright", pal = pal, values = ~prev, title = "Farms tested positive for Ec") %>%
  addScaleBar(position = c("bottomleft"))

d$prev <- as.character(d$prev)

#Plotting the locations of samples using leaflet
sample_locations <- ggplot(mapping=aes(x=long,y=lat))+
  geom_point(data = d , aes(x = Longitude, y = Latitude, col = prev), size = 4, 
             shape = 20,)+
  scale_color_discrete(breaks = c("0", "1"))+
  scale_color_manual(values = c("0" = "green", "1" = "red"), labels = c("Negative farm", "Positive farm"), name = "Farm location")+
  geom_polygon(data=it_map,mapping=aes(group=group),
               fill=NA,color='black')+
  coord_map()+
  theme_bw()+
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Farm Locations")+
  theme(legend.position = c(0.8,0.8), legend.key.size = unit(1.4, 'cm'), legend.key.height = unit(1, 'cm'),
        legend.title = element_text(size=16), legend.text = element_text(size=15),
        plot.title = element_text(hjust = 0.5, face = "bold", size = (25)), axis.title=element_text(size=18))

png("Buffalo sample locations.png", width = 1200, height = 1200,)
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

save.image("INLAData_buffalo.RData")

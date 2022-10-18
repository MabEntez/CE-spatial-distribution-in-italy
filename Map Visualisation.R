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



#### Combining 4 maps from different models ####
library(grid)
library(gridExtra)

png("4 maps.png", width = 3500, height = 4000)
grid.arrange(
  grobs = list(m.s, m.g, m.c, m.b ),
  name = "test",
  nrow = 2,
  ncol = 2,
  layout_matrix = rbind(c(1, 2),
                        c(3, 4)
  )
)
dev.off()

png("upper and lower limit.png", width = 3400, height = 2000)
grid.arrange(
  grobs = list(mll, mul),
  name = "test",
  ncol = 2,
  layout_matrix = rbind(c(1, 2)
  )
)
dev.off()

png("Figure 2.png", width = 2900, height = 2000)
grid.arrange(
  grobs = list(m, mll, mul),
  name = "test",
  nrow = 2,
  ncol = 3,
  layout_matrix = rbind(c(1, 1, 3),
                        c(1, 1, 2)
  )
)
dev.off()

png("Figure 2.png", width = 4500, height = 2000)
par(fig=c(0, 0.5, 0, 1))
m
par(fig=c(0.5, 1, 0, 0.5), new = TRUE)
mll
par(fig=c(0,5, 1, 0.5, 1), new = TRUE)
mul
dev.off()

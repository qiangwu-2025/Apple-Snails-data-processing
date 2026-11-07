setwd("C:\\Users\\31998\\Desktop\\环境DNA\\论文写作\\福寿螺\\数据处理")

library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(sf)
library(ggspatial)
library(ggnewscale)
library(patchwork)
library(vegan)
library(ape)
library(performance)
library(car)
library(plotly)
library(terra)
library(exactextractr)
library(mgcv)
library(gratia)
library(ggrepel)

########## Fig 1A sampling point 

sample_map1<- st_read("Map\\Zhejiang_Province.shp")
target_map1 <- st_read("Map\\Wenzhou_wgs84_Area.shp", crs = 4326)
target_map2 <- st_read("Map\\Wenzhou_wgs84_Line.shp", crs = 4326)

fig1A <- ggplot() +
  geom_sf(
    data = sample_map1, 
    fill = "white",
    color = "black",
    linewidth = 0.3) +
  geom_sf_text(
    data = sample_map1, 
    aes(label = Name), 
    size = 5, 
    color = "black") +
  annotation_north_arrow(
    location = "tr",
    style = north_arrow_minimal()) +
  annotation_scale(
    location = "br",
    width_hint = 0.3) +
  theme_bw() +
  theme(
    text = element_text(size = 24),
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) 


############### Fig 2A/B COI Haplotype network plot 
### Fig 2A Pomacea canaliculata

Pomacea_canaliculata <- read.FASTA("Sequence\\Pomacea canaliculata HAP.fasta")
class(Pomacea_canaliculata)

h<-pegas::haplotype(Pomacea_canaliculata,strict=FALSE,trailingGapsAsN=TRUE)
h
hname<-paste("H",1:nrow(h),sep="")
hname
rownames(h)<-hname
h
net<-pegas::haploNet(h,d=NULL,getProb = TRUE)
net
ind.hap<-with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, individuals=names(Pomacea_canaliculata))
)
ind.hap

Fig_2A<-plot(net, size=attr(net, "freq"),
             scale.ratio = 2, 
             cex = 1, 
             labels=F, 
             pie = ind.hap, 
             show.mutation=1, 
             font=2, fast=TRUE
)



### Fig 2B Pomacea occulta

Pomacea_occulta<- read.FASTA("Sequence\\Pomacea occulta HAP.fasta")
class(Pomacea_occulta)

h<-pegas::haplotype(Pomacea_occulta,strict=FALSE,trailingGapsAsN=TRUE)
h
hname<-paste("H",1:nrow(h),sep="")
hname
rownames(h)<-hname
h
net<-pegas::haploNet(h,d=NULL,getProb = TRUE)
net
ind.hap<-with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, individuals=names(Pomacea_occulta))
)
ind.hap

Fig_2B<-plot(net, size=attr(net, "freq"),
             scale.ratio = 2, 
             cex = 1.5, 
             labels=F, 
             pie = ind.hap, 
             show.mutation=1, 
             font=2, fast=TRUE
)




########## Fig 2C Haplotype Accumulation Curves（单倍型物种累积曲线）

species_table<-read_csv("Sequence\\Haplotype table.csv") %>% 
  column_to_rownames("sample")

dat_spacc <- specaccum(species_table, method = "random", permutations = 5000) %>% 
  {data.frame(sample=.[["sites"]],richness=.[["richness"]],sd=.[["sd"]])}

Fig_2C<-ggplot(data=dat_spacc,aes(x=sample,y=richness))+
  geom_line(linewidth = 1)+
  geom_ribbon(
    aes(ymin = richness-1.96*sd, 
        ymax = richness+1.96*sd),
    alpha=0.3)+
  theme_bw()+
  theme(text=element_text(size=25),
        panel.grid = element_blank(),
        axis.title = element_text(family = "serif"))+
  xlab("Individuals")+
  ylab("Haplotype")


########## 95% detection limit

### qPCR LOD

data_LOD <- read_csv("95% detection limit\\QPCR-LOD.csv") %>% 
  select(1,2) %>% 
  rename(copy = names(.)[1],
         detection = names(.)[2]) %>% 
  filter(copy>0)

mod_lod <- glm(data=data_LOD,
               detection~log(copy), 
               family = binomial(link = "logit")
)

summary(mod_lod)

qPCR_LOD_n1 <- exp((log(0.95/(1-0.95))-coef(mod_lod)[1])/coef(mod_lod)[2])
cat("LoD_0.95 = ",qPCR_LOD_n1, "≈", ceiling(qPCR_LOD_n1),"copies/reaction")

qPCR_LOD_n3 <- exp((log( (1 - 0.05^(1/3)) / (0.05^(1/3)) ) - coef(mod_lod)[1]) / coef(mod_lod)[2])
cat("LoD_3rep = ",qPCR_LOD_n3, "≈", ceiling(qPCR_LOD_n3),"copies/reaction\n")

dat_pred <- tibble(copy = seq(0.1,25,length.out = 200)) %>% 
  cbind(
    predict(mod_lod, newdata = ., type = "link", se.fit = T) %>% 
      .[c("fit","se.fit")] %>% 
      bind_cols()
  ) %>% 
  mutate(prob = plogis(fit),
         prob_min = plogis(fit - 1.96*se.fit),
         prob_max = plogis(fit + 1.96*se.fit))

Fig.3A <- ggplot(dat_pred,
                 aes(x = copy,
                     y = prob))+
  geom_line()+
  geom_ribbon(aes(ymin = prob_min,
                  ymax = prob_max), 
              alpha = 0.2)+
  geom_jitter(data = data_LOD,
              aes(x = copy, y = detection),
              width = 0,
              height = 0.02,
              alpha = 0.25,
              size=5)+
  geom_hline(yintercept = 0.95,
             linetype=2)+
  geom_vline(aes(xintercept = qPCR_LOD_n1, color="1 replicate"), linewidth=1)+
  geom_vline(aes(xintercept = qPCR_LOD_n3, color="3 replicates"), linewidth=1)+
  scale_color_manual(name="reaction replicate",
                     values = c("1 replicate" = "red",
                                "3 replicates" = "blue")) +
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = c(0.05, 0.80),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = NA, color = NA)
  ) +
  labs(x="copies/Reaction",
       y = "Probability") +
  scale_x_log10()



### Sample LOD（copies/L）

data_Sample_LOD<-tibble(copy=seq(0.1,200,0.1))

p<-cbind(data_Sample_LOD,
         p=predict(mod_lod, 
                   newdata =data_Sample_LOD/50,
                   type="response" )) %>% 
  mutate(p1=1-(1-p)^3)

p2 <- tibble(copy = seq(0.1,200,0.1)) %>% 
  cbind(
    predict(mod_lod,
            newdata = ./50, 
            type = "link", 
            se.fit = T) %>% 
      .[c("fit","se.fit")] %>% 
      bind_cols()) %>% 
  mutate(prob = plogis(fit),
         prob_min = plogis(fit - 1.96*se.fit),
         prob_max = plogis(fit + 1.96*se.fit),
         prob2=1-(1-prob)^3,
         prob_min2=1-(1-prob_min)^3,
         prob_max2=1-(1-prob_max)^3)

Fig.3B<-ggplot(p2,aes(x = copy, y = prob2))+
  geom_line()+
  geom_ribbon(aes(ymin = prob_min2, ymax = prob_max2), alpha = 0.2)+
  geom_hline(yintercept = 0.95, linetype=2)+
  geom_vline(xintercept = 117)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 20)
  ) +
  labs(x="copies/L", 
       y = "Probability"
  )


########## Environmental DNA degradation and shedding rates

###### Environmental DNA degradation rates in different water bodies (/s)

data_tank1 <- read_csv("Tank experiment\\QPCR-mydata1.csv") %>% 
  mutate(
    tank=factor(tank,levels=c("A","B","C","D","E","F")),
    water=factor(water),
    conc = (q1 + q2 + q3)/3 * total_DNA / template / filtered_water
  )


### filed water model
## Exponential decay equation of eDNA
model1.1 <- lmer(log(conc) ~ time+ (1|tank) , data = data_tank1 %>% filter(water=="filed"))
r2(model1.1)

## Nonlinear degradation of environmental DNA
model1.2 <- lmer(log(conc) ~ time +I(time^2)+ (1|tank), data = data_tank1 %>% filter(water=="filed"))
summary(model1.2)
r2(model1.2)



### tap water model 
## Exponential decay equation of eDNA
model2.1 <- lmer(log(conc) ~ time+ (1|tank), data = data_tank1 %>% filter(water=="tap"))
r2(model2.1)

## Nonlinear degradation of environmental DNA
model2.2 <- lmer(log(conc) ~ time+I(time^2) + (1|tank), data = data_tank1 %>% filter(water=="tap"))
summary(model2.2)
r2(model2.2)


# field water K1
k1 <- fixef(model1.2)["time"] %>% abs()

# tap water K2
k2 <- fixef(model2.2)["time"] %>% abs()


dat_plot1 <- expand.grid(
  tank=data_tank1 %>% 
    filter(water=="filed") %>% 
    pull(tank) %>% unique,
  time=seq(min(data_tank1$time),
           max(data_tank1$time),1)) %>% 
  mutate(
    tank=factor(tank),
    water="filed",) %>% 
  mutate(fit=predict(model1.2,newdata=.))

dat_plot2 <- expand.grid(
  tank=data_tank1 %>% filter(water=="tap") %>%
    pull(tank) %>% unique,
  time=seq(min(data_tank1$time),
           max(data_tank1$time),1)) %>% 
  mutate(
    tank=factor(tank),
    water="tap",
  ) %>% 
  mutate(fit=predict(model2.2,newdata=.))

dat_plot <- bind_rows(dat_plot1,dat_plot2)


Fig.3C <- ggplot() +
  geom_point(
    data = data_tank1,
    aes(x = time, y = conc, col = water,group=tank),
    size=5,alpha=0.5
  ) +
  geom_line(
    data = dat_plot,
    aes(x = time, y = exp(fit), col = water,group=tank)
  ) +
  scale_y_log10() +
  labs(
    x = "Time (h)",
    y = "eDNA Concentration (copies/L)"
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    text = element_text(size = 20,family = "serif"),
    legend.position = c(0.95, 0.95), 
    legend.justification = c(1, 1),   
    legend.background = element_blank()
  )




######## Shedding rate a （copies/s）
data_tank2<-read_csv("Tank experiment\\QPCR-mydata2.csv") %>% 
  mutate(conc=(q1+q2+q3)/3*total_DNA/template/filtered_water,
         ind_L=ind/total_water,
         bio_L=bio/total_water,
         a=conc * k2/3600 * total_water)   # Use formulas: a = x * k2/3600 * V

### ind_L
model3.1<-lm(data=data_tank2,log(a)~log(ind_L))
summary(model3.1)

###bio_L
model3.2<-lm(data=data_tank2,log(a)~log(bio_L))
summary(model3.2)
β0 <- coef(model3.2)["(Intercept)"]
β1<- coef(model3.2)["log(bio_L)"]


fig3_label <- paste0(
  str_c(
    "y=",
    round(β1,2),
    "*x+",
    round(β0,2)
  ),
  "\n",
  "R² = ", sprintf("%.3f", summary(model3.2)$adj.r.squared)
)
x_pos <- 10^(min(range(log10(data_tank2$bio_L)))*1) %>% log()
y_pos <- 10^(max(log10(data_tank2$a)) *1) %>% log()

Fig.3D <- ggplot(data_tank2, aes(x = log(bio_L), y = log(a))) +
  geom_point(size = 5,
             alpha = 0.8) +  
  geom_smooth(method = "lm", 
              se = TRUE,formula = "y~x") +
  labs(
    x = "log(biomass density) (g/L)", 
    y = "log(eDNA shedding rate) (copies/s)"
  ) +
  annotate("text", 
           x = x_pos, 
           y = y_pos, 
           label = fig3_label, 
           parse = FALSE, 
           hjust = 0, 
           vjust = 1,
           size = 7, 
           lineheight = 0.9,
           family = "serif", 
           color = "black") +
  theme_minimal() + 
  theme(
    panel.border = element_rect(color = "black",
                                fill = NA,
                                linewidth = 0.5),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    text = element_text(family = "serif", size = 20),
  )


fig3 <- (Fig.3A + Fig.3B) / (Fig.3C + Fig.3D) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "",
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  )



########## Field biomass density model  BIO(g/m^3)

data_filed_eDNA<- read_csv("Filed sample\\QPCR-mydata3.csv") %>%
  mutate(conc=(q1+q2+q3)/3*total_DNA/template/filtered_water,
         conc_s=conc*k1/3600*1,
         BIO=1000*exp((log(conc_s)-β0)/β1),
         season=factor(season,levels=c("spring","summer","autumn","winter"))
  )


########## Map of average biomass density in SY and WR

###### SY map
data_sy <- data_filed_eDNA %>% 
  filter(region == "SY") %>%
  group_by(region,site,longitude,latitude) %>% 
  summarise(BIO=mean(BIO))

longitude_range_sy <- range(data_sy$longitude)
latitude_range_sy <- range(data_sy$latitude)
longitude_breaks_sy <- seq(
  from = floor(longitude_range_sy[1] * 50) / 50,
  to = ceiling(longitude_range_sy[2] * 50) / 50,
  by = 0.04
)
latitude_breaks_sy <- seq(
  from = floor(latitude_range_sy[1] * 50) / 50,
  to   = ceiling(latitude_range_sy[2] * 50) / 50,
  by = 0.02   
)


fig1B <- ggplot() +
  geom_sf(data = target_map1, fill = "skyblue", alpha = 0.4) +
  geom_sf(data = target_map2, fill = "skyblue", alpha = 0.4) +
  geom_point(data = data_sy, 
             aes(longitude, latitude, 
                 col = log(BIO)),
             size = 5) +
  geom_text(data = data_sy, 
            aes(x = longitude,
                y = latitude,
                label = site ),
            nudge_y=-0.002,
            size = 6, color="black")+
  scale_color_gradient(low = "green",
                       high = "red", 
                       name = "log(BIO)") +
  scale_x_continuous(breaks = longitude_breaks_sy) +
  scale_y_continuous(breaks = latitude_breaks_sy) +
  coord_sf(xlim = range(data_sy$longitude),
           ylim = range(data_sy$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br", 
                   width_hint = 0.3) +
  theme_bw() +
  theme( 
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )

###### WR map
data_wr <-data_filed_eDNA  %>%
  filter(region == "WR") %>%
  group_by(region, site, longitude, latitude) %>% 
  summarise(BIO = mean(BIO))

longitude_range_wr <- range(data_wr$longitude)
latitude_range_wr  <- range(data_wr$latitude)

longitude_breaks_wr <- seq(
  from = floor(longitude_range_wr[1] * 50) / 50,
  to   = ceiling(longitude_range_wr[2] * 50) / 50,
  by   = 0.04
)

latitude_breaks_wr <- seq(
  from = floor(latitude_range_wr[1] * 50) / 50,
  to   = ceiling(latitude_range_wr[2] * 50) / 50,
  by   = 0.04
)

fig1C <- ggplot() +
  geom_sf(data = target_map1, fill = "skyblue", alpha = 0.4) +
  geom_sf(data = target_map2, fill = "skyblue", alpha = 0.4) +
  geom_point(data = data_wr, 
             aes(longitude, latitude, 
                 col = log(BIO)),
             size = 5) +
  geom_text(data = data_wr, 
            aes(x = longitude,
                y = latitude,
                label = site),
            nudge_y = -0.003,
            size = 6, color = "black") +
  scale_color_gradient(low = "green",
                       high = "red", 
                       name = "log(BIO)") +
  scale_x_continuous(breaks = longitude_breaks_wr) +
  scale_y_continuous(breaks = latitude_breaks_wr) +
  coord_sf(xlim = range(data_wr$longitude),
           ylim = range(data_wr$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br", 
                   width_hint = 0.3) +
  theme_bw() +
  theme( 
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )

Fig1BC <- fig1B / fig1C 




########## Seasonal biomass variations at SY and WR

seasons <- c("spring", "summer", "autumn", "winter")
regions <- c("SY", "WR")

nudge_map  <- c(SY = -0.002, WR = -0.003)
lat_by_map <- c(SY = 0.02,  WR = 0.04)

season_plots <- list()
for (reg in regions) {
  for (seas in seasons) {
    dat <-data_filed_eDNA %>%
      filter(region == reg, season == seas) %>%
      group_by(region, site, BIO, longitude, latitude) %>%
      summarise(.groups = "drop")
    dat_pos <- dat %>% filter(BIO > 0)
    dat_na  <- dat %>% filter(BIO == 0)
    lon_range <- range(dat$longitude)
    lat_range <- range(dat$latitude)
    
    lon_breaks <- seq(
      from = floor(lon_range[1] * 50) / 50,
      to   = ceiling(lon_range[2] * 50) / 50,
      by   = 0.04
    )
    
    lat_breaks <- seq(
      from = floor(lat_range[1] * 50) / 50,
      to   = ceiling(lat_range[2] * 50) / 50,
      by   = lat_by_map[[reg]]
    )
    p <- ggplot() +
      geom_sf(data = target_map1, fill = "skyblue", alpha = 0.4) +
      geom_sf(data = target_map2, fill = "skyblue", alpha = 0.4) +
      geom_point(data = dat_pos,
                 aes(longitude, latitude, col = log(BIO)),
                 size = 4) +
      scale_color_gradient(low = "green",
                           high = "red",
                           name = "log(BIO)") +
      new_scale_color() +
      geom_point(data = dat_na,
                 aes(longitude, latitude, col = "NA"),
                 size = 4) +
      scale_color_manual(name = NULL,
                         values = c("NA" = "grey50"),
                         labels = c("NA")) +
      geom_text(data = dat,
                aes(x = longitude, y = latitude, label = site),
                nudge_y = nudge_map[[reg]],
                size = 4, color = "black") +
      
      scale_x_continuous(breaks = lon_breaks) +
      scale_y_continuous(breaks = lat_breaks) +
      coord_sf(xlim = range(dat$longitude),
               ylim = range(dat$latitude),
               expand = TRUE) +
      annotation_north_arrow(location = "tr",
                             style = north_arrow_minimal()) +
      annotation_scale(location = "br",
                       width_hint = 0.3) +
      theme_bw() +
      theme(
        text = element_text(size = 20),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.position = "right",        
        legend.box = "vertical",
        legend.spacing.y = unit(0.2, "cm")
      )
    
    season_plots[[paste(reg, seas, sep = "_")]] <- p
  }
}

Figs3 <- (
  (season_plots[["SY_spring"]] | season_plots[["SY_summer"]]) /
    (season_plots[["SY_autumn"]] | season_plots[["SY_winter"]]) /
    (season_plots[["WR_spring"]] | season_plots[["WR_summer"]]) /
    (season_plots[["WR_autumn"]] | season_plots[["WR_winter"]])
) +
  plot_annotation (
    tag_levels = "A",
    tag_suffix = "",
    theme = theme (
      plot.tag = element_text (size = 15, face = "bold")
    )
  )



########## Water surface area (km²)

###SY

sy_bbox <- st_bbox(c(
  xmin = longitude_range_sy[1], 
  ymin = latitude_range_sy[1], 
  xmax = longitude_range_sy[2], 
  ymax = latitude_range_sy[2]
), crs = st_crs(target_map1))

sy_polygon <- st_as_sfc(sy_bbox)

water_in_sy1 <- st_intersection(target_map1, sy_polygon)
water_in_sy2 <- st_intersection(target_map2, sy_polygon)
water_in_sy <- bind_rows(water_in_sy1, water_in_sy2)

if(nrow(water_in_sy) > 0) {
  water_area_sy <- sum(st_area(water_in_sy))
  water_area_sy_km2 <- as.numeric(water_area_sy) / 1e6
} else {
  water_area_sy_km2 <- 0
}

cat("SY area:", round(water_area_sy_km2, 4), "km²\n")

SY_mean_BIO<- mean(data_sy$BIO)
SY_water_volume <- water_area_sy_km2 * 3* 1e6
SY_total_BIO <- SY_mean_BIO * SY_water_volume

datasy <- data_filed_eDNA  %>%
  filter(region == "SY")
hist(datasy$BIO %>% log())
datasy$BIO %>% mean
replicate(5000,mean(sample(datasy$BIO,replace = T))) %>% quantile(c(0.025,0.5,0.975))


###WR

wr_bbox <- st_bbox(c(
  xmin = longitude_range_wr[1], 
  ymin = latitude_range_wr[1], 
  xmax = longitude_range_wr[2], 
  ymax = latitude_range_wr[2]
), crs = st_crs(target_map1))

wr_polygon <- st_as_sfc(wr_bbox)
water_in_wr1 <- st_intersection(target_map1, wr_polygon)
water_in_wr2 <- st_intersection(target_map2, wr_polygon)
water_in_wr <- bind_rows(water_in_wr1, water_in_wr2)

if(nrow(water_in_wr) > 0) {
  water_area_wr <- sum(st_area(water_in_wr))
  water_area_wr_km2 <- as.numeric(water_area_wr) / 1e6
} else {
  water_area_wr_km2 <- 0
}

cat("WR area:", round(water_area_wr_km2, 4), "km²\n")

WR_mean_BIO<- mean(data_wr$BIO)
WR_water_volume <- water_area_wr_km2 * 3* 1e6
WR_total_BIO <- WR_mean_BIO * WR_water_volume

datawr<- data_filed_eDNA  %>%
  filter(region == "WR")
hist(datawr$BIO %>% log())
datawr$BIO %>% mean
replicate(5000,mean(sample(datawr$BIO,replace = T))) %>% quantile(c(0.025,0.5,0.975))


########## Cultivated land data processing
# Read and cut the grid data of the study area

data_grid <- rast("Map\\CLCD_v01_2025_.tif")
study_grid <- as.polygons(ext(120, 121, 27, 29), crs = "EPSG:4326")  %>% 
  project(crs(data_grid))
grid_utm <- crop(data_grid, study_grid) %>% 
  project("EPSG:32651", method = "near")

# calculation of buffers
buffers <-data_filed_eDNA  %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(32651)  %>% 
  st_buffer(500)

# Extracting the area of each class in the buffer zone
clcd_labels <- c(
  "1" = "Cropland",
  "2" = "Forest", 
  "3" = "Shrub",
  "4" = "Grassland",
  "5" = "Water", 
  "6" = "Snow/Ice",
  "7" = "Barren", 
  "8" = "Impervious",
  "9" = "Wetland"
)

area_long <- exact_extract(grid_utm, buffers, coverage_area = TRUE) %>%
  bind_rows(.id = "buffer_id") %>%
  mutate(
    sample_id = mydata5$sample_id[as.numeric(buffer_id)],
    land_type = factor(value, levels = 1:9, labels = clcd_labels)
  ) %>%
  group_by(sample_id, land_type) %>%
  summarise(area = sum(coverage_area), .groups = "drop")

# Cropland area extracted
cropland_area <- area_long  %>% 
  filter(land_type == "Cropland")  %>% 
  select(sample_id, cropland_area = area) %>% 
  mutate(cropland_area = cropland_area / 1e6) #Km²

data_filed_eDNA <-data_filed_eDNA  %>% 
  left_join(cropland_area, by = "sample_id")  %>% 
  mutate(cropland_area = replace_na(cropland_area, 0))


########## Distribution model of Pomacea
model4 <- glm(
  data = data_filed_eDNA %>%
    mutate(BIO=if_else(BIO>0,1,0)),
  BIO~log(ph)+log(do)+log(tur)+log(tm)+log(chla)+cropland_area, 
  family = binomial(link="logit")
)
summary(model4)
vif(model4)

plot_dat1 <- expand.grid(
  ph = mean(data_filed_eDNA$ph, na.rm = TRUE),
  do = mean(data_filed_eDNA$do, na.rm = TRUE),
  tur = mean(data_filed_eDNA$tur, na.rm = TRUE),
  chla = mean(data_filed_eDNA$chla, na.rm = TRUE),
  tm = seq(
    min(data_filed_eDNA$tm[data_filed_eDNA$tm > 0], na.rm = TRUE),
    max(data_filed_eDNA$tm, na.rm = TRUE),
    length.out = 100
  ),
  cropland_area = seq(
    min(data_filed_eDNA$cropland_area, na.rm = TRUE),
    max(data_filed_eDNA$cropland_area, na.rm = TRUE),
    length.out = 100
  )
)

pred <- predict(
  model4,
  newdata = plot_dat1,
  type = "link",
  se.fit = TRUE
)
plot_dat1$probability <- plogis(pred$fit)

z_matrix <- matrix(
  plot_dat1$probability,
  nrow = length(unique(plot_dat1$tm)),
  ncol = length(unique(plot_dat1$cropland_area)),
  byrow = FALSE
)

Fig4 <- plotly::plot_ly(
  x = unique(plot_dat1$tm),
  y = unique(plot_dat1$cropland_area),
  z = z_matrix,
  type = "surface",
  colorscale = "Viridis",
  colorbar = list(
    title = list(
      text = "Probability",
      font = list(size = 15)
    ),
    tickfont = list(size = 12)
  ),
  contours = list(
    x = list(
      show = TRUE,
      color = "black",      
      width = 0.5,             
      start = min(unique(plot_dat1$tm)),
      end = max(unique(plot_dat1$tm)),
      size = diff(range(unique(plot_dat1$tm))) / 25,  
      highlight = FALSE      
    ),
    y = list(
      show = TRUE,
      color = "black",      
      width = 0.5,
      start = min(unique(plot_dat1$cropland_area)),
      end = max(unique(plot_dat1$cropland_area)),
      size = diff(range(unique(plot_dat1$cropland_area))) / 25,
      highlight = FALSE
    )
  )
) %>%
  plotly::layout(
    scene = list(
      xaxis = list(
        title = list(text = "Water Temperature(°C)", font = list(size = 16)),
        tickfont = list(size = 13)
      ),
      yaxis = list(
        title = list(text = "Cropland area(km²)", font = list(size = 16)),
        tickfont = list(size = 13)
      ),
      zaxis = list(
        title = list(text = "Probability", font = list(size = 16)),
        tickfont = list(size = 13),
        range = c(0, 1)
      ),
      camera = list(
        eye = list(x = -1.25, y = -1.25, z = 1.25),
        center = list(x = 0, y = 0, z = -0.2)
      )
    )
  )



########## Effect of biomass density of apple snail on plankton
# Calculation of α diversity index of plankton
# phytoplankton
data_plant<-read_csv("plankton\\data_plant.csv") %>% 
  column_to_rownames("sample")*1e6 #density（×10^6cells/L）
ai_plant = abs(scale(data_plant))
avd_plant = rowSums(ai_plant) / ncol(data_plant)  # average variation degree，AVD
index_plant <- data.frame(
  Abundance_plant = rowSums(data_plant, na.rm = TRUE), #abundance
  S_plant = specnumber(data_plant), #richness
  H_plant = vegan::diversity(data_plant, index = "shannon") ) %>% # shannon-wiener
  mutate(
    J_plant = H_plant / log(S_plant),
    avd_plant = avd_plant) %>%
  rownames_to_column("sample")

# zooplankton
data_zoo<-read_csv("plankton\\data_zoo.csv") %>% 
  column_to_rownames("sample")
data_spe_inf<-read_csv("plankton\\data_species information.csv")
ai_zoo <- abs(scale(data_zoo))
avd_zoo <- rowSums(ai_zoo) / ncol(data_zoo)
index_zoo <- data.frame(
  Abundance_zoo = rowSums(data_zoo, na.rm = TRUE),
  S_zoo = specnumber(data_zoo),
  H_zoo = vegan::diversity(data_zoo, index = "shannon")
) %>% 
  mutate(
    J_zoo = H_zoo / log(S_zoo),
    avd_zoo = avd_zoo
  )%>%
  rownames_to_column("sample")

data_filed_eDNA <-data_filed_eDNA %>%
  left_join(index_plant, by = c("sample_id" = "sample")) %>%
  left_join(index_zoo,   by = c("sample_id" = "sample")) %>%
  filter(tolower(season) != "winter") %>%  # Preserve sample data for spring, summer and autumn
  mutate(site=as.factor(site))


########## Effect of biomass density of apple snail on phytoplankton
# Establish response variables
plant_response <- c(
  "log(Abundance_plant)", #Total abundance of phytoplankton
  "log(S_plant)", #Species richness of phytoplankton
  "H_plant", #Shannon index of phytoplankton
  "J_plant", #Pielou evenness index of phytoplankton
  "avd_plant" #Average variation degree of phytoplankton
)
plant_model_names <- paste0("gam", 1:5)

# smooth function
plant_rhs <- "~ s(BIO) + te(longitude, latitude) + s(tm) + s(ph) + s(site, bs = 're') + s(tur)"

# model fitting (all data)
plant_gam_list <- list()
for (i in seq_along(plant_response)) {
  formula_full <- as.formula(paste(plant_response[i], plant_rhs))
  plant_gam_list[[i]] <- gam(
    formula_full, 
    data = data_filed_eDNA, 
    method = "REML", 
    select = TRUE
  )
  names(plant_gam_list)[i] <- plant_model_names[i]
  
  cat("\n", "==========", plant_model_names[i], "|", plant_response[i], "==========\n")
  print(summary(plant_gam_list[[i]]))   # Model summary
  gam.check(plant_gam_list[[i]])        # model diagnosis
  # plot(plant_gam_list[[i]], residuals = TRUE, main = plant_model_names[i]) 
}


########## Effect of biomass density of apple snail on zooplankton
# Establish response variables
zoo_response <- c(
  "log(Abundance_zoo)", #Total abundance of zooplankton
  "log(S_zoo)",         #Species richness of zooplankton
  "H_zoo",              #Shannon index of zooplankton
  "J_zoo",              #Pielou index of zooplankton
  "avd_zoo"             #Average variation degree of zooplankton
)
zoo_model_names <- paste0("gam", 6:10)

# smooth function
zoo_rhs <- "~ s(BIO) + s(log(Abundance_plant)) + te(longitude, latitude) + s(tm) + s(ph) + s(site, bs = 're') + s(tur)"

# model fitting (all data)
zoo_gam_list <- list()
for (i in seq_along(zoo_response)) {
  formula_full <- as.formula(paste(zoo_response[i], zoo_rhs))
  zoo_gam_list[[i]] <- gam(
    formula_full, 
    data = data_filed_eDNA,
    method = "REML", 
    select = TRUE
  )
  names(zoo_gam_list)[i] <- zoo_model_names[i]
  
  cat("\n", "==========", zoo_model_names[i], "|", zoo_response[i], "==========\n")
  print(summary(zoo_gam_list[[i]]))   # Model summary
  gam.check(zoo_gam_list[[i]])        # model diagnosis
}


######### Partial effect plot
p_to_label <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  }
  if (p < 0.001) {
    "p < 0.001***"
  } else if (p < 0.01) {
    "p < 0.01**"
  } else if (p < 0.05) {
    "p < 0.05*"
  } else {
    sprintf("p = %.3f", p)
  }
}

# plant
plant_labels <- c("Abundance", "S", "H", "J", "AVD")
plant_terms <- list(
  gam1 = c("BIO", "tm", "ph", "tur"),
  gam2 = c("BIO", "tm", "ph", "tur"),
  gam3 = c("BIO", "tm", "ph", "tur"),
  gam4 = c("BIO", "tm", "ph", "tur"),
  gam5 = c("BIO", "tm", "ph", "tur")
)
plant_variable_order <- c("BIO", "tm", "ph", "tur")
plant_variable_labels <- c(
  BIO = "BIO",
  tm  = "WT",
  ph  = "pH",
  tur = "Turbidity"
)

plant_extract_smooth_data <- function(model, model_name, response_label, variables) {
  # smooth effect
  smooth_data <- map_dfr(variables, function(var) {
    smooth_term <- paste0("s(", var, ")")
    est <- smooth_estimates(model, select = smooth_term)
    tibble(
      response = response_label,
      model    = model_name,
      variable = var,
      x        = est[[var]],
      estimate = est$.estimate,
      se       = est$.se
    )
  })
  
  # partial residuals
  resid_all <- add_partial_residuals(model, data = model$model, type = "working")
  residual_data <- map_dfr(variables, function(var) {
    smooth_term <- paste0("s(", var, ")")
    tibble(
      response      = response_label,
      model         = model_name,
      variable      = var,
      x             = resid_all[[var]],
      partial_resid = resid_all[[smooth_term]]
    ) %>%
      filter(!is.na(x), !is.na(partial_resid))
  })
  
  smooth_table <- summary(model)$s.table
  p_data <- map_dfr(variables, function(var) {
    smooth_term <- paste0("s(", var, ")")
    if (smooth_term %in% rownames(smooth_table)) {
      p_value <- smooth_table[smooth_term, "p-value"]
    } else {
      p_value <- NA_real_
    }
    tibble(
      response = response_label,
      model    = model_name,
      variable = var,
      p_value  = p_value,
      p_label  = p_to_label(p_value)
    )
  })
  
  list(smooth = smooth_data, residual = residual_data, p = p_data)
}

plot_plant_partial <- function(gam_list) {
  plant_effect_list <- map2(gam_list, plant_model_names, function(model, model_name) {
    response_label <- plant_labels[match(model_name, plant_model_names)]
    plant_extract_smooth_data(
      model          = model,
      model_name     = model_name,
      response_label = response_label,
      variables      = plant_terms[[model_name]]
    )
  })
  
  plant_smooth_data   <- map_dfr(plant_effect_list, "smooth")
  plant_residual_data <- map_dfr(plant_effect_list, "residual")
  plant_p_data        <- map_dfr(plant_effect_list, "p")
  
  plant_smooth_data <- plant_smooth_data %>% mutate(
    response = factor(response, levels = plant_labels),
    variable = factor(variable, levels = plant_variable_order, labels = plant_variable_labels)
  )
  plant_residual_data <- plant_residual_data %>% mutate(
    response = factor(response, levels = plant_labels),
    variable = factor(variable, levels = plant_variable_order, labels = plant_variable_labels)
  )
  plant_p_data <- plant_p_data %>% mutate(
    response = factor(response, levels = plant_labels),
    variable = factor(variable, levels = plant_variable_order, labels = plant_variable_labels)
  )
  fig_plant <- ggplot() +
    geom_ribbon(
      data = plant_smooth_data,
      aes(x = x, ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
      fill = "grey50", alpha = 0.20
    ) +
    geom_line(
      data = plant_smooth_data,
      aes(x = x, y = estimate),
      colour = "black", linewidth = 0.8
    ) +
    geom_point(
      data = plant_residual_data,
      aes(x = x, y = partial_resid),
      colour = "#2c7fb8", alpha = 0.25, size = 1, shape = 16
    ) +
    geom_hline(
      yintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey40"
    ) +
    geom_text(
      data = plant_p_data,
      aes(x = Inf, y = Inf, label = p_label),
      hjust = 1.05, vjust = 1.20, 
      size = 5,               
      family = "serif",
      fontface = "italic"  
    ) +
    facet_grid(
      response ~ variable,
      scales = "free",
      drop = FALSE,
      labeller = labeller(plant_variable_labels)
    ) +
    labs(x = NULL, y = "Partial effect") +
    theme_bw() +
    theme(
      text = element_text(family = "serif", size = 18),
      strip.text.y = element_text( face = "italic")
    )
  
  return(fig_plant)
}

fig_plant_effect <- plot_plant_partial(plant_gam_list)
fig_plant_effect


# zoo
zoo_labels <- c("Abundance", "S", "H", "J", "AVD")
zoo_terms <- list(
  gam6  = c("BIO", "log(Abundance_plant)", "tm", "ph", "tur"),
  gam7  = c("BIO", "log(Abundance_plant)", "tm", "ph", "tur"),
  gam8  = c("BIO", "log(Abundance_plant)", "tm", "ph", "tur"),
  gam9  = c("BIO", "log(Abundance_plant)", "tm", "ph", "tur"),
  gam10 = c("BIO", "log(Abundance_plant)", "tm", "ph", "tur")
)
zoo_variable_order <- c("BIO", "log(Abundance_plant)", "tm", "ph", "tur")
zoo_variable_labels <- c(
  BIO                    = "BIO",
  'log(Abundance_plant)' = "PA",
  tm                     = "WT",
  ph                     = "pH",
  tur                    = "Turbidity"
)

zoo_extract_smooth_data <- function(model, model_name, response_label, variables) {
  smooth_data <- map_dfr(variables, function(var) {
    smooth_term <- paste0("s(", var, ")")
    est <- smooth_estimates(model, select = smooth_term)
    tibble(
      response = response_label,
      model    = model_name,
      variable = var,
      x        = est[[var]],
      estimate = est$.estimate,
      se       = est$.se
    )
  })
  
  resid_all <- add_partial_residuals(model, data = model$model, type = "working")
  residual_data <- map_dfr(variables, function(var) {
    smooth_term <- paste0("s(", var, ")")
    tibble(
      response      = response_label,
      model         = model_name,
      variable      = var,
      x             = resid_all[[var]],
      partial_resid = resid_all[[smooth_term]]
    ) %>%
      filter(!is.na(x), !is.na(partial_resid))
  })
  
  smooth_table <- summary(model)$s.table
  p_data <- map_dfr(variables, function(var) {
    smooth_term <- paste0("s(", var, ")")
    if (smooth_term %in% rownames(smooth_table)) {
      p_value <- smooth_table[smooth_term, "p-value"]
    } else {
      p_value <- NA_real_
    }
    tibble(
      response = response_label,
      model    = model_name,
      variable = var,
      p_value  = p_value,
      p_label  = p_to_label(p_value)
    )
  })
  
  list(smooth = smooth_data, residual = residual_data, p = p_data)
}

plot_zoo_partial <- function(gam_list) {
  zoo_effect_list <- map2(gam_list, zoo_model_names, function(model, model_name) {
    response_label <- zoo_labels[match(model_name, zoo_model_names)]
    zoo_extract_smooth_data(
      model          = model,
      model_name     = model_name,
      response_label = response_label,
      variables      = zoo_terms[[model_name]]
    )
  })
  
  zoo_smooth_data   <- map_dfr(zoo_effect_list, "smooth")
  zoo_residual_data <- map_dfr(zoo_effect_list, "residual")
  zoo_p_data        <- map_dfr(zoo_effect_list, "p")
  
  zoo_smooth_data <- zoo_smooth_data %>% mutate(
    response = factor(response, levels = zoo_labels),
    variable = factor(variable, levels = zoo_variable_order, labels = zoo_variable_labels)
  )
  zoo_residual_data <- zoo_residual_data %>% mutate(
    response = factor(response, levels = zoo_labels),
    variable = factor(variable, levels = zoo_variable_order, labels = zoo_variable_labels)
  )
  zoo_p_data <- zoo_p_data %>% mutate(
    response = factor(response, levels = zoo_labels),
    variable = factor(variable, levels = zoo_variable_order, labels = zoo_variable_labels)
  )
  fig_zoo <- ggplot() +
    geom_ribbon(
      data = zoo_smooth_data,
      aes(x = x, ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
      fill = "grey50", alpha = 0.20
    ) +
    geom_line(
      data = zoo_smooth_data,
      aes(x = x, y = estimate),
      colour = "black", linewidth = 0.8
    ) +
    geom_point(
      data = zoo_residual_data,
      aes(x = x, y = partial_resid),
      colour = "#2c7fb8", alpha = 0.25, size = 1, shape = 16
    ) +
    geom_hline(
      yintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey40"
    ) +
    geom_text(
      data = zoo_p_data,
      aes(x = Inf, y = Inf, label = p_label),
      hjust = 1.05, vjust = 1.20, 
      size = 5,               
      family = "serif",
      fontface = "italic"  
    ) +
    facet_grid(
      response ~ variable,
      scales = "free",
      drop = FALSE,
      labeller = labeller(zoo_variable_labels)
    ) +
    labs(x = NULL, y = "Partial effect") +
    theme_bw() +
    theme(
      text = element_text(family = "serif", size = 18),
      strip.text.y = element_text( face = "italic")
    )
  return(fig_zoo)
}

fig_zoo_effect <- plot_zoo_partial(zoo_gam_list)
fig_zoo_effect


######## Model diagnosis drawing
extract_diag_data <- function(model_list, response_names, bio_group_label) {
  map2_dfr(model_list, response_names,
           function(model, resp_name) {
             tibble(
               resi = residuals(model, type = "deviance"),
               ftvl = fitted(model),
               resp = model$y,
               response = resp_name,
               bio_group = bio_group_label
             )
           })
}

plant_resi_raw <- extract_diag_data(plant_gam_list, plant_response, bio_group_label = "all_BIO")
zoo_resi_raw   <- extract_diag_data(zoo_gam_list,   zoo_response,   bio_group_label = "all_BIO")

plant_order <- c(
  "log(Abundance_plant)",
  "log(S_plant)",
  "H_plant",
  "J_plant",
  "avd_plant"
)
zoo_order <- c(
  "log(Abundance_zoo)",
  "log(S_zoo)",
  "H_zoo",
  "J_zoo",
  "avd_zoo"
)

plant_resi <- plant_resi_raw %>%
  mutate(response = factor(response, levels = plant_order, labels = plant_labels))
zoo_resi <- zoo_resi_raw %>%
  mutate(response = factor(response, levels = zoo_order, labels = zoo_labels))

plot_diagnostics <- function(resi_df) {
  # QQ plot
  p_qq <- ggplot(resi_df, aes(sample = resi)) +
    facet_wrap(~ response, scales = "free", ncol = 1) +
    geom_qq(shape = 1, alpha = 0.25, color = "#2c7fb8") +
    geom_qq_line(color = "black") +
    theme_bw() +
    theme(
      text = element_text(size = 15, family = "serif"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = "none"
    ) +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
  
  # residual histogram
  p_hist <- ggplot(resi_df, aes(x = resi)) +
    facet_wrap(~ response, scales = "free", ncol = 1) +
    geom_histogram(bins = 12, col = "black", fill = "#2c7fb8", alpha = 0.4) +
    theme_bw() +
    theme(
      text = element_text(size = 15, family = "serif"),
      strip.background = element_blank(),
      strip.text = element_blank(),
    ) +
    labs(x = "Residuals", y = "Frequency")
  
  # residuals vs fitted values
  p_res_fit <- ggplot(resi_df, aes(x = ftvl, y = resi)) +
    facet_wrap(~ response, scales = "free", ncol = 1) +
    geom_point(alpha = 0.25, color = "#2c7fb8") +
    geom_hline(yintercept = 0, linetype = 2) +
    theme_bw() +
    theme(
      text = element_text(size = 15, family = "serif"),
      strip.background = element_blank(),
      strip.text = element_blank(),
    ) +
    labs(x = "Fitted values", y = "Residuals")
  
  # response vs fitted values
  p_resp_fit <- ggplot(resi_df, aes(x = ftvl, y = resp)) +
    facet_wrap(
      ~ response,
      scales = "free",
      ncol = 1,
      strip.position = "right"
    ) +
    geom_point(alpha = 0.25, color = "#2c7fb8") +
    geom_abline(slope = 1, linetype = 2) +
    coord_cartesian(clip = "on") +
    theme_bw() +
    theme(
      text = element_text(size = 15, family = "serif"),
      strip.background = element_rect(
        fill = "grey85",
        colour = "grey50",
        linewidth = 0.5
      ),
      strip.text.y = element_text(
        face = "italic",
        size = 12,
        angle = -90),
      panel.spacing = unit(0.8, "lines")
    ) +
    labs(x = "Fitted values", y = "Response")
  
  combined <- p_qq + p_hist + p_res_fit + p_resp_fit +
    plot_layout(guides = "collect", nrow = 1)
  return(combined)
}

fig_plant_diag <- plot_diagnostics(plant_resi)
fig_zoo_diag   <- plot_diagnostics(zoo_resi)



########### Community impact analysis
#  Bray-Curtis distance matrix 
bray_dist_plant <- data_plant %>% vegdist("bray")  # phytoplankton 
bray_dist_zoo   <- data_zoo   %>% vegdist("bray")  # zooplankton 

# 2. VIF collinearity test 
# phytoplankton
env_plant <- data_filed_eDNA[match(labels(bray_dist_plant), data_filed_eDNA$sample_id), ]
lm_plant  <- lm(BIO ~ tm + ph + tur + longitude + latitude, data = env_plant)
vif_plant <- car::vif(lm_plant)
vif_plant

# zooplankton
env_zoo <- data_filed_eDNA[match(labels(bray_dist_zoo), data_filed_eDNA$sample_id), ]
lm_zoo  <- lm(BIO ~ tm + ph + tur + longitude + latitude + Abundance_plant, data = env_zoo)
vif_zoo <- car::vif(lm_zoo)
vif_zoo

# PERMANOVA
# phytoplankton
per_plant <- adonis2(
  bray_dist_plant ~ tm + ph + tur + longitude * latitude + BIO,
  data = data_filed_eDNA[match(labels(bray_dist_plant), data_filed_eDNA$sample_id), ],
  permutations = 999, by = "margin", add = "cailliez"
)
per_plant

# zooplankton
per_zoo <- adonis2(
  bray_dist_zoo ~ tm + ph + tur + longitude * latitude + BIO + log(Abundance_plant),
  data = data_filed_eDNA[match(labels(bray_dist_zoo), data_filed_eDNA$sample_id), ],
  permutations = 999, by = "margin", add = "cailliez"
)
per_zoo

# Community dominance calculation 
mf_cal_domi <- function(x){
  {colSums(x)/sum(x)*colSums(x>0)/nrow(x)} %>% 
    enframe(name = "species", value = "domi")
}

# phytoplankton dominance
domi_plant <- data_plant %>%
  mf_cal_domi() %>%
  left_join(data_spe_inf, by = join_by(species == Latin_name)) %>%
  mutate(species = str_replace(species, "_", " ")) %>%
  arrange(-domi) 
write_csv(domi_plant, "domi_plant.csv")


#  Phytoplankton PCoA + BIO ordisurf plot 
#  PCoA
pcoa_plant <- wcmdscale(bray_dist_plant, eig = TRUE, add = "cailliez")

# sample scores
score_sample_plant <- pcoa_plant$point %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  select(sample_id, Dim1, Dim2) %>%
  left_join(data_filed_eDNA, by = "sample_id")

# top 10 dominant species scores
score_spe_plant <- wascores(pcoa_plant$points[,1:2],
                            data_plant, expand = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  mutate(species = str_replace_all(species, "_", " ")) %>%
  filter(species %in% domi_plant$species[1:10])

#  PCoA axis explanation
eig_plant     <- pcoa_plant$eig
var_per_plant <- round(eig_plant / sum(eig_plant) * 100, 2)
xlab_plant    <- paste0("PCoA 1 (", var_per_plant[1], "%)")
ylab_plant    <- paste0("PCoA 2 (", var_per_plant[2], "%)")

# BIO GAM contour (ordisurf)
gam_plant      <- ordisurf(pcoa_plant ~ BIO, data = score_sample_plant,
                           plot = FALSE, knots = 10)
gam_grid_plant <- gam_plant$grid
gam_data_plant <- expand.grid(x = gam_grid_plant$x, y = gam_grid_plant$y)
gam_data_plant$z <- as.vector(gam_grid_plant$z)
gam_data_plant <- na.omit(gam_data_plant)


Fig7 <- ggplot() +
  metR::geom_contour2(
    data = gam_data_plant,
    aes(x = x, y = y, z = z, color = after_stat(level)),
    linewidth = 0.8
  ) +
  metR::geom_text_contour(
    data = gam_data_plant,
    aes(x = x, y = y, z = z),
    stroke = 0.2, min.size = 5, rotate = FALSE, check_overlap = TRUE,
    color = "black",
    family = "serif"
  ) +
  geom_segment(
    data = score_spe_plant,
    aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
    linewidth = 0.4, arrow = arrow(length = unit(0.015, "npc")),
    alpha = 0.7, color = "black"
  ) +
  geom_text_repel(
    data = score_spe_plant,
    aes(x = Dim1, y = Dim2, label = species),
    color = "black", size = 4, family = "serif", fontface = "italic"
  ) +
  geom_point(
    data = score_sample_plant,
    aes(x = Dim1, y = Dim2),
    shape = 21, size = 2, fill = "grey40", color = "white", stroke = 0.5
  ) +
  scale_color_distiller(name = "BIO\n(g/m³)", palette = "Spectral") +
  labs(x = xlab_plant, y = ylab_plant) +
  theme_bw() +
  theme(
    axis.title   = element_text(size = 15, family = "serif"),
    axis.text    = element_text(size = 15, family = "serif"),
    legend.title = element_text(size = 15, family = "serif"),
    legend.text  = element_text(size = 15, family = "serif"),
    panel.grid   = element_blank()
  ) +
  guides(color = guide_colorbar(barwidth = 0.8, barheight = 10))






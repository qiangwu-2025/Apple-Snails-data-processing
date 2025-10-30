
setwd("path")
fig_type <- "pdf"

library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(sf)
library(ggspatial)
library(patchwork)
library(vegan)
library(ape)
library(performance)


########## Fig 1A sampling point 

sample_map1<- st_read("zhejiang.shp")
target_map1 <- st_read("zenzhou_area.shp", crs = 4326)
target_map2 <- st_read("wenzhou_line.shp", crs = 4326)

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

 fig1A


############### Fig 2A/B COI Haplotype network plot 
### Fig 2A Pomacea canaliculata

Pomacea_canaliculata <- read.FASTA("Pomacea canaliculata HAP1.fasta")
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

Pomacea_occulta<- read.FASTA("Pomacea occulta HAP2.fasta")
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




########## Fig 2C Haplotype Accumulation Curves

data_HAC<-species_table<-read_csv("hap.csv") %>% 
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

Fig_2C





########## Fig 3  95% detection limit

### Fig 3A  qPCR LOD

data_qPCR <- read_csv("QPCR-LOD.csv") %>% 
  select(1,2) %>% 
  rename(copy = names(.)[1],
         detection = names(.)[2]) %>% 
  filter(copy>0)

mod_lod <- glm(data=data_qPCR,
               detection~log(copy), 
               family = binomial(link = "logit")
)

summary(mod_lod)

qPCR_LOD <- exp((log(0.95/(1-0.95))-coef(mod_lod)[1])/coef(mod_lod)[2])
cat("LoD_0.95 = ",qPCR_LOD, "≈", ceiling(qPCR_LOD),"copies/reaction")

dat_pred <- tibble(copy = seq(0.1,25,length.out = 200)) %>% 
  cbind(
    predict(mod_lod, newdata = ., type = "link", se.fit = T) %>% 
      .[c("fit","se.fit")] %>% 
      bind_cols()
  ) %>% 
  mutate(prob = plogis(fit),
         prob_min = plogis(fit - 1.96*se.fit),
         prob_max = plogis(fit + 1.96*se.fit))

Fig.3A<-ggplot(dat_pred,
               aes(x = copy,
                   y = prob))+
  geom_line()+
  geom_ribbon(aes(ymin = prob_min,
                  ymax = prob_max), 
              alpha = 0.2)+
  geom_jitter(data = data_qPCR,
              aes(x = copy, y = detection),
              width = 0,
              height = 0.02,
              alpha = 0.25,
              size=5)+
  geom_hline(yintercept = 0.95,
             linetype=2)+
  geom_vline(xintercept = qPCR_LOD)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 25)) +
  labs(x="copies/Reaction",
       y = "Probability",
       title="A")+scale_x_log10()


### Fig 3B  Sample LOD（copies/L）

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
    text = element_text(size = 25)
  ) +
  labs(x="copies/L", 
       y = "Probability",
       title="B"
  )

Fig3<-Fig.3A/Fig.3B

Fig3






########## Environmental DNA degradation and shedding rates

###### Environmental DNA degradation rates in different water bodies (/s)

mydata1 <- read_csv("QPCR-mydata1.csv") %>% 
  mutate(
    tank=factor(tank,levels=c("A","B","C","D","E","F")),
    water=factor(water),
    conc = (q1 + q2 + q3)/3 * total_DNA / template / filtered_water
  )


### field water model
## Exponential decay equation of eDNA
model1.1 <- lmer(log(conc) ~ time+ (1|tank) , data = mydata1 %>% filter(water=="field"))
r2(model1.1)

## Nonlinear degradation of environmental DNA
model1 <- lmer(log(conc) ~ time +I(time^2)+ (1|tank), data = mydata1 %>% filter(water=="field"))
summary(model1)
r2(model1)



### tap water model 
## Exponential decay equation of eDNA
model2.1 <- lmer(log(conc) ~ time+ (1|tank), data = mydata1 %>% filter(water=="tap"))
r2(model2.1)

## Nonlinear degradation of environmental DNA
model2 <- lmer(log(conc) ~ time+I(time^2) + (1|tank), data = mydata1 %>% filter(water=="tap"))
summary(model2)
r2(model2)


# field water K1
k1 <- fixef(model1)["time"] %>% abs()

# tap water K2
k2 <- fixef(model2)["time"] %>% abs()


dat_plot1 <- expand.grid(
  tank=mydata1 %>% 
    filter(water=="field") %>% 
    pull(tank) %>% unique,
  time=seq(min(mydata1$time),
           max(mydata1$time),1)) %>% 
  mutate(
    tank=factor(tank),
    water="field",) %>% 
  mutate(fit=predict(model1,newdata=.))

dat_plot2 <- expand.grid(
  tank=mydata1 %>% filter(water=="tap") %>%
    pull(tank) %>% unique,
  time=seq(min(mydata1$time),
           max(mydata1$time),1)) %>% 
  mutate(
    tank=factor(tank),
    water="tap",
  ) %>% 
  mutate(fit=predict(model2,newdata=.))

dat_plot <- bind_rows(dat_plot1,dat_plot2)


Fig.4A <- ggplot() +
  geom_point(
    data = mydata1,
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
    y = "eDNA Concentration (copies/L)",
    title="A") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    text = element_text(size = 25,family = "serif")
  )




######## Shedding rate a （copies/s）
mydata2<-read_csv("QPCR-mydata2.csv") %>% 
  mutate(conc=(q1+q2+q3)/3*total_DNA/template/filtered_water,
         ind_L=ind/total_water,
         bio_L=bio/total_water,
         a=conc * k2/3600 * total_water)###### Use formulas: a = x * k2/3600 * V

### ind_L
model3.1<-lm(data=mydata2,log(a)~log(ind_L))
summary(model3.1)

###bio_L
model3<-lm(data=mydata2,log(a)~log(bio_L))
summary(model3)
β0 <- coef(model3)["(Intercept)"]
β1<- coef(model3)["log(bio_L)"]


fig_label <- paste0(
  str_c(
    "y=",
    round(β1,2),
    "*x+",
    round(β0,2)
  ),
  "\n",
  "R² = ", sprintf("%.3f", summary(model3)$adj.r.squared)
)
x_pos <- 10^(min(range(log10(mydata2$bio_L)))*1) %>% log()
y_pos <- 10^(max(log10(mydata2$a)) *1) %>% log()

Fig.4B <- ggplot(mydata2, aes(x = log(bio_L), y = log(a))) +
  geom_point(size = 5,
             alpha = 0.8) +  
  geom_smooth(method = "lm", 
              se = TRUE,formula = "y~x") +
  labs(
    x = "log(biomass density) (g/L)", 
    y = "log(eDNA shedding rate) (copies/s)",
    title="B") +
  annotate("text", 
           x = x_pos, 
           y = y_pos, 
           label = fig_label, 
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
    text = element_text(family = "serif", size = 25),
  )

fig4<-Fig.4A/Fig.4B
 
fig4







########## Field biomass density model  BIO(g/m^3)

mydata3<- read_csv("QPCR-mydata3.csv") %>%
  mutate(conc=(q1+q2+q3)/3*total_DNA/template/filtered_water,
         conc_s=conc*k1/3600*1,
         BIO=1000*exp((log(conc_s)-β0)/β1),
         season=factor(season,levels=c("spring","summer","autumn","winter"))
  )






########## Map of average biomass density in SY and WR

mydata3 %>%
  group_by(region,type) %>% 
  summarise(BIO1=quantile(BIO,0.25),
            BIO3=quantile(BIO,0.75))


###### SY map
mydata4_sy <- mydata3 %>%
  filter(region == "SY") %>%
  group_by(region,type,longitude,latitude) %>% 
  summarise(BIO=mean(BIO))


longitude_range_sy <- range(mydata4_sy$longitude)
longitude_breaks_sy <- seq(
  from = floor(longitude_range_sy[1] * 50) / 50,
  to = ceiling(longitude_range_sy[2] * 50) / 50,
  by = 0.04
)

fig1B <- ggplot() +
  geom_sf(data = target_map1, fill = "skyblue", alpha = 0.4) +
  geom_sf(data = target_map2, fill = "skyblue", alpha = 0.4) +
  geom_point(data = mydata4_sy, 
             aes(longitude, latitude, 
                 col = log(BIO), 
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green",
                       high = "red", 
                       name = "log(BIO)") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_sy) +
  coord_sf(xlim = range(mydata4_sy$longitude),
           ylim = range(mydata4_sy$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br", 
                   width_hint = 0.3) +
  theme_bw() +
  theme( 
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )

###### WR map
mydata4_wr <- mydata3 %>%
  filter(region == "WR") %>%
  group_by(region,type,longitude,latitude) %>% 
  summarise(BIO=mean(BIO))


longitude_range_wr <- range(mydata4_wr$longitude)
longitude_breaks_wr <- seq(
  from = floor(longitude_range_wr[1] * 50) / 50,
  to = ceiling(longitude_range_wr[2] * 50) / 50,
  by = 0.04
)

fig1C <- ggplot() +
  geom_sf(data = target_map1, fill = "skyblue", alpha = 0.4) +
  geom_sf(data = target_map2, fill = "skyblue", alpha = 0.4) +
  geom_point(data = mydata4_wr, 
             aes(longitude, latitude, 
                 col = log(BIO), 
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green", 
                       high = "red", 
                       name = "log(BIO)",
                       guide = "none") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_wr) +
  coord_sf(xlim = range(mydata4_wr$longitude),
           ylim = range(mydata4_wr$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br",
                   width_hint = 0.3) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )

Fig1BC <- fig1B / fig1C +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.box.margin = margin(l = 20)  
  )

Fig1BC




########## Seasonal biomass variations at SY and WR

### SY
mydata4_sy <-mydata3 %>%
  filter(region == "SY", season == "spring") %>%
  group_by(region, type, longitude, latitude)

longitude_range_sy <- range(mydata4_sy$longitude)
longitude_breaks_sy <- seq(
  from = floor(longitude_range_sy[1] * 50) / 50,
  to = ceiling(longitude_range_sy[2] * 50) / 50,
  by = 0.04
)

SY_spring <- ggplot() +
  geom_sf(data = target_map1, 
          fill = "skyblue",
          alpha = 0.4) +
  geom_sf(data = target_map2, 
          fill = "skyblue", 
          alpha = 0.4) +
  geom_point(data = mydata4_sy,
             aes(longitude, 
                 latitude,
                 col = log(BIO), 
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green", 
                       high = "red",
                       name = "log(BIO)") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_sy) +
  coord_sf(xlim = range(mydata4_sy$longitude),
           ylim = range(mydata4_sy$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br",
                   width_hint = 0.3) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )


mydata4_sy <-mydata3 %>%
  filter(region == "SY", season == "summer") %>%
  group_by(region, type, longitude, latitude)

longitude_range_sy <- range(mydata4_sy$longitude)
longitude_breaks_sy <- seq(
  from = floor(longitude_range_sy[1] * 50) / 50,
  to = ceiling(longitude_range_sy[2] * 50) / 50,
  by = 0.04
)

SY_summer <- ggplot() +
  geom_sf(data = target_map1,
          fill = "skyblue", 
          alpha = 0.4) +
  geom_sf(data = target_map2, 
          fill = "skyblue", 
          alpha = 0.4) +
  geom_point(data = mydata4_sy,
             aes(longitude, 
                 latitude, 
                 col = log(BIO), 
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green",
                       high = "red",
                       name = "log(BIO)") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_sy) +
  coord_sf(xlim = range(mydata4_sy$longitude),
           ylim = range(mydata4_sy$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br",
                   width_hint = 0.3) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )


mydata4_sy <-mydata3 %>%
  filter(region == "SY", season == "autumn") %>%
  group_by(region, type, longitude, latitude)

longitude_range_sy <- range(mydata4_sy$longitude)
longitude_breaks_sy <- seq(
  from = floor(longitude_range_sy[1] * 50) / 50,
  to = ceiling(longitude_range_sy[2] * 50) / 50,
  by = 0.04
)

SY_autumn <- ggplot() +
  geom_sf(data = target_map1, 
          fill = "skyblue",
          alpha = 0.4) +
  geom_sf(data = target_map2,
          fill = "skyblue", 
          alpha = 0.4) +
  geom_point(data = mydata4_sy,
             aes(longitude, latitude,
                 col = log(BIO), 
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green",
                       high = "red",
                       name = "log(BIO)") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_sy) +
  coord_sf(xlim = range(mydata4_sy$longitude),
           ylim = range(mydata4_sy$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br", 
                   width_hint = 0.3) +
  theme_bw() +
  guides(
    shape = guide_legend(order = 1),  
    color = guide_colorbar(order = 2)  
  ) +
  theme(
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )


mydata4_sy <-mydata3 %>%
  filter(region == "SY", season == "winter") %>%
  group_by(region, type, longitude, latitude)

longitude_range_sy <- range(mydata4_sy$longitude)
longitude_breaks_sy <- seq(
  from = floor(longitude_range_sy[1] * 50) / 50,
  to = ceiling(longitude_range_sy[2] * 50) / 50,
  by = 0.04
)

SY_winter <- ggplot() +
  geom_sf(data = target_map1,
          fill = "skyblue", 
          alpha = 0.4) +
  geom_sf(data = target_map2,
          fill = "skyblue", 
          alpha = 0.4) +
  geom_point(data = mydata4_sy,
             aes(longitude, latitude, 
                 col = log(BIO), 
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green",
                       high = "red",
                       name = "log(BIO)") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_sy) +
  coord_sf(xlim = range(mydata4_sy$longitude),
           ylim = range(mydata4_sy$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br",
                   width_hint = 0.3) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )

SY_season <- (SY_spring | SY_summer) / (SY_autumn | SY_winter) +
  theme(
    legend.position = "right",
    legend.box.margin = margin(l = 20)
  )
SY_season 




### WR

mydata4_wr <-mydata3 %>%
  filter(region == "WR", season == "spring") %>%
  group_by(region, type, longitude, latitude)

longitude_range_wr <- range(mydata4_wr$longitude)
longitude_breaks_wr <- seq(
  from = floor(longitude_range_wr[1] * 50) / 50,
  to = ceiling(longitude_range_wr[2] * 50) / 50,
  by = 0.04
)

WR_spring <- ggplot() +
  geom_sf(data = target_map1, 
          fill = "skyblue",
          alpha = 0.4) +
  geom_sf(data = target_map2,
          fill = "skyblue", 
          alpha = 0.4) +
  geom_point(data = mydata4_wr,
             aes(longitude, 
                 latitude, 
                 col = log(BIO), 
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green", high = "red",
                       name = "log(BIO)") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_wr) +
  coord_sf(xlim = range(mydata4_wr$longitude),
           ylim = range(mydata4_wr$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br",
                   width_hint = 0.3) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )


mydata4_wr <-mydata3 %>%
  filter(region == "WR", season == "summer") %>%
  group_by(region, type, longitude, latitude)

longitude_range_wr <- range(mydata4_wr$longitude)
longitude_breaks_wr <- seq(
  from = floor(longitude_range_wr[1] * 50) / 50,
  to = ceiling(longitude_range_wr[2] * 50) / 50,
  by = 0.04
)

WR_summer <- ggplot() +
  geom_sf(data = target_map1,
          fill = "skyblue",
          alpha = 0.4) +
  geom_sf(data = target_map2, 
          fill = "skyblue", 
          alpha = 0.4) +
  geom_point(data = mydata4_wr,
             aes(longitude,
                 latitude, 
                 col = log(BIO), 
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green",
                       high = "red",
                       name = "log(BIO)") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_wr) +
  coord_sf(xlim = range(mydata4_wr$longitude),
           ylim = range(mydata4_wr$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br",
                   width_hint = 0.3) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )


mydata4_wr <-mydata3 %>%
  filter(region == "WR", season == "autumn") %>%
  group_by(region, type, longitude, latitude)

longitude_range_wr <- range(mydata4_wr$longitude)
longitude_breaks_wr <- seq(
  from = floor(longitude_range_wr[1] * 50) / 50,
  to = ceiling(longitude_range_wr[2] * 50) / 50,
  by = 0.04
)

WR_autumn <- ggplot() +
  geom_sf(data = target_map1,
          fill = "skyblue",
          alpha = 0.4) +
  geom_sf(data = target_map2,
          fill = "skyblue", 
          alpha = 0.4) +
  geom_point(data = mydata4_wr,
             aes(longitude, latitude,
                 col = log(BIO), 
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green", 
                       high = "red",
                       name = "log(BIO)") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_wr) +
  coord_sf(xlim = range(mydata4_wr$longitude),
           ylim = range(mydata4_wr$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br",
                   width_hint = 0.3) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )


mydata4_wr <-mydata3 %>%
  filter(region == "WR", season == "winter") %>%
  group_by(region, type, longitude, latitude)

longitude_range_wr <- range(mydata4_wr$longitude)
longitude_breaks_wr <- seq(
  from = floor(longitude_range_wr[1] * 50) / 50,
  to = ceiling(longitude_range_wr[2] * 50) / 50,
  by = 0.04
)

WR_winter <- ggplot() +
  geom_sf(data = target_map1,
          fill = "skyblue", 
          alpha = 0.4) +
  geom_sf(data = target_map2, 
          fill = "skyblue", 
          alpha = 0.4) +
  geom_point(data = mydata4_wr,
             aes(longitude, 
                 latitude, 
                 col = log(BIO),
                 shape = type),
             size = 4) +
  scale_color_gradient(low = "green", 
                       high = "red",
                       name = "log(BIO)") +
  scale_shape_discrete(name = "type") +
  scale_x_continuous(breaks = longitude_breaks_wr) +
  coord_sf(xlim = range(mydata4_wr$longitude),
           ylim = range(mydata4_wr$latitude),
           expand = T) +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br", 
                   width_hint = 0.3) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank()
  )


WR_season <- (WR_spring | WR_summer) / (WR_autumn | WR_winter) +
  theme(
    legend.position = "right",
    legend.box.margin = margin(l = 20)
  )
WR_season




########## Biomass density across different land types

mydata4 <-mydata3 %>% 
  mutate(type = factor(type, levels = c("agriculture", "exurban", "urban")))


fig5<-ggplot(data = mydata4, 
             aes(x = type, 
                 y = BIO, 
                 fill = type)) +
  facet_grid(~region)+
  geom_boxplot(position = position_dodge(width = 0.8), 
               width = 0.7, 
               na.rm = TRUE) +  
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 20, 
                        family = "serif"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
  labs(fill = "type",
       x = "type",
       y = bquote(BIO~(g/m^3)))+
  ylim(0,10)

fig5





########## Water surface area (km²)

###SY

mydata4_sy <-mydata3 %>%
  filter(region == "SY") %>%
  group_by(region, type, longitude, latitude) %>% 
  summarise(BIO = mean(BIO))

sy_xlim <- range(mydata4_sy$longitude)
sy_ylim <- range(mydata4_sy$latitude)

sy_bbox <- st_bbox(c(
  xmin = sy_xlim[1], 
  ymin = sy_ylim[1], 
  xmax = sy_xlim[2], 
  ymax = sy_ylim[2]
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

SY_mean_BIO<- mean(mydata4_sy$BIO)
SY_water_volume <- water_area_sy_km2 * 3* 1e6
SY_total_BIO <- SY_mean_BIO * SY_water_volume




###WR

mydata4_wr <-mydata3 %>%
  filter(region == "WR") %>%
  group_by(region, type, longitude, latitude) %>% 
  summarise(BIO = mean(BIO))


wr_xlim <- range(mydata4_wr$longitude)
wr_ylim <- range(mydata4_wr$latitude)
wr_bbox <- st_bbox(c(
  xmin = wr_xlim[1], 
  ymin = wr_ylim[1], 
  xmax = wr_xlim[2], 
  ymax = wr_ylim[2]
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

WR_mean_BIO<- mean(mydata4_wr$BIO)
WR_water_volume <- water_area_wr_km2 * 3* 1e6
WR_total_BIO <- WR_mean_BIO * WR_water_volume






########## Distribution model of Pomacea

model4 <- glm(
  data = mydata3 %>%
    mutate(BIO=if_else(BIO>0,1,0)),
  BIO~log(ph)+log(do)+log(tur)+log(tm)+log(chla)+type+region,
  family = binomial(link="logit")
)
summary(model4)


plot_dat1 <- expand.grid(
  ph=mean(mydata3$ph),
  do=mean(mydata3$do),
  tur=mean(mydata3$tur),
  chla=mean(mydata3$chla),
  tm=seq(1,max(mydata3$tm),0.1),
  region=c("SY","WR"),
  type=c("agriculture","urban","exurban")
)
plot_dat1 <- predict(model4,newdata = plot_dat1,se=T) %>% 
  .[1:2] %>% 
  bind_cols() %>% 
  bind_cols(plot_dat1) %>% 
  mutate(
    fit.m=plogis(fit),
    fit.max=plogis(fit+1.96*se.fit),
    fit.min=plogis(fit-1.96*se.fit)
  )


fig6 <- ggplot() +
  geom_point(data = mydata3, 
             aes(x = tm, 
                 y = if_else(BIO > 0, 1, 0)), 
             size = 3, 
             alpha = 0.1) +
  scale_x_log10() +
  geom_line(data = plot_dat1, 
            aes(x = tm, 
                y = fit.m,
                color = type,
                linetype = region), 
            linewidth = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "serif",
                        size = 20)
  ) +
  xlab("Water temperature (°C)") +
  ylab("Probability") +
  labs(color = "type", 
       linetype = "region", fill = "type")

fig6



########## Linear model between biomass and influencing factors
model4.1 <- lm(
  data = mydata3,
  BIO~log(ph)+log(do)+log(tur)+log(tm)+log(chla)+type+region
)
summary(model4.1)



########## Linear model of biomass in two agricultural regions
model4.2 <- lm(
  data = mydata3 %>% filter(type=="agriculture"),
  BIO~region
)
summary(model4.2)


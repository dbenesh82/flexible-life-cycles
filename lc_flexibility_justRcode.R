## ---- message=FALSE, warning=FALSE---------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

options(stringsAsFactors = FALSE)

#set theme for plots
theme.o <- theme_update(axis.text = element_text(colour="black", size = 15),
                        axis.title = element_text(colour="black", size = 18, face = "bold", lineheight=0.25),
                        axis.ticks = element_line(colour="black"),
                        panel.border = element_rect(colour = "black",fill=NA),
                        panel.grid.minor=element_blank(),
                        panel.grid.major=element_line(color="gray",linetype = "dotted"),
                        panel.background= element_rect(fill = NA))

# import data
dataH <- read.csv(file="data/CLC_database_hosts.csv", header = TRUE, sep=",")
dataL <- read.csv(file="data/CLC_database_lifehistory.csv", header = TRUE, sep=",")

## ---- message=FALSE, warning=FALSE---------------------------------------
maxLCL <- group_by(dataH, Parasite.species)%>%
  summarize(maxLCL = max(Host.no))
minLCL <- filter(dataH, Facultative == "no")%>%
  group_by(Parasite.species)%>%
  summarize(minLCL = length(unique(Host.no)))

dataH <- left_join(dataH, maxLCL)
dataH <- left_join(dataH, minLCL)

## ---- message=FALSE, warning=FALSE---------------------------------------
# reduce to species level for plot
plot.dat <- select(dataH, Parasite.species, minLCL, maxLCL)%>%
  distinct()

plot.dat <- mutate(plot.dat, Flex.cycle = 
                        if_else(minLCL != maxLCL, "Some hosts\nfacultative", "All hosts\nobligate"),
                      maxLCL.fac = if_else(maxLCL > 3, "4", as.character(maxLCL)))%>%
  mutate(maxLCL.fac = factor(maxLCL.fac, labels = c("1", "2", "3", ">3")))

## ---- message=FALSE, warning=FALSE---------------------------------------
mypalette <- mypalette <- brewer.pal(3, "Accent") # colors to use in plot

ggplot(plot.dat, aes(x = maxLCL.fac, fill = Flex.cycle)) + 
  geom_bar() +
  scale_fill_manual(values = mypalette[1:2]) +
  labs(x = "Life cycle length", y = "Number of species")

## ---- message=FALSE, warning=FALSE---------------------------------------
mos.prop <- group_by(plot.dat, maxLCL.fac, Flex.cycle)%>%
  summarize(n = n())
t <- tapply(mos.prop$n, INDEX = mos.prop$maxLCL.fac, sum)

mos.prop$prop <- mos.prop$n/
  t[match(mos.prop$maxLCL.fac, names(t))]
rm(t)
mos.prop

## ---- message=FALSE, warning=FALSE---------------------------------------
ss <- group_by(mos.prop, maxLCL.fac)%>%summarize(n = sum(n)) # sample size for fig

outfig <- ggplot(mos.prop, 
                 aes(x = maxLCL.fac, y = prop, fill = Flex.cycle)) + 
  geom_bar(stat = "identity") +
  labs(x = "\nLife cycle length", y = "Proportion of species\n") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = mypalette[1:2]) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
  annotate("text", x = c(1, 2, 3, 4), y = 0.03, label = as.character(ss$n))
outfig

# want to keep so save to file
ggsave(filename = "figs/prop_fac_vs_lcl.png", width = 5, height = 4, units = "in")

## ---- message=FALSE, warning=FALSE---------------------------------------
flexs <- left_join(select(dataH, Parasite.species, Host.no, Stage, Def.int, Facultative),
                   plot.dat)
flexs <- distinct(flexs)%>%
  mutate(Fac_dummy = if_else(Facultative != "no", "yes", "no"))

## ---- message=FALSE, warning=FALSE---------------------------------------
ggplot(filter(flexs, Flex.cycle == "Some hosts\nfacultative"),
       aes(x = Host.no, fill = Fac_dummy)) + 
  geom_bar() + 
  scale_fill_manual(values = mypalette[1:2]) +
  facet_grid(~maxLCL.fac) +
  labs(x = "Host in cycle", y = "Count", fill = "Facultative?")

## ---- message=FALSE, warning=FALSE---------------------------------------
rig <- filter(flexs, Fac_dummy == "yes")%>%
  select(Parasite.species, Def.int, Fac_dummy)%>%
  distinct()

lef <- filter(flexs, Flex.cycle == "Some hosts\nfacultative")%>%
  select(Parasite.species, Def.int, maxLCL.fac)%>%
  distinct()

flexs2 <- left_join(lef, rig)
flexs2 <- mutate(flexs2, Fac_dummy = if_else(is.na(Fac_dummy), 'no', 'yes'))
rm(lef, rig)

## ---- message=FALSE, warning=FALSE---------------------------------------
flexs2 <- mutate(flexs2, Def.int = factor(Def.int, levels = c("int", "def")))%>%
                   mutate(Def.int = factor(Def.int, labels = c("Intermediate", "Definitive") ))
ggplot(flexs2,
       aes(x = Def.int, fill = Fac_dummy)) + 
  geom_bar() + 
  scale_fill_manual(values = mypalette[1:2]) +
  facet_grid(~maxLCL.fac) +
  labs(x = "Host", y = "Count", fill = "Facultative?") +
  theme(axis.text.x = element_text(size = 11))

## ---- message=FALSE, warning=FALSE---------------------------------------
dataL <- mutate(dataL, biovolume = 
                  if_else(Shape %in% c("cylinder", "thread-like", "whip"), 
                          pi * (Width/2)^2 * Length, # calculate volume as a cylinder
                          if_else(Shape %in% c("coiled", "sphere", "ellipsoid"),
                                  4/3 * pi * Length/2 * Width/4, # calculate volume as a ellipsoid
                                  Length * Width # calculate volume as area for remaining ribbon, leaf shapes
                                  )),
                biovolume = biovolume * 1.1) # covert to biomass with assumed 1.1. g/cm3 tissue density 

## ---- message=FALSE, warning=FALSE---------------------------------------
dataL <- filter(dataL, is.na(Asexual))%>% # remove data for asexual species
  filter( !(Stage == 'adult' & Sex == 'm') ) # remove adult males

## ---- message=FALSE, warning=FALSE---------------------------------------
# id species that hatch or not
eggos <- filter(dataL, Host.no == 0)%>%
  select(Parasite.species, Egg.hatch)%>%
  mutate(propagule_selector = if_else(Egg.hatch != "eaten", "free larva", "egg"))%>%
  select(-Egg.hatch)%>%
  na.omit%>%distinct()

# determine whether there is a size measurement for embryo or egg stages
eggos2 <- filter(dataL, Host.no == 0)%>%
  select(Parasite.species, Stage, biovolume)%>%
  group_by(Parasite.species, Stage)%>%
  summarize(x = sum(!is.na(biovolume)))

# combine and spread these two tables
eggos2 <- left_join(eggos, eggos2)
eggos2 <- spread(na.omit(eggos2), Stage, x)

# identify the stage where growth starts for each species
eggos2 <- mutate(eggos2, propagule_selector = if_else(propagule_selector == 'free larva', 'free larva',
                                                       if_else(embryo > 0, 'embryo', 'egg')))

# add selector variable to main life history table
eggos2 <- select(eggos2, Parasite.species, propagule_selector)
dataL <- left_join(dataL, eggos2)
rm(eggos, eggos2)

## ---- message=FALSE, warning=FALSE---------------------------------------
dataL <- filter(dataL, !(Host.no == 0 & Stage != propagule_selector))

## ---- message=FALSE, warning=FALSE---------------------------------------
dataL.sp <- group_by(dataL, Parasite.species, Host.no, Stage)%>%
  summarize(biovolume = mean(biovolume, na.rm=T))

## ---- message=FALSE, warning=FALSE---------------------------------------
dataL.sp <- arrange( ungroup(dataL.sp), Parasite.species, Host.no)%>% # arrange by species and host.no
  mutate(biov = lag(x = biovolume, 1))%>% # make a variable representing size in previous stage
  mutate(abs_diff = biovolume - biov, # absolute size diff
         rel_diff = log10(biovolume) - log10(biov)) # relative size diff

# remove growth values for egg stages; arise when calculating 'species B egg' - 'species A adult'
dataL.sp$abs_diff[which(dataL.sp$Host.no == 0)] <- NA 
dataL.sp$rel_diff[which(dataL.sp$Host.no == 0)] <- NA

# reduce data to just rows where a growth could be calculated
dataL.sp <- filter(dataL.sp, !is.na(abs_diff))%>%
  select(Parasite.species, Host.no, Stage, abs_diff, rel_diff)

## ---- message=FALSE, warning=FALSE---------------------------------------
lc_growth <- left_join(flexs, dataL.sp)

## ---- message=FALSE, warning=FALSE---------------------------------------
ggplot(lc_growth, 
       aes(x = Fac_dummy, y = rel_diff)) +
  geom_boxplot(outlier.color = "white") +
  geom_jitter(width = 0.2, height = 0, color = "red", alpha = 0.15) +
  labs(x = "Facultative?", y = "Orders of magnitude increase") + 
  theme(panel.grid.major.x = element_blank())

## ---- message=FALSE, warning=FALSE---------------------------------------
group_by(lc_growth, Fac_dummy)%>%
  summarize(fold_size_inc = mean(rel_diff, na.rm=T))%>%
  mutate(fold_size_inc = 10^fold_size_inc)

## ---- message=FALSE, warning=FALSE---------------------------------------
dataL.sp <- group_by(dataL, Parasite.species, Host.no, Stage)%>%
  summarize(biovolume = mean(biovolume, na.rm=T))

## ---- message=FALSE, warning=FALSE---------------------------------------
slope.plot <- full_join(lc_growth, dataL.sp)%>%arrange(Parasite.species, Host.no)

# add variables for plots
slope.plot <- mutate(slope.plot, 
                   Host.nofac = factor(Host.no, labels = c("propagule", "1st", "2nd", "3rd", "4th", "5th")),
                   trans.to.fac = "no")

slope.plot$trans.to.fac[which(slope.plot$Facultative != "no") - 1] <- "yes" # create variable for coloring the lines

## ---- message=FALSE, warning=FALSE---------------------------------------
outfig <- ggplot(data = slope.plot,
       aes(x = Host.nofac, y = log10(biovolume), 
           group = Parasite.species)) +
  geom_point(size = 3, alpha = 0.1, color = "gray") +
  geom_line(aes(color = trans.to.fac, alpha = trans.to.fac)) + 
  labs(x="\nHost",y="Log(Biomass)\n") +
  guides(color = FALSE, alpha = FALSE) +
  scale_alpha_discrete(range = c(0.15, 1)) +
  scale_color_manual(values = c("gray", "red")) +
  scale_x_discrete(expand=c(0.05,0.05)) +
  annotate("text", x = 4.5, y = -4.5, label = "Transmission to facultative host", 
           size = 7, color = "red")
outfig

# want to keep so save to file
ggsave(filename = "figs/lc_flex_slopegraph.png", width = 7, height = 5, units = "in")


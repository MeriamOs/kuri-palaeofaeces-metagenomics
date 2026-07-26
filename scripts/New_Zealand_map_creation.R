# NEW ZEALAND MAP CREATION
#Script written by Meriam van Os

library(sf)
library(ggplot2)
library(ggspatial)
library(rnaturalearth)

nz <- ne_countries(country = "New Zealand",
  scale = 10, returnclass = "sf")

ggplot(nz) +
  geom_sf(fill = "darkolivegreen3", colour = "black", linewidth = 0.2) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "br",     
                   width_hint = 0.25,
                   text_cex = 1.4) +
  coord_sf(xlim = c(155, 180),
            ylim = c(-48, -34)) +
  theme_bw() +
  theme(axis.text = element_text(size = 18))

ggsave("pub_map_NZ_rnaturalearth.png", width = 9, height = 7.5)

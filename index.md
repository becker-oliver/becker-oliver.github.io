---
title: "Random sample of Berlin areas"
author: "Oliver Becker"
date: "2020-05-20"
output: html_document
---



# Package osmdata and bounding boxes/polygons

The R package **osmdata** provides a functions that delivers either the bounding box (rectangular) or the bounding polygon of an object in OpenStreetMaps.


```r
library(osmdata); library(dplyr); library(sf)
# Retrieve outer border of Berlin via osmdata-package
berlin_border <- getbb ("berlin germany", format_out = "polygon") %>% data.frame()
# Convert it to a simple feature object
berlin_border <- berlin_border %>% st_as_sf(coords = c("X1","X2")) %>% st_set_crs(4326)
```

Att this point we have a bunch of points. Now we combine them to a polygon:

```r
berlin_polygon <- berlin_border %>% summarise(geometry = st_combine(geometry)) %>% 
  st_cast("POLYGON")
# Plot
plot(berlin_polygon, col = "cyan")
```

<img src="figure/unnamed-chunk-2-1.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

# Create raster for Berlin and sample

The function `sp::bbox` returns the bounding box of a polygon.


```r
library(raster)
r <- raster(xmn = bbox(st_coordinates(berlin_polygon))[1],   # set minimum x coordinate
                      xmx = bbox(st_coordinates(berlin_polygon))[3],    # set maximum x coordinate
                      ymn = bbox(st_coordinates(berlin_polygon))[2],     # set minimum y coordinate
                      ymx = bbox(st_coordinates(berlin_polygon))[4],     # set maximum y coordinate
                      res = c(0.1, 0.05)) # resolution in c(x,y) direction
r
```

```
## class      : RasterLayer 
## dimensions : 7, 7, 49  (nrow, ncol, ncell)
## resolution : 0.1, 0.05  (x, y)
## extent     : 13.08835, 13.78834, 52.32551, 52.67551  (xmin, xmax, ymin, ymax)
## crs        : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0
```

```r
# Assign values to raster cells
r[] <- 1:ncell(r)
plot(r)
```

<img src="figure/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

Here we save the raster in order to convert it to a sf-object via stars-package:

```r
writeRaster(r,'~/Desktop/test.tif', overwrite = TRUE)
library(stars)
tif=read_stars("~/Desktop/test.tif")
sf=st_as_sf(tif)
# Use st_intersects to filter cells that intersect with the Berlin polygon
sf_berlin <- sf %>% filter(row_number() %in% st_intersects(berlin_polygon, sf)[[1]])

# Use fasterize to convert sf object back to raster object
library(fasterize)
r <- fasterize(sf_berlin, r, field = "test.tif", fun="sum")
plot(r)
```

<img src="figure/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```r
x <- sampleRandom(r, ncell(r)*.3, asRaster=TRUE)
plot(x)
```

<img src="figure/unnamed-chunk-4-2.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

To underlie the the OSM map, we can use the leaflet-package based on the sf sample of raster-polygons:

```r
library(leaflet)
sf_berlin_random <- sample_n(sf_berlin, round(nrow(sf_berlin)*0.3))
leaflet(sf_berlin_random) %>% addTiles() %>% addPolygons()
```

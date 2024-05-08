#### The function PROC_DISTRIB provides workflow 
#### to compute geomorphic processes distribution 
#### over the flow path distance from divide to outlet

### CREATOR - Sergey Kharchenko, LomonosovMSU ###

proc_distrib <- function(dem, # digital elevation model file
                         proc, # vector file (e.g. *.shp) with the processes contours 
                         proc_col, # name of column with the processes names
                         eq_areas = T, # if TRUE - equi-areal distribution, if FALSE - equi-distant distribution
                         bands_n = 20, # number of bands/levels
                         acc = F, # TRUE - left2right is down-top, if FALSE - left2right is top-down
                         c2bb = T, # crop the territory to the biggest basin 
                         plot_progress = F, # plot intermediate raster plots
                         areaplot = F, # plot stacked chart just after finish
                         save_metrics = F, # save physical distance too
                         save_data = F, # save intermediate raster and vector files
                         ...)
{ 
  install.packages(c("terra", "whitebox")[!(c("terra", "whitebox") %in% installed.packages()[,1])])
  
  require(terra)
  require(whitebox)
  if (!is.character(wbt_version())) install_whitebox()
  
  proc_initial <- vect(proc)
  border <- aggregate(proc_initial)
  writeVector(border, "border.shp", overwrite=TRUE)
  wbt_remove_polygon_holes("border.shp", "border.shp", verbose_mode = F)
  proc <- vect("border.shp")
  
  writeRaster(crop(rast(dem), proc, mask = T), "dem_crop.tif", overwrite=TRUE)
  
  wbt_fill_depressions(dem = "dem_crop.tif", output = "dem_filled.tif", verbose_mode = F)
  #wbt_max_upslope_flowpath_length(dem = "dem_filled.tif", output = "fpl_upslope.tif")
  wbt_d8_pointer(dem = "dem_filled.tif", output = "d8p.tif", verbose_mode = F)
  
  if (c2bb) {
    wbt_basins(d8_pntr = "d8p.tif", output = "basins.tif")
    bas <- as.polygons(rast("basins.tif"))
    bas <- bas[which.max(expanse(bas)),]
    writeRaster(crop(rast("d8p.tif"), bas, mask = T), "d8p.tif", overwrite=TRUE)
  }
  
  wbt_downslope_flowpath_length(d8_pntr = "d8p.tif", output = "fpl_dslope.tif", verbose_mode = F)
  
  fpl_crop <- rast("fpl_dslope.tif")
  plot(fpl_crop, main = "Downslope Flow Path Length")
  
  if (eq_areas) {
    qs <- global(fpl_crop, function(x) quantile(x, seq(0, 1, 1/bands_n), na.rm = T))
  } else {
    qs <- seq(0, global(fpl_crop, function(x) max(x, na.rm = T))[[1]], len = bands_n + 1)
  }
  
  if (!acc) {
    qs <- rev(qs)
  }
  
  pr_data <- matrix(NA, ncol = (length(qs)-1), nrow = nrow(unique(as.data.frame(proc_initial)[proc_col])))
  rownames(pr_data) <- unlist(unique(as.data.frame(proc_initial)[proc_col]))
  
  for (i in 1:(length(qs)-1)) {
    fpl_temp <- fpl_crop
    if (acc) {
      fpl_temp[fpl_temp < qs[i][[1]] | fpl_temp > qs[i + 1][[1]]] <- NA
    } else {
      fpl_temp[fpl_temp > qs[i][[1]] | fpl_temp < qs[i + 1][[1]]] <- NA
    }
    
    
    if (plot_progress) {
      plot(fpl_temp, main = paste0(i, " of ", bands_n))
      plot(proc_initial, add = T)
    } else {
      print(paste0(i, " of ", bands_n))
    }
    
    fpl_temp[!is.na(fpl_temp)] <- 1
    fpl_temp <- as.polygons(fpl_temp)
    
    proc_crop <- crop(proc_initial, fpl_temp)
    eval(parse(text = paste0("agr_table <- aggregate(expanse(proc_crop) ~ ",proc_col,", as.data.frame(proc_crop), FUN = function(x) sum(x, na.rm = T))")))
    
    pr_data[,i] <- agr_table[match(rownames(pr_data), agr_table[,1]),2]
  }
  pr_data[is.na(pr_data)] <- 0
  pr_data <- t(pr_data)
  
  if (areaplot) {
    install.packages(c("areaplot")[!(c("areaplot") %in% installed.packages()[,1])])
    require(areaplot)

    cols <- hcl.colors(ncol(pr_data), palette = "viridis", alpha = 0.8)
    
    areaplot(1:bands_n, pr_data, prop = TRUE, col = cols,
             legend = TRUE, xlab = paste0("Percentiles by ",100/bands_n, " %"),
             ylab = "Part of the area",
             args.legend = list(x = "topleft", cex = 0.65,
                                bg = "white", bty = "o"))
  }
  
  if (!save_data) {unlink(c("dem_crop.tif", "dem_filled.tif", "d8p.tif", "fpl_dslope.tif", "border.shp"))}
  
  if (save_metrics) {
    return(list(qs, pr_data))
  } else {
    return(pr_data)
  }
}


#longl <- proc_distrib(dem = "D:/d_mosaic2.tif", 
#                      proc = "C:/rates/processes_en.shp", 
#                      proc_col = "process_en", 
#                      c2bb = T,
#                      eq_areas = T, 
#                      bands_n = 20,
#                      plot_progress = T,
#                      areaplot = T,
#                      acc = F,
#                      save_data = F)
#longl <- proc_distrib(dem = "D:/d_mosaic2.tif", 
#                      proc = "C:/rates/processes_en.shp", 
#                      proc_col = "process_en")
#colnames(longl)
#colnames(longl) <- c("Glacier", "Rock falls", "Rock creep", "Soil creep", "Stable areas", "Gully erosion","Stream erosion", "Sheet wash")
#longl <- longl[,c("Glacier", "Rock falls", "Rock creep", "Soil creep", "Sheet wash", "Gully erosion", "Stream erosion", "Stable areas")]
#cols <- hcl.colors(ncol(longl), palette = "viridis", alpha = 0.8)
#
#longl <- proc_distrib(dem = "H:/FP/demdonguz.tif", 
#                      proc = "C:/rates/processes_dif_erosion.shp", 
#                      proc_col = "Process", 
#                      c2bb = T,
#                      eq_areas = T, 
#                      bands_n = 20,
#                      plot_progress = T,
#                      areaplot = T,
#                      acc = F,
#                      save_data = F)
#colnames(longl)
#colnames(longl) <- c("Glacier", "Stream erosion", "Lake sed.", "Rock falls", "Soil creep", "Sheet wash",  "Stable areas", "Rock creep", "Gully erosion")
#longl <- longl[,c("Glacier", "Rock falls", "Rock creep", "Soil creep", "Sheet wash", "Gully erosion", "Stream erosion", "Stable areas", "Lake sed.")]
#cols <- c(hcl.colors(ncol(longl)-1, palette = "viridis", alpha = 0.8), "blue")
#
#
#areaplot(1:nrow(longl), longl, prop = TRUE, col = cols,
#         legend = TRUE, xlab = paste0("Percentiles by ",100/bands_n, " %"),
#         ylab = "Part of the area",
#         args.legend = list(x = "topleft", cex = 0.8,
#                            bg = "white", bty = "o"))
dm <- downslope_movement()

downslope_movement <- function(proc_vect = "C:/rates/processes_dif_erosion.shp", 
                               proc_col = "Process", 
                               proc_of_int = c("Крип толщи грубообломочного материала"), # , "Крип и солифлюкция толщи мелкозема"
                               dem = "H:/FP/demdonguz2.tif", 
                               bufsize = 5) {
  proc <- vect(proc_vect)
  poi <- proc[as.data.frame(proc)[,proc_col] %in% proc_of_int]
  poi$Process <- "poi"
  poi <- aggregate(poi)
  poi <- fillHoles(poi)
  poi <- disagg(poi)
  
  writeVector(poi, "process_of_interest.shp", overwrite = T)
  
  other_proc <- proc[!as.data.frame(proc)[,proc_col] %in% proc_of_int]
  
  list_results <- list(len = length(poi))
  
  dem <- rast(dem)
  
  if (!as.logical(grep('metre|meter', crs(dem)))) {
    extdem <- project(as.polygons(ext(dem), crs = crs(dem)), "EPSG:4326")
    extdem <- crds(centroids(extdem))[1,1]
    utmepsg <- paste0("EPSG:326",30 + extdem %/% 6 + ifelse(extdem > 0, 1, -1))
    dem <- project(dem, utmepsg)
  }
  
  aspect <- terrain(dem, "aspect")
  
  if (is.null(bufsize)) {
    bufsize <- 5*(res(dem)[1])
  }
  
  for (i in 1:length(poi)) {
    poi_sample <- poi[i,]
    plot(poi_sample)
    crdss <- crds(poi_sample)
    crdss <- data.frame(head(crdss, -1), tail(crdss, -1), rbind(tail(crdss, 2)[1,], head(crdss, -2)))
    dir1 <- 180 * atan2((crdss[,3] - crdss[,1]), (crdss[,4] - crdss[,2])) / pi
    dir1 <- ifelse(dir1 > 0, dir1, dir1 + 360)
    dir2 <- 180 + 180 * atan2((crdss[,5] - crdss[,1]), (crdss[,6] - crdss[,2])) / pi
    
    mean_dir <- (dir1 + dir2) / 2
    mean_dir <- mean_dir - 90
    mean_dir <- ifelse(mean_dir < 0, mean_dir + 360, mean_dir)
    
    mean_segment <- (((crdss[,3] - crdss[,1])^2 + (crdss[,4] - crdss[,2])^2)^0.5 + ((crdss[,5] - crdss[,1])^2 + (crdss[,6] - crdss[,2])^2)^0.5)/2
    
    nodes <- as.points(poi_sample)
    bufs <- buffer(nodes, bufsize)
    
    aspect_dir <- extract(aspect, bufs, mean)[,2]
    
    print(paste0("There are ", nrow(bufs), " nodes in the ",i, " polygon"))
    down_dir <- sapply(1:nrow(bufs), function(x) {
      #plot(crop(dem, bufs[x,], mask = T))
      #plot(proc, add = T)
      dat_elevs <- extract(crop(dem, bufs[x,], mask = T), proc, function(a) mean(a, na.rm = T))
      dat_elevs <- na.omit(dat_elevs)
      is.poi <- as.data.frame(proc[dat_elevs$ID,])[which.max(dat_elevs[,2]),proc_col] %in% proc_of_int
      
      #print(x)
      if (nrow(dat_elevs) == 1) {
        list(ifelse(dat_elevs[,2] > global(crop(dem, bufs[x,], mask = T), function(x) mean(x, na.rm = T))[[1]], TRUE, FALSE), "Unknown")
      } else if (is.poi) {
        cross_point <- c(crdss[x,1] + 0.5*bufsize * sin(pi * aspect_dir[x] / 180), 
                         crdss[x,2] + 0.5*bufsize * cos(pi * aspect_dir[x] / 180))
        cross_point <- matrix(cross_point, nrow = 1)
        whereto <- as.data.frame(other_proc)[relate(other_proc, vect(cross_point, type = "points", crs = crs(other_proc)), "contains"),proc_col]
        if (length(whereto) == 0) {
          if (length(crop(other_proc, bufs[x,])) == 0) {
            list(is.poi, "Unknown")
          } else {
            whereto <- as.data.frame(crop(other_proc, bufs[x,])[which.max(expanse(crop(other_proc, bufs[x,]))),])[,proc_col]
            list(is.poi, whereto)
          }
        } else {
          list(is.poi, whereto)
        }
        
      } else {
        list(is.poi, "Unknown")
      }
    })
    if (is.null(dim(down_dir))) {
      down_dir <- do.call(rbind, down_dir)
    } else if (dim(down_dir)[2] != 2 & dim(down_dir)[2] == nrow(bufs) & dim(down_dir)[1] == 2) {
      down_dir <- t(down_dir)
    }
    
    #down_dir2 <- down_dir[down_dir[,1] == TRUE,]
    #View(down_dir)
    #which(unlist(down_dir[,1]) == TRUE & lapply(down_dir[,2], function(x) ifelse(is.null)) == NULL)
    
    if (is.null(nrow(down_dir))) {
      index2down <- unlist(down_dir[1])
    } else {
      index2down <- unlist(down_dir[,1])
    }
    #table(unlist(down_dir[,2]))
    
    if (sum(index2down)) {
      wei <- abs(cos(pi * (aspect_dir - mean_dir) / 180))[index2down] / sum(abs(cos(pi * (aspect_dir - mean_dir) / 180))[index2down])
      
      proc_diff <- aggregate(wei ~ procs, data.frame(wei, procs = unlist(down_dir[index2down,2])), FUN = sum)
      list_results[[i]] <- list(sum(mean_segment[index2down]),
                                sum(mean_segment[index2down] * abs(cos(pi * (aspect_dir - mean_dir) / 180))[index2down]),
                                proc_diff)
    } else {
      list_results[[i]] <- list(NULL, NULL, NULL)
    }
    
    print(paste0(i, " object from ", length(poi), " was done"))
    timestamp()
  }
  return(list_results)
}

setwd("F:/Donguz/")
library(terra)
pde <- vect("Processes/processes_different_erosion.shp")
pde$Area <- expanse(pde)
aggregate(round(Area/10000, 2) ~ Process, as.data.frame(pde), sum)

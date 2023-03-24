library(raster)
library(ncdf4)
source("cfg.r")


dir_in = 'data/carboncci/'

dir_our = 'outputs/carboncci/'

files = list.files(dir_in)

ncuts = 72
nboots = 50

egRaster = raster('data/GFED4s_burnt_area.nc')
egRaster[] = 0.0
egArea = raster::area(egRaster)
tfile0 = 'temp/make_cci_carbon_map-cut'

write_temp_text_file <- function(tfile, txt) {
    fileConn = file(tfile)
    writeLines(c(txt), fileConn)
    close(fileConn)
}

regridFile <- function(file) {
    file_in = paste0(dir_in, file)
    fileID =  substr(file, 1, nchar(file)- 3)
    if (grepl("CHANGE", file)) return()
    nc = nc_open(file_in)
    variables = names(nc$var)[1:2]
    
    mu = raster(file_in, varname = variables[1])
    sigma = raster(file_in, varname = variables[2])

    extent = extent(mu)

    cut_extents <- function(i, j) {
        out = seq(extent[i], extent[j], length.out = ncuts+1)
        return(head(out, -1))
    }
    xs = cut_extents(1, 2); ys = cut_extents(3, 4)
    
    xdiff = diff(xs[1:2]); ydiff = diff(ys[1:2])
    regrid4cut <- function(x, y) {
        
        tfile = paste(tfile0, ncuts, nboots, x, y, xdiff, ydiff,
                      fileID, '.txt', sep = '-')
        print(file_in)
        print(x)
        print(y)
        tfile_nc = paste(tfile0, 'nc', ncuts, nboots, x, y, xdiff, ydiff,
                         fileID, '.nc', sep = '-')
        if (file.exists(tfile)) {
            info = read.table(tfile, stringsAsFactors=FALSE)
            if (info[1,1] == "empty") {
                return()        
            } else {
                if (file.exists(tfile_nc)) return(brick(tfile_nc))
            }
        }
        
        cut_extent = c(x, x+xdiff, y, y+ydiff)  
        
        mu = terra::crop(mu, cut_extent)
        if (is.na(maxValue(mu)) || maxValue(mu) == 0) {
             write_temp_text_file(tfile, "empty")
            return()
        }
            
        print("yay")
        sigma = terra::crop(sigma, cut_extent)
        
        mu_boot = mu
        
        egRasteri = raster::crop(egRaster, extent(mu) + c(-1, 1, -1, 1))
        mu_agg = terra::resample(mu_boot, egRasteri) 
        
        doBoot <- function(id) {
            print(id)
            mu_boot = sigma* rnorm(ncell(sigma), 0, 1)            
            boot = mu_agg + terra::resample(mu_boot, egRasteri) 
            boot[boot<0] = 0
            return(boot)
        }
        out = layer.apply(1:nboots, doBoot) * egArea * 100
        
        out = writeRaster(out, file = tfile_nc, overwrite = TRUE)
        write_temp_text_file(tfile, tfile_nc)
        return(out)        
    }

    
    combine_all_ens <- function(outs, resam = TRUE) {
        combine_ens <- function(i) {
            print(i)
            
            selectResample <- function(out) {
                out = out[[min(nlayers(out), i)]]
                if (resam) out = raster::resample(out, egRaster)
                out
            }
            ens = layer.apply(outs, selectResample)
            
            sum(ens/egArea, na.rm = TRUE)
        }
        layer.apply(1:max(sapply(outs, nlayers)), combine_ens)
    }

    nn <<- 0
    applyOverY <- function(y) {
        
        tfile = paste(tfile0, ncuts, nboots,'ymerge', y, fileID, '.nc', sep = '-')
        print(tfile)
        
        if (file.exists(tfile)) return( brick(tfile))
        #nn <<- nn + 1
        #print(nn)
        #if (nn < 15) return() 
        return()
        lock_file = paste(tfile0, 'lock', y, fileID, '.txt', sep = '-')
        
        if (file.exists(lock_file)) return('locked')
        file.create(lock_file)
       
        out = lapply(xs, regrid4cut, y)
        #return()
        out = unlist(out)
        if (is.null(out))  out = egRaster
            else out = combine_all_ens(out)
        #return(out)
        writeRaster(out, file = tfile, overwrite = TRUE)
        
        return(out)
    }
    
    outs = lapply(rev(ys), applyOverY)
    
    if (any(sapply(outs, function(i) is.character(i)))) stop("some locked files")
    outs = outs0 = outs[!sapply(outs, is.null)]
    
    outs = combine_all_ens(outs, FALSE)
    browser()    
}

lapply(files, regridFile)

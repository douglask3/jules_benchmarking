source("cfg.r")

paths = c(GFED5 = "/data/dynamic/dkelley/fireMIPbenchmarking/data/benchmarkData/ISIMIP3a_obs/GFED5_Burned_Percentage.nc",
        JULES = "/scratch/hadea/isimip3a/u-cc669_isimip3a_fire/20CRv3-W5E5_obsclim/jules-fire-vn6p3_20crv3-w5e5_obsclim_histsoc_default_burntarea-total_global_monthly_1901_2019.nc")

mask = raster('data/remapcon_ce_new.nc')

months = c(5, 9)
year_range = c(2002, 2020)
extent = c(-65, -40,-30, 0)

reds9   = c("#FFF7FB", "#F7DEEB", "#EFC6DB", "#E19ECA", "#D66BAE",
            "#C64292", "#B52171", "#9C0851", "#6B0830")

width = 1
centre = seq(-1, 10, 0.1) + width/2


mask = raster::crop(mask, extent)
forDat <- function(path, name = '') {
    dat = brick(path)
    dat = raster::crop(dat, extent)
    dat[is.na(mask)] = NaN
    dat[dat > 9E9] = NaN

    dates = names(dat)
    years = as.numeric(substr(dates, 2, 5))
    dat = dat[[which((years > year_range[1])  & (years < year_range[2]))]]

    getMnth <- function(mn) {
        mn = seq(mn, nlayers(dat), by = 12)
        dat[[mn]][]
    }
    datv = sapply(months, getMnth)
    datv = datv[apply(datv, 1, function(i) !any(is.na(i))),]

    colnames(datv) = c('may', 'sep')
    #datv = log10(datv+0.1)
    

    cols = reds9
    cols = cols[unlist(mapply(rep, 3:9, 9 + (3:9)^2))]
    cols = densCols(datv[,1],datv[,2], colramp = colorRampPalette(cols), bandwidth = 0.1)
    plot(datv[,2]~datv[,1], pch = 20, col = cols, cex = 0.5, xlab = 'May', ylab = 'September', 
         axes = F, xlim = c(0, 8.5), ylim = c(0, 30))
    mtext(name, line = -1.5, adj = 0.5, side = 3)
    #labels = c(0, 0.1, 0.2, 0.5, 1,2, 5, 10, 20, 50, 100)
    #at = labels#log10(labels + 0.1)

    axis(1)#, labels = labels, at = at)
    axis(2)# labels = labels, at = at)


    add_rq <- function(tau) {
        #fit = rq(sep ~ may, tau = tau, data = data.frame(datv))
        #sep = predict(fit, newdata = data.frame(may = may))
        #lines(may, sep, lty = 2)
        forBin <- function(x) {
            test = ((x-width) < datv[,'may']) & (datv[,'may'] < (x+width))
            if (sum(test) < 1/(1-tau)) return(NaN)
            quantile(datv[test, 'sep'], tau)
        }
        sep = sapply(centre, forBin)
        lines(centre, sep, lty = 2)
        pos = tail(which(!is.na(sep)), 1)
        text(x = centre[pos], y = sep[pos], tau, xpd = NA, adj = c(1,-0.5))    
    }
    lapply(c(0.1, 0.5, 0.9, 0.95, 0.99) , add_rq)
}
png("figs/May_Sep_burnt_area.png", height = 7, width = 5, res = 300, units = 'in')
par(mfrow = c(2,1), mar = c(1, 0, 1, 0), oma = c(3, 3, 0, 3))
mapply(forDat, paths, names(paths))
dev.off()





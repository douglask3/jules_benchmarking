source("cfg.r")

load("outputs/regional_process_ts.Rd")

models = c("GFDL-ESM2M", "HADGEM2-ES", "IPSL-CM5A-LR", "MIROC5")

metricComparison <- function(var, name, FUN = NME, Null = null.NME, obs_sc = 1) {
    obss = var[[1]]
    mods = var[[2]]
    obs_start = var[[3]]
    mod_start = var[[4]]
    regionAndAver <- function(dat, start = 1, end = NULL) {
        out = sapply(2:15, function(rn) dat[[2]][[rn]])
        if (!is.null(end)) out = out[start:end,]
        out = apply(out, 2, mean)
        out / regionArea
    }
    compareObs<- function(ON) {
        start = obs_start[ON]-mod_start+1
        end = length(obss[[ON]][[2]][[1]])
        
        
        x = regionAndAver(obss[[ON]])
        compareMod <- function(SN) {
            y = regionAndAver(mods[[SN]], start, end)
            score(FUN(x, y * obs_sc, w = regionArea))
        }
        nullms = Null(x, w = regionArea)
        nullms = summary(nullms)[c(1:3, 3)] + c(0, 0, c(-1, 1) *summary(nullms)[4])
        scores = sapply(1:length(mods), compareMod)
        nullms =  t(matrix(nullms, nrow = 4, ncol = nrow(scores)))
        cbind(nullms, scores)
    }
    
    
    out = lapply(1:length(obss), compareObs)
    out = do.call(rbind, out)
    obs_names = names(obss)

    rownames(out) =  paste(rep(obs_names, each = nrow(out)/length(obs_names)),
                           rownames(out), sep = ' - ')

    colnames(out) = c('Median Null', 'Mean Null', 
                      'Randomly Resampled lower', 'Randonly Resampled upper', models)
    
    if (length(obss) > 1) {
        combined_obs = sapply(2:15, function(rgn)
            range(sapply(obss, function(obs) mean(obs[[2]][[rgn]])))/regionArea[rgn-1])
        NME_maxmin <- function(SN) {
            y = regionAndAver(mods[[SN]])
            top = apply(abs(sweep(combined_obs, 2, y, '-')), 2, min) * regionArea
            bot = apply(abs(combined_obs - mean(combined_obs)), 2, min) * regionArea
            score = sum(top)/sum(bot)
        }
        scores_combined = sapply(1:length(mods), NME_maxmin)
        median = sum(apply(abs(combined_obs - median(combined_obs)), 2, min) * regionArea)/
                 sum(apply(abs(combined_obs - mean(combined_obs)), 2, min) * regionArea)
        combined = c(median, 1, NaN, NaN, scores_combined)
        out = rbind(out, combined)
    }
    
    write.csv(out, file = paste0('outputs/benchmarking_by_regions-', name, '.csv'))
    return(out)
}

metricComparison(burnt_area, 'burnt_area')

metricComparison(tallTrees_on, 'tallTrees_on')
metricComparison(tallTrees_off, 'tallTrees_off')

metricComparison(fire_emissions, 'fire_emissions')

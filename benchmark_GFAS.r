source("cfg.r")

dir = '../jules_outputs/'

mods = c('u-bi607_HIST_new_annual')
varName = 'veg_c_fire_emission_gb;burnt_carbon_dpm;burnt_carbon_rpm'
varName = "fire_em_CO2_gb"
modScale = 1#60*60*24*360


extent = c(-90, -30, -60, 15)
#extent = c(-180, 180, -90, 90)

obss = paste0("../fireMIPbenchmarking/data/benchmarkData/",
              c("GFAS.nc"))
names(obss) = c("GFAS")
years = list(2001:2013)
obsLayers = list(11:166)
obsScale = 12

limits_aa = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)                
cols_aa = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
            '#fc4e2a','#e31a1c','#bd0026','#800026')

dlimits_aa = c(-20, -10, -5, -1, -0.5, -0.1, 0.1, 0.5, 1, 5, 10, 20)                  
dcols_aa = rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

limits_tr = c(-0.1, -0.01, -0.005, -0.001, -0.0001, 0.0001, 0.001, 0.005, 0.01, 0.1)#c(-2, -1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1, 2)/10         
cols_tr = rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

dlimits_tr = limits_tr*10            
dcols_tr = cols_tr

cols_phase = c('#313695', '#a50026', '#ffff00','#313695')
limits_phase = 0.5:11.5 

dcols_phase = c("#003300", "#542788", "#f7f7f7", "#b35806", "#003300")
dlimits_phase = (-5.5:5.5)


cols_conc = c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4',
                                       '#1d91c0','#225ea8','#253494','#081d58')
limits_conc =  seq(0, 0.9, 0.1)


dcols_conc = c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30')
dlimits_conc =  c(-0.8, -0.6, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.6, 0.8)

cols_modal = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a')#c('#34eeba','#d95f02','#7570b3')
limits_modal = c(1, 1.1, 1.2, 1.5, 2)


dcols_modal = c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')
dlimits_modal = c(-1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1)

score_aa = runComparison(mods, obss, 'figs/fire_emissions_aa.png', 
                               cols_aa, limits_aa, dcols_aa, dlimits_aa, TRUE)

score_tC = runComparison(mods, obss, 'figs/fire_emission_tC.png',
                              cols_tr, limits_tr, dcols_tr, dlimits_tr, 
                              FALSE, TRUE, trendTransFUN = NULL, extend_min = TRUE)

out = mapply(function(i, j) cbind(j, rownames(i), round(i,2)),
             list(score_aa, score_tC),
             c("Annual average", "Trend"), SIMPLIFY = FALSE)
out = do.call(rbind, out)
write.csv(out, "outputs/fire_emissions_comparison.csv")

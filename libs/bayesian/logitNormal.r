logitNormal <- function(obs, sim) {
    
    obs0 = obs == 0
    obs1 = !obs0
    obs = logit(obs)
    simt = logit(sim)
    sim2 = (1-sim)^2
    FUN <- function(params) {
        P0 = sim2 * (1 - params[1])
        P1 =  (1 - P0[obs1]) * (dnorm(obs[obs1], simt[obs1], params[2]))
        
        -mean(log(P0[obs0])) - sum(log(P1))
    }
    -optim(c(0.1, 1), FUN)$value
}


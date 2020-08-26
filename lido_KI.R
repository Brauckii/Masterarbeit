### Hier wird der Simulationsdatensatz für die Konfidenzintervallanalyse erzeugt

### Konfidenzintervallanalyse

source(file = "Functions.R")

reg = makeRegistry(file.dir = NA)
#reg = makeRegistry(file.dir = NA, work.dir = "/work/smrnbrau", seed = 2506)


show_calc_PAR_withKI <- function(n, M, pE, rho, IA = 0){
  
  realPAR <- bind_cols("n" = n, extract_RealPar(M, pE, rho, IA))
  
  simdat <- dat_sim(n, M, pE, rho, interaction = IA)
  
  #Boot_perc <- PAR_rbn(simdat, Var = "boot", CI = "perc", B = 500) %>% dplyr::select(c("conf.low","conf.high")) %>% setNames(paste0("Boot_perc:",colnames(.)))
  
  Boot_rbn <- PAR_rbn(simdat, Var = "boot", CI = "normal", B = 500) %>% dplyr::select(c("conf.low","conf.high")) %>% unite_(col = "Boot_rbn",from = c("conf.low", "conf.high"),sep = "_")
  
  Delta_rbn <- PAR_rbn(simdat, Var = "delta", CI = "normal") %>% dplyr::select(c("conf.low","conf.high")) %>% unite_(col = "Delta_rbn",from = c("conf.low", "conf.high"),sep = "_")
  
  Boot_partial <- PAR_partial(simdat, Var = "boot", CI = "normal", B = 500) %>% dplyr::select(c("conf.low","conf.high")) %>% unite_(col = "Boot_partial",from = c("conf.low", "conf.high"),sep = "_")
  
  Delta_partial <- PAR_partial(simdat, Var = "delta", CI = "normal") %>% dplyr::select(c("conf.low","conf.high")) %>% unite_(col = "Delta_partial",from = c("conf.low", "conf.high"),sep = "_")
  
  # Es kann passieren, dass NANs bzw NAS produziert werden, in diesem Fall werden Tibbles mit 0 Zeilen ausgegeben
  modoutput <- function(TIB){
    TIB <- as_tibble(matrix(rep(NaN,M*ncol(TIB)), nrow = M)) %>% setNames(colnames(TIB))
    return(TIB)
  }
  
  if(nrow(Boot_rbn) == 0) Boot_rbn <- modoutput(Boot_rbn)
  if(nrow(Delta_rbn) == 0) Delta_rbn <- modoutput(Delta_rbn)
  
  result <- bind_cols(realPAR, Boot_rbn, Delta_rbn, Boot_partial, Delta_partial)
  
  return(result)
  
}

repls = 200
git_KI <- as.matrix(expand.grid("n" = c(50,100,1000) ,"M" = 3:5, "pE" = c(0.1,0.25,0.5), "rho" = seq(0,0.5,0.25), repl = seq_len(repls)))

batchMap(show_calc_PAR_withKI, n = git_KI$n, M = git_KI$M, pE = git_KI$pE, rho = git_KI$rho)

#res = list(measure.memory = TRUE, partition =  "short", walltime = 3600)
submitJobs(reg = reg)

waitForJobs(reg = reg)

Result <- reduceResultsList(fun = as.data.frame, reg = reg)
KI_result <- do.call(rbind, Result)


# Speichern der Daten
save(KI_result, file = "KI_result.Rdata")

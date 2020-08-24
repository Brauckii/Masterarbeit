### Hier wird der Datensatz für die Laufzeitanalyse erzeugt

### Funktion basierend auf der sys.time

source(file = "Functions.R")

reg = makeRegistry(file.dir = NA)

show_runtime_sys.time <- function(n,M){
  
  simdat <- dat_sim(n, M, pE = 0.25, rho = 0)
  
  # fmla aus der Datensimulation ist fuer den Algorithmus aus dem pARtial Paket unbrauchbar
  fmla = paste("D ~ ", paste(colnames(simdat$E), collapse = "+"))
  
  runtime_PAR_rbn <- system.time(PAR_rbn(simdat))[3]
  
  runtime_PAR_partial <- system.time(pARtial::PartialAR(D = simdat$D, x = simdat$E, model = TRUE, fmla = fmla))[3]
  
  runtime_PAR_rbn_mf <- system.time(PAR_rbn(simdat, model = FALSE))[3]
  
  runtime_PAR_partial_mf <- system.time(pARtial::PartialAR(D = simdat$D, x = simdat$E, model = NULL))[3]
  
  runtime <- c(n,M,runtime_PAR_rbn, runtime_PAR_partial, runtime_PAR_rbn_mf, runtime_PAR_partial_mf)
  
  return(runtime)
  
}


repls = 100

git_1 <- expand.grid("n" = 100 , "M" = c(3,4,6,8,10,12,14), repl = seq_len(repls))
git_2 <- expand.grid("n" = c(50,100,500,1000,5000,10000,1E5), "M" = 6, repl = seq_len(repls))

git_runtime <- rbind(git_1, git_2)
rm(git_1)
rm(git_2)

batchMap(show_runtime_sys.time, n = git_runtime$n, M = git_runtime$M)

# res = list(measure.memory = TRUE, partition =  "short", walltime = 7200)
submitJobs(reg = reg)

waitForJobs()

# An der Stelle weiß ich nicht was richtig ist!!!
Result <- reduceResultsList(fun = as.data.frame, reg = reg)
Runtime_result <- t(do.call(cbind, Result))

Runtime_result <- set_names(as_tibble(Runtime_result), c("n","M", "runtime_PAR_rbn", "runtime_PAR_partial", "runtime_PAR_rbn_mf", "runtime_PAR_partial_mf"))


# Speichern der Daten
save(Runtime_result, file = "Runtime_result.Rdata")

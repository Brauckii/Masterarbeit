##### R-Code mit allen wichtigen Funktionen; gedacht fuer den Durchlauf auf dem LiDO3

inst_packages = installed.packages()[,1]
inst_if_missing = function(x, ...) {
  if (!(x %in% inst_packages)) install.packages(x, ...)
}
### Pakete laden
inst_if_missing("e1071")
inst_if_missing("mvtnorm")
inst_if_missing("tidyverse")
inst_if_missing("rockchalk")
inst_if_missing("boot")
inst_if_missing("broom")
#install.packages("~/Partial Attributable RIsk/Alte Pakete/pARtial_0.1.1.tar.gz", repos = NULL, type = "source")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/pARtial/pARtial_0.1.1.tar.gz"	
if(length(which(installed.packages()[,1] == "pARtial")) == 0) install.packages(packageurl, repos=NULL, type="source")
inst_if_missing("parallel")
inst_if_missing("numDeriv")
inst_if_missing("forcats")
inst_if_missing("batchtools")
#inst_if_missing("microbenchmark")
#inst_if_missing("rlist")


library(e1071)
library(tidyverse)
library(mvtnorm) 
library(rockchalk)
library(boot)
library(broom)
library(pARtial)
library(parallel)
library(numDeriv)
library(forcats)
library(batchtools)
#library(microbenchmark)
#library(rlist)

# Vorarbeit 

logit <- function(p) log(p/(1-p))
logit.inv <- function(y) 1/(1+exp(-y))
# (= exp(y)/1+exp(y))

### Simulieren der Daten

# Daten werden in Abhängigkeit von 
# n (# Beobachtungen), M (# Expositionen), pE (Auftrittsw'keit der Expositionen) 
# rho (Korrelation zwischen den Expositionen), IA (Interaktionseffekte zwischen den Expos) gezogen

dat_sim <- function(n, M, pE, rho = 0, interaction = 0){
  
  # Erstelle Expositionen, 0 oder 1
  
  ### generate random expositions from a multivariate normal distribution
  # Varianzmatrix: M x M (Varianz 1, alle Kovarianzen 0.49)
  # Confounding Struktur
  #Sigma <- matrix(data = c(rep(c(1, rep(0.49, M)), M-1), 1), ncol = M)
  # Expositionsmatrix n x M wird aus multivariater Normalverteilung gezogen
  # Falls Werte größer 0 -> TRUE -> *1 = 1, FALSE -> 0 
  # M Werte -> final 0.5, das heißt 50% Quantil -> mu = 0 qnorm(0.5,0,1) = 0
  # W'keit das Exposition vorliegt liegt bei maximal 50%
  #E     <- (mvrnorm(n = n, mu = qnorm(seq(0.2, 0.5, (0.5-0.2)/(M-1))), Sigma = Sigma) > 0)*1
  
  Sigma <- lazyCov(Rho = rho, Sd = 1, d = M)
  
  # mu = qnorm(pE) gewährleistet, dass im Mittel pE * n die Exposition auftritt (da Standardabweichung = 1)
  E <- (rmvnorm(n = n, mean = rep(qnorm(pE),M), sigma = Sigma) > 0)*1
  
  # interaction = 0 means no interaction
  if(interaction > 0){  
    V <- combn(M,2)
    lI <- ncol(V)
    # Betrachte jede Spalte
    I <- apply(V, 2, function(x) E[,x[1]] * E[,x[2]])
    colnames(I) <- paste0("I", 1:lI)
    if(interaction == 1) E_new <- cbind(E,I) else E_new <- cbind(I,E)
    # Wie viele Effekte haben einen Einfluss auf die Erkrankung?
    lE <- M + lI
  } 
  else{
    E_new <- E
    lE <- M
  }
  
  ### derive D from a logistic model
  
  beta <- logit(c(0.2, rev(logit.inv(2*logit(0.2)/-(sum(1:lE))*1:lE))))
  # sum(beta) = 0 -> 1. beta = Intercept ist negativ; danach alle beta größer 0 -> alle Expositionen haben negativen Einfluss auf Krankheit
  # bei den weiteren beta wird ein linearer Anstieg simuliert
  # für M = 1 -> 0.2 und 0.8 berechnet, für M sehr groß logit 0.2 und M mal logit 0.5
  # d.h. auch: keine Exposition liegt vor, W'keit, dass Kranheit vorliegt = 0.2
  
  # D (bzw. Zielvariable) erstellen; Designmatrix (mit Intercept) 
  lincomb <- cbind(rep(1, n), E_new) %*% beta 
  
  D <- logit.inv(lincomb)
  
  # Dichtomisiere D, vorherige D fungieren als W'keit für "Münzwurf"
  D <- rbinom(length(D),1,D)
  
  colnames(E) <- paste0("E",1:M)
  
  fml <- as.formula(paste("D ~", paste(colnames(E), collapse = "+")))
  
  MAT <- combn(paste0("E",1:M), 2)
  intterm <- apply(MAT, 2, function(i) paste0(i, collapse = ":"))
  fml2 <- as.formula(paste("D ~", paste(c(colnames(E),intterm), collapse = "+")))
  
  dat<-list(beta = beta, D=D, E=E, fmla=fml, fmla_IA = fml2)
  
  return(dat)
  
}

### Funktion zur Berechnung des PARs

PAR_rbn <- function(dat, model = TRUE, fmla = dat$fmla ,Var = c("none","boot","delta"), CI = c("none", "normal" ,"basic", "perc", "bca"), alpha = 0.05, B = 500, parallel = FALSE, n_cpus = 7){
  
  # Anzahl an Expositionen
  M <- ncol(dat$E)
  
  rawdata <- cbind.data.frame(D = dat$D, dat$E)
  
  ## delete rows with NA, Beobachtungen mit fehlenden Werten
  del <- which(is.na(rowSums(rawdata)))
  if (length(del) > 0 ) rawdata <- rawdata[-del, ]
  
  # Zielvariable und Exposition müssen binär sein
  if(any(!(rawdata == 1 | rawdata == 0))) stop("Error: Only dichotomized (0/1) variables are allowed!")
  
  # Nehme ersten Eintrag bei Varianz und Konfidenzintervall -> Default: "none"
  # keine Varianzschaetzung und kein Konfidenzintervall schaetzen
  Var <- Var[1]
  CI <- CI[1]
  
  par_boot <- function(rawdata, sv){
    
    rawdata <- rawdata[sv,]
    
    D <- rawdata[,1]
    E <- rawdata[,-1]
    
    # Vorbereiten einer Matrix, die alle relevanten Informationen aus dem Datensatz enthaelt
    
    # Anzahl an Expositionen  
    M<-ncol(E)
    
    #praevalenz (P(D))
    p.d = mean(D)
    # W enthält alle möglichen Kombinationen von Expositionen,   
    W<-bincombinations(M)
    
    colnames(W)<-colnames(E)
    
    zp<-2^((M-1):0)
    
    # Vektor von 0 bis 2^M - 1
    W1<-W%*%zp #Zahl=pb*Linearkombi for all
    
    #Prävalenzen in Data, Häufigkeit der verschiedenen Expositionszusammensetzungen
    D3<-as.matrix(E)%*%zp #für alle
    # summe pij aufsummiert = 1
    pij<-sapply(W1, function(x) sum(D3==x)/nrow(D3)) #unter allen -->relevant
    
    # modellbasierter Ansatz: prediction basiert auf logistischer Regression
    if (model){
      mod <- glm(fmla, family = binomial(link="logit"), data = rawdata)
      
      # Vorhersage fuer die jeweile Exposition
      predict<-predict(mod, data.frame(W), "response")
    } 
    else{
      
      index <- lapply(W1, function(x) which(D3 ==x))
      predict <- unlist(lapply(index, function(x) mean(D[x])))
      # ist hier möglich, da ein NAN nur entsteht, falls pij ebenfalls 0 ist
      predict[is.nan(predict)] <- 0
      predict <- predict
      #predict <- unlist(lapply(lapply(W1, function(y) which(DE ==y)), function(x) mean(D[x])))
    }
    
    # W_new 
    # Verwende eine Matrix zur schnelleren Berechnung (alle Variablen numerisch!)
    W_new <- cbind(W,pij,predict)
    # W_new<-data.frame(W,W1,pij, predict)
    
    if(Var == "delta"){
      
      # Funktion die theoretisch die predictions berechnet!
      
      # Delta Methode basiert ausschließlich auf der logistischen Regression ohne Interaktionen
      
      if(model){
        
        f <- function(beta, Wv){
          XBET <- Wv %*% beta
          out <- logit.inv(XBET)
          return(out)
        }
        
        # MAT <- cbind(rep(1,dim(W)[1]),W)
        MAT <- model.matrix(as.formula(paste0("~", str_split(fmla, "~")[[3]])), as.data.frame(W))
        
        # Berechne Ableitung für jede relevante Expositionskombination
        varMAT <- t(apply(MAT, 1, function(i) grad(f, x = as.vector(mod$coefficients), Wv = i)))
        
      }
      
      else{
        
        # Ableitungen ergeben eine Diagonalmatrix
        varMAT <- diag(nrow(W))
        
        # Kovarianzmatrix basiert auf der Kovarianz einer multibinomialverteilung
        
        covMAT <- W_new[,"predict"]*(1-W_new[,"predict"])
        covMAT_mf <- diag(covMAT/(W_new[,"pij"]*nrow(rawdata)))
        
      }
      
    }
    
    # Hilfsfunktion zur Berechnung des ARs, s sei eine beliebige Sequenz
    ARadj_com<-function(s){
      
      j<-which(rowSums(W_new[,s, drop=FALSE]==0)==ncol(W_new[,s, drop=FALSE]))
      #j<-which(W_new[,s]==0)  (versagt, falls s Vektor und keine Zahl ist!)
      
      ARadj<-1-(((t(W_new[j,"pij"])%*%W_new[j,"predict"])/sum(W_new[j,"pij"]))/p.d)
      
      return(ARadj)
      
    }
    
    SAR<-function(S,i){
      #i=which Exposure
      
      q<-which(S==i) #Position of Exposure in Sequenz
      
      
      s1<-S[1:q]  # Betrachte bis zur Stelle, an dem die Exposition in der Sequenz steht
      s2<-S[1:(q-1)] # eine Stelle weniger 
      
      sar<- ARadj_com(s1)-dplyr::if_else(q==1,0,ARadj_com(s2)) # Verwende bekannte Formel! 
      
      return(sar)      
    }
    
    # Ableitung Nenner (Ableitung der Praevalenz) (delta(p(D)))
    if(Var == "delta"){
      delta_p.d_mat <- W_new[,"pij"] * varMAT
      delta_p.d <- colSums(delta_p.d_mat)
    }
    
    delta_SAR <- function(S,i){
      
      q<-which(S==i) #Position of Exposure in Sequenz
      
      s1<-S[1:q]  # Betrachte bis zur Stelle, an dem die Exposition in der Sequenz steht
      s2<-S[1:(q-1)] # eine Stelle weniger 
      
      # Ableitung Nenner (Ableitung der Praevalenz) (delta(p(D)))
      #delta_p.d_mat <- W_new[,"pij"] * varMAT
      #delta_p.d <- colSums(delta_p.d_mat)
      
      delta_AR<-function(s,...){
        
        j<-which(rowSums(W_new[,s, drop=FALSE]==0)==ncol(W_new[,s, drop=FALSE]))
        #j<-which(W_new[,s]==0)  (versagt, falls s Vektor und keine Zahl ist!)
        
        #ARadj<-1-(((t(W_new[j,"pij"])%*%W_new[j,"predict"])/sum(W_new[j,"pij"]))/p.d)
        
        # Verwendet wird die Quotientenregel
        
        # Ableitung Zaehler
        DD <- matrix(delta_p.d_mat[j,], nrow = length(j))
        DD <- colSums(DD)
        
        nen <- sum(W_new[j,"pij"])*(p.d^2)
        #delta_AR<-1-((colSums(DD)/sum(W_new[j,"pij"]))/p.d)
        #delta_AR <- - DD/(sum(W_new[j,"pij"])*p.d)
        delta_AR <- - (DD*p.d - as.numeric(t(W_new[j,"pij"])%*%W_new[j,"predict"])*delta_p.d)/nen
        
        return(delta_AR)
        
      }
      
      delta_s2 <- if(q==1) 0 else  delta_AR(s2) # Verwende bekannte Formel! 
      delta_sar <- delta_AR(s1) - delta_s2
      
      return(delta_sar)      
    }
    
    
    # Funktion verhindert, dass unnötig alle Sequenzen berechnet werden
    # Bsp: SAR(1,(1,2,3)) = SAR(1,(1,3,2))
    reduce_seq.fun<-function(M,i){
      
      if(M == 2){
        ALLKOM <- list(seq_len(M)[-i])
      }
      else{
        ALLKOM<-unlist(lapply(seq_len(M-1), function(x) combn(seq_len(M)[-i],x, simplify = FALSE)), recursive=FALSE)
      }
      # combn(seq_len(M)[-i],x, simplify = FALSE)) Liste mit verschiedenen Möglichkeiten, ziehe x Zahlen aus 1 bis M (ohne i)
      # M (- i), die vor i stehen können, was fehlt ist NULL (i ganz vorne)
      L<-unlist(lapply(ALLKOM, length))
      # Vektor mit der Vektorlänge
      
      Lv <- c(0,L)
      Anzahl <- factorial(Lv)*factorial(rev(Lv))
      
      
      ALLKOM<-rbind(c(i,setdiff(seq_len(M),i)),
                    matrix(unlist(lapply(ALLKOM, function(x) c(x,i,setdiff(seq_len(M)[-i],x)))), ncol=M, byrow=TRUE))
      
      seq.sort<-cbind(ALLKOM,Anzahl)
      
      return(seq.sort) 
      
    }  
    
    # Berechne das PAR
    
    par_loop <-function(i){
      
      # Anzahl ist irrelevant
      SAR_vor <- reduce_seq.fun(M,i)
      seq.red <- SAR_vor[,1:M]
      anzahl <- SAR_vor[,M+1]
      
      # Position von i in jeder Zeile
      SAR_all<- sapply(seq_len(nrow(seq.red)), function(x) SAR(seq.red[x,],i))
      
      # alle benötigten SAR's werden hier berechnet
      
      par<-(t(anzahl)%*%SAR_all)/factorial(M)
      return(par)
      
    }
    
    delta_par_loop <- function(i){
      
      SAR_vor <- reduce_seq.fun(M,i)
      seq.red <- SAR_vor[,1:M]
      anzahl <- SAR_vor[,M+1]
      
      var_SAR_all<- do.call(rbind, lapply(seq_len(nrow(seq.red)), function(x) delta_SAR(seq.red[x,],i)))
      
      stdpar<-(colSums(anzahl*var_SAR_all))/factorial(M)
      
      # möglicherweise anpassen
      #stdpar <- sqrt(as.numeric(t(stdpar) %*% vcov(mod) %*% stdpar)) 
      
      if(model){
        
        stdpar <- sqrt(as.numeric(t(stdpar) %*% vcov(mod) %*% stdpar)) 
        
      }
      
      else{
        
        stdpar <- sqrt(abs(as.numeric(t(stdpar) %*% covMAT_mf %*% stdpar)))
        #stdpar <- as.numeric(t(stdpar) %*% covMAT_mf %*% stdpar)
        
      }
      # Anzahlen werden relevant, in der Formel werden alle SARs betrachtet (aufsummiert), hier Anzahl mal das SAR
      
      return(stdpar)
      
    }
    
    if(Var == "delta"){
      par.est <- sapply(1:M, par_loop)
      std.par <- sapply(1:M, delta_par_loop)
      #varPar.est <- varMAT
      TIB <- tibble(PAR = par.est, std.error = std.par) 
      if(CI == "none"){
        TIB <- TIB %>% mutate(variance = std.error^2)
      }
      else if(CI == "normal"){
        TIB <- TIB %>% mutate(conf.low = par.est - qnorm(1- alpha/2)*std.error, conf.high = par.est + qnorm(1-alpha/2)*std.error, variance = std.error^2)
      }
      else{
        message("Please choose CI = none or normal")
      }
      # TIB reicht vermutlich aus
      return(TIB)
    }#names(par.est) <- paste0("PAR(E",1:M,")")
    else{
      par.est <- sapply(1:M, par_loop)
      return(par.est)
    }
    
    
  }
  
  if(Var == "none"){
    PAR <- par_boot(rawdata)
    names(PAR) <- paste0("PAR(E",1:M,")")
    result <- PAR
  }
  else if(Var == "delta"){
    PAR <- par_boot(rawdata)
    result <- PAR
  }
  else{
    
    if(Var == "boot"){
      parallel <- dplyr::if_else(parallel == TRUE, "multicore","no")
      boot.mod <- boot::boot(data = rawdata, statistic = par_boot, R = B, parallel = parallel, ncpus = n_cpus)
      #normaler boot
      #PAR <- tidy(boot.mod) 
      #rownames(PAR) <- colnames(dat$E)
      
      if(CI == "none"){
        result <- tidy(boot.mod) %>% mutate(variance = std.error^2)
      }
      else if(CI == "normal"){
        KI <- tidy(boot.mod) %>% mutate(conf.low = statistic - bias - std.error * qnorm((1 + (1-alpha))/2)) %>% 
          mutate(conf.high = statistic - bias  + std.error * qnorm((1 + (1-alpha))/2)) %>% mutate(variance = std.error^2)   
        result <- KI 
      }
      else{
        #KI <- boot.ci(boot.mod, conf = 1-alpha, type = CI)
        KI <- tidy(boot.mod, conf.int=TRUE, conf.level= 1-alpha ,conf.method= CI) %>% mutate(variance = std.error^2)
        result <- KI
      }
      result <- result %>% rename(PAR = statistic)
    }  
  }
  
  return(result)
  
}

### Funktion zur Berechnung des PARs (basierend auf dem pARtial- Paket)

PAR_partial <- function(data, model = TRUE, Var = c("none","delta","boot","bayes","jackknife"), CI = c("none","normal","percentile","BCa"), alpha = 0.05, B=500){
  
  #Var <- match.arg(Var)
  #CI <- match.arg(CI)
  
  fmla = dplyr::if_else(model == TRUE, paste("D ~ ", paste(colnames(data$E), collapse = "+")), "")
  result <- pARtial::PartialAR(D = data$D, x = data$E, model = model, fmla = fmla, Var = Var, CI = CI)
  
  Var <- Var[1]
  CI <- CI[1]
  
  if(Var != "none"){
    if(CI == "none"){
      result <-  tibble(PAR = result$PAR, variance = result$VarPAR)
    }
    else{
      result <- tibble(PAR = result$PAR, conf.low = result$CIPAR[1,], conf.high = result$CIPAR[2,], variance = result$VarPAR)
    }
  }
  
  return(result)
  
}

# Laden eines Datensatz, der für die verschiedenen Kombinationen die wahren PARs enthält

load("Result_RealPAR.Rdata")

# Funktion welches das wahre PAR fuer die gegebene Datensituation auswirft

extract_RealPar <- function(M,pE,rho,IA){
  
  ll <- c(M,pE,rho,IA)
  real_PAR <- Result_RealPAR %>% filter(M == ll[1], pE == ll[2], rho == ll[3], IA == ll[4])
  
  return(real_PAR)
  
}



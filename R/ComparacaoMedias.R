#'Comparacao de medias
#'@description "This function returns the multiple
#'  comparison of means by the t, t with Bonferroni protection, Duncan, SNK,
#'  Tukey and Scott-Knott tests".
#'  Esta funcao retorna a comparacao multipla de medias obtidas pelos
#'  testes t, t com protecao de Bonferroni, Duncan, SNK, Tukey e Scott-Knott
#'  (Funcao adaptada do expDes.pt).
#'@usage ComparacaoMedias(y, trt, DFerror, MSerror, alpha = 0.05, group =
#'  TRUE,main = NULL)
#'@param         y Vetor numerico contendo a variavel resposta.
#'@param         trt Vetor numerico ou complexo contendo os tratamentos.
#'@param         DFerror  Grau de liberdade do residuo.
#'@param         MSerror Quadrado medio do residuo.
#'@param         alpha Nivel nominal de significancia.
#'@param         group TRUE ou FALSE
#'@param         main Titulo.
#'@return Retorna a comparacao multipla de medias obtida por varios testes.
#'@references BANZATTO, D. A.; KRONKA, S. N. Experimentacao Agricola. 4 ed.
#'  Jaboticabal: Funep. 2006. 237 p. ISBN: 85-87632-71-X
#'
#'  GOMES, F. P. Curso de Estatistica Experimental. 10a ed. Piracicaba:
#'  ESALQ/USP. 1982. 430.
#'


ComparacaoMedias=function(y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL) {
######################################
lastC=function (x)
{
  y <- sub(" +$", "", x)
  p1 <- nchar(y)
  cc <- substr(y, p1, p1)
  return(cc)
}
#####################################
tapply.stat=function (y, x, stat = "mean")
{
  cx <- deparse(substitute(x))
  cy <- deparse(substitute(y))
  x <- data.frame(c1 = 1, x)
  y <- data.frame(v1 = 1, y)
  nx <- ncol(x)
  ny <- ncol(y)
  namex <- names(x)
  namey <- names(y)
  if (nx == 2)
    namex <- c("c1", cx)
  if (ny == 2)
    namey <- c("v1", cy)
  namexy <- c(namex, namey)
  for (i in 1:nx) {
    x[, i] <- as.character(x[, i])
  }
  z <- NULL
  for (i in 1:nx) {
    z <- paste(z, x[, i], sep = "&")
  }
  w <- NULL
  for (i in 1:ny) {
    m <- tapply(y[, i], z, stat)
    m <- as.matrix(m)
    w <- cbind(w, m)
  }
  nw <- nrow(w)
  c <- rownames(w)
  v <- rep("", nw * nx)
  dim(v) <- c(nw, nx)
  for (i in 1:nw) {
    for (j in 1:nx) {
      v[i, j] <- strsplit(c[i], "&")[[1]][j + 1]
    }
  }
  rownames(w) <- NULL
  junto <- data.frame(v[, -1], w)
  junto <- junto[, -nx]
  names(junto) <- namexy[c(-1, -(nx + 1))]
  return(junto)
}


##################################
#-> order.groups
#####
order.group= function (trt, means, N, MSerror, Tprob, std.err, parameter = 1)
{
  N <- rep(1/mean(1/N), length(N))
  n <- length(means)
  letras <- letters
  if (n > 26) {
    l <- floor(n/26)
    for (i in 1:l) letras <- c(letras, paste(letters, i,
                                             sep = ""))
  }
  z <- data.frame(trt, means, N, std.err)
  w <- z[order(z[, 2], decreasing = TRUE), ]
  M <- rep("", n)
  k <- 1
  j <- 1
  k <- 1
  cambio <- n
  cambio1 <- 0
  chequeo = 0
  M[1] <- letras[k]
  while (j < n) {
    chequeo <- chequeo + 1
    if (chequeo > n)
      break
    for (i in j:n) {
      minimo <- Tprob * sqrt(parameter * MSerror * (1/N[i] +
                                                      1/N[j]))
      s <- abs(w[i, 2] - w[j, 2]) <= minimo
      if (s) {
        if (lastC(M[i]) != letras[k])
          M[i] <- paste(M[i], letras[k], sep = "")
      }
      else {
        k <- k + 1
        cambio <- i
        cambio1 <- 0
        ja <- j
        for (jj in cambio:n) M[jj] <- paste(M[jj], " ",
                                            sep = "")
        M[cambio] <- paste(M[cambio], letras[k], sep = "")
        for (v in ja:cambio) {
          if (abs(w[v, 2] - w[cambio, 2]) > minimo) {
            j <- j + 1
            cambio1 <- 1
          }
          else break
        }
        break
      }
    }
    if (cambio1 == 0)
      j <- j + 1
  }
  w <- data.frame(w, stat = M)
  trt <- as.character(w$trt)
  means <- as.numeric(w$means)
  N <- as.numeric(w$N)
  std.err <- as.numeric(w$std.err)

  output <- data.frame(Tratamentos=trt, Medias=means, M, N, std.err)
  return(output)
}

##########################################
order.stat.SNK=function (treatment, means, minimum)
{
  n <- length(means)
  z <- data.frame(treatment, means)
  w <- z[order(z[, 2], decreasing = TRUE), ]
  M <- rep("", n)
  k <- 1
  k1 <- 0
  j <- 1
  i <- 1
  r <- 1
  cambio <- n
  cambio1 <- 0
  chequeo = 0
  M[1] <- letters[k]
  while (j < n) {
    chequeo <- chequeo + 1
    if (chequeo > n)
      break
    for (i in j:n) {
      if (abs(j - i) == 0) {
        r <- 1
      }
      else {
        r <- abs(j - i)
      }
      s <- abs(w[i, 2] - w[j, 2]) <= minimum[r]
      if (s) {
        if (lastC(M[i]) != letters[k])
          M[i] <- paste(M[i], letters[k], sep = "")
      }
      else {
        k <- k + 1
        cambio <- i
        cambio1 <- 0
        ja <- j
        for (jj in cambio:n) M[jj] <- paste(M[jj], " ",
                                            sep = "")
        M[cambio] <- paste(M[cambio], letters[k], sep = "")
        for (v in ja:cambio) {
          if (abs(v - cambio) == 0) {
            r <- 1
          }
          else {
            r <- abs(v - cambio)
          }
          if (abs(w[v, 2] - w[cambio, 2]) > minimum[r]) {
            j <- j + 1
            cambio1 <- 1
          }
          else break
        }
        break
      }
    }
    if (cambio1 == 0)
      j <- j + 1
  }
  w <- data.frame(w, stat = M)
  trt <- as.character(w$treatment)
  means <- as.numeric(w$means)
    output <- data.frame(trt, means, M)
  return(output)
}

####################################################
##########-> Teste t
###################################################
TesteT=function (y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL){
  SSerror <- MSerror*DFerror
  name.y <- paste(deparse(substitute(y)))
  name.t <- paste(deparse(substitute(trt)))
  junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
  means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
  sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
  nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
  means <- data.frame(means, std.err = sds[, 2]/sqrt(nn[, 2]),
                      replication = nn[, 2])
  names(means)[1:2] <- c(name.t, name.y)
  ntr <- nrow(means)
  Tprob <- qt(1 - (alpha/2), DFerror) * sqrt(2)
  nr <- unique(nn[, 2])
  nfila <- c("Alpha", "Error Degrees of Freedom", "Error Mean Square",
             "Critical Value of Studentized Range")
  nvalor <- c(alpha, DFerror, MSerror, Tprob)
  xtabla <- data.frame(...... = nvalor)
  row.names(xtabla) <- nfila
  if (group) {
    if (length(nr) == 1) {
      HSD <- Tprob * sqrt(MSerror/nr)
    }
    else {
      nr1 <- 1/mean(1/nn[, 2])
      HSD <- Tprob * sqrt(MSerror/nr1)
    }
    #cat("\nTeste t (LSD)\n------------------------------------------------------------------------")
    #cat("\nGrupos  Tratamentos  Medias\n")
    output= order.group(means[, 1], means[, 2], means[,
                                                        4], MSerror, Tprob, means[, 3], parameter = 0.5)
    #cat("------------------------------------------------------------------------\n")
  }
  if (!group) {
    comb <- combn(ntr, 2)
    nn <- ncol(comb)
    dif <- rep(0, nn)
    pvalue <- rep(0, nn)
    for (k in 1:nn) {
      i <- comb[1, k]
      j <- comb[2, k]
      dif[k] <- abs(means[i, 2] - means[j, 2])
      sdtdif <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,
                                                        4]))
      pvalue[k] <- round(1 - ptukey(dif[k] * sqrt(2)/sdtdif,
                                    ntr, DFerror), 4)
    }
    tr.i <- comb[1, ]
    tr.j <- comb[2, ]
    output <- data.frame(trt = means[, 1], means = means[,
                                                         2], M = "", N = means[, 4], std.err = means[, 3])
  }
  return(output[,1:3])
}

#######################################################
####################################################
##########-> Teste t protegido
###################################################
TesteTProtegido=function (y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL){
  SSerror <- MSerror*DFerror
  name.y <- paste(deparse(substitute(y)))
  name.t <- paste(deparse(substitute(trt)))
  junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
  means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
  sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
  nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
  means <- data.frame(means, std.err = sds[, 2]/sqrt(nn[, 2]),
                      replication = nn[, 2])
  names(means)[1:2] <- c(name.t, name.y)
  ntr <- nrow(means)
  alphap = (2 * alpha)/(ntr * (ntr - 1))
  Tprob <- qt(1 - (alphap/2), DFerror) * sqrt(2)
  nr <- unique(nn[, 2])
  nfila <- c("Alpha", "Error Degrees of Freedom", "Error Mean Square",
             "Critical Value of Studentized Range")
  nvalor <- c(alpha, DFerror, MSerror, Tprob)
  xtabla <- data.frame(...... = nvalor)
  row.names(xtabla) <- nfila
  if (group) {
    if (length(nr) == 1) {
      HSD <- Tprob * sqrt(MSerror/nr)
    }
    else {
      nr1 <- 1/mean(1/nn[, 2])
      HSD <- Tprob * sqrt(MSerror/nr1)
    }
    #cat("\nTeste t de Bonferroni (LSD protegido)\n------------------------------------------------------------------------")
    #cat("\nGrupos  Tratamentos  Medias\n")
    output <- order.group(means[, 1], means[, 2], means[,
                                                        4], MSerror, Tprob, means[, 3], parameter = 0.5)
    #cat("------------------------------------------------------------------------\n")
  }
  if (!group) {
    comb <- combn(ntr, 2)
    nn <- ncol(comb)
    dif <- rep(0, nn)
    pvalue <- rep(0, nn)
    for (k in 1:nn) {
      i <- comb[1, k]
      j <- comb[2, k]
      dif[k] <- abs(means[i, 2] - means[j, 2])
      sdtdif <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,
                                                        4]))
      pvalue[k] <- round(1 - ptukey(dif[k] * sqrt(2)/sdtdif,
                                    ntr, DFerror), 4)
    }
    tr.i <- comb[1, ]
    tr.j <- comb[2, ]
    output <- data.frame(trt = means[, 1], means = means[,
                                                         2], M = "", N = means[, 4], std.err = means[, 3])
  }
  return(output[,1:3])
}

###############################################################
####################################################
##########-> Duncan
###################################################
Duncan=function (y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL){
  SSerror <- MSerror*DFerror
  name.y <- paste(deparse(substitute(y)))
  name.t <- paste(deparse(substitute(trt)))
  junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
  means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
  sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
  nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
  means <- data.frame(means, std.err = sds[, 2]/sqrt(nn[, 2]),
                      replication = nn[, 2])
  names(means)[1:2] <- c(name.t, name.y)
  ntr <- nrow(means)
  k.snk <- ntr - 1
  Tprob <- vector(mode = "integer", k.snk)
  kk <- 1
  for (kk in 1:k.snk) {
    alphap = 1 - (1 - alpha)^((kk + 1) - 1)
    xxxx=qtukey(1 - alphap, kk + 1, DFerror)
    if(is.na(xxxx)){xxxx=4}
    Tprob[kk] <- xxxx
  }
  p.nan <- as.vector(na.action(na.omit(Tprob)))[1]
  ult <- p.nan - 1
  if (ntr == 50)
    Tprob[p.nan:length(Tprob)] <- seq(Tprob[ult], 3.61, length = length(Tprob) -
                                        ult)
  if (ntr == 100)
    Tprob[p.nan:length(Tprob)] <- seq(Tprob[ult], 3.67, length = length(Tprob) -
                                        ult)
  nr <- unique(nn[, 2])
  nfila <- c("Alpha", "Error Degrees of Freedom", "Error Mean Square")
  nfila1 <- c("Distances between averages", "Critical Value of Studentized Range")
  nvalor <- c(alpha, DFerror, MSerror)
  nvalor1 <- rbind(t(seq(2, ntr)), t(Tprob))
  xtabla <- data.frame(...... = nvalor)
  xtabla1 <- data.frame(...... = nvalor1)
  row.names(xtabla) <- nfila
  row.names(xtabla1) <- nfila1
  HSD <- vector(mode = "integer", k.snk)
  if (group) {
    if (length(nr) == 1) {
      kk <- 1
      for (kk in 1:k.snk) {
        HSD[kk] <- Tprob[kk] * sqrt(MSerror/nr)
      }
    }
    else {
      nr1 <- 1/mean(1/nn[, 2])
      kk <- 1
      for (kk in 1:k.snk) {
        HSD[kk] <- Tprob[kk] * sqrt(MSerror/nr1)
      }
    }
    #cat("\nTeste de Duncan \n------------------------------------------------------------------------")
    #cat("\nGrupos  Tratamentos  Medias\n")
    output <- order.stat.SNK(means[, 1], means[, 2], HSD)
    #cat("------------------------------------------------------------------------\n")
  }
  return(output)
}

####################################################
##########-> SNK
###################################################
SNK=function (y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL){
  SSerror <- MSerror*DFerror
  name.y <- paste(deparse(substitute(y)))
  name.t <- paste(deparse(substitute(trt)))
  junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
  means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
  sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
  nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
  means <- data.frame(means, std.err = sds[, 2]/sqrt(nn[, 2]),
                      replication = nn[, 2])
  names(means)[1:2] <- c(name.t, name.y)
  ntr <- nrow(means)
  k.snk <- ntr - 1
  Tprob <- vector(mode = "integer", k.snk)
  kk <- 1
  for (kk in 1:k.snk) {
    Tprob[kk] <- qtukey(1 - alpha, kk + 1, DFerror)
  }
  nr <- unique(nn[, 2])
  nfila <- c("Alpha", "Error Degrees of Freedom", "Error Mean Square")
  nfila1 <- c("Distances between averages", "Critical Value of Studentized Range")
  nvalor <- c(alpha, DFerror, MSerror)
  nvalor1 <- rbind(t(seq(2, ntr)), t(Tprob))
  xtabla <- data.frame(...... = nvalor)
  xtabla1 <- data.frame(...... = nvalor1)
  row.names(xtabla) <- nfila
  row.names(xtabla1) <- nfila1
  HSD <- vector(mode = "integer", k.snk)
  if (group) {
    if (length(nr) == 1) {
      kk <- 1
      for (kk in 1:k.snk) {
        HSD[kk] <- Tprob[kk] * sqrt(MSerror/nr)
      }
    }
    else {
      nr1 <- 1/mean(1/nn[, 2])
      kk <- 1
      for (kk in 1:k.snk) {
        HSD[kk] <- Tprob[kk] * sqrt(MSerror/nr1)
      }
    }
    #cat("\nTeste de Student-Newman-Keuls (SNK)\n------------------------------------------------------------------------")
    #cat("\nGrupos  Tratamentos  Medias\n")
    output <- order.stat.SNK(means[, 1], means[, 2], HSD)
    #cat("------------------------------------------------------------------------\n")
  }
  return(output[,1:3])
}

####################################################
##########-> Tukey
###################################################
tukey=function (y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL){
  SSerror <- MSerror*DFerror
  name.y <- paste(deparse(substitute(y)))
  name.t <- paste(deparse(substitute(trt)))
  junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
  means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
  sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
  nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
  means <- data.frame(means, std.err = sds[, 2]/sqrt(nn[, 2]),
                      replication = nn[, 2])
  names(means)[1:2] <- c(name.t, name.y)
  ntr <- nrow(means)
  Tprob <- qtukey(1 - alpha, ntr, DFerror)
  nr <- unique(nn[, 2])
  nfila <- c("Alpha", "Error Degrees of Freedom", "Error Mean Square",
             "Critical Value of Studentized Range")
  nvalor <- c(alpha, DFerror, MSerror, Tprob)
  xtabla <- data.frame(...... = nvalor)
  row.names(xtabla) <- nfila
  if (group) {
    if (length(nr) == 1) {
      HSD <- Tprob * sqrt(MSerror/nr)
    }
    else {
      nr1 <- 1/mean(1/nn[, 2])
      HSD <- Tprob * sqrt(MSerror/nr1)
    }
    #cat("\nTeste de Tukey\n------------------------------------------------------------------------")
    #cat("\nGrupos Tratamentos Medias\n")
    output <- order.group(means[, 1], means[, 2], means[,
                                                        4], MSerror, Tprob, means[, 3], parameter = 0.5)
    #cat("------------------------------------------------------------------------\n")
  }
  if (!group) {
    comb <- combn(ntr, 2)
    nn <- ncol(comb)
    dif <- rep(0, nn)
    pvalue <- rep(0, nn)
    for (k in 1:nn) {
      i <- comb[1, k]
      j <- comb[2, k]
      dif[k] <- abs(means[i, 2] - means[j, 2])
      sdtdif <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,
                                                        4]))
      pvalue[k] <- round(1 - ptukey(dif[k] * sqrt(2)/sdtdif,
                                    ntr, DFerror), 4)
    }
    tr.i <- comb[1, ]
    tr.j <- comb[2, ]
    output <- data.frame(trt = means[, 1], means = means[,
                                                         2], M = "", N = means[, 4], std.err = means[, 3])
  }
  return(output[,1:3])
}

#########################################
##########->ScottKnott
#########################################
ScottKnott= function (y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,
                      main = NULL)
{
  SSerror=MSerror*DFerror
  sk <- function(medias, s2, dfr, prob) {
    bo <- 0
    si2 <- s2
    defr <- dfr
    parou <- 1
    np <- length(medias) - 1
    for (i in 1:np) {
      g1 <- medias[1:i]
      g2 <- medias[(i + 1):length(medias)]
      B0 <- sum(g1)^2/length(g1) + sum(g2)^2/length(g2) -
        (sum(g1) + sum(g2))^2/length(c(g1, g2))
      if (B0 > bo) {
        bo <- B0
        parou <- i
      }
    }
    g1 <- medias[1:parou]
    g2 <- medias[(parou + 1):length(medias)]
    teste <- c(g1, g2)
    sigm2 <- (sum(teste^2) - sum(teste)^2/length(teste) +
                defr * si2)/(length(teste) + defr)
    lamb <- pi * bo/(2 * sigm2 * (pi - 2))
    v0 <- length(teste)/(pi - 2)
    p <- pchisq(lamb, v0, lower.tail = FALSE)
    if (p < prob) {
      for (i in 1:length(g1)) {
        cat(names(g1[i]), "\n", file = "skresult", append = TRUE)
      }
      cat("*", "\n", file = "skresult", append = TRUE)
    }
    if (length(g1) > 1) {
      sk(g1, s2, dfr, prob)
    }
    if (length(g2) > 1) {
      sk(g2, s2, dfr, prob)
    }
  }
  medias <- sort(tapply(y, trt, mean), decreasing = TRUE)
  dfr <- DFerror
  rep <- tapply(y, trt, length)
  s0 <- MSerror <- SSerror/DFerror
  s2 <- s0/rep[1]
  # cat("\nTeste de Scott-Knott\n------------------------------------------------------------------------\n")
  prob <- alpha
  sk(medias, s2, dfr, prob)
  f <- names(medias)
  names(medias) <- 1:length(medias)
  resultado <- data.frame(r = 0, f = f, m = medias)
  if (file.exists("skresult") == FALSE) {
    stop
  }
  else {
    xx <- read.table("skresult")
    file.remove("skresult")
    x <- xx[[1]]
    x <- as.vector(x)
    z <- 1
    for (j in 1:length(x)) {
      if (x[j] == "*") {
        z <- z + 1
      }
      for (i in 1:length(resultado$f)) {
        if (resultado$f[i] == x[j]) {
          resultado$r[i] <- z
        }
      }
    }
  }
  letras <- letters
  if (length(resultado$r) > 26) {
    l <- floor(length(resultado$r)/26)
    for (i in 1:l) letras <- c(letras, paste(letters, i,
                                             sep = ""))
  }
  res <- 1
  for (i in 1:(length(resultado$r) - 1)) {
    if (resultado$r[i] != resultado$r[i + 1]) {
      resultado$r[i] <- letras[res]
      res <- res + 1
      if (i == (length(resultado$r) - 1)) {
        resultado$r[i + 1] <- letras[res]
      }
    }
    else {
      resultado$r[i] <- letras[res]
      if (i == (length(resultado$r) - 1)) {
        resultado$r[i + 1] <- letras[res]
      }
    }
  }
  names(resultado) <- c("Grupos", "Tratamentos", "Medias")
  return(resultado[,c(2,3,1)])
  #cat("------------------------------------------------------------------------\n")
}

Resultado=cbind(
TesteT(y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL)[,1:2],
TesteT=TesteT(y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL)[,3],
TesteT_Bonferroni=TesteTProtegido(y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL)[,3],
Duncan=Duncan(y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE,main = NULL)[,3],
SNK=SNK(y, trt, DFerror, MSerror,alpha = 0.05, group = TRUE,main = NULL)[,3],
Tukey=tukey(y, trt, DFerror, MSerror,alpha = 0.05, group = TRUE,main = NULL)[,3],
ScottKnott=ScottKnott(y, trt, DFerror, MSerror, alpha = 0.05)[,3] )

return(Resultado)}


RegressaoPolinomial=function (resp, trat, glres, SQres, gltrat, SQtrat,verbose=T) {

  ginv<-function(X, tol = sqrt(.Machine$double.eps))
  {
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
      stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X))
      X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X))
      Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive))
      Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive))
      array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
                                                 t(Xsvd$u[, Positive, drop = FALSE]))
  }

  tapply.stat=function (y, x, stat = "mean") {
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
  mean.table <- tapply.stat(resp, trat, mean)
  colnames(mean.table) <- c("  Niveis", "   Medias Observadas")
  if(verbose){print(mean.table)}

 QMres <- SQres/glres
 if(verbose) {cat("\nAjuste de modelos polinomiais de regressao\n------------------------------------------------------------------------\n")}
  X <- matrix(1, length(trat), 4)
  X[, 2] <- trat
  X[, 3] <- trat^2
  X[, 4] <- trat^3

  b = ginv(t(X[, 1:2]) %*% X[, 1:2], tol = .Machine$double.eps) %*%
    t(X[, 1:2]) %*% resp
  ep = sqrt(diag(ginv(t(X[, 1:2]) %*% X[, 1:2], tol = .Machine$double.eps) *
                   QMres))
  tc = b/ep
  pv = 2 * pt(abs(tc), glres, lower.tail = FALSE)
  tm1 <- data.frame(Estimativa = round(b, 8), `Erro padrao` = round(ep,
                                                                    5), tc = round(tc, 5), `valor-p` = round(pv, 5))
  rownames(tm1) <- c("b0", "b1")
  aov.m1 <- anova(lm(resp ~ trat))
  if (dim(mean.table)[1] == 2) {
    r2m1 <- 1
  }
  if (dim(mean.table)[1] > 2) {
    r2m1 <- aov.m1[1, 2]/SQtrat
  }
  nomes1 <- c("Efeito linear", "Desvios de Regressao",
              "Residuos")
  anava1 <- data.frame(GL = c(1, (gld = c(gltrat - 1)), (glr = glres)),
                       SQ = c(round(c(aov.m1[[2]][1], (sqd = c(SQtrat - aov.m1[[2]][1])),
                                      SQres), 5)), QM = c(round(c(aov.m1[[3]][1], (if (gld ==
                                                                                       0) {
                                        qmd = 0
                                      } else {
                                        qmd = (sqd/gld)
                                      }), QMres), 5)), Fc = c(round(c((fcl = aov.m1[[3]][1]/QMres),
                                                                      (fc = qmd/QMres)), 2), ""), `valor-p` = c(round(c(pf(fcl,
                                                                                                                           1, glr, lower.tail = FALSE), (if (gld == 0) {
                                                                                                                             pv = 1
                                                                                                                           } else {
                                                                                                                             pv = 1 - pf(fc, gld, glr)
                                                                                                                           })), 5), ""))
  rownames(anava1) <- nomes1
  if(verbose){print("Modelo linear --------------------------------------------")}
    if(verbose){print(tm1)}
      if(verbose){print(" ")}
        if(verbose){print(paste("R2 do modelo linear =", r2m1))}
          if(verbose){print("Analise de variancia da regressao linear ")}

            if(verbose){print(anava1)}


  if(verbose){print("------------------------------------------------------------------------")}
  if (dim(mean.table)[1] > 2) {
    b2 = ginv(t(X[, 1:3]) %*% X[, 1:3], tol = .Machine$double.eps) %*%
      t(X[, 1:3]) %*% resp
    ep2 = sqrt(diag(ginv(t(X[, 1:3]) %*% X[, 1:3], tol = .Machine$double.eps) *
                      QMres))
    tc2 = b2/ep2
    pv2 = 2 * pt(abs(tc2), glres, lower.tail = FALSE)
    tm2 <- data.frame(Estimativa = round(b2, 8), `Erro padrao` = round(ep2,
                                                                       5), tc = round(tc2, 5), `valor-p` = round(pv2,
                                                                                                                 5))
    rownames(tm2) <- c("b0", "b1", "b2")
    t2 <- trat^2
    aov.m2 <- anova(lm(resp ~ trat + t2))
    if (dim(mean.table)[1] == 3) {
      r2m2 <- 1
    }
    if (dim(mean.table)[1] > 3) {
      r2m2 <- (aov.m2[1, 2] + aov.m2[2, 2])/SQtrat
    }
    nomes2 <- c("Efeito linear", "Efeito quadratico",
                "Desvios de Regressao", "Residuos")
    anava2 <- data.frame(GL = c(aov.m2[[1]][1:2], (gld = c(gltrat -
                                                             2)), (glr = glres)), SQ = c(round(c(aov.m2[[2]][1:2],
                                                                                                 (sqd = c(SQtrat - sum(aov.m2[[2]][1:2]))), SQres),
                                                                                               5)), QM = c(round(c(aov.m2[[3]][1:2], (if (gld ==
                                                                                                                                          0) {
                                                                                                 qmd = 0
                                                                                               } else {
                                                                                                 qmd = (sqd/gld)
                                                                                               }), QMres), 5)), Fc = c(round(c((fcl = aov.m2[[3]][1:2]/QMres),
                                                                                                                               (fc = qmd/QMres)), 2), ""), `valor-p` = c(round(c(pf(fcl,
                                                                                                                                                                                    1, glr, lower.tail = FALSE), (if (gld == 0) {
                                                                                                                                                                                      pv = 1
                                                                                                                                                                                    } else {
                                                                                                                                                                                      pv = 1 - pf(fc, gld, glr)
                                                                                                                                                                                    })), 5), ""))
    rownames(anava2) <- nomes2
    if(verbose){print("Modelo quadratico --------------------------------------------")}
      if(verbose){print(tm2)}
        if(verbose){ print(" ")}
          if(verbose){ print(paste("R2 do modelo quadratico =", round(r2m2,6))
        , `Analise de variancia do modelo quadratico` = anava2)}


    if(verbose){ print("Analise de variancia da regressao  quadratica")}
      if(verbose){print(anava2)}
        if(verbose){print("------------------------------------------------------------------------")}
  }
  if (dim(mean.table)[1] > 3) {
    b3 = ginv(t(X[, 1:4]) %*% X[, 1:4], tol = .Machine$double.eps) %*%
      t(X[, 1:4]) %*% resp
    ep3 = sqrt(diag(ginv(t(X[, 1:4]) %*% X[, 1:4], tol = .Machine$double.eps) *
                      QMres))
    tc3 = b3/ep3
    pv3 = 2 * pt(abs(tc3), glres, lower.tail = FALSE)
    tm3 <- data.frame(Estimativa = round(b3, 8), `Erro padrao` = round(ep3,
                                                                       5), tc = round(tc3, 5), `valor-p` = round(pv3,
                                                                                                                 5))
    rownames(tm3) <- c("b0", "b1", "b2",
                       "b3")
    t3 <- trat^3
    aov.m3 <- anova(lm(resp ~ trat + t2 + t3))
    if (dim(mean.table)[1] == 4) {
      r2m3 <- 1
    }
    if (dim(mean.table)[1] > 4) {
      r2m3 <- (aov.m3[1, 2] + aov.m3[2, 2] + aov.m3[3,
                                                    2])/SQtrat
    }
    nomes3 <- c("Efeito linear", "Efeito quadratico",
                "Efeito cubico", "Desvios de Regressao",
                "Residuos")
    anava3 <- data.frame(GL = c(aov.m3[[1]][1:3], (gld = c(gltrat -
                                                             3)), (glr = glres)), SQ = c(round(c(aov.m3[[2]][1:3],
                                                                                                 (sqd = c(SQtrat - sum(aov.m3[[2]][1:3]))), SQres),
                                                                                               5)), QM = c(round(c(aov.m3[[3]][1:3], (if (gld ==
                                                                                                                                          0) {
                                                                                                 qmd = 0
                                                                                               } else {
                                                                                                 qmd = (sqd/gld)
                                                                                               }), QMres), 5)), Fc = c(round(c((fcl = aov.m3[[3]][1:3]/QMres),
                                                                                                                               (fc = qmd/QMres)), 2), ""), `valor-p` = c(round(c(pf(fcl,
                                                                                                                                                                                    1, glr, lower.tail = FALSE), (if (gld == 0) {
                                                                                                                                                                                      pv = 1
                                                                                                                                                                                    } else {
                                                                                                                                                                                      pv = 1 - pf(fc, gld, glr)
                                                                                                                                                                                    })), 5), ""))
    rownames(anava3) <- nomes3
    if(verbose){print("Modelo cubico --------------------------------------------")}
    if(verbose){print(tm3)}
    if(verbose){print(" ")}
    if(verbose){print("`R2 do modelo cubico` =")}
    if(verbose){print(r2m3)}
    if(verbose){print("Analise de variancia da regressao cubica")}
    if(verbose){print(anava3)}

      if(verbose){print("-----------------------------------------------------------------")}
  }
  if (dim(mean.table)[1] > 3) {
    return(list(`Quadro de medias` = mean.table, `Coeficientes reta` = b,
                `R2 reta` = r2m1, `Coeficientes parabola` = b2,
                `R2 parabola` = r2m2, `Coeficientes cubica` = b3,
                `R2 cubica` = r2m3))
  }
  if (dim(mean.table)[1] == 3) {
    return(list(`Quadro de medias` = mean.table, `Coeficientes reta` = b,
                `R2 reta` = r2m1, `Coeficientes parabola` = b2,
                `R2 parabola` = r2m2))
  }
  if (dim(mean.table)[1] < 3) {
    return(list(`Quadro de medias` = mean.table, `Coeficientes reta` = b,
                `R2 reta` = r2m1))
  }
  if(verbose){print("------------------------------------------------------------------------")}
}


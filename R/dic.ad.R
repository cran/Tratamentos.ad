#' Analise de experimento  conduzido no delineamento inteiramente casualizado
#' com testemunhas adicionais
#' @description Esta funcao retorna a comparacao multipla de medias (obtidas
#'  pelos testes t, t com protecao de Bonferroni, Duncan, Dunnet, SNK, Tukey e
#'  Scott-Knott) se os tratamentos for qualitativos. Ou a analise de regressao
#'  se os tratamentos forem quantitativos. Para comparar a testemunha adicional
#'  com os demais e utilizado o teste Dunnet. Esta funcao considera o
#'  delineamento inteiramente casualizado. "This function returns the multiple
#'  comparison tests (t, t tests with Bonferroni protection, Duncan, Dunnet,
#'  SNK, Tukey and Scott-Knott ) if the treatments are qualitative. Or
#'  regression analysis if treatments are quantitative. To compare the
#'  additional control with the others treatments, the Dunnet test is used. This
#'  function considers the completely randomized design."
#' @usage dic.ad(Dados,alfa=0.05,quali=TRUE,verbose=TRUE,plot=2)
#' @param  Dados Matriz contendo na primeira coluna a identificacao dos
#'   tratamentos testemunhas (tratamentos comuns deve ter valor zero ou NA). A
#'   segunda coluna deve ter a identificacao de todos os tratamentos.  A
#'   terceira coluna a identificacao das repeticoes. Na quarta coluna os
#'   resultados experimentais.
#' @param quali Valor logico (TRUE/FALSE). TRUE indica que o Tratamento e
#'   qualitativo, realizando-se o teste de medias. FALSE indica que o fator e
#'   quantitativo, sendo feita a analise de regressao.
#'@param  verbose Valor logico (TRUE/FALSE). TRUE apresenta os resultados da analise.
#'@param  plot Valor numerico indicando o grafico desejado para analise dos residuos:
#'   \itemize{
#'    \item 1: Residuals vs Fitted
#'   \item 2:  QQ-plot
#'    \item 3:  Scale-Location
#'     \item 4:  Cook's distance
#'      \item 5: Histogram
#'      }
#' @importFrom stats anova lm na.action na.omit pchisq  pf pt ptukey qt qtukey
#' @importFrom utils combn read.table
#' @importFrom stats aov bartlett.test residuals shapiro.test
#' @importFrom graphics hist
#' @param  alfa valor indicando o nivel de significancia. Deve ser
#'   obrigatoriamente "0.001", "0.01", "0.05" ou "0.10" (default = 0.05).
#' @return Retorna a comparacao multipla de medias obtida por varios testes.
#' @author Alcinei Mistico Azevedo, \email{alcineimistico@@hotmail.com}
#'@references
#' <https://www.youtube.com/playlist?list=PLvth1ZcREyK4wSzwg-IxvrzaNzSLLrXEB>
#'
#'  BANZATTO, D. A.; KRONKA, S. N. Experimentacao Agricola. 4 ed.
#'  Jaboticabal: Funep. 2006. 237 p.
#'
#'  GOMES, F. P. Curso de Estatistica Experimental. 10a ed. Piracicaba:
#'  ESALQ/USP. 1982. 430.
#' @examples
#'  ######
#'  #Exemplo de um experimento em DIC com tratamentos qualitativos e uma
#'  #testemunha adicional
#'  data(Dados1)
#'  dic.ad(Dados = Dados1,alfa = 0.05,quali =TRUE)
#'  #Exemplo de um experimento em DIC com tratamentos quantitativos e tres
#'  #testemunhas adicionais
#' data(Dados2)
#' dic.ad(Dados = Dados2,alfa = 0.05,quali =FALSE)
#' @export
dic.ad=function (Dados, alfa=0.05,quali=TRUE,verbose=TRUE,plot=2){
  pp=plot
  Dados[,1]=as.character(Dados[,1])
  Dados[is.na(Dados[,1]),1]=0
  Dados[,2]=as.character(Dados[,2])
  Dados[is.na(Dados[,2]),2]=0
  DadosTrat3=Dados

    #require(DunnettTests)
  DadosTrat = Dados[Dados[,1]==0,-1]
  DadosTest = Dados[Dados[,1]!=0,-1]

  DadosTest[,1]=Dados[Dados[,1]!=0,1]

  sq = function(Fator, Y) {
        X = as.factor(Fator)
        NumTrat = length(unique(as.factor(Fator)))
        NomeTrat = unique(as.factor(Fator))
        NumParc = length(X)
        Matrizout = matrix(0, ncol = NumTrat, nrow = NumParc)
        for (i in 1:NumParc) {
            for (j in 1:NumTrat) {
                if (X[i] == NomeTrat[j]) {
                  Matrizout[i, j] = 1
                }
            }
        }
        X = Matrizout
        return(c(GL = NumTrat - 1, SQ = sum((X %*% solve(t(X) %*%
            X) %*% t(X) %*% Y - mean(Y))^2)))
    }
    sqI2 = function(Fator1, Fator2, Y) {
        Fator1 = as.factor(Fator1)
        Fator2 = as.factor(Fator2)
        X = paste(Fator1, Fator2)
        sqA = sq(Fator1, Y)
        sqB = sq(Fator2, Y)
        sqA.B = sq(X, Y)
        return(sqA.B - sqA - sqB)
    }


    i=3
        DadosTrat2 = cbind(DadosTrat[, 1:2], Y = DadosTrat[,
            i])
        DadosTest2 = cbind(DadosTest[, 1:2], Y = DadosTest[,
            i])
        if (length(unique(DadosTest[, 1])) == 1) {
            anova1 = function(DadosTrat2, DadosTest2) {
                DadosTrat = DadosTrat2
                DadosTest = DadosTest2
                DadosTest[,1]=paste("T",DadosTest[,1],sep="_")
                DT = rbind(DadosTrat, DadosTest)
                Trat = sq(DT[, 1], DT[, 3])
                Comum = sq(DadosTrat[, 1], DadosTrat[, 3])
                DI = rbind(cbind(I = 1, DadosTrat), cbind(I = 2,
                  DadosTest))
                Comum_vs_Test = sq(DI[, 1], DI[, 4])
                #Bloco = sq(DT[, 2], DT[, 3])
                Total = sq(1:nrow(DT), DT[, 3])
                Residuo = Total - Trat# - Bloco
                anova = rbind(Trat, Comum, Comum_vs_Test,
                  Residuo, Total)
                QM = anova[, 2]/anova[, 1]
                Fc = QM/QM[length(QM) - 1]
                QMR = QM[length(QM) - 1]
                GLR = anova[length(QM) - 1, 1]
                pValor = round(1 - pf(Fc, anova[, 1], anova[length(QM) -
                  1, 1]), 4)
                pValor[Fc < 1] = 1
                QM[length(QM)] = ""
                Fc[(length(Fc) - 1):length(Fc)] = ""
                pValor[(length(pValor) - 1):length(pValor)] = ""
                anova = cbind(anova, QM, Fc, pValor)

                if(verbose){print("Analise de variancia")}
                  if(verbose){ print(anova)}
       CV=100 * sqrt(QMR)/mean(DT[, 3])
       if(verbose){print(paste("CV= ",round(CV,4)) )}

       TT=paste(Dados[,1],Dados[,2])
       m=aov(Dados[,4]~TT)

       if(pp<5){(plot(m,pp))}
       if(pp==5){(hist(residuals(m)))}
       if(verbose){print("")}
       if(verbose){print("Teste de normalidade")}
       ST=shapiro.test(residuals(m))
       if(verbose){print(ST)}
       BT=bartlett.test(residuals(m),TT)
       if(verbose){print("")}
       if(verbose){print("Teste de homogeneidade de variancias")}
       if(verbose){print(BT)}





         if(verbose){print("")}
           if(verbose){print("Testes")}
                if(quali==T){
                  TesteM=ComparacaoMedias(DadosTrat2[,3],trt =DadosTrat2[,1],
                          DFerror = GLR,MSerror = QMR,alpha = as.numeric(alfa) )
                  if(verbose){print(TesteM)}
                }

                if(quali==F){
                  TesteM=RegressaoPolinomial(resp = DadosTrat2[,3],trat = as.numeric(as.character(DadosTrat2[,1])),
                                             glres =GLR,SQres =GLR*QMR,gltrat = Comum[1],
                                             SQtrat = Comum[2],verbose=verbose )
                }

                ntratComum = length(unique(DadosTrat[, 1]))
                nrep = length(unique(DadosTrat[, 2]))


                #q1=dunnet1[Residuo[1],ntratComum]
                q1=DunnetCritic(Residuo[1],ntratComum,alfa)
                DMS1=sqrt(QMR/nrep+QMR/nrep)*q1


                Test = tapply(DadosTest[, 3], as.factor(as.matrix(DadosTest[,
                  1])), mean, na = T)
                Comum = tapply(DadosTrat[, 3], as.factor(as.matrix(DadosTrat[,
                  1])), mean, na = T)
                Dunnet = NULL
                for (it in 1:length(unique(DadosTest[, 1]))) {
                  Comum2 = Comum

                  if(sum(abs(Comum - Test[it]) > DMS1)>0){
                  Comum2[abs(Comum - Test[it]) > DMS1] = paste(Comum2[abs(Comum -
                    Test[it]) > DMS1], "*", sep = "")
                  }

                  if(sum(abs(Comum - Test[it]) < DMS1)>0){
                    Comum2[abs(Comum - Test[it]) < DMS1] = paste(Comum2[abs(Comum -
                                                                              Test[it]) < DMS1], "ns", sep = "")
                  }
                  Comum2 = c(Comum2, Testemunha = Test[it])
                  NomeTest = names(Test)[it]
                  Dunnet = c(Dunnet, list(TesteDunnett = Comum2))
                }
                #Teste Dunnet
                if(verbose){print(Dunnet)}
                #par = c(MediaGeral = mean(DT[, 3]), MediaTratComum = mean(DadosTrat[,
                 # 3]), MediaTest = mean(DadosTest[, 3]), QM = QMR,
                 # CV = 100 * sqrt(QMR)/mean(DT[, 3]), `DMS` = DMS1)
                Resultado = list(Anova = anova, Testes=TesteM, Dunnett = Dunnet)

                return(Resultado)
            }
            Resultado = anova1(DadosTrat2, DadosTest2)
        }




        if (length(unique(DadosTest[, 1])) > 1) {
            anova2 = function(DadosTrat2, DadosTest2) {
                DadosTrat = DadosTrat2
                DadosTest = DadosTest2
                DT = rbind(DadosTrat, DadosTest)
                Trat = sq(DT[, 1], DT[, 3])
                Comum = sq(DadosTrat[, 1], DadosTrat[, 3])
                Testemunha = sq(DadosTest[, 1], DadosTest[, 3])
                DI = rbind(cbind(I = 1, DadosTrat), cbind(I = 2,
                  DadosTest))
                Comum_vs_Test = sq(DI[, 1], DI[, 4])
                #Bloco = sq(DT[, 2], DT[, 3])
                Total = sq(1:nrow(DT), DT[, 3])
                Residuo = Total - Trat# - Bloco
                anova = rbind(Trat, Comum, Testemunha, Comum_vs_Test,
                   Residuo, Total)
                QM = anova[, 2]/anova[, 1]
                Fc = QM/QM[length(QM) - 1]
                QMR = QM[length(QM) - 1]
                GLR = anova[length(QM) - 1, 1]
                pValor = round(1 - pf(Fc, anova[, 1], anova[length(QM) -
                  1, 1]), 5)
                pValor[Fc < 1] = 1
                QM[length(QM)] = ""
                Fc[(length(Fc) - 1):length(Fc)] = ""
                pValor[(length(pValor) - 1):length(pValor)] = ""
                anova = cbind(anova, QM, Fc, pValor)

                if(verbose){print("Analise de variancia")}
                if(verbose){print(anova)}
                CV=100 * sqrt(QMR)/mean(DT[, 3])
                if(verbose){print(paste("CV= ",round(CV,4)) )}

                TT=paste(Dados[,1],Dados[,2])
                m=aov(Dados[,4]~TT)
                if(pp<5){(plot(m,pp))}
                if(pp==5){(hist(residuals(m)))}
                if(verbose){print("")}
                if(verbose){print("Teste de normalidade")}
                ST=shapiro.test(residuals(m))
                if(verbose){print(ST)}
                BT=bartlett.test(residuals(m),TT)
                if(verbose){print("")}
                if(verbose){print("Teste de homogeneidade de variancias")}
                if(verbose){print(BT)}





                  if(verbose){print("")}
                    if(verbose){print("Testes")}
                if(quali==T){
                  TesteM=ComparacaoMedias(DadosTrat2[,3],trt =DadosTrat2[,1],
                                          DFerror = GLR,MSerror = QMR,alpha = as.numeric(alfa) )
                  if(verbose){print(TesteM)}
                }

                if(quali==F){
                  TesteM=RegressaoPolinomial(resp = DadosTrat2[,3],trat = as.numeric(as.character(DadosTrat2[,1])),
                                             glres =GLR,SQres =GLR*QMR,gltrat = Comum[1],
                                             SQtrat = Comum[2], verbose=verbose)
                }






                ntratComum = length(unique(DadosTrat[, 1]))
                nrep = length(unique(DadosTrat[, 2]))






                q1=DunnetCritic(Residuo[1],ntratComum,alfa)
                DMS5=sqrt(QMR/nrep+QMR/nrep)*q1

                Test = tapply(DadosTest[, 3], as.factor(as.matrix(DadosTest[,
                  1])), mean, na = T)
                Comum = tapply(DadosTrat[, 3], as.factor(as.matrix(DadosTrat[,
                  1])), mean, na = T)
                Dunnet = NULL
                for (it in 1:length(unique(DadosTest[, 1]))) {
                  Comum2 = Comum

                  if(sum(abs(Comum - Test[it]) > DMS5)>0){
                    Comum2[abs(Comum - Test[it]) > DMS5] = paste(Comum2[abs(Comum -
                                                                              Test[it]) > DMS5], "*", sep = "")
                  }

                  if(sum(abs(Comum - Test[it]) < DMS5)>0){
                    Comum2[abs(Comum - Test[it]) < DMS5] = paste(Comum2[abs(Comum -
                                                                              Test[it]) < DMS5], "ns", sep = "")
                  }




                  Comum2 = c(Comum2, Testemunha = Test[it])
                  NomeTest = names(Test)[it]
                  Dunnet = c(Dunnet, list(TesteDunnett = Comum2))
                }

                if(verbose){print(Dunnet)}
                #par = c(MediaGeral = mean(DT[, 3]), MediaTratComum = mean(DadosTrat[,
                 # 3]), MediaTest = mean(DadosTest[, 3]), QM = QMR,
                 # CV = 100 * sqrt(QMR)/mean(DT[, 3]),
                 # `DMS` = DMS5)
                Resultado = list(Anova = anova,  Dunnett = Dunnet)

                return(Resultado)
            }
            Resultado = anova2(DadosTrat2, DadosTest2)



        }
       # return(Resultado)
    }









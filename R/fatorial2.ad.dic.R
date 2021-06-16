#'Analise de experimento em esquema fatorial duplo com testemunhas em DIC
#'@description Esta funcao considera experimentos em esquema fatorial no
#'  delineamento inteiramente casualizado.Esta funcao retorna a comparacao
#'  multipla de medias (obtidas pelos testes t, t com protecao de Bonferroni,
#'  Duncan, Dunnet, SNK, Tukey e Scott-Knott) se o fator for qualitativo. Ou a
#'  analise de regressao se o fator for quantitativo. Para comparar a testemunha
#'  adicional com os demais e utilizado o teste Dunnet.  "This function
#'  considers experiments in factorial scheme for the completely randomized
#'  design. This function returns the multiple comparison tests ( t, t with
#'  Bonferroni protection, Duncan, Dunnet, SNK, Tukey and Scott-Knott) if the
#'  factor is qualitative. Or regression analysis if the factor is quantitative.
#'  To compare the additionals controls with the others treatments, the Dunnet
#'  test is used."
#'@usage fatorial2.ad.dic(Dados,Protegido=FALSE,alfa=0.05,quali=c(TRUE,TRUE),verbose=TRUE)
#'@param Dados Matriz contendo na primeira coluna a identificacao das
#'  testemunhas (tratamentos comuns deve ter valor "NA" ou zero). A segunda
#'  coluna deve ter a identificacao dos niveis do primeiro fator. A terceira
#'  coluna a identificacao dos niveis do segundo fator. A quarta coluna a
#'  identificacao das repeticoes. Na quinta coluna os resultados experimentais.
#'@param Protegido Valor logico. Se for FALSE os testes de medias serao feitos
#'  indendente das significancias pelo teste F. Se for TRUE os testes de medias
#'  serao feitas apenas quando haver efeito significativo pelo teste F.
#'@param alfa Valor indicando o nivel de significancia deve ser obrigatoriamente
#'  "0.001", "0.01", "0.05" ou "0.10" (default = 0.05).
#'@param quali Vetor com dois valores logicos (TRUE/FALSE). TRUE indica que o
#'  respectivo fator e qualitativo, realizando-se o teste de medias. FALSE
#'  indica que o fator e quantitativo, sendo feita a analise de regressao.
#'@param  verbose Valor logico (TRUE/FALSE). TRUE apresenta os resultados da analise.


#'@references BANZATTO, D. A.; KRONKA, S. N. Experimentacao Agricola. 4 ed.
#'  Jaboticabal: Funep. 2006. 237 p. ISBN: 85-87632-71-X
#'
#'  GOMES, F. P. Curso de Estatistica Experimental. 10a ed. Piracicaba:
#'  ESALQ/USP. 1982. 430.
#'@return Retorna a comparacao multipla de medias obtida por varios testes para
#'  tratamentos qualitativos e regressao para testes quantitativos. O teste
#'  Dunnet e feito para comparar os tratamentos testemunhas com os demais.
#' @examples

#' ##Exemplo de experimento com duas testemunhas adicionais e fatores qualitativos
#'
#' data(Dados3)
#' fatorial2.ad.dic(Dados3,Protegido=FALSE,alfa = 0.05,quali = c(TRUE,TRUE))
#'
#' ##Exemplo de experimento com uma testemunha adicional e um fator qualitativo
#' data(Dados4)
#' fatorial2.ad.dic(Dados4,Protegido=TRUE,alfa = 0.05,quali = c(FALSE,TRUE))
#'
#' ##Exemplo com tres testemunhas adicionais e um fator qualitativo
#' data(Dados5)
#' fatorial2.ad.dic(Dados5,Protegido=TRUE,alfa = 0.05,quali = c(TRUE,FALSE))
#' @export



fatorial2.ad.dic=function(Dados,Protegido=FALSE,alfa=0.05,quali=c(TRUE,TRUE),verbose=TRUE) {
######################################

sq=function(Fac,x){
  GL=length(unique(Fac))-1
  y=tapply(as.numeric(x),as.factor(Fac),mean,na=TRUE)
  r=mean(table(Fac))
  SQ=sum((y-mean(y))^2)*r
  QM=SQ/GL
  return(c(GL=GL,SQ=SQ,QM=QM))
}
sq2=function(FacA,FacB,x){
  total=sq(paste(FacA,FacB),x)
  FatorA=sq(FacA,x)
  FatorB=sq(FacB,x)
  GL=FatorA[1]*FatorB[1]
  SQ=total[2]-FatorA[2]-FatorB[2]
  QM=SQ/GL
  return(c(GL=GL,SQ=SQ,QM=QM))
}
var=Dados[,5]
for(i in 1:4){Dados[,i]=as.character(Dados[,i])}
Dados[is.na(Dados[,1]),1]=0
Dados[is.na(Dados[,2]),2]=0
Dados[is.na(Dados[,3]),3]=0

Resultado=NULL

id=(Dados[,1])==0
Bloco=sq(c(Dados[,4]),var)
FatorA=sq(c(Dados[id,2]),var[id])
FatorB=sq(c(Dados[id,3]),var[id])
AxB=sq2(c(Dados[id,2]),c(Dados[id,3]),var[id])
Parte1=rbind(FatorA,FatorB,AxB)

# Caso 1
if (length(unique(Dados[,1]))==2){
Test_Vs_Comuns=sq(paste(Dados[,1],Dados[,2],Dados[,3]),var)-colSums(Parte1)
Test_Vs_Comuns[3]=Test_Vs_Comuns[2]/Test_Vs_Comuns[1]
Parte2=rbind(Test_Vs_Comuns)
}

# Caso 2
if (length(unique(Dados[,1]))>2){
id2=Dados[,1]!=0
Testemunha=sq(c(Dados[id2,1]),Dados[id2,5])
Test_Vs_Comuns=sq(paste(Dados[,1],Dados[,2],Dados[,3]),var)-colSums(rbind(Parte1,Testemunha))
Test_Vs_Comuns[3]=Test_Vs_Comuns[2]/Test_Vs_Comuns[1]
Parte2=rbind(Testemunha,Test_Vs_Comuns)
}

Total=sq(1:nrow(Dados),var)
Residuo=Total-colSums(rbind(Parte1,Parte2))
Residuo[3]=Residuo[2]/Residuo[1]

ANOVA=data.frame(rbind(Parte1,Parte2,Residuo,Total))
ANOVA=cbind(ANOVA,Fc=ANOVA[,3]/Residuo[3])
ANOVA=round(cbind(ANOVA,pValor=1-pf(ANOVA$Fc,ANOVA$GL,Residuo[1])),5)
n=nrow(ANOVA)
ANOVA[n,3:5]=""
ANOVA[n-1,4:5]=""

if(verbose){print("Analise de variancia")}
if(verbose){print(ANOVA)}
if(verbose){print(paste("CV%=", round(100*sqrt(Residuo[3])/mean(var),4)))}

FAsig=as.numeric(ANOVA[1,5])<alfa
FBsig=as.numeric(ANOVA[2,5])<alfa
INTsig=as.numeric(ANOVA[3,5])<alfa

if(Protegido==F){INTsig=F;FAsig=T;FBsig=T}


if(INTsig==F){

  if(FAsig==T){
if(quali[1]==T){
if(verbose){print("Comparacao de medias dos niveis do fator A")}
if(verbose){print(ComparacaoMedias((var[id]), as.character(Dados[id,2]), Residuo[1], Residuo[3], alpha = alfa, group = TRUE,main = NULL))}
if(verbose){print("______________________________________________________")}
}
if(quali[1]==F){
  if(verbose){print("Analise de regressao para os niveis do fator A")}
  #   Dados[,2]=as.numeric(as.factor((Dados[,2])))

 # ComparacaoMedias((var[id]), as.character(Dados[id,2]), Residuo[1], Residuo[3], alpha = alfa, group = TRUE,main = NULL)
  Reg=RegressaoPolinomial((var[id]), as.numeric(as.character((Dados[id,2]))), glres=Residuo[1], Residuo[2], FatorA[1], FatorA[2],
                          verbose=verbose)
  if(verbose){print("______________________________________________________")}
}



  if(verbose){print("Teste Dunnet para comparar as testemunhas com os niveis de A")}
  ########################
  ntratComum=length(unique(Dados[id,2]))
  nrep1=mean(table(Dados[id,2]))
  nrep2=mean(table(Dados[id==0,1]))
  QMR=Residuo[3]
  NL=Residuo[1]


  #q1=  qNCDun(1-alfa,nu =NL,rho=0.5, delta=rep(0,times=ntratComum), two.sided=T)
  q1=DunnetCritic(gln = ntratComum,gld = NL,sig = alfa)

  DMS1=sqrt(QMR/nrep1+QMR/nrep2)*q1


  Test=round(tapply(Dados[id==0,5], as.factor(as.matrix(Dados[id==0,1])),mean,na=T),5)
  Comum=round(tapply(var[id], as.factor(as.matrix(Dados[id,2])),mean,na=T),5)


  Dunnet=NULL
  for (it in 1:length(Test)){
    Comum2=Comum
    Comum2[abs(Comum-Test[it])>DMS1]=paste(Comum2[abs(Comum-Test[it])>DMS1], "*",sep="")
    Comum2[abs(Comum-Test[it])<DMS1]=paste(Comum2[abs(Comum-Test[it])<DMS1], "ns",sep="")
    Comum2=c(Comum2,Testemunha=Test[it])
    NomeTest=names(Test)[it]
    Dunnet=c(Dunnet,list(TesteDunnett=Comum2))
  }
  if(verbose){print(Dunnet)}
  if(verbose){print("____________________________________________________________")}

  }

  if(FBsig==T){

if(quali[2]==T){
if(verbose){print("Comparacao de medias dos niveis do fator B")}
if(verbose){print(ComparacaoMedias((var[id]), as.character(Dados[id,3]), Residuo[1], Residuo[3], alpha = alfa, group = TRUE,main = NULL))}
if(verbose){print("______________________________________________________")}
}
if(quali[2]==F){
  if(verbose){print("Comparacao de medias dos niveis do fator B")}
  #   Dados[,3]=as.numeric(as.factor((Dados[,3])))
 # ComparacaoMedias((var[id]), as.character(Dados[id,3]), Residuo[1], Residuo[3], alpha = alfa, group = TRUE,main = NULL)
  reg=RegressaoPolinomial((var[id]), as.numeric(as.character(Dados[id,3])),
                          glres=Residuo[1], SQres=Residuo[2], gltrat =FatorB[1],
                          SQtrat=FatorB[2],verbose=verbose)
  if(verbose){print("______________________________________________________")}
}


  if(verbose){print("Teste Dunnet para comparar as testemunhas com os niveis de B")}
  ########################
  ntratComum=length(unique(Dados[id,3]))
  nrep1=mean(table(Dados[id,3]))
  nrep2=mean(table(Dados[id==0,1]))
  QMR=Residuo[3]
  NL=Residuo[1]
  #q1=  qNCDun(1-alfa,nu =NL,rho=0.5, delta=rep(0,times=ntratComum), two.sided=T)
  q1=DunnetCritic(gln = ntratComum,gld = NL,sig = alfa)

  DMS1=sqrt(QMR/nrep1+QMR/nrep2)*q1

  Test=round(tapply(Dados[id==0,5], as.factor(as.matrix(Dados[id==0,1])),mean,na=T),5)
  Comum=round(tapply(var[id], as.factor(as.matrix(Dados[id,3])),mean,na=T),5)


  Dunnet=NULL
  for (it in 1:length(Test)){
    Comum2=Comum
    Comum2[abs(Comum-Test[it])>DMS1]=paste(Comum2[abs(Comum-Test[it])>DMS1], "*",sep="")
    Comum2[abs(Comum-Test[it])<DMS1]=paste(Comum2[abs(Comum-Test[it])<DMS1], "ns",sep="")
    Comum2=c(Comum2,Testemunha=Test[it])
    NomeTest=names(Test)[it]
    Dunnet=c(Dunnet,list(TesteDunnett=Comum2))
  }
  if(verbose){print(Dunnet)}


}
}

if(Protegido==F){INTsig=T}
if(INTsig==T){
if(verbose){print("Desdobramento dos niveis de A dentro do Fator B")}
for(i in (unique(as.character(Dados[id,3])))){
  if(verbose){print(paste("Desdobramento dos niveis de A dentro do nivel",
                                i,"do fator B"))}
   id3=(id*(Dados[,3]==i))==1

   if(quali[1]==T){
   if(verbose){print(ComparacaoMedias(Dados[id3,5], as.character(Dados[id3,2]),
                Residuo[1], Residuo[3], alpha = alfa, group = TRUE,main = NULL))}
}
   if(quali[1]==F){
     A.B=sq((as.character((Dados[id3,2]))),Dados[id3,5])
     Reg=RegressaoPolinomial(Dados[id3,5], as.numeric(as.character((Dados[id3,2]))),
                          glres=Residuo[1], Residuo[2], A.B[1], A.B[2],verbose=verbose)

   }
    if(verbose){print("______________________________________________________")}
  }


if(verbose){print("Desdobramento dos niveis de B dentro do Fator A")}
for(i in (unique(as.character(Dados[id,2])))){
  if(verbose){print(paste("Desdobramento dos niveis de B dentro do nivel",i,"do fator A"))}
  id3=id*(Dados[,2]==i)==1

  if(quali[2]==T){
  if(verbose){print(ComparacaoMedias(Dados[id3,5], as.character(Dados[id3,3]),
                  Residuo[1], Residuo[3], alpha = alfa, group = TRUE,main = NULL))}
  }
  if(quali[2]==F){
    A.B=sq((as.character((Dados[id3,3]))),Dados[id3,5])
    Reg=RegressaoPolinomial(Dados[id3,5], as.numeric(as.character(Dados[id3,3])),
                            glres=Residuo[1], Residuo[2],  A.B[1], A.B[2],verbose=verbose)

  }


    if(verbose){print("______________________________________________________")}
}





if(verbose){print("Teste Dunnet para comparar as testemunhas todos os tratamentos da combinacao A x B")
########################
ntratComum=length(unique(paste(Dados[id,2],Dados[id,3])))
nrep1=mean(table(paste(Dados[id,2],Dados[id,3])))
nrep2=mean(table(Dados[id==0,1]))
QMR=Residuo[3]
NL=Residuo[1]
#q1=  qNCDun(1-alfa,nu =NL,rho=0.5, delta=rep(0,times=ntratComum), two.sided=T)
q1=DunnetCritic(gln = ntratComum,gld = NL,sig = alfa)

DMS1=sqrt(QMR/nrep1+QMR/nrep2)*q1


Test=round(tapply(Dados[id==0,5], as.factor(as.matrix(Dados[id==0,1])),mean,na=T),5)
Comum=round(tapply(var[id], as.factor(paste(Dados[id,2],Dados[id,3],sep="_")),mean,na=T),5)


Dunnet=NULL
for (it in 1:length(Test)){
  Comum2=Comum
  Comum2[abs(Comum-Test[it])>DMS1]=paste(Comum2[abs(Comum-Test[it])>DMS1], "*",sep="")
  Comum2[abs(Comum-Test[it])<DMS1]=paste(Comum2[abs(Comum-Test[it])<DMS1], "ns",sep="")
  Comum2=c(Comum2,Testemunha=Test[it])
  NomeTest=names(Test)[it]
  Dunnet=c(Dunnet,list(TesteDunnett=Comum2))
}
if(verbose){print(Dunnet)}

}


}
}




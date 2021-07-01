

Interpolacao=function(Mat,GLtrat,GLres,gln,gld){

  if((sum(GLtrat==gln)+sum(GLres==gld))==2){
    return(Mat[GLres==gld,GLtrat==gln])
  }

  if((sum(GLtrat==gln)+sum(GLres==gld))==1){
    #Interpolacao para trat
    if(sum(GLtrat==gln)==0){

      GLtrat2=GLtrat
      GLtrat2[GLtrat>gln]=10000
      ordem=order(abs(GLtrat2-gln))
      ID1=abs(GLtrat2-gln)==abs(GLtrat2-gln)[ordem[1]]
      GLtrat2=GLtrat
      GLtrat2[GLtrat<gln]=10000
      ordem=order(abs(GLtrat2-gln))
      ID2=abs(GLtrat2-gln)==abs(GLtrat2-gln)[ordem[1]]
      idmin=GLtrat[ID1]
      idmax=GLtrat[ID2]
      Vmin=Mat[GLres==gld,ID1]
      Vmax=Mat[GLres==gld,ID2]

      #    (idmax-idmin) --- (Vmax-Vmin)
      #    (idmax-gln)    --- x

      return(Vmax-(idmax-gln)*(Vmax-Vmin)/(idmax-idmin))
    }

    #Interpolacao para res
    if(sum(GLres==gld)==0){

      GLres2=GLres
      GLres2[GLres>gld]=10000
      ordem=order(abs(GLres2-gld))
      ID1=abs(GLres2-gld)==abs(GLres2-gld)[ordem[1]]
      GLres2=GLres
      GLres2[GLres<gld]=10000
      ordem=order(abs(GLres2-gld))
      ID2=abs(GLres2-gld)==abs(GLres2-gld)[ordem[1]]
      idmin=GLres[ID1]
      idmax=GLres[ID2]
      Vmin=Mat[ID1,GLtrat==gln]
      Vmax=Mat[ID2,GLtrat==gln]

      #    (idmax-idmin) --- (Vmax-Vmin)
      #    (idmax-gln)    --- x

      return(Vmax-(idmax-gld)*(Vmax-Vmin)/(idmax-idmin))
    }
  }

  if((sum(GLtrat==gln)+sum(GLres==gld))==0){

    GLtrat2=GLtrat
    GLtrat2[GLtrat>gln]=10000
    ordem=order(abs(GLtrat2-gln))
    ID1=abs(GLtrat2-gln)==abs(GLtrat2-gln)[ordem[1]]
    GLtrat2=GLtrat
    GLtrat2[GLtrat<gln]=10000
    ordem=order(abs(GLtrat2-gln))
    ID2=abs(GLtrat2-gln)==abs(GLtrat2-gln)[ordem[1]]
    gln2=GLtrat[ID1]
    gln1=GLtrat[ID2]

    ##################
    GLres2=GLres
    GLres2[GLres>gld]=10000
    ordem=order(abs(GLres2-gld))
    ID1=abs(GLres2-gld)==abs(GLres2-gld)[ordem[1]]
    GLres2=GLres
    GLres2[GLres<gld]=10000
    ordem=order(abs(GLres2-gld))
    ID2=abs(GLres2-gld)==abs(GLres2-gld)[ordem[1]]
    idmin=GLres[ID1]
    idmax=GLres[ID2]
    Vmin=Mat[ID1,GLtrat==gln1]
    Vmax=Mat[ID2,GLtrat==gln1]
    TabR1=(Vmax-(idmax-gld)*(Vmax-Vmin)/(idmax-idmin))

    Vmin=Mat[ID1,GLtrat==gln2]
    Vmax=Mat[ID2,GLtrat==gln2]
    TabR2=(Vmax-(idmax-gld)*(Vmax-Vmin)/(idmax-idmin))


    idmin=gln1
    idmax=gln2
    Vmin=TabR1
    Vmax=TabR2

    #    (idmax-idmin) --- (Vmax-Vmin)
    #    (idmax-gln)    --- x

    return(Vmax-(idmax-gln)*(Vmax-Vmin)/(idmax-idmin))




  }


}

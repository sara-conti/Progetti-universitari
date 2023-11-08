dati <- read.csv("dati_R.csv", sep=";"
                 , dec = ",",  
                 stringsAsFactors=TRUE, na.strings=c("NA","NaN",""))


library(metafor)
library(meta)
dati
dati$LN.OR.<-as.numeric(dati$LN.OR.)
dati$var<-as.numeric(dati$SE.LN.OR..**2)
#dati1$var<-as.numeric((dati1$Se**2))
dati_medio<-dati[which(dati$CAT_ESP == "medio"),]
dati_alto<-dati[which(dati$CAT_ESP == "alto"),]
fe_medio<-rma(dati_medio$LN.OR.,dati_medio$var,data=dati_medio,method="FE")
fe_alto<-rma(dati_alto$LN.OR.,dati_alto$var,data=dati_alto,method="FE")
summary(fe_medio)
summary(fe_alto)

#FE effetti fissi
#DL effetti casuali
#transf=exp da aggiungere
forest(fe_medio, xlim=c(-1,2), alim=c(0,2),transf=exp, xlab="OR", mlab="FE Estimate", cex=1, 
       main="FOREST PLOT CALSSE MEDIA", cex.main=1.5, digits=2, font.lab=2,  
       family="serif", refline=1)
forest(fe_alto, xlim=c(-1,2), alim=c(0,4), xlab="OR",transf = exp, mlab="FE Estimate", cex=1, 
       main="FOREST PLOT CLASSE ALTA", cex.main=1.5, digits=2, font.lab=2,  
       family="serif", refline=1)


funnel(fe_medio, main="Funnel plot - Medium vs Low consumption", xlab="log OR", ylab="invse")#, studlab = TRUE)
funnel(fe_alto, main="Funnel plot - High vs Low consumption", xlab="log OR", ylab="invse")#, slab = TRUE)

#test per l'asimmetria del funnel plot
regtest(fe_medio) 
regtest(fe_alto)
























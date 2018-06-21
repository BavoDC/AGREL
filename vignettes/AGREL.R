## ---- results = 'asis', echo = F, message = F, warnings = F--------------
library(knitr)
library(AGREL)
library(htmlTable)
options(table_counter = TRUE)
data("Agreement_deVetArticle")
tmp = Agreement_deVetArticle[,c(1,2)]
for(i in 1:2) tmp[,i] = factor(tmp[,i], level=c("Satisfied", "Dissatisfied"), ordered=T)
TableExample = table(tmp, dnn = c("Rater A", "Rater B"))
TableExample = addmargins(TableExample)
rownames(TableExample) <- sprintf('<b>%s</b>', rownames(TableExample))
print(htmlTable(TableExample, rowlabel = "Rater A", cgroup = c("Rater B",""), n.cgroup = c(2,1),
                align = "lccc", caption = "Two raters, dichotomous variable"))

## ------------------------------------------------------------------------
data("Agreement_deVetArticle")
tmp = Agreement_deVetArticle[,c(1,2)]
for(i in 1:2) tmp[,i] = factor(tmp[,i], level=c("Satisfied", "Dissatisfied"), ordered=T)
TableExample = table(tmp, dnn = c("Rater A", "Rater B"))
TableExample = addmargins(TableExample)
TableExample

## ------------------------------------------------------------------------
head(tmp)
Results1 = DiagnosticAgreement.deVet(tmp)
Results1

## ------------------------------------------------------------------------
DiagnosticAgreement.deVet(tmp, correction = "Fleiss")

## ------------------------------------------------------------------------
DiagnosticAgreement.deVet(tmp, correction = "bootstrap", NrBoot = 500)

## ---- results='asis', echo = F, message = F, warnings = F----------------
tmp = Agreement_deVetArticle[,c(1,4)]
for(i in 1:2) tmp[,i] = factor(tmp[,i], level=c("Satisfied", "Dissatisfied"), ordered=T)
TableExample2 = table(tmp)
TableExample2 = addmargins(TableExample2)
rownames(TableExample2) <- sprintf('<b>%s</b>', rownames(TableExample2))
print(htmlTable(TableExample2, rowlabel = "Rater A", cgroup = c("Rater C",""), n.cgroup = c(2,1),
                align = "lccc", caption = "> 2 raters, dichotomous variable"))

tmp = Agreement_deVetArticle[,c(2,4)]
for(i in 1:2) tmp[,i] = factor(tmp[,i], level=c("Satisfied", "Dissatisfied"), ordered=T)
TableExample3 = table(tmp)
TableExample3 = addmargins(TableExample3)
rownames(TableExample3) <- sprintf('<b>%s</b>', rownames(TableExample3))
print(htmlTable(TableExample3, rowlabel = "Rater B", cgroup = c("Rater C",""), n.cgroup = c(2,1),
                align = "lccc", caption = "> 2 raters, dichotomous variable"))

## ------------------------------------------------------------------------
# 3 raters

# A vs C
tmp = Agreement_deVetArticle[,c(1,4)]
for(i in 1:2) tmp[,i] = factor(tmp[,i], level=c("Satisfied", "Dissatisfied"), ordered=T)
TableExample2 = table(tmp)
TableExample2 = addmargins(TableExample2)
TableExample2

# B vs C
tmp = Agreement_deVetArticle[,c(2,4)]
for(i in 1:2) tmp[,i] = factor(tmp[,i], level=c("Satisfied", "Dissatisfied"), ordered=T)
TableExample3 = table(tmp)
TableExample3 = addmargins(TableExample3)
TableExample3

## ------------------------------------------------------------------------
DiagnosticAgreement.deVet(Agreement_deVetArticle[,c(1:2,4)])

## ------------------------------------------------------------------------
data("PsychMorbid")
AgreemGeneralizedVector(PsychMorbid)

## ------------------------------------------------------------------------
Df = PsychMorbid[,1:2]
CohenK(Df)

## ------------------------------------------------------------------------
?CohenK
ResultsCohenK = CohenK(Df)
str(ResultsCohenK)

## ------------------------------------------------------------------------
ConfintJack(ResultsCohenK)

## ------------------------------------------------------------------------
data(Agreement_deVet)
CohenK(Agreement_deVet[,2:3], weight = "squared")

## ------------------------------------------------------------------------
data(PsychMorbid)
FleissK(PsychMorbid)

## ------------------------------------------------------------------------
ResultsFleissK = FleissK(PsychMorbid)
str(ResultsFleissK)
ConfintJack(ResultsFleissK)

## ------------------------------------------------------------------------
KMatrix = Kappa.Matrix(PsychMorbid)
round(KMatrix$KappaMatrix, 3)

## ------------------------------------------------------------------------
# Minimum values possible 
round(KMatrix$MinLambda, 3)

# Maximum values possible
round(KMatrix$MaxLambda, 3)

## ---- warning = F--------------------------------------------------------
data(DentalStudy)
par(cex = 0.5, cex.lab = 0.75, cex.axis = 0.7)
ModifBAplot = BAplotMultipleR(dentist, patient, DMFS, DentalStudy)
par(cex = 1, cex.lab = 1, cex.axis = 1)

## ------------------------------------------------------------------------
ModifBAplot

## ------------------------------------------------------------------------
Df = cbind.data.frame(ID = sort(rep(1:10, 4)), var = sample(letters[1:2], 40, T), raters = rep(paste("rater",1:4), 10))
head(Df)

## ------------------------------------------------------------------------
NewDf = TransfData(ID, var, raters, Df)
FleissK(NewDf)


# rm(list=ls())
# gc()
# setwd("~/Download")

if (!require("ed90")) {
  install.packages("ed90")
}
library("ed90")
if (!require("boot")) {
  install.packages("boot")
}
library("boot")
patient.data<-read.csv("~/Downloads/PIEB.CSV")
pieb<-read.csv("~/Downloads/PIEB.CSV",1,encoding='UTF-8')  #这里导入数据

pavaData <- preparePava(pieb)  #这里是pavaData

bootResult <- boot(data = pieb,
                   statistic = bootIsotonicRegression,
                   R = 100,
                   sim = 'parametric',
                   ran.gen = bootIsotonicResample,
                   mle = list(baselinePava = pavaData,
                              firstDose = 0,
                              PROBABILITY.GAMMA = 0.9),
                   baselinePava = pavaData,
                   PROBABILITY.GAMMA = 0.9)

result=bootBC.ci(tObserved = bootResult$t0[3],
           tBoot = bootResult$t[, 3],
          conf = 0.95)

#结果都在result里面


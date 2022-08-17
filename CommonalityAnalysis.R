# commonality analysis

#----- load tools -----#
rm(list = ls())
library(foreign)
library(R.matlab)
library(gridExtra)
library(ggplot2)
library(yhat)

#----- load data -----#
MRPETmemory<-read.spss('/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_memoryTable_withLC-CR.sav', to.data.frame=TRUE)

## Predict paragraph comprehension based on three verbal
## tests: general info, sentence comprehension, & word
## classification

## Use HS dataset in MBESS 
#if (require ("MBESS")){
#  data(HS)
  
#  ## All-possible-subsets regression
#  apsOut=aps(HS,"t6_paragraph_comprehension",
#             list("t5_general_information", "t7_sentence","t8_word_classification"))
#  
#  ## Commonality analysis
#  commonality(apsOut)
#}

## Use HS dataset in MBESS 
if (require("MBESS")){
  data(HS)
  
  ## Regression
  lm.out<-lm(t6_paragraph_comprehension~
               t5_general_information+t7_sentence+t8_word_classification,data=HS)
  
  ## Regression Indices
  regr.out<-calc.yhat(lm.out)
}


# ---------- run analysis ----------- #
#         * data *     * Y *
apsOut=aps(MRPETmemory,"LC_CR_mean",
           list("Dprime_highDA_delayed","Dprime_lowDA_delayed","Dprime_rew_highDA_delayed","Dprime_rew_lowDA_delayed"))

commonality(apsOut)




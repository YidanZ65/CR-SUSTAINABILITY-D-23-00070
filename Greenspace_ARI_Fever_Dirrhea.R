##################################################################################
####       Greenspace and burden of acute respiratory infection, fever,       ####
####    and diarrhea, among children in 49 low- and middle-income countries   ####     
##################################################################################


library(survey)

################# MAIN ANALYSIS ###############################
###Dataset using KRbind, including the following variables:
## "ARI","Fever","Diarrhea",""cluster,"survey","NDVIaverage","NDVImax","Greenspace",
##"age","sex","parent_education","wealth_index","URBAN_RURA","Continent"
## "PM25", "Temperature", "SPEI","diet_diversity"



#### Main model ####
greenness <- c("NDVIaverage","NDVImax","Greenspace")
outcome <- c("ARI","Fever","Diarrhea")

formu <- list()
for (y in 1:3) {
  for (x in 1:3) {
    base1 <- as.formula(paste0(outcome[y], 
                               " ~ age+sex+parent_education+wealth_index+URBAN_RURA+Continent",
                               greenness[x]))
    base <- list(base1)
    formu <- c(formu, base)
  }
}


functionmain<-function(datan){
  data3<-data.frame()
  for (i in 1:9) {
    dclus1<-svydesign(id=~cluster,weights=~new_weight,data=datan)
    assign('model',svyglm(update(formu[[i]], . ~ .),
                          family=quasibinomial(), design = dclus1))
    estimate<-as.data.frame(coef(summary(model)))
    estimate_x<-estimate[7, which(colnames(estimate)=="Estimate")]
    estimate_sd<-estimate[7,which(colnames(estimate)=="Std. Error")]
    estimate_p<-estimate[7,which(colnames(estimate)=="Pr(>|t|)")]
    ci<-data.frame(confint(model))
    ci_x<-ci[7,]
    ci_x$estimate<-estimate_x
    ci_x$sd<-estimate_sd
    ci_x$p<-estimate_p
    ci_x$rown<-rownames(ci_x)
    ci_x$total<-paste0(sprintf("%.2f",ci_x$estimate)," (",sprintf("%.2f",ci_x$X2.5..),", ",sprintf("%.2f",ci_x$X97.5..),")")
    lower<-ci_x[,1]
    upper<-ci_x[,2]
    ORl<-exp(lower)
    ORu<-exp(upper)
    ORe<-exp(estimate_x)
    ci_x$ORll<-ORl
    ci_x$ORuu<-ORu
    ci_x$ORee<-ORe
    ci_x$ORR<-paste0(sprintf("%.2f",ci_x$ORee)," (",sprintf("%.2f",ci_x$ORll),", ",sprintf("%.2f",ci_x$ORuu),")")
    ci_x$numb<-rep(dim(datan)[1],nrow(ci_x))
    data3<-rbind(data3,ci_x)
  }
  return(data3)
}

main_greenspace <- functionmain(KRbind)


#### Startified analysis ####

greenness <- c("NDVIaverage","NDVImax","Greenspace")
outcome <- c("ARI","Fever","Diarrhea")


sub1 <- list()
for (y in 1:3) {
  for (x in 1:3) {
    base1 <- as.formula(paste0(outcome[y], 
                               " ~ age+wealth_index+parent_edu+Continent+URBAN_RURA +",
                               greenness[x]))
    base <- list(base1)
    sub1 <- c(sub1, base)
  }
}


sub2 <- list()
for (y in 1:3) {
  for (x in 1:3) {
    base1 <- as.formula(paste0(outcome[y], 
                               " ~ age+sex+parent_edu+Continent+URBAN_RURA +",
                               greenness[x]))
    base <- list(base1)
    sub2 <- c(sub2, base)
  }
}


sub3 <- list()
for (y in 1:3) {
  for (x in 1:3) {
    base1 <- as.formula(paste0(outcome[y], 
                               " ~ age_cat+sex+parent_edu+Continent+wealth_index +",
                               greenness[x]))
    base <- list(base1)
    sub3 <- c(sub3, base)
  }
}


functionsub<-function(sub,datan){
  data3<-data.frame()
  for (i in 1:9) {
    dclus1<-svydesign(id=~cluster,weights=~new_weight,data=datan)
    assign('model',svyglm(update(sub[[i]], . ~ .),
                          family=quasibinomial(), design = dclus1))
    estimate<-as.data.frame(coef(summary(model)))
    estimate_x<-estimate[8, which(colnames(estimate)=="Estimate")]
    estimate_sd<-estimate[8,which(colnames(estimate)=="Std. Error")]
    estimate_p<-estimate[8,which(colnames(estimate)=="Pr(>|t|)")]
    ci<-data.frame(confint(model))
    ci_x<-ci[8,]
    ci_x$estimate<-estimate_x
    ci_x$sd<-estimate_sd
    ci_x$p<-estimate_p
    ci_x$rown<-rownames(ci_x)
    ci_x$total<-paste0(sprintf("%.2f",ci_x$estimate)," (",sprintf("%.2f",ci_x$X2.5..),", ",sprintf("%.2f",ci_x$X97.5..),")")
    lower<-ci_x[,1]
    upper<-ci_x[,2]
    ORl<-exp(lower)
    ORu<-exp(upper)
    ORe<-exp(estimate_x)
    ci_x$ORll<-ORl
    ci_x$ORuu<-ORu
    ci_x$ORee<-ORe
    ci_x$ORR<-paste0(sprintf("%.2f",ci_x$ORee)," (",sprintf("%.2f",ci_x$ORll),", ",sprintf("%.2f",ci_x$ORuu),")")
    ci_x$numb<-rep(dim(datan)[1],nrow(ci_x))
    data3<-rbind(data3,ci_x)
  }
  return(data3)
}

sex1 <- functionsub(sub1, KRbind[which(KRbind$sex == 1),])
sex2 <- functionsub(sub1, KRbind[which(KRbind$sex == 2),])

wealth1 <- functionsub2(sub2, KRbind[which(KRbind$wealth_binary == 1),])
wealth2 <- functionsub2(sub2, KRbind[which(KRbind$wealth_binary == 2),])

UR1 <- functionsub3(sub3, KRbind[which(KRbind$URBAN_RURA == "U"),])
UR2 <- functionsub3(sub3, KRbind[which(KRbind$URBAN_RURA == "R"),])

sex <- cbind(sex1[,c(3,4)], sex2[,c(3,4)])
sex$p_difference <- 2*pnorm(-abs((sex1$estimate-sex2$estimate)/sqrt(sex1$sd^2+sex2$sd^2)))

wealth <- cbind(wealth1[,c(3,4)], wealth2[,c(3,4)])
wealth$p_difference <- 2*pnorm(-abs((wealth1$estimate-wealth2$estimate)/sqrt(wealth1$sd^2+wealth2$sd^2)))

UR <- cbind(UR1[,c(3,4)], UR2[,c(3,4)])
UR$p_difference <- 2*pnorm(-abs((UR1$estimate-UR2$estimate)/sqrt(UR1$sd^2+UR2$sd^2)))


#### Mediation analysis ####

#####PM25, Temperature and SPEI as Mediator
Exposure <- "Greenspace"
Mediator <- c("PM25", "Temperature", "SPEI")
Outcome <- c("ARI","Fever","Diarrhea")

####calculate estimates
get_stat<-function(data,covar2,mediator,exposure,interactionterm,outcome,M=NULL,O=NULL){
  
  data <- filter(data, !is.na(get(mediator)))
  
  # Pick only the covariates
  df_cov <- data[ , covar2, drop = FALSE]
  
  # Loop over every column and get the mean
  # mcv<-sapply(df_cov, mean, na.rm = T)
  
  #form_M <- paste(mediator, "~", exposure, "+", paste(covar1, collapse = " + "))
  form_M <- paste("LST_c", "~", exposure, "+", paste(covar2, collapse = " + "))
  
  if(M==1){
    outm.fit<-glm(as.formula(form_M), family=binomial(link=logit), data = data) 
    ssm<-NULL  
  }
  
  if(M==0){
    outm.fit<-lm(as.formula(form_M), data = data)
    ssm<-(summary(outm.fit)$sigma)**2
  }
  
  
  #bcv<-coef(summary(outm.fit))[c(3:9),"Estimate"]
  b0<-coef(summary(outm.fit))["(Intercept)","Estimate"]
  b1<-coef(summary(outm.fit))[exposure,"Estimate"]
  bcc<-0
  
  
  #browser()
  form_Y <- paste( outcome, "~", exposure, "+", mediator, "+", interactionterm, "+", paste(covar2, collapse = " + "))
  if(O==1){
    
    dclus1<-svydesign(id=~cluster,weights=~new_weight,data=data)
    
    outy.fit <- svyglm(as.formula(form_Y),
                       family=quasibinomial(), design = dclus1)
    
  }
  
  if(O==0){
    outy.fit<-glm(as.formula(form_Y), family=gaussian(link=identity), data = data)
    
  }
  
  
  
  t1<-coef(summary(outy.fit))[exposure,"Estimate"]
  t2<-coef(summary(outy.fit))[mediator, "Estimate"]
  t3<-coef(summary(outy.fit))[paste(exposure,mediator,sep=':'), "Estimate"]
  
  return(c( bcc,b0,t1,t2,t3,b1,ssm ))
  
}


####calculate four components
boot.cMbO<- function(data, indices){

  data<-data[indices, ]
  
  outcome<- Y
  exposure<- A
  mediator<- M
  #  covar1<-COVAR1
  covar2<-COVAR
  
  interactionterm<- paste(mediator,exposure,sep='*')
  
  
  statistics<-get_stat(data,covar2,mediator,exposure,interactionterm, outcome,M=0,O=1)
  
  bcc<-statistics[1]
  b0<-statistics[2]
  t1<-statistics[3]
  t2<-statistics[4]
  t3<-statistics[5]
  b1<-statistics[6]
  ssm<-statistics[7]
  
  
  
  total<-exp(t1 + t2*b1 + (t3*(b0 + b1*astar + b1*a + bcc + t2*ssm)*(a - astar)) + (0.5*(t3*t3)*ssm*((a*a)-(astar*astar))))
  #browser()
  cde_comp1<-t1*(a - astar) + t2*mstar + t3*a*mstar
  cde_comp2<-(t2 + t3*astar)*(b0 + b1*astar + bcc)
  cde_comp3<-0.5*(t2 + t3*astar)*(t2 + t3*astar)*ssm
  cde_comp4<-t2*mstar + t3*astar*mstar
  cde_comp5<-(t2 + t3*astar)*(b0 + b1*astar + bcc)
  cde_comp6<-0.5*(t2 + t3*astar)*(t2 + t3*astar)*ssm
  cde_comp<-exp(cde_comp1 - cde_comp2 - cde_comp3) - exp(cde_comp4 - cde_comp5 - cde_comp6)
  
  
  
  intref_comp1<-(t1 + t3*(b0 + b1*astar + bcc + t2*ssm)*(a - astar)) + (0.5*(t3*t3)*ssm*((a*a)-(astar*astar)))
  intref_comp2<-(t1*(a - astar) + t2*mstar + t3*a*mstar) - ((t2 + t3*astar)*(b0 + b1*astar + bcc)) - (0.5*(t2 + t3*astar)*(t2 + t3*astar)*ssm)
  intref_comp3<-(t2*mstar + t3*astar*mstar) - ((t2 + t3*astar)*(b0 + b1*astar + bcc)) - (0.5*(t2 + t3*astar)*(t2 + t3*astar)*ssm)
  intref_comp<-exp(intref_comp1) - 1  - exp(intref_comp2) + exp(intref_comp3)
  
  intmed_comp1<-(t1 + t2*b1 + t3*(b0 + b1*astar + b1*a + bcc + t2*ssm))*(a - astar) + (0.5*(t3*t3)*ssm*((a*a)-(astar*astar)))
  intmed_comp2<-(t2*b1 + t3*b1*astar)*(a - astar)
  intmed_comp3<-(t1 + t3*(b0 + b1*astar + bcc + t2*ssm))*(a - astar) + (0.5*(t3*t3)*ssm*((a*a)-(astar*astar)))
  intmed_comp<-exp(intmed_comp1) - exp(intmed_comp2) - exp(intmed_comp3) + 1
  
  pie_comp<-exp(((t2*b1 + t3*b1*astar)*(a - astar))) - 1
  
  terr<-cde_comp + intref_comp + intmed_comp + pie_comp
  
  excess_rr_total<-total - 1
  excess_rr_cde<-cde_comp*(total - 1)/terr
  excess_rr_intref<-intref_comp*(total - 1)/terr
  excess_rr_intmed<-intmed_comp*(total - 1)/terr
  excess_rr_pie<-pie_comp*(total - 1)/terr
  
  prop_cde<-cde_comp/terr
  prop_intref<-intref_comp/terr
  prop_intmed<-intmed_comp/terr
  prop_pie<-pie_comp/terr
  prop_med<-(pie_comp + intmed_comp)/terr
  prop_int<-(intmed_comp + intref_comp)/terr
  prop_elm<-(pie_comp + intmed_comp + intref_comp)/terr
  
  
  
  ests<- c(total, excess_rr_total, excess_rr_cde, excess_rr_intref, excess_rr_intmed, excess_rr_pie, 
           prop_cde, prop_intref, prop_intmed, prop_pie, prop_med, prop_int, prop_elm)
  return(ests)
}

####save results
save_results<-function(boot_function=NULL,N=NULL){

  system.time(results<- boot(data = data, statistic = boot.cMbO, R = N_r))
  
  #collecting results so that they can be output into a single file
  
  df_list<-list()
  
  if (outcome==1) {
    
    ESTIMAND<-c(
      "Total Effect Risk Ratio", "Total Excess Relative Risk", "Excess Relative Risk due to CDE", "Excess Relative Risk due to INTref",
      "Excess Relative Risk due to INTmed","Excess Relative Risk due to PIE", "Proportion CDE", "Proportion INTref",
      "Proportion INTmed", "Proportion PIE", "Overall Proportion Mediated",
      "Overall Proportion Attributable to Interaction","Overall Proportion Eliminated"
    )
    
    index=13
    
  }
  else if (outcome==0) {
    
    ESTIMAND<-c(
      "Total Effect", "CDE", "INTref", "INTmed",
      "PIE","Proportion CDE","Proportion INTref",
      "Proportion INTmed", "Proportion PIE", "Overall Proportion Mediated",
      "Overall Proportion Attributable to Interaction","Overall Proportion Eliminated"
    )
    index=12
    
  }
  
  for (i  in 1:index){
    r<-boot.ci(results, type = "basic", index=i) 
    estq<-r[[2]]
    cis<- r[[4]]
    colnames(cis)<-NULL
    lcl<-cis[1,4]
    ucl<-cis[1,5]
    estimand<-ESTIMAND[i]
    rdf<-data.frame(estimand, estq, lcl, ucl)
    names(rdf)<-c("Estimand", "Estimate", "LCL", "UCL")
    rdf$X <- A
    rdf$M <- M
    rdf$Y <- Y
    df_list[[i]]<-rdf
  }
  
  # Putting all this together
  allresults <- rbindlist(df_list)
  
  return(allresults)
  
}


for (Y in Outcome) {
  for (A in Exposure) {
    for (M in Mediator) {
      
      COVAR<<-c('age','sex','parent_edu',"wealth_index","URBAN_RURA","Continent")
      
      
      #1=binary 0=continuous  
      outcome=1
      mediator=0
      
      #Assign levels for the exposure that are being compared; 
      #for mstar it is the level at which to compute the CDE and the remainder of the decomposition  
      a<<-1     
      astar<<-0 
      mstar<<-0 
      
      #Boostrap number of iterations
      N_r=1000
      
      results <- save_results(boot_function=boot.cMbO, N=N_r)
      
      med_result <- rbind(results, med_result)
    }
  }
}


#####Dietary diversity as Mediator
Exposure <- "Greenspace"
Mediator <- c("Diet_diversity")
Outcome <- c("ARI","Fever","Diarrhea")


boot.bMbO<- function(data, indices){
  
  
  outcome<- Y
  exposure<- A
  mediator<- M
  covar<-COVAR
  
  interactionterm<- paste(mediator,exposure,sep='*')
  
  data<-data[indices, ]     
  
  statistics<-get_stat(data,covar,mediator,exposure,interactionterm,outcome,M=1,O=1)
  
  bcc<-statistics[1]
  b0<-statistics[2]
  t1<-statistics[3]
  t2<-statistics[4]
  t3<-statistics[5]
  b1<-statistics[6]
  
  total<- exp(t1*a)*(1+exp(b0+b1*astar+bcc))*(1+exp(b0+b1*a+bcc+t2+t3*a)) / ( exp(t1*astar)*(1+exp(b0+b1*a+bcc))* (1+ exp(b0+b1*astar+bcc+t2+t3*astar)))
  
  cde_comp1<- exp(t1*(a - astar) + t2*mstar + t3*a*mstar)
  cde_comp2<-(1+exp(b0+b1*astar+bcc))
  cde_comp3<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  cde_comp4<- exp(t2*mstar+t3*astar*mstar)
  cde_comp5<-(1+exp(b0+b1*astar+bcc))
  cde_comp6<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  cde_comp<- cde_comp1*cde_comp2/cde_comp3-cde_comp4*cde_comp5/cde_comp6
  
  intref_comp1<- exp(t1*(a-astar))
  intref_comp2<-(1+exp(b0+b1*astar+bcc+t2+t3*a))
  intref_comp3<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intref_comp4<- exp(t1*(a-astar)+t2*mstar+t3*a*mstar)
  intref_comp5<-(1+exp(b0+b1*astar+bcc))
  intref_comp6<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intref_comp7<- exp(t2*mstar+t3*astar*mstar)
  intref_comp8<-(1+exp(b0+b1*astar+bcc))
  intref_comp9<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intref_comp<-intref_comp1*intref_comp2/intref_comp3-1-intref_comp4*intref_comp5/intref_comp6 + intref_comp7*intref_comp8/intref_comp9
  
  intmed_comp1<- exp(t1*(a-astar))
  intmed_comp2<-(1+exp(b0+b1*a+bcc+t2+t3*a))
  intmed_comp3<-(1+exp(b0+b1*astar+bcc))
  intmed_comp4<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intmed_comp5<-(1+exp(b0+b1*a+bcc))
  intmed_comp6<-(1+exp(b0+b1*a+bcc+t2+t3*astar))
  intmed_comp7<-(1+exp(b0+b1*astar+bcc)) 
  intmed_comp8<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intmed_comp9<-(1+exp(b0+b1*a+bcc))
  intmed_comp10<- exp(t1*(a-astar))
  intmed_comp11<-(1+exp(b0+b1*astar+bcc+t2+t3*a))
  intmed_comp12<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intmed_comp<- intmed_comp1*intmed_comp2*intmed_comp3/(intmed_comp4*intmed_comp5)-intmed_comp6*intmed_comp7/(intmed_comp8*intmed_comp9)-intmed_comp10*intmed_comp11/intmed_comp12 + 1
  
  pie_comp1<-(1+exp(b0+b1*astar+bcc))
  pie_comp2<-(1+exp(b0+b1*a+bcc+t2+t3*astar))
  pie_comp3<-(1+exp(b0+b1*a+bcc))
  pie_comp4<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  pie_comp<-(pie_comp1*pie_comp2)/(pie_comp3*pie_comp4) - 1
  
  terr<- cde_comp + intref_comp + intmed_comp + pie_comp
  
  excess_rr_total<- total - 1
  excess_rr_cde<- cde_comp*(total - 1)/terr
  excess_rr_intref<- intref_comp*(total - 1)/terr
  excess_rr_intmed<- intmed_comp*(total - 1)/terr
  excess_rr_pie<- pie_comp*(total - 1)/terr
  
  prop_cde<- cde_comp/terr
  prop_intmed<- intmed_comp/terr
  prop_intref<- intref_comp/terr
  prop_pie<- pie_comp/terr
  prop_med<- (pie_comp + intmed_comp)/terr
  prop_int<- (intmed_comp + intref_comp)/terr
  prop_elm<- (pie_comp + intmed_comp + intref_comp)/terr
  
  ests<- c(total, excess_rr_total, excess_rr_cde, excess_rr_intref, excess_rr_intmed, excess_rr_pie, 
           prop_cde, prop_intref, prop_intmed, prop_pie, prop_med, prop_int, prop_elm)
  
  return(ests)
}

save_results<-function(boot_function=NULL,N=NULL){
  
  system.time(results<- boot(data = data, statistic = boot.bMbO, R = N_r))
  
  #collecting results so that they can be output into a single file
  
  df_list<-list()
  
  if (outcome==1) {
    
    ESTIMAND<-c(
      "Total Effect Risk Ratio", "Total Excess Relative Risk", "Excess Relative Risk due to CDE", "Excess Relative Risk due to INTref",
      "Excess Relative Risk due to INTmed","Excess Relative Risk due to PIE", "Proportion CDE", "Proportion INTref",
      "Proportion INTmed", "Proportion PIE", "Overall Proportion Mediated",
      "Overall Proportion Attributable to Interaction","Overall Proportion Eliminated"
    )
    
    index=13
    
  }
  else if (outcome==0) {
    
    ESTIMAND<-c(
      "Total Effect", "CDE", "INTref", "INTmed",
      "PIE","Proportion CDE","Proportion INTref",
      "Proportion INTmed", "Proportion PIE", "Overall Proportion Mediated",
      "Overall Proportion Attributable to Interaction","Overall Proportion Eliminated"
    )
    index=12
    
  }
  
  for (i  in 1:index){
    r<-boot.ci(results, type = "basic", index=i) 
    estq<-r[[2]]
    cis<- r[[4]]
    colnames(cis)<-NULL
    lcl<-cis[1,4]
    ucl<-cis[1,5]
    estimand<-ESTIMAND[i]
    rdf<-data.frame(estimand, estq, lcl, ucl)
    names(rdf)<-c("Estimand", "Estimate", "LCL", "UCL")
    rdf$X <- A
    rdf$M <- M
    rdf$Y <- Y
    df_list[[i]]<-rdf
  }
  
  # Putting all this together
  allresults <- rbindlist(df_list)
  
  return(allresults)
  
}


med_result <- data.frame()

for (Y in Outcome) {
  for (A in Exposure) {
    for (M in Mediator) {
      
      COVAR<<-c('age_','sex','parent_edu',"wealth_index","URBAN_RURA","Continent")
      
      
      #1=binary 0=continuous  
      outcome=1
      mediator=1
      
      #Assign levels for the exposure that are being compared; 
      #for mstar it is the level at which to compute the CDE and the remainder of the decomposition  
      a<<-1     
      astar<<-0 
      mstar<<-0 
      
      #Boostrap number of iterations
      N_r=1000
      
      results <- save_results(boot_function=boot.cMbO, N=N_r)
      
      med_result <- rbind(results, med_result)
    }
  }
}

library(tidyverse)
library(pROC)

load('dpp_manuscript_data.RData')

#### Typhoid versus culture negative
mod.comb = glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                  data=data)
roc.comb = roc(mod.comb$y,mod.comb$fitted.values,ci=T)
roc.comb

#plot
textsize=10
fig_roc_all<-ggroc(roc.comb) + theme_classic() +
  theme(axis.text.x = element_text(size = textsize),axis.text.y=element_text(size=textsize), 
        legend.position='none',axis.title = element_text(size=textsize)) +
  scale_x_reverse(expand=c(0,0)) + scale_y_continuous(expand = c(0, 0)) +
  annotate("text",x=0.20,y=0.10,label="AUC: 0.920 (95% CI: 0.897-0.944)") 
  
#### Typhoid versus alternative etiologies

mod.comb.alt<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
               data=data %>% filter(typhoid.d==1 | alter_etiol==1))
roc.comb.alt = roc(mod.comb.alt$y,mod.comb.alt$fitted.values,ci=T)
roc.comb.alt

fig_roc_alt <-ggroc(roc.comb.alt) +  theme_classic()  +
  theme(axis.text.x = element_text(size = textsize),axis.text.y=element_text(size=textsize), 
        legend.position='none',axis.title = element_text(size=textsize)) +
  scale_x_reverse(expand=c(0,0)) + scale_y_continuous(expand = c(0, 0))  +
  annotate("text",x=0.20,y=0.10,label="AUC: 0.969 (95% CI: 0.943-0.994)") 

#### Age
#under 5
mod.comb.alt.under5<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                   data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(age<5))
roc.comb.alt.under5 = roc(mod.comb.alt.under5$y,mod.comb.alt.under5$fitted.values,ci=T)
roc.comb.alt.under5

#over 5
mod.comb.alt.over5<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                          data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(age>=5))
roc.comb.alt.over5 = roc(mod.comb.alt.over5$y,mod.comb.alt.over5$fitted.values,ci=T)
roc.comb.alt.over5

roc.test(roc.comb.alt.under5,roc.comb.alt.over5,method="delong")

#### Fever duration

mod.comb.alt.feverunder5<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                   data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(fever_duration<5))
roc.comb.alt.feverunder5 = roc(mod.comb.alt.feverunder5$y,mod.comb.alt.feverunder5$fitted.values,ci=T)
roc.comb.alt.feverunder5

mod.comb.alt.feverover5<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                               data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(fever_duration>=5))
roc.comb.alt.feverover5 = roc(mod.comb.alt.feverover5$y,mod.comb.alt.feverover5$fitted.values,ci=T)
roc.comb.alt.feverover5

roc.test(roc.comb.alt.feverunder5,roc.comb.alt.feverover5,method="delong")

#### Antibiotic use

mod.comb.alt.noabx<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                   data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(abx==0))
roc.comb.alt.noabx = roc(mod.comb.alt.noabx$y,mod.comb.alt.noabx$fitted.values,ci=T)
roc.comb.alt.noabx

mod.comb.alt.abx<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                         data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(abx==1))
roc.comb.alt.abx = roc(mod.comb.alt.abx$y,mod.comb.alt.abx$fitted.values,ci=T)
roc.comb.alt.abx

roc.test(roc.comb.alt.noabx,roc.comb.alt.abx,method="delong")

#### Typhoid vs paratyphoid

mod.comb.alt.typhoid<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                   data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(Paratyphi.d != 1))
roc.comb.alt.typhoid = roc(mod.comb.alt.typhoid$y,mod.comb.alt.typhoid$fitted.values,ci=T)
roc.comb.alt.typhoid

mod.comb.alt.paratyphoid<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                           data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(Typhi.d != 1))
roc.comb.alt.paratyphoid = roc(mod.comb.alt.paratyphoid$y,mod.comb.alt.paratyphoid$fitted.values,ci=T)
roc.comb.alt.paratyphoid

roc.test(roc.comb.alt.typhoid,roc.comb.alt.paratyphoid,method="delong")



#### Capillary and venous blood

mod.comb.capillary <-  glm(typhoid.d ~ dpp_lps_3_capillary+ dpp_hlye_3_capillary, family=binomial(logit),
                           data=data %>% filter(typhoid.d==1 | alter_etiol==1))
roc.comb.capillary = roc(mod.comb.capillary$y,mod.comb.capillary$fitted.values,ci=T)
roc.comb.capillary                      
                           
mod.comb.venous <- glm(typhoid.d ~ dpp_lps_venous + dpp_hlye_venous, family=binomial(logit),
                       data=data %>% filter(typhoid.d==1 | alter_etiol==1))
roc.comb.venous = roc(mod.comb.venous$y,mod.comb.venous$fitted.values,ci=T)
roc.comb.venous                      




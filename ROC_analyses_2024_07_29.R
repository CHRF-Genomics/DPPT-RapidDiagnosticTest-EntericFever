
library(tidyverse)
library(pROC)
library(stringr)
library(scales)

load('dpp_manuscript_data.RData')

#### Typhoid versus culture negative
mod.comb = glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                  data=data)
roc.comb = roc(mod.comb$y,mod.comb$fitted.values,ci=T)
roc.comb
comb.ci = sprintf('%.3f',roc.comb$ci)

roc.hlye = roc(typhoid.d ~ dpp_hlye_combined,ci=T,data=data)
roc.hlye
hlye.ci = sprintf('%.3f',roc.hlye$ci)

roc.lps = roc(typhoid.d ~ dpp_lps_combined,ci=T,data=data)
roc.lps
lps.alt.ci = sprintf('%.3f',roc.lps$ci)


#plot
textsize=10

fig_roc <-ggroc(list(roc.comb,roc.hlye,roc.lps)) +  theme_classic()  +
  theme(axis.text.x = element_text(size = textsize),axis.text.y=element_text(size=textsize), 
        legend.position='none',axis.title = element_text(size=textsize)) +
  scale_x_reverse(expand=c(0,0)) + scale_y_continuous(expand = c(0, 0))  +
  annotate("text",x=0.30,y=0.15,
           color=hue_pal()(3)[1],
           label=str_glue('Combined, AUC: {comb.ci[2]} (95% CI: {comb.ci[1]}-{comb.ci[3]})')) +
  annotate("text",x=0.33,y=0.11,
           color=hue_pal()(3)[2],
           label=str_glue('HlyE, AUC: {hlye.ci[2]} (95% CI: {hlye.ci[1]}-{hlye.ci[3]})')) +  
  annotate("text",x=0.33,y=0.07,
           color=hue_pal()(3)[3],
           label=str_glue('LPS, AUC: {lps.ci[2]} (95% CI: {lps.ci[1]}-{lps.ci[3]})'))


#### Typhoid versus alternative etiologies

mod.comb.alt<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
               data=data %>% filter(typhoid.d==1 | alter_etiol==1))
roc.comb.alt = roc(mod.comb.alt$y,mod.comb.alt$fitted.values,ci=T)
roc.comb.alt
comb.alt.ci = sprintf('%.3f',roc.comb.alt$ci)

roc.hlye.alt = roc(typhoid.d ~ dpp_hlye_combined,ci=T,data=data %>% filter(typhoid.d==1 | alter_etiol==1))
roc.hlye.alt
hlye.alt.ci = sprintf('%.3f',roc.hlye.alt$ci)

roc.lps.alt = roc(typhoid.d ~ dpp_lps_combined,ci=T,data=data %>% filter(typhoid.d==1 | alter_etiol==1))
roc.lps.alt
lps.alt.ci = sprintf('%.3f',roc.lps.alt$ci)


fig_roc_alt <-ggroc(roc.comb.alt) +  theme_classic()  +
  theme(axis.text.x = element_text(size = textsize),axis.text.y=element_text(size=textsize), 
        legend.position='none',axis.title = element_text(size=textsize)) +
  scale_x_reverse(expand=c(0,0)) + scale_y_continuous(expand = c(0, 0))  +
  annotate("text",x=0.20,y=0.10,label=str_glue("AUC: {comb.alt.ci[2]} (95% CI: {comb.alt.ci[1]}-{comb.alt.ci[3]})"))

show_col(hue_pal()(3))

fig_roc_alt <-ggroc(list(roc.comb.alt,roc.hlye.alt,roc.lps.alt)) +  theme_classic()  +
  theme(axis.text.x = element_text(size = textsize),axis.text.y=element_text(size=textsize), 
        legend.position='none',axis.title = element_text(size=textsize)) +
  scale_x_reverse(expand=c(0,0)) + scale_y_continuous(expand = c(0, 0))  +
  annotate("text",x=0.30,y=0.15,
           color=hue_pal()(3)[1],
           label=str_glue('Combined, AUC: {comb.alt.ci[2]} (95% CI: {comb.alt.ci[1]}-{comb.alt.ci[3]})')) +
  annotate("text",x=0.33,y=0.11,
           color=hue_pal()(3)[2],
           label=str_glue('HlyE, AUC: {hlye.alt.ci[2]} (95% CI: {hlye.alt.ci[1]}-{hlye.alt.ci[3]})')) +  
  annotate("text",x=0.33,y=0.07,
           color=hue_pal()(3)[3],
           label=str_glue('LPS, AUC: {lps.alt.ci[2]} (95% CI: {lps.alt.ci[1]}-{lps.alt.ci[3]})'))






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


##### Sex
mod.comb.alt.male<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                          data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(sex==1))
roc.comb.alt.male = roc(mod.comb.alt.male$y,mod.comb.alt.male$fitted.values,ci=T)
roc.comb.alt.male

mod.comb.alt.female<- glm(typhoid.d ~ dpp_hlye_combined+ dpp_lps_combined, family=binomial(logit),
                        data=data %>% filter(typhoid.d==1 | alter_etiol==1) %>% filter(sex==2))
roc.comb.alt.female = roc(mod.comb.alt.female$y,mod.comb.alt.female$fitted.values,ci=T)
roc.comb.alt.female

roc.test(roc.comb.alt.male,roc.comb.alt.female,method="delong")





library(tidyverse)
library(pROC)
library(gridExtra)
library(ggbeeswarm)
library(GGally)

load('dpp_manuscript_data.RData')

#### Figure 2

textsize=10
p1<-data %>%
  mutate(all_etiol = factor(all_etiol,levels=c("Dengue","RSV","RSV and Dengue","Rickettsia","Influenza","None","Typhoid/Paratyphoid"))) %>%
  ggplot(aes(x=all_etiol,y=dpp_hlye_3_capillary,color=all_etiol)) + geom_quasirandom(varwidth = TRUE) + 
  theme_bw() + ylab("Capillary anti-HlyE IgA (DPP)") +  theme(legend.position = "none", axis.text.x= element_text(angle=45,hjust=1,size=textsize),
                                                         axis.text.y= element_text(size=textsize)) +
  xlab('') +
  scale_y_continuous(trans='log10')


p2<- data %>%
  mutate(all_etiol = factor(all_etiol,levels=c("Dengue","RSV","RSV and Dengue","Rickettsia","Influenza","None","Typhoid/Paratyphoid"))) %>%
  ggplot(aes(x=all_etiol,y=dpp_lps_3_capillary,color=all_etiol)) + geom_quasirandom(varwidth = TRUE) + 
  scale_y_continuous(trans='log10') +
  theme_bw() + ylab("Capillary anti-LPS IgA (DPP)") +  theme(legend.position = "none", axis.text.x= element_text(angle=45,hjust=1,size=textsize),
                                                        axis.text.y= element_text(size=textsize),axis.title.y = element_text(size = textsize)) +
  xlab('') 

fig3<- grid.arrange(p1,p2,nrow=1)


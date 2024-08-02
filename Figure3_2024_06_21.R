library(tidyverse)
library(GGally)
library(gridExtra)

load('dpp_manuscript_data.RData')

cor_func <- function(data, mapping, method, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor(x, y, method=method, use='complete.obs')
  
  
  ggally_text(
    label = as.character(round(corr, 2)), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...
  )
}

fig_correlation_lps<-data %>% select(dpp_lps_venous,dpp_lps_plasma,dpp_lps_3_capillary) %>% 
  ggpairs(columnLabels = c("Venous", "Plasma", "Capillary"),upper = list(continuous = wrap(cor_func,method = 'pearson')),
          diag = list(continuous = "blankDiag")) +
  theme_bw() +
  scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') + 
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(tag='A',title = "Anti-LPS IgA")

fig_correlation_hlye<-data %>% select(dpp_hlye_venous,dpp_hlye_2_plasma,dpp_hlye_3_capillary) %>% 
  ggpairs(columnLabels = c("Venous", "Plasma", "Capillary"),upper = list(continuous = wrap(cor_func,method = 'pearson')),
          diag = list(continuous = "blankDiag")) +
  theme_bw() +
  scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') + 
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(tag='B', title="Anti-HlyE IgA")





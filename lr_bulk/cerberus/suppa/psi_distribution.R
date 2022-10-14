library(ggsci)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(tidyverse)

metadata <- read.table('swan_metadata.txt', sep = '\t', header = T) # metatdata of biosa,ples
extract_PSIvalue <- function(psi_file, metadata, type){
  df <- read.table(psi_file) %>% rownames_to_column('event') %>% gather(dataset, value, -event) %>% merge(metadata, by = 'dataset')
  df <- df %>%  na.omit()  %>% mutate(type = type)
  return(df)
}
A3 <- extract_PSIvalue('psi/cerberus_A3.psi', metadata, 'A3')
A5 <- extract_PSIvalue('psi/cerberus_A5.psi', metadata, 'A5')
AF <- extract_PSIvalue('psi/cerberus_AF.psi', metadata, 'AF')
AL <- extract_PSIvalue('psi/cerberus_AL.psi', metadata, 'AL')
MX <- extract_PSIvalue('psi/cerberus_MX.psi', metadata, 'MX')
RI <- extract_PSIvalue('psi/cerberus_RI.psi', metadata, 'RI')
SE <- extract_PSIvalue('psi/cerberus_SE.psi', metadata, 'SE')

p <- rbind(A3, A5, AF, AL, MX, RI, SE) %>%
  ggplot(aes(x = sample, y = value, fill = type)) + geom_violin() + labs(x = '', y = 'PSI value') + facet_wrap(~type, ncol = 1) + scale_fill_npg() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +  theme(legend.title=element_blank())
ggsave('psi.pdf', p, height = 10, width = 8)
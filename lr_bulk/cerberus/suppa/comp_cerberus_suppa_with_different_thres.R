library(dplyr)
library(ggplot2)
library(ggsci)
library(forcats)
## different PSI threshold (0.05/0.95, 0.1/0.9, 0.25/0.75)
AF_combined <- read.table('cerberus_suppa/0.1_0.9/tss_combined.tsv', sep = '\t', header = T) %>% na.omit() %>%
  dplyr::filter(suppa == TRUE | cerberus == TRUE) %>% mutate(type = ifelse((suppa == TRUE & cerberus == TRUE), 'Consensus', ifelse(suppa == TRUE, 'Only suppa', 'Only cerberus')))

AL_combined <- read.table('cerberus_suppa/0.1_0.9/tes_combined.tsv', sep = '\t', header = T) %>% na.omit() %>%
  dplyr::filter(suppa == TRUE | cerberus == TRUE) %>% mutate(type = ifelse((suppa == TRUE & cerberus == TRUE), 'Consensus', ifelse(suppa == TRUE, 'Only suppa', 'Only cerberus')))

ggplot(AF_combined, aes(x = fct_infreq(sample), fill=type)) + geom_bar() + scale_fill_npg() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + labs(x = '', y = '', title = 'TSS') + theme(legend.title=element_blank(), legend.position="bottom")
ggsave('cerberus_suppa/0.1_0.9//TSS.pdf', width = 10, height = 5)

ggplot(AL_combined, aes(x = fct_infreq(sample), fill=type)) + geom_bar() + scale_fill_npg() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + labs(x = '', y = '', title = 'TES') + theme(legend.title=element_blank(), legend.position="bottom")
ggsave('cerberus_suppa/0.1_0.9/TES.pdf', width = 10, height = 5)
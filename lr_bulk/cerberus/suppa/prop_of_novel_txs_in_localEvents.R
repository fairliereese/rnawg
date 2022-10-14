library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(ggsci)
library(patchwork)
## using new transcript table
transcripts <- read.table('cerberus_transcript_novelty.tsv', sep = '\t', header = T, check.names = F, comment.char = '')
transcripts$transcript_id <- gsub(',','_', transcripts$transcript_id)

AF <- read.table('localEvents/cerberus.events_AF_strict.ioe', header = T)
events <- unlist(lapply(AF$alternative_transcripts, function(x) str_split(x, ',')[[1]]))
AF <- subset(transcripts, transcripts$transcript_id %in% events) %>% mutate(event = 'AF')

A3 <- read.table('localEvents/cerberus.events_A3_strict.ioe', header = T)
events <- unlist(lapply(A3$alternative_transcripts, function(x) str_split(x, ',')[[1]]))
A3 <- subset(transcripts, transcripts$transcript_id %in% events) %>% mutate(event = 'A3')

A5 <- read.table('localEvents/cerberus.events_A5_strict.ioe', header = T)
events <- unlist(lapply(A5$alternative_transcripts, function(x) str_split(x, ',')[[1]]))
A5 <- subset(transcripts, transcripts$transcript_id %in% events) %>% mutate(event = 'A5')

AL <- read.table('localEvents/cerberus.events_AL_strict.ioe', header = T)
events <- unlist(lapply(AL$alternative_transcripts, function(x) str_split(x, ',')[[1]]))
AL <- subset(transcripts, transcripts$transcript_id %in% events) %>% mutate(event = 'AL')

MX <- read.table('localEvents/cerberus.events_MX_strict.ioe', header = T)
events <- unlist(lapply(MX$alternative_transcripts, function(x) str_split(x, ',')[[1]]))
MX <- subset(transcripts, transcripts$transcript_id %in% events) %>% mutate(event = 'MX')

RI <- read.table('localEvents/cerberus.events_RI_strict.ioe', header = T)
events <- unlist(lapply(RI$alternative_transcripts, function(x) str_split(x, ',')[[1]]))
RI <- subset(transcripts, transcripts$transcript_id %in% events) %>% mutate(event = 'RI')

SE <- read.table('localEvents/cerberus.events_SE_strict.ioe', header = T)
events <- unlist(lapply(SE$alternative_transcripts, function(x) str_split(x, ',')[[1]]))
SE <- subset(transcripts, transcripts$transcript_id %in% events) %>% mutate(event = 'SE')

tss <- rbind(A3, A5, AF, AL, MX, RI, SE) %>% 
  select(transcript_id, type = tss_novelty, event) %>%
  distinct(transcript_id, type, event,.keep_all = TRUE) %>%
  ggplot(aes(x = fct_infreq(event), fill = type)) + geom_bar(position = 'fill') + scale_fill_npg() +
  coord_flip() + labs(x = 'Local events', y = 'Proportion of novel transcript events defined by TSS') + theme_classic() +
  theme(legend.title=element_blank())

tes <- rbind(A3, A5, AF, AL, MX, RI, SE) %>% 
  select(transcript_id, type = tes_novelty, event) %>%
  distinct(transcript_id, type, event,.keep_all = TRUE) %>%
  ggplot(aes(x = fct_infreq(event), fill = type)) + geom_bar(position = 'fill') + scale_fill_npg() +
  coord_flip() + labs(x = 'Local events', y = 'Proportion of novel transcript events defined by TES') + theme_classic() +
  theme(legend.title=element_blank())

ic <- rbind(A3, A5, AF, AL, MX, RI, SE) %>% 
  select(transcript_id, ic_novelty, event) %>% mutate(type = if_else(ic_novelty == 'Known', 'Known','Novel')) %>%
  distinct(transcript_id, type, event,.keep_all = TRUE) %>%
  ggplot(aes(x = fct_infreq(event), fill = type)) + geom_bar(position = 'fill') + scale_fill_npg() +
  coord_flip() + labs(x = 'Local events', y = 'Proportion of novel transcript events defined by IC') + theme_classic() +
  theme(legend.title=element_blank())


tss + tes + ic+plot_layout(guides = 'collect')
ggsave('prop.pdf', width = 15, height = 3)
## total number of local events in cerberus
count <- rbind(A3, A5, AF, AL, MX, RI, SE) %>% 
  select(transcript_id, type = ic_novelty, event) %>%
  distinct(transcript_id, type, event,.keep_all = TRUE) %>%
  group_by(event) %>% summarise(n = n())
write.table(count, file = 'tot_number_of_localEvents.csv', quote = F, row.names = F, col.names = F, sep = '\t')
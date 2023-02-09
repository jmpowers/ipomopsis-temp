library(tidyverse)
library(reshape2)
library(lubridate)

# Read chromatograms ------------------------------------------------------

source("read_shimadzu.R")
batches <- read_tsv("data/volatiles/CWu Ipomopsis Scent Metadata - batches.tsv") %>% 
  distinct(batch) %>% pull(batch)
#former Shimadzu output directories:
#~/Downloads/ipoieldVols21/batches #has same versions as in RMBL Batches
#/media/C/GCMSsolution/Data/heather Ipo/Temperature expt CWu 2018/ #copied to RMBL Batches from backup drive
datapath <- "~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/"
# slow steps, uncomment to reload Shimadzu data
#ipo.data <- map_dfr(set_names(paste0(datapath, "RMBL Batches/", batches), batches), read.shimadzu, .id="batch")
#save(ipo.data, file="./data/volatiles/ipotemp181921.Rdata")

load("data/volatiles/ipotemp181921.Rdata") 
ipo.data %>% count(batch, Filename) %>% count(batch)

ipo.all <- dcast(ipo.data, Filename~Name, sum, value.var="Area")
rownames(ipo.all) <- ipo.all[,1]
ipo.all[,1] <- NULL
ipo.cut <- ipo.all[,colSums(ipo.all)>5e8]#arbitrary cutoff

# k-means  --------------------------------------------------------------------
k <- 40
set.seed(1)
km <- kmeans(decostand(ipo.cut, method="log"), k, nstart=3)
#save(km, file="km30.Rdata")

ipo.km <- tibble(FileName=rownames(ipo.all)) %>% 
  mutate(rowSum = rowSums(ipo.all),
         Project = str_extract(FileName, "Blank|[aA]ir|Aeven|Conditioning|Corydalis|PH|SE|SP|Sterile|Yeast|Veg|Ambient|GNA|SC|Lupinus") %>% replace_na("sample"),
         Type = fct_collapse(Project, blank="Blank", air=c("air","Air"), other_level = "sample"),
         nameBlank = Type=="blank",
         runYear = str_extract(FileName, "2018|2019|2020|2021") %>% replace_na("2018") %>% factor,
         Cluster = km$cluster) %>% # Figure out which k-means clusters are the blanks
  mutate(kBlank = Cluster %in% (count(., nameBlank, Cluster) %>% filter(nameBlank, n>2) %>% pull(Cluster)),
         Mixup = nameBlank != kBlank)

with(ipo.km, table(kBlank, nameBlank))

# Join Ipo samples with Markes inventory -----------------------------------------------------------
# sequ.summary is the Markes inventory fuzzy joined by time with the directory listing by markes_sequence.R
sequ.summary <- read_csv(paste0(datapath,"Inventory/sequfile_181921.csv"))
# the fuzzy joining introduces 2+6+0 erroneous duplicates, 14 of Lucas' 2019 samples, and 42+18+18 files with no Markes match
sequ.summary %>% count(fuzzy_n, year(eithertime)) #large count is result=NA from no Markes match, eithertime=NA are skips

# get all the Markes batches involved in files for temperature experiment
ipo.batchids <- sequ.summary %>% filter(FileName %in% ipo.km$FileName) %>% pull(id) %>% unique() %>% na.omit()

ipogc <- sequ.summary %>% filter(id %in% ipo.batchids | is.na(id)) %>%  #keep in files with no Markes match
  left_join(ipo.km %>% select(FileName, nameBlank, Mixup, kBlank, Cluster)) %>% 
  mutate(verdict="", sample="", index=row_number()) %>% 
  #get batch file names and peak counts from Shimadzu output - one blank filename is not unique
  left_join(ipo.data %>% rename(FileName=Filename) %>% group_by(FileName,batch) %>% tally(name="n_peaks"), multiple="first") %>% 
  select(c("index", "sequence.start", "batch", "Desorb.Start.Time", "CreationTime", "eithertime", "status", 
           "Tube", "markes_n", "GC_n", "either_n", "markes_GC", "create_desorb", "desorb.Start.diff", 
           "Mixup", "nameBlank", "kBlank", "Cluster", "n_peaks", "verdict", "FileName", "sample", "user", "FullName", "id", "fuzzy_n")) %>% 
  write_csv("data/volatiles/ipogc_181921.csv")

#Manually wrote in verdicts and correct skip names for FileName in sample
read_tsv("data/volatiles/RMBL GC-MS Data Inventory - ipogc181921.tsv") %>% 
  filter(is.na(verdict) | verdict != "mismatched")

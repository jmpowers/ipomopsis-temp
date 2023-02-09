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
#datapath <- "~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/"
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
sequ.summary <- read_csv("data/volatiles/sequfile_181921.csv")
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
ipogc <- read_tsv("data/volatiles/RMBL GC-MS Data Inventory - ipogc181921.tsv")


# # Overview - no metadata
# ## Inventory
# ```{r inventory, eval=FALSE}
# #TODO merge F, G, H, I, J into GNA for 2018
# with(ipos, table(Year, Type))
# with(ipos, table(Type, DNLeaf))
# ```
# 
# ## NMDS
# ```{r nmds, eval=FALSE}
# ipo.cut <- ipo.all[,colSums(ipo.all)>1e9]
# set.seed(1)
# nmds.ipo <- metaMDS(decostand(ipo.cut, "hellinger"), dist="bray", autotransform = FALSE, try=5, trymax=5)
# 
# par(bg="grey40")
# ordiplot(nmds.ipo, type = "n")
# ipos$rs <- rowSums(ipo.all)
# ipos$nameBlank <- ipos$Type %in% c("Blank","Bake")
# ipos$nameAir <- ipos$Type == "Air"
# thisyear <- ipos$Year == 2019
# #with(ipos, points(nmds.ipo, display="sites", col=viridis(200)[round(200*sqrt(rs/max(rs)))], pch=as.integer(isblank)+1))
# with(ipos, points(nmds.ipo, display="sites", col=rainbow(nlevels(ipos$Type))[ipos$Type], pch=as.integer(nameAir*2+nameBlank)+1))
# #with(ipos, points(nmds.ipo, display="sites", col=rainbow(nlevels(ipos$Type))[ipos$Type], pch=ipos$Year-2017))
# legend("topright", legend=levels(ipos$Type), fill=rainbow(nlevels(ipos$Type)), cex=0.8)
# ```
# 
# ## Kmeans clustering
# ```{r kmeans, eval=FALSE}
# k <- 8
# set.seed(1)
# km <- kmeans(decostand(ipo.cut, "log"), k, nstart=3)
# ipos$Cluster <- km$cluster
# ipos$kBlank <- kblank <- km$cluster %in% c(4)
# ipos$Mixup <- ipos$nameBlank != ipos$kBlank
# with(ipos, table(kBlank, nameBlank))
# 
# ordiplot(nmds.ipo, type = "n")
# points(nmds.ipo, display="sites", col=ifelse(kblank, "black", rainbow(k)[km$cluster]), pch=with(ipos,as.integer(nameAir*2+nameBlank))+1) 
# text(nmds.ipo, display="species",col="grey80",labels=ifelse(colSums(ipo.cut)>4e8, colnames(ipo.cut), ""))
# text(nmds.ipo, display="sites", col=with(ipos, as.integer(thisyear)*4+as.integer(nameBlank)*2+as.integer(nameAir)+1), labels=as.character(km$cluster))
# ```
# 
# 
# ## CAP - blanks and years
# ```{r cap_blankyear, eval=FALSE}
# ipo.cap <- capscale(ipo.cut ~ as.factor(ipos$kBlank) * as.factor(thisyear), distance="bray", metaMDSdist = F)
# plot(ipo.cap, type="n")
# points(ipo.cap, display="sites", col=ifelse(ipos$kBlank, "black", rainbow(k)[ipos$Cluster]), pch=as.integer(ipos$nameAir*2+ipos$nameBlank)+1) 
# #View(ipo.cap$CCA$v)
# ```
# 
# ## CAP - DNLeaf
# ```{r cap_DNLeaf, eval=FALSE}
# dnl <- ipos$DNLeaf!="" & ipos$Type !="Air"
# ipo.cap.dnl <- capscale(ipo.cut[dnl,] ~ ipos[dnl,"DNLeaf"], distance="bray", metaMDSdist = F)
# plot(ipo.cap.dnl, type="n")
# points(ipo.cap.dnl, display="sites", pch=as.integer(ipos[dnl,"DNLeaf"])-1, col=as.integer(factor(ipos[dnl,"Type"]))) 
# legend("bottomleft", levels(factor(ipos[dnl,"Type"])), fill=1:nlevels(ipos$Type), cex=0.8)
# legend("left", levels(ipos[dnl,"DNLeaf"]), pch=1:nlevels(ipos[dnl,"DNLeaf"])-1, cex=0.8)
# #View(ipo.cap.dnl$CCA$v)
# ```

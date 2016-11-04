# Doing some quick analysis of the spread in the initial conditions ensemble (ICE)
# Load the libraries
library(ggplot2)
library(abind)

# Set up the file path
dir.ice <- in.base  <- "~/Desktop/init_test_data/lat42.75lon-77.25/SAS_out"

# Get list of base file names
files.base <- dir(dir.ice, ".css")
files.base <- substr(files.base, 1, nchar(files.base)-4)

# Extract the data, summarize it, & put it into a data frame
# plant.summary <- data.frame(IC=1:length(files.base))
# site.summary  <- data.frame(IC=1:length(files.base))
for(i in 1:length(files.base)){
  css <- read.table(file.path(dir.ice, paste0(files.base[i], ".css")), header=T, sep=" ")
  pss <- read.table(file.path(dir.ice, paste0(files.base[i], ".pss")), header=T, sep=" ")
  
  # Calculate site-level statistics on the initialization
  css <- merge(css, pss[,c("patch", "area")])
  css$biomass <- rowSums(css[,c("bdead", "balive")]) * css$n * css$area 
  css$basal.area <- pi*(css$dbh/2)^2 * css$n * css$area 
  css$density    <- css$n * css$area

  css.summary <- aggregate(css[,c("basal.area", "density", "biomass")], 
                           by=css[,c("pft", "year")],
                           FUN=sum)
  css.summary$IC <- i
  
  # Calculate the patch summaries
  pss.summary <- data.frame(total.C = sum(pss[1,c("fsc", "stsc", "ssc")]),
                            total.N = sum(pss[1,c("msn", "fsn")]))
  pss.summary$IC <- i
  
  if(i==1){
    css.all <- css.summary
    pss.all <- pss.summary
  } else {
    css.all <- rbind(css.all, css.summary)
    pss.all <- rbind(pss.all, pss.summary)
  }
}

css.all$pft <- as.factor(css.all$pft)
summary(css.all)
summary(pss.all)

summary(css.all[css.all$pft==5,])

ggplot(data=css.all[css.all$pft!=5,]) +
  facet_wrap(~pft, scales="free") +
  geom_density(aes(x=basal.area, fill=as.factor(pft)), alpha=0.5)+
  scale_x_continuous(expand=c(0,0))

ggplot(data=css.all[css.all$pft!=5,]) +
  facet_wrap(~pft, scales="free") +
  geom_freqpoly(aes(x=basal.area, fill=as.factor(pft)), alpha=0.5)+
  scale_x_continuous(expand=c(0,0))


ggplot(data=css.all[css.all$pft!=5,]) +
  geom_density(aes(x=basal.area, fill=as.factor(pft)), alpha=0.5)+
  scale_x_continuous(expand=c(0,0))

ggplot(data=css.all[css.all$pft!=5,]) +
  geom_density(aes(x=density, fill=as.factor(pft)), alpha=0.5)+
  scale_x_continuous(expand=c(0,0))

ggplot(data=css.all[css.all$pft!=5,]) +
  facet_wrap(~pft, scales="free") +
  geom_density(aes(x=biomass, fill=as.factor(pft)), alpha=0.5)+
  scale_x_continuous(expand=c(0,0))

ggplot(data=pss.all) +
  geom_density(aes(x=total.C), fill="blue", alpha=0.5) +
  scale_x_continuous(expand=c(0,0))

ggplot(data=pss.all) +
  geom_density(aes(x=total.N), fill="blue", alpha=0.5) +
  scale_x_continuous(expand=c(0,0))
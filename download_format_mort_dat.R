####be careful with the downloading section of this code because the connection sometimes times out (or has a related issue), causing certain datasets to only partially download

library(readr)
library(dplyr)
library(here)
library(data.table)
library(XML)
library(dplyr)
library(sqldf)
library(RCurl)
library(here)

setwd(here())

######download data
#links
file_links <- c(paste(rep("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHIS_", 33), 1986:2018, rep("_MORT_2019_PUBLIC.dat", 33), sep=""))
                                                                                                                     
#set destinations
destinations <- paste(1986:2018, ".dat", sep="")

#download files
for(i in seq_along(file_links)){
  download.file(url=file_links[i], destfile=destinations[i])
}

######format data

files <- paste(1986:2018, ".dat", sep="")
df_comb <- NULL

for(i in 1:33){
  print(i)
  df <- read_fwf(file=files[i],
                 col_types = "ciiiiiiidd",
                 fwf_cols(publicid = c(1,14),
                          eligstat = c(15,15),
                          mortstat = c(16,16),
                          ucod_leading = c(17,19),
                          diabetes = c(20,20),
                          hyperten = c(21,21),
                          dodqtr = c(22,22),
                          dodyear = c(23,26),
                          wgt_new = c(27,34),
                          sa_wgt_new = c(35,42)
                 ),
                 na = c("", "."))
  df$year <- 1985 + i
  df_comb <- rbind(df_comb, df)
}

#####save data
write.csv(df_comb, "mort_dat_1986_2018_2019.csv")

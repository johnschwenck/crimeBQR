# Google Drive package to access data files
library(googledrive)
library(googlesheets4)

# Authenticate with Google Drive
drive_auth()
drive_user()

data_path <- 'https://drive.google.com/drive/u/1/folders/1qFN7uwvR0vx3-UyjKCtu_g_GmxhY6uNW'
files <- drive_ls(data_path)
files$name

gs4_auth()
datasets <- c("court_summons_hist","court_summons_ytd","shootings_hist",
              "shootings_ytd","mvc_crashes","mvc_person","mvc_vehicles",
              "arrests_hist","arrests_ytd","complaints_hist","complaints_ytd")

for(i in 1: length(datasets)){
  
  tmp = gs4_find(datasets[i])
  tmp$drive_resource[1][[1]]$exportLinks$`text/csv`
  
}


pre_path <- 'https://docs.google.com/spreadsheets/d/'
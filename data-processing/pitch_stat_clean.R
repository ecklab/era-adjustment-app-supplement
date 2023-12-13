rm(list=ls())
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)
library(readr)
library(EnvStats)
library(dplyr)
ncores <- detectCores() -1

setwd("~/Desktop/PhD_UIUC/pythonlearning")

###################
## pitchers

## pitching stats from baseball reference

years <- c(1871:2022)
pitching_stats <- do.call(rbind, lapply(years, function(x) {
  m <- read_csv(paste("pitching", x, ".csv", sep = ""))
  m <- m %>% mutate(name = gsub('[^[:alnum:] ]',' ',Name)) %>%
    mutate(name = ifelse(str_sub(name, -1) == ' ', gsub('.{1}$', '',name), name))
  cbind(m[-nrow(m),], yearID = x)
}))

pitching_stats[is.na(pitching_stats$Lg),]$Lg <- 'NL'

pitching_stats <- pitching_stats %>% filter(Lg %in% c('NL', 'AL')) %>% 
  filter(Tm != "TOT")


## park factor adjustment for Runs, Hits and Home Runs
setwd("~/Documents/GitHub/retrosplits/daybyday")

## park factor adjustment from 1901 - 2022

years <- c(1901:2022)
team_log <- do.call(rbind, lapply(years, function(x) 
  cbind(read_csv(paste("teams-", x, ".csv", sep = "")), yearID = x)))

team_log <- team_log %>% select(c(team.alignment, team.key, opponent.key, 
                                  B_R, B_H, B_HR, P_R, P_H, P_HR, yearID))

park_factor_year <- function(year){
  team_log_year <- team_log %>% filter(yearID == year)
  teamID <- unique(team_log_year$team.key)
  teamID <- teamID[!teamID %in% c('ALS', 'NLS', 'NL1', 'NL2', 'AL1', 'AL2')]
  do.call(rbind, mclapply(teamID, mc.cores = ncores, FUN = function(xx){
    home_RS <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 1) %>%
                     select(B_R))
    home_HS <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 1) %>%
                     select(B_H))
    home_HRS <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 1) %>%
                      select(B_HR))
    
    home_RA <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 1) %>%
                     select(P_R))
    home_HA <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 1) %>%
                     select(P_H))
    home_HRA <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 1) %>%
                      select(P_HR))
    
    away_RS <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 0) %>%
                     select(B_R))
    away_HS <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 0) %>%
                     select(B_H))
    away_HRS <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 0) %>%
                      select(B_HR))
    
    away_RA <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 0) %>%
                     select(P_R))
    away_HA <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 0) %>%
                     select(P_H))
    away_HRA <- sum(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 0) %>%
                      select(P_HR))
    
    homeG <- nrow(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 1))
    awayG <- nrow(team_log_year %>% filter(team.key == xx) %>% filter(team.alignment == 0))
    
    n <- nrow(unique(team_log_year %>% select(team.key)))
    
    park_index_R <- ((home_RS+home_RA)/(homeG)) / ((away_RS+away_RA)/(awayG))
    park_index_H <- ((home_HS+home_HA)/(homeG)) / ((away_HS+away_HA)/(awayG))
    park_index_HR <- ((home_HRS+home_HRA)/(homeG)) / ((away_HRS+away_HRA)/(awayG))
    
    park_factor_R <- ((park_index_R +1)/2)/((park_index_R+n-1)/n)
    park_factor_H <- ((park_index_H +1)/2)/((park_index_H+n-1)/n)
    park_factor_HR <- ((park_index_HR +1)/2)/((park_index_HR+n-1)/n)
    
    data.frame(teamID = xx, yearID = year, park_factor_R = park_factor_R, 
               park_factor_H = park_factor_H, park_factor_HR = park_factor_HR)
  }))
}

park_factor <- do.call(rbind, mclapply(years, mc.cores = ncores, FUN = function(xx){
  park_factor_year(xx)
}))

## add baseball reference teamID to chadwick teamID

setwd("~/Documents/GitHub/baseballdatabank/core")

teams_chadwick <- read.csv('Teams.csv')

teams_chadwick_1901 <- teams_chadwick %>% filter(yearID >= 1901) %>% 
  filter(lgID %in% c('NL', 'AL')) %>% 
  select(yearID, teamID, teamIDBR)

## teamID in 2022 has not been updated on Chadwick
## manually add 2022 teamID

teams_chadwick_1901 <- rbind(teams_chadwick_1901, 
                             cbind(yearID = rep(2022,30), 
                                   teams_chadwick_1901[2557:2586,2:3]))

m <- park_factor
m$teamID <- as.character(m$teamID)
m$teamID[m$yearID >= 2005] <- gsub('ANA', 'LAA', m$teamID[m$yearID >= 2005])
m$teamID[m$yearID >= 1953 & m$yearID <= 1965] <- gsub('MLN', 'ML1', m$teamID[m$yearID >= 1953 & m$yearID <= 1965])
m$teamID[m$yearID >= 1970 & m$yearID <= 1997] <- gsub('MIL', 'ML4', m$teamID[m$yearID >= 1970 & m$yearID <= 1997])
park_factor <- m


## merge park factor from Chadwick with teamID from baseball reference
park_factor_1901 <- merge(park_factor, teams_chadwick_1901, 
                          by = c('yearID', 'teamID'))

## moving window for park factor. 

moving_window_mean <- function(x){
  n <- length(x)
  y <- rep(0,n)
  if (n == 1) {
    y[1] = x[1]
  }
  if (n == 2) {
    y[1:2] = mean(x)
  }
  if (n >= 3) {
    y[1] <- (x[1] + x[2])/2
    y[2:(n-1)] <- (x[1:(n-2)] + x[2:(n-1)] + x[3:n])/3
    y[n] <- (x[n-1] + x[n])/2
  }
  y
}
park_factor_averaged <- do.call(rbind, mclapply(unique(park_factor_1901$teamID), mc.cores = ncores, FUN = function(xx){
  if (xx == "LAA") {
    m <- park_factor_1901 %>% filter(teamID == xx)
    m$park_factor_R[1:4] <- moving_window_mean(m$park_factor_R[1:4])
    m$park_factor_H[1:4] <- moving_window_mean(m$park_factor_H[1:4])
    m$park_factor_HR[1:4] <- moving_window_mean(m$park_factor_HR[1:4])
    m$park_factor_R[5:21] <- moving_window_mean(m$park_factor_R[5:21])
    m$park_factor_H[5:21] <- moving_window_mean(m$park_factor_H[5:21])
    m$park_factor_HR[5:21] <- moving_window_mean(m$park_factor_HR[5:21])
    m
  } else{
    m <- park_factor_1901 %>% filter(teamID == xx)
    m$park_factor_R <- moving_window_mean(m$park_factor_R)
    m$park_factor_H <- moving_window_mean(m$park_factor_H)
    m$park_factor_HR <- moving_window_mean(m$park_factor_HR)
  }
  m
}))

pitchers <- pitching_stats[,-c(1,2)]
colnames(pitchers)[c(2,34)] <- c('teamID', 'playerID')

park_factor_br <- park_factor_averaged[,-2]
colnames(park_factor_br)[5] <- 'teamID'

pitchers_park_factor <- merge(pitchers %>% filter(yearID >= 1901), park_factor_br, 
                             by = c("teamID", "yearID"))

## No data for the teams before 1901
## No park factor for them

pre01_pitchers <- pitchers %>% filter(yearID <= 1900) 
pitchers_adj <- rbind(cbind(pre01_pitchers, park_factor_R = 1, park_factor_H = 1, 
                           park_factor_HR = 1), pitchers_park_factor)

pitchers_adj <- pitchers_adj %>% mutate(H_PK =  H / park_factor_H) %>% 
  mutate(HR_PK =  HR / park_factor_HR) %>% 
  mutate(ER_PK =  ER / park_factor_R)

## merge with population
## years to consider
years <- 1871:2022

setwd("~/Documents/GitHub/era_adjustment_resources/eligible_pop/")

MLB_pop <- read.csv('datMLBpop.csv') 

## add weighted populations to data
pops_data <- rbind(
  data.frame(pops = MLB_pop$NL_pop, lgID = "NL", yearID = 1870 + 1:(15*10 + 2)),
  data.frame(pops = MLB_pop$AL_pop, lgID = "AL", yearID = 1870 + 1:(15*10 + 2)))

colnames(pitchers_adj)[3] <- "lgID"
pop_pitchers <- merge(pitchers_adj, pops_data, by = c("yearID", "lgID"))

pop_pitchers$HBP[is.na(pop_pitchers$HBP)] <- 0

pitchers_stats <- pop_pitchers %>% group_by(playerID, yearID) %>% 
  summarise(name = unique(name), 
            age = unique(Age),
            G = sum(G), 
            lgID = lgID[which.max(IP)], 
            teamID = teamID[which.max(IP)],
            IP = sum(IP), 
            obs_ER = sum(ER), 
            obs_HR = sum(HR), 
            obs_H = sum(H), 
            SO = sum(SO), 
            ER = sum(ER_PK), 
            HR_PF = sum(HR_PK), 
            H_PF = sum(H_PK), 
            BB = sum(BB), 
            HBP = sum(HBP),
            pops = pops[which.max(IP)])


## batting bwar from baseball reference
pitching_bwar = readr::read_csv("https://www.baseball-reference.com/data/war_daily_pitch.txt", na = "NULL")
pitching_bwar <- pitching_bwar[!is.na(pitching_bwar$WAR),]
pitching_bWAR <- pitching_bwar %>% filter(lg_ID %in% c('NL', 'AL', 'NA')) %>% 
  select(player_ID, year_ID, WAR) %>% group_by(player_ID, year_ID) %>% 
  summarise(bWAR = sum(WAR))
colnames(pitching_bWAR)[c(1,2)] <- c('playerID', 'yearID')

## all batting stats except fWAR are merged together
pitching_stats_bWAR <- merge(pitchers_stats, pitching_bWAR, 
                            by = c('playerID', 'yearID'))

## fWAR
setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/expansion/updated_data/pitchers/WAR/data")
pitching_fWAR_2021 <- read.csv('pitchers_fWAR.csv')[,-1]
pitching_fWAR_2022 <- rbind(cbind(read.csv('pitchers_fWAR_2022_AL.csv'), 
                                 yearID = 2022, lgID = 'AL'), 
                           cbind(read.csv('pitchers_fWAR_2022_NL.csv'), 
                                 yearID = 2022, lgID = 'NL'))

pitching_fWAR <- rbind(pitching_fWAR_2021, pitching_fWAR_2022)

pitchers_fWAR <- pitching_fWAR %>% group_by(playerid, yearID) %>% summarise(
  name = unique(Name), fWAR = sum(WAR)
)

## merge fangraph ID with baseball reference ID

setwd("~/Documents/GitHub/register/data/")
people <- read.csv('people.csv') %>% select(key_bbref, key_fangraphs)

bbref_fangraphs_ID <- people[!is.na(people$key_fangraphs),]
colnames(bbref_fangraphs_ID)[2] <- 'playerid'

pitchers_fWAR_bWAR <- merge(pitchers_fWAR, bbref_fangraphs_ID, by = 'playerid')

## combine fWAR with batting stats and bWAR

colnames(pitchers_fWAR_bWAR)[5] <- 'playerID'

pitchers <- merge(pitching_stats_bWAR, 
                 pitchers_fWAR_bWAR %>% select(playerID, yearID, fWAR), 
                 by = c('playerID', 'yearID'))

setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/expansion/updated_data/pitchers")
write.csv(pitchers, 'pitchers.csv')

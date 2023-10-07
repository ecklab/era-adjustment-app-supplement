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
## batters


## batting stats from baseball reference

years <- c(1871:2022)
batting_stats <- do.call(rbind, lapply(years, function(x) {
  m <- read_csv(paste("batting", x, ".csv", sep = ""))
  m <- m %>% mutate(name = gsub('[^[:alnum:] ]',' ',Name)) %>%
    mutate(name = ifelse(str_sub(name, -1) == ' ', gsub('.{1}$', '',name), name))
  cbind(m[-nrow(m),], yearID = x)
}))

batting_stats[is.na(batting_stats$Lg),]$Lg <- 'NL'

batting_stats <- batting_stats %>% filter(Lg %in% c('NL', 'AL')) %>% 
  filter(Tm != "TOT")

c <- pre01_batter %>% group_by(yearID, teamID) %>% summarise(n = n())
table(c$yearID)


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

batters <- batting_stats[,-c(1,2,30)]
colnames(batters)[c(2,28)] <- c('teamID', 'playerID')

park_factor_br <- park_factor_averaged[,-2]
colnames(park_factor_br)[5] <- 'teamID'

batters_park_factor <- merge(batters %>% filter(yearID >= 1901), park_factor_br, 
                             by = c("teamID", "yearID"))

## No data for the teams before 1901
## No park factor for them

pre01_batter <- batters %>% filter(yearID <= 1900) 
batters_adj <- rbind(cbind(pre01_batter, park_factor_R = 1, park_factor_H = 1, 
                           park_factor_HR = 1), batters_park_factor)

batters_adj <- batters_adj %>% mutate(H_PK =  H / park_factor_H) %>% 
  mutate(HR_PK =  HR / park_factor_HR) %>% 
  mutate(R_PK =  R / park_factor_R)

## merge with population
## years to consider
years <- 1871:2022

setwd("~/Documents/GitHub/era_adjustment_resources/eligible_pop/")

MLB_pop <- read.csv('datMLBpop.csv') 

## add weighted populations to data
pops_data <- rbind(
  data.frame(pops = MLB_pop$NL_pop, lgID = "NL", yearID = 1870 + 1:(15*10 + 2)),
  data.frame(pops = MLB_pop$AL_pop, lgID = "AL", yearID = 1870 + 1:(15*10 + 2)))

colnames(batters_adj)[3] <- "lgID"
pop_batters <- merge(batters_adj, pops_data, by = c("yearID", "lgID"))

pop_batters$BA[is.na(pop_batters$BA)] <- 0
pop_batters$OBP[is.na(pop_batters$OBP)] <- 0
pop_batters$SH[is.na(pop_batters$SH)] <- 0
pop_batters$SF[is.na(pop_batters$SF)] <- 0
pop_batters$HBP[is.na(pop_batters$HBP)] <- 0

batters_stats <- pop_batters %>% group_by(playerID, yearID) %>% 
  summarise(name = unique(name), 
            age = unique(Age),
            G = sum(G), 
            lgID = lgID[which.max(PA)], 
            PA = sum(PA), 
            AB = sum(AB), 
            BB = sum(BB), 
            obs_hits = sum(H), 
            obs_HR = sum(HR), 
            HBP = sum(HBP),
            SH = sum(SH), 
            SF = sum(SF), 
            H = sum(H_PK), 
            HR = sum(HR_PK), 
            pops = pops[which.max(PA)])


## batting bwar from baseball reference
batting_bwar = readr::read_csv("https://www.baseball-reference.com/data/war_daily_bat.txt", na = "NULL")
batting_bwar <- batting_bwar[!is.na(batting_bwar$WAR),]
batters_bWAR <- batting_bwar %>% filter(lg_ID %in% c('NL', 'AL', 'NA')) %>% 
  select(player_ID, year_ID, WAR) %>% group_by(player_ID, year_ID) %>% 
  summarise(bWAR = sum(WAR))
colnames(batters_bWAR)[c(1,2)] <- c('playerID', 'yearID')

## all batting stats except fWAR are merged together
batters_stats_bWAR <- merge(batters_stats, batters_bWAR, 
                            by = c('playerID', 'yearID'))

## fWAR
setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/expansion/updated_data/batters/WAR/data")
batting_fWAR_2021 <- read.csv('batters_fWAR.csv')[,-1]
batting_fWAR_2022 <- rbind(cbind(read.csv('batters_fWAR_2022_AL.csv'), 
                                 yearID = 2022, lgID = 'AL'), 
                           cbind(read.csv('batters_fWAR_2022_NL.csv'), 
                                 yearID = 2022, lgID = 'NL'))

batting_fWAR <- rbind(batters_fWAR_2021, batters_fWAR_2022)

batters_fWAR <- batting_fWAR %>% group_by(playerid, yearID) %>% summarise(
  name = unique(Name), fWAR = sum(WAR)
)

## merge fangraph ID with baseball reference ID

setwd("~/Documents/GitHub/register/data/")
people <- read.csv('people.csv') %>% select(key_bbref, key_fangraphs)

bbref_fangraphs_ID <- people[!is.na(people$key_fangraphs),]
colnames(bbref_fangraphs_ID)[2] <- 'playerid'

batters_fWAR_bWAR <- merge(batters_fWAR, bbref_fangraphs_ID, by = 'playerid')

## combine fWAR with batting stats and bWAR

colnames(batters_fWAR_bWAR)[5] <- 'playerID'

batters <- merge(batters_stats_bWAR, 
                 batters_fWAR_bWAR %>% select(playerID, yearID, fWAR), 
                 by = c('playerID', 'yearID'))

setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/expansion/updated_data/batters")
write.csv(batters, 'batters.csv')

source("helper_functions.R")
cr_check_packages()

finalYear <- 2022

## NOTE: this script is provided so as to illustrate the importation
## and broad processing steps taken in preparing the data for
## modelling. Access to the primary database and datasets is not
## provided within this repository. Access to the primary data can be
## requested by emailing the senior author directly.

## Manta
 {
   ## get data   
   {
     ## SQL for extracting manta tow data from oracle
     writeLines("
SELECT V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME,V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.REEF_LAT,
  V_RM_SAMPLE.REPORT_YEAR, RM_MANTA.TOW_SEQ_NO, RM_MANTA.LIVE_CORAL, V_RM_SAMPLE.SAMPLE_CLASS
FROM RM_MANTA INNER JOIN V_RM_SAMPLE ON RM_MANTA.SAMPLE_ID = V_RM_SAMPLE.SAMPLE_ID
WHERE (((V_RM_SAMPLE.SAMPLE_CLASS) In ("K","C","G","Z") Or (V_RM_SAMPLE.SAMPLE_CLASS) Is Null))
ORDER BY V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LAT, V_RM_SAMPLE.REPORT_YEAR,
RM_MANTA.TOW_SEQ_NO", paste0(DATA_PATH, "primary/manta.sql"))

     ## if (goto_database_manta) system(paste0("java -jar ", DATA_PATH, "scripts/dbExport.jar ",
     system(paste0("java -jar dbExport.jar ",
       DATA_PATH, "primary/manta.sql ",
       DATA_PATH, "primary/manta.csv reef reefmon"))

     manta <- read.csv(paste0(DATA_PATH, "primary/manta.csv"),strip.white = TRUE)

     save(manta, file = paste0(DATA_PATH, "primary/manta.RData"))
   }
   ## Process data
   {
     pcode.mod <- read.csv(paste0(DATA_PATH, "primary/P_CODE_MOD.csv"),
       strip.white = TRUE
     )
     manta <- manta |>
       left_join(pcode.mod) |>
       droplevels()

     #####################################################################
     ## Data are collected per Manta Tow.                               ##
     ## The most appropriate unit for these analyses is the reef level. ##
     ## - Only include data collected after 1985                        ##
     ## - Convert data into percent cover                               ##
     ## - Summarize the data to the reef/year level                     ##
     #####################################################################
     ## ---- processManta_towlevel
     manta.tow <- manta |>
       filter(REPORT_YEAR > 1985, REPORT_YEAR < (finalYear + 1)) |>
       droplevels() |>
       mutate(Cover = CoralTrends_calcPercent(LIVE_CORAL)) |>
       mutate(Latitude = REEF_LAT, Longitude = REEF_LONG) |>
       CoralTrends_calc3ZoneLocation() |>
       filter(!is.na(Region)) |>
       droplevels() |>
       mutate(Region = factor(Region,
         levels = c("Northern GBR", "Central GBR", "Southern GBR")),
         Zone = Region) |>
       mutate(Year = factor(REPORT_YEAR)) |>
       mutate(Cover = ifelse(Cover == 0, 0.01, Cover)) |>
       mutate(P_CODE.mod = factor(ifelse(is.na(P_CODE.mod), "Other", P_CODE.mod))) |>
       as.data.frame()
     ## ----end
     ## ---- processManta_aggregateToReef
     manta.sum <-
       manta |>
       filter(REPORT_YEAR > 1985,
         REPORT_YEAR < (finalYear + 1)) |>
       droplevels() |>
       mutate(Cover = CoralTrends_calcPercent(LIVE_CORAL)) |>
       group_by(P_CODE.mod, A_SECTOR, SHELF, REEF_NAME, REEF_ID, REPORT_YEAR) |>
       summarise(Cover = mean(Cover, na.rm = TRUE),
         CoverCat.median = median_cover_cat(LIVE_CORAL),
         Tows = length(unique(TOW_SEQ_NO)),
         Latitude = mean(REEF_LAT, na.rm = TRUE),
         Longitude = mean(REEF_LONG, na.rm = TRUE)) |>
       mutate(Cover_from_cat = CoralTrends_calcPercent(CoverCat.median)) |>
       ungroup() |>
       group_by(REEF_NAME) |>
       ungroup() |>
       CoralTrends_calc3ZoneLocation() |>
       filter(!is.na(Region)) |>
       droplevels() |>
       mutate(Region = factor(Region, levels = c("Northern GBR", "Central GBR", "Southern GBR")),
         Zone = Region) |>
       as.data.frame()
     ## ----end

     save(manta.sum, file = "../data/processed/manta.sum.RData")
     save(manta.tow, file = "../data/processed/manta.tow.RData")

   }
 }

load('../data/processed/manta.sum.RData')

## Bleaching
{
## Reef lookup
  {
    lookup.reefs <-
      manta.sum |>
      dplyr::select(REEF_ID, REEF_NAME, Latitude, Longitude, Zone) |>
      distinct()
  }
   ## Bleaching from Mike (LTMP + COE)
  {
    bleaching <- read_csv("../data/primary/bleach.csv") |>
      ## The following forces the listed reefs to be in central (which is what their lat/long would do)
      ## the data do not have lat/long, Mike already put them into Regions - however his method was not
      ## consistent with De"ath
      mutate(Region = ifelse(Reef %in% c("BOULDER REEF","15077S","EDGELL REEFS (NO 5)",
        "15047S", "EGRET REEF", "IRENE REEF", "LENA REEF",
        "RIBBON NO 1 REEF", "RIBBON NO 3 REEF",
        "ROSSER REEF"),
        "Central", Region)) |>
      mutate(Region = ifelse(Full.Reef.ID %in% c("20113S"), "Central", Region)) |>
      mutate(Zone = Region,
        Location = Region,
        Severe = ifelse(mid.point >= 0.6, 1, 0),
        BleachingCAT = COE_BLEACHINGcategories(mid.point)) |>
      filter(!is.na(mid.point),
        REPORT_YEAR < (finalYear + 1)) |> droplevels() |>
      group_by(Zone, REPORT_YEAR, Reef, Full.Reef.ID) |>
      mutate(n = n()) |>
      ## NEW....  Retain the highest of COE and LTMP
      mutate(flag = ifelse(any(Project == "Aerial"), 1,0),
        flag2 = ifelse(Project == "LTMP", 1, 0),
        flag3 = ifelse(n>1, 1, 0),
        flag4 = ifelse(flag & flag2 & flag3, 1, 0),
        Date = as.Date(SAMPLE_DATE, "%d-%b-%y"),
        Season = ifelse(months(Date) %in% c("December","January","February",
          "March"), "Bleaching", "Non-bleaching")) |>
      filter(flag4 == 0) |>
      dplyr::select(-starts_with("flag"))
  }
  ## Process data
  {
    bleaching <-
      bleaching |>
      mutate(REEF_ID = stringr::str_sub(as.character(Full.Reef.ID),
        start = 1, end = 5)) |> 
      mutate(BleachingCAT = ifelse(!is.na(BleachingCAT),B leachingCAT,
        ifelse(!is.na(RAYcat), RAYcat, NA))) |>
      mutate(Project = ifelse(is.na(Project), 'RB', Project)) |>
      mutate(Reef = ifelse(is.na(Reef), REEF_ID, Reef))
  }
  ## Bleaching full
  {
    bleaching.full_3Zone <- bleaching |>
      ungroup() |>
      dplyr::select(REEF_ID, # = Full.Reef.ID,
        REPORT_YEAR,
        REEF_NAME = Reef,
        BLEACHINGcat = BleachingCAT,
        Zone,
        Date,
        Season,
        Project)
    save(bleaching.full_3Zone, file="../data/modelled/bleaching.full_3Zone.RData")
  }
  ## Bleaching Sum
  {
    bleaching.sum.all <-
      bleaching |>
      group_by(Zone, REPORT_YEAR, BleachingCAT) |>
      summarise(BLEACHING = n(),
        BLEACHING.type2 = sum((Project == "LTMP" & Season == "Bleaching") |
                                Project %in% c("Aerial", "RB")), na.rm = TRUE,
        BLEACHING.type2 = ifelse(is.na(BLEACHING.type2), 0, BLEACHING.type2)) |>
      group_by(Zone, REPORT_YEAR) |>
      mutate(N = sum(BLEACHING),
        BLEACHING.p = BLEACHING/N,
        N.type2 = sum(BLEACHING.type2),
        BLEACHING.type2.p = BLEACHING.type2/N.type2,
        BLEACHING.type2.p = ifelse(is.na(BLEACHING.type2.p), 0, BLEACHING.type2.p))|>
      mutate(Zone = factor(Zone, levels = c("Northern", "Central", "Southern"),
        labels = c("Northern GBR", "Central GBR", "Southern GBR"))) |>
      arrange(Zone, REPORT_YEAR, BleachingCAT) |>
      dplyr::select(Zone,REPORT_YEAR, BLEACHINGcat = BleachingCAT,
        BLEACHING, N, BLEACHING.p,
        BLEACHING.type2, N.type2, BLEACHING.type2.p) |>
      mutate(Location = Zone) |>
      mutate(BLEACHINGcat = factor(BLEACHINGcat),
        BLEACHING.type2 = factor(BLEACHING.type2))
    save(bleaching.sum.all, file="../data/modelled/bleaching.sum.all_3Zone.RData")    
  }
}
## COTS
{
  ## Get data
  {
    writeLines("
SELECT V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_LAT,
V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.OPENORCLOSED, V_RM_SAMPLE.OPENORCLOSED_AFTER2004, V_RM_SAMPLE.REPORT_YEAR, V_RM_SAMPLE.VISIT_NO,
Avg(RM_MEDIAN.MEAN_LIVE) AS AvgOfMEAN_LIVE, Avg(RM_MEDIAN.MEAN_COTS) AS AvgOfMEAN_COTS
FROM V_RM_SAMPLE INNER JOIN RM_MEDIAN ON V_RM_SAMPLE.SAMPLE_ID = RM_MEDIAN.SAMPLE_ID
GROUP BY V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_LAT,
V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.OPENORCLOSED, V_RM_SAMPLE.OPENORCLOSED_AFTER2004, V_RM_SAMPLE.REPORT_YEAR, V_RM_SAMPLE.VISIT_NO
HAVING (((V_RM_SAMPLE.P_CODE) Not Like 'TS' And (V_RM_SAMPLE.P_CODE) Not Like 'WA_NI'))
ORDER BY V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REPORT_YEAR", "../data/primary/cots.sql")

    system("java -jar dbExport.jar ../data/primary/cots.sql ../data/primary/cots.csv reef reefmon")

    cots <- read.csv("../data/primary/cots.csv", strip.white = TRUE) |>
      filter(REPORT_YEAR < (finalYear + 1)) |>
      droplevels()
    summary(cots)
  }
  ## Process data
  {
    cots.sum <-
      cots |>
      filter(REEF_NAME %in% as.character(unique(manta.sum$REEF_NAME))) |>
      mutate(COTScat = COTScategories(AVGOFMEAN_COTS),
        COTScat = factor(COTScat, levels = c("Zero", "NO", "IO", "AO"))) |>
      left_join(manta.sum |>
                  dplyr:::select(REEF_NAME, Zone)) |>
      group_by(Zone, REPORT_YEAR) |>
      mutate(N = n()) |>
      group_by(Zone, REPORT_YEAR, COTScat) |>
      summarize(COTS = n(), N = mean(N), COTS.p = 100 * COTS/N)

    cots.full <- cots |>
      filter(REEF_NAME %in% as.character(unique(manta.sum$REEF_NAME))) |>
      mutate(COTScat = COTScategories(AVGOFMEAN_COTS),
        COTScat = factor(COTScat, levels = c("Zero", "NO", "IO", "AO"))) |>
      left_join(manta.sum |>
                  dplyr:::select(REEF_NAME, Zone)) |>
      dplyr:::select(REEF_ID, REEF_NAME, Zone, REPORT_YEAR, COTScat)
    save(cots.full, file = "../data/modelled/cots.full_3Zone.RData")

    cots.all <-
      cots |>
      filter(REEF_NAME %in% as.character(unique(manta.sum$REEF_NAME))) |>
      mutate(COTScat = COTScategories(AVGOFMEAN_COTS),
        COTScat = factor(COTScat, levels = c("Zero", "NO", "IO", "AO"))) |>
      left_join(manta.sum |>
                  dplyr:::select(REEF_NAME, Zone)) |>
      group_by(REPORT_YEAR) |>
        mutate(N = n()) |>
        group_by(REPORT_YEAR, COTScat) |>
      summarize(COTS = n(), N = mean(N), COTS.p = 100 * COTS / N) |>
      mutate(Zone = "All")

    cots.sum.all <- rbind(cots.sum, cots.all)
    cots.sum.all$Location <- factor(cots.sum.all$Zone,
      levels = c("All", "Northern GBR", "Central GBR", "Southern GBR"),
      labels = c("Great Barrier Reef", "Northern GBR", "Central GBR", "Southern GBR"))

    save(cots.sum.all, file = "../data/modelled/cots.sum.all_3Zone.RData")
  }
}
## Cyclones
{
  ## Reef lookup
  {
    all.reefs <-
      manta.sum |>
      dplyr:::select(REEF_NAME, REEF_ID, Latitude, Longitude, Zone) |>
      group_by(REEF_NAME, REEF_ID, Zone) |>
      summarize_at(vars(Latitude, Longitude), funs(mean)) |>
      as.data.frame()
  }
  ## Data
  {
    cyclones.all <-
      read.csv("../data/primary/20221114 Cyclone wave data from Marji.csv",
        strip.white = TRUE) |>
      dplyr::select(-max)

    cyclones <-
      cyclones.all |>
      dplyr::rename(REEF_NAME = Reef,A_SECTOR = TMP_sector, SHELF = Shelf,
        REEF_LAT = lat, REEF_LONG = long_) |>
      filter(!REEF_NAME == "") |> 
      dplyr::select(-Project, -gbrmpa_sector, -full_reef_id,
        -gazetted_name, -GBRMPA_ID) |> 
      gather(key = REPORT_YEAR, value = S, -A_SECTOR:-REEF_LONG) |>
      mutate(REPORT_YEAR = as.numeric(as.character(gsub("X", "", REPORT_YEAR)))) |>
      mutate(CYCLONEcat = S)
  }
## Processing
  {
    cyclones <-
      all.reefs  |>
    left_join(cyclones |> dplyr::select(-A_SECTOR,-SHELF)) |>
    filter(!is.na(CYCLONEcat))

    cyclones.full <-
      cyclones |>
      dplyr::select(REEF_ID, REEF_NAME, REPORT_YEAR, CYCLONEcat)
    save(cyclones.full, file = "../data/modelled/cyclones.full_3Zone.RData")

    cyclones$Location <- factor(cyclones$Zone,
      levels = c("Northern GBR", "Central GBR", "Southern GBR"),
      labels = c("Northern GBR", "Central GBR", "Southern GBR")
    )

## Fill in the gaps
cyclones.lookup <-
  expand.grid(
    Location = c("Northern GBR", "Central GBR", "Southern GBR"),
    REPORT_YEAR = seq.int(1985, 2022, by = 1)
  )

    cyclones <-
      cyclones |>
      full_join(cyclones.lookup) |>
      arrange(REEF_NAME, REPORT_YEAR)

    cyclones <-
      cyclones |>
      mutate(Zone = ifelse(Location == "Northern GBR", "Northern GBR",
        ifelse(Location == "Central GBR", "Central GBR",
          ifelse(Location == "Southern GBR", "Southern GBR", "Great Barrier Reef")
        )
      ))
  }
  {
    ## Generate a summary that calculates the number and percentage of reefs in each zone per year
    ## that are impacted by each level of severity category
    cyclones.sum <-
      cyclones |>
      group_by(REPORT_YEAR,Zone) |>
      mutate(N = n()) |>
      group_by(Zone, REPORT_YEAR, CYCLONEcat) |>
      summarize(CYCLONE = n(),
        N = mean(N),
        CYCLONE.p = 100 * CYCLONE/N) |>
      mutate(CYCLONE = ifelse(is.na(CYCLONEcat), NA, CYCLONE),
        N = ifelse(is.na(CYCLONEcat), NA, N),
        CYCLONE.p = ifelse(is.na(CYCLONEcat), 0, CYCLONE.p))
    ## as above but only for all the whole GBR
    cyclones.all <-
      cyclones |>
      ungroup() |>
      droplevels() |>
      group_by(REPORT_YEAR) |>
      mutate(N = n()) |>
      group_by(REPORT_YEAR, CYCLONEcat) |>
      summarize(CYCLONE = n(), N = mean(N), CYCLONE.p = 100 * CYCLONE / N) |>
      mutate(
        CYCLONE = ifelse(is.na(CYCLONEcat), NA, CYCLONE),
        N = ifelse(is.na(CYCLONEcat), NA, N),
        CYCLONE.p = ifelse(is.na(CYCLONEcat), 0, CYCLONE.p)
      ) |>
      mutate(Zone = "All")

    cyclones.sum.all <- rbind(cyclones.sum, cyclones.all)
    cyclones.sum.all <-
      cyclones.sum.all |>
      ungroup() |>
      mutate(
        REPORT_YEAR = as.integer(REPORT_YEAR),
        CYCLONEcat = factor(CYCLONEcat)
      )
    cyclones.sum.all$Location <- factor(cyclones.sum.all$Zone,
      levels = c("All", "Northern GBR", "Central GBR", "Southern GBR"),
      labels = c("Great Barrier Reef", "Northern GBR", "Central GBR", "Southern GBR")
    )

    ## Fill in the gaps
    cyclones.sum.all <-
      cyclones.sum.all |>
      mutate(CYCLONEcat = ifelse(is.na(CYCLONEcat), "0", as.character(CYCLONEcat)),
        CYCLONE.p = ifelse(is.na(CYCLONE.p), 0, CYCLONE.p)) |>
      droplevels()
    save(cyclones.sum.all, file = "../data/modelled/cyclones.sum.all_3Zone.RData")
  }

}

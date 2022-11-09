read_hobo_csv2 <- function (csv_file, units_out = c("as.is", "metric", "imperial")) 
{
  units_out <- match.arg(units_out)
  con <- file(csv_file, encoding = "UTF-8")
  header <- readLines(con = con, n = 2)
  close(con)
  header_bits <- unlist(strsplit(header[2], "\",\\\""))
  SNs <- stringr::str_extract(header_bits, "(?<=S\\/N:\\s)[0-9]+")
  SN <- unique(SNs[!is.na(SNs)])
  if (length(SN) > 1) 
    stop("multiple serial numbers in file header")
  tz <- stringr::str_extract(header_bits[grep("Date Time", 
                                              header_bits)], "GMT[+-][0-9][0-9]")
  olsonsign <- chartr("+-", "-+", substr(tz,4,4)) #POSIX Etc/GMT have the sign reversed, see ?timezones
  if (substr(tz, 5, 5) == 0) {
    olsontz <- paste("Etc/", substr(tz, 1, 3), olsonsign, substr(tz, 6, 6), sep = "")
  }
  else {
    olsontz <- paste("Etc/", substr(tz, 1, 3), olsonsign, substr(tz, 5, 5), sep = "")
  }
  hobofile <- read.csv(csv_file, skip = 2, header = FALSE, 
                       stringsAsFactors = FALSE, na.strings = "")
  hobofile$timestamp <- lubridate::mdy_hms(hobofile[, 2], 
                                           tz = olsontz)
  df_out <- data.frame(Timestamp = hobofile$timestamp, Logger.SN = rep(SN, 
                                                                       nrow(hobofile)))
  hobofile <- tidyr::separate(hobofile, timestamp, c("Year", 
                                                     "Month", "Day", "Hour", "Minute", "Second"), remove = FALSE, 
                              convert = TRUE)
  hobofile$tz <- rep(tz, nrow(hobofile))
  tempregex <- "Avg: (Hourly )?Temperature" #average hourly temp, not min or max for day or hour
  #TODO modify regex for percivals in 109 (column just says Temp). Don't match max/min
  if (length(grep(tempregex , header_bits)) > 0) {
    temp <- hobofile[, grep(tempregex, header_bits)]
    # units::units(temp) <- dplyr::case_when(any(grepl("Temp, .F", 
    #                                           header_bits)) ~ "fahrenheit", any(grepl("Temp, .C", 
    #                                                                                   header_bits)) ~ "celsius")
    # units::units(temp) <- with(units::ud_units, dplyr::case_when(units_out == 
    #                                                         "metric" ~ "celsius", units_out == "imperial" ~ 
    #                                                         "fahrenheit", TRUE ~ toString(units(temp))))
    df_out$Temp <- temp 
  }
  if (length(grep("RH", header_bits)) > 0) {
    rh <- hobofile[, grep("RH", header_bits)]
    df_out$RH <- rh
  }
  if (length(grep("Intensity", header_bits)) > 0) {
    illum <- hobofile[, grep("Intensity", header_bits)]
    # units::units(illum) <- dplyr::case_when(any(grepl("Intensity, lum/ft", 
    #                                            header_bits)) ~ "lumen/ft^2", any(grepl("Intensity, Lux", 
    #                                                                                    header_bits)) ~ "lux")
    # units::units(illum) <- with(units::ud_units, dplyr::case_when(units_out == 
    #                                                          "metric" ~ "lux", units_out == "imperial" ~ "lumen/ft^2", 
    #                                                        TRUE ~ toString(units(illum))))
    df_out$Illum <- illum
  }
  df_out <- cbind(subset(hobofile, select = c("Year", "Month", 
                                              "Day", "Hour", "Minute", "Second", "tz")), df_out)
  df_env <- df_out#[complete.cases(df_out), ]
  logger_events <- numeric(4)
  HOBO_names <- c("Host Connected", "Coupler Detached", "Coupler Attached", 
                  "End Of File", "Stopped")
  df_names <- stringr::str_replace_all(HOBO_names, " ", "")
  for (i in 1:5) {
    if (length(grep(HOBO_names[i], header_bits)) > 0) {
      names(hobofile)[grep(HOBO_names[i], header_bits)] <- df_names[i]
      logger_events[i] <- grep(HOBO_names[i], header_bits)
    }
  }
  if (sum(logger_events) > 0) {
    df_logger <- tidyr::gather(hobofile, logger, event, 
                               logger_events[logger_events > 0], factor_key = TRUE)
    df_logger <- subset(df_logger[!is.na(df_logger$event), 
    ], select = c("timestamp", "logger"))
    df_logger$Logger.SN = rep(SN, nrow(df_logger))
  }
  else {
    df_logger <- NULL
  }
  df_units_base <- data.frame(variable = c("Temp", "RH", "Illum"),
                              unit = c(ifelse(exists("temp"), "temp",
                                              NA), ifelse(exists("rh"), "percent (%)", NA), ifelse(exists("illum"),
                                                                                                   "illum", NA)), stringsAsFactors = FALSE)
  #Encoding(df_units_base$unit) <- rep("UTF-8", nrow(df_units_base))
  df_units <- df_units_base[complete.cases(df_units_base), ]
  return(structure(list(df_env = df_env, df_logger = df_logger, 
                        df_units = NULL), class = "microclim"))
}

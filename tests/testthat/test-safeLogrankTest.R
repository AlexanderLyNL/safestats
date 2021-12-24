test_that("Computation logrank z statistic left truncation works", {
  handResult <- c(-0.81649658093, -1.48820571771, -0.77208803219, -0.75021409365)

  # Example with left trunctation due to Judith ter Schure
  designObj <- designSafeLogrank(hrMin=1/2)

  enrollment <- 10     # 5 treatment, 5 placebo
  lambdaC <- 0.03943723
  hr1 <- 0.5           # hazard ratio between treatment en placebo group
  fup <- 40            # folow up of 40 days
  data <- generateSurvData(nP = 5, nT = 5, lambdaP = lambdaC, lambdaT = hr1*lambdaC,
                           endTime = fup, seed = 2006)

  # Add different time of randomisation
  dateRandStart <- as.Date("2020-05-04")
  dateRandEnd <- as.Date("2020-05-15")

  set.seed(2005)
  data$"dateRand" <- sample(seq.Date(from = dateRandStart, to = dateRandEnd, by = "day"),
                            size = enrollment, replace = TRUE)
  data$"dateEvent/LastFup" <- as.Date(data$dateRand + data$time)
  data$"dateLastFup" <- as.Date("2020-06-15")
  data$"participantID" <- 1:nrow(data)
  data$"participantID"[order(data$"dateRand")] <- 1:nrow(data)
  data <- data[order(data$"dateRand"), ]
  data$time <- data$"dateEvent/LastFup" - dateRandStart

  # Add additional complication with multiple events at the same time
  data$"dateEvent/LastFup"[data$participantID == 3] <-
    data$"dateEvent/LastFup"[data$participantID == 5]
  data$time[data$participantID == 3] <-
    data$"dateEvent/LastFup"[data$participantID == 3] - dateRandStart
  data$dateRand[data$participantID == 3] <-
    data$"dateEvent/LastFup"[data$participantID == 3] -
    data$time[data$participantID == 3]

  # Interim analyses events 1 to 4

  result <- numeric(length(handResult))

  for (i in seq_along(result)) {
    calDate <- sort(data$"dateEvent/LastFup")[i]

    dataSoFar <- data[data$dateRand < calDate, ]
    dataSoFar$dateLastFup <- calDate
    dataSoFar$time <- pmin(dataSoFar$time, calDate - dateRandStart)
    dataSoFar$status[dataSoFar$"dateEvent/LastFup" > calDate] <- 1
    survObj <- survival::Surv(time = dataSoFar$dateRand - dateRandStart,
                              time2 = dataSoFar$time,
                              event = dataSoFar$status,
                              type = "counting")

    interimResult <- safeLogrankTest(survObj ~ dataSoFar$group,
                                     designObj = designObj, exact=FALSE)

    result[i] <- interimResult$statistic
  }

  expect_equal("object"=result, "expected"=handResult)

  # Test that Gaussian confidence interval (based on summary statistics) are the same as the sequential one
  sumStatResult <- safeLogrankTestStat(z=interimResult$sumStats$z,
                                       nEvents=interimResult$sumStats$nEvents,
                                       designObj=designObj)

  expect_equal("object"=interimResult$confSeq, "expected"=sumStatResult$confSeq)
})

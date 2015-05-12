#inFolder <- "~/Dropbox/Aaron_work/cuy_migration/models/"
#inFolder <- "~/src/"
inFolder <- getwd() # store the current working directory

# Function to set the important parameters as global variables.
initialize <- function(cuysamplesizeIn = 10, vectorsamplesizeIn = 1000,
                       dogsamplesizeIn = 2, humansamplesizeIn = 4,
                       chickensamplesizeIn = 3,
                       ncuyinfectedIn = 1, nbuginfectedIn = 0,
                       feedIntervalIn = 14,
                       incubationIn = 45, BtoCInfectProbIn = 0.00068,
                       BtoDInfectProbIn = 0.00068, # prob of infecting a dog
                       BtoHInfectProbIn = 0.00058, # prob of infecting a human
                       cuyLifetimeIn = 120, dogLifetimeIn = 1825, # 5 years
                       bugLifetimeIn = 200,
                       humanLifetimeIn = 18250, # 50 years
                       probSuperCuyIn = 0,
                       cuyPersistTimeIn = 56, # days before infection clears
                       cuyPrefFactorIn = 2.5, dogPrefFactorIn = 2.5,
                       chickenPrefFactorIn = 4.8,
                       bugBottleStartTimeIn = -1,
                       bugBottlePopSizeIn = 0,
                       bugBottleEndTimeIn = 0,
                       bugBottlePopRegainTimeIn = 0) {

  # Store all inputs as global variables
  cuysamplesize <<- cuysamplesizeIn
  vectorsamplesize <<- vectorsamplesizeIn
  dogsamplesize <<- dogsamplesizeIn
  humansamplesize <<- humansamplesizeIn
  chickensamplesize <<- chickensamplesizeIn
  ncuyinfected <<- ncuyinfectedIn
  nbuginfected <<- nbuginfectedIn
  feedInterval <<- feedIntervalIn
  incubation <<- incubationIn
  BtoCInfectProb <<- BtoCInfectProbIn
  BtoDInfectProb <<- BtoDInfectProbIn
  BtoHInfectProb <<- BtoHInfectProbIn
  cuyLifetime <<- cuyLifetimeIn
  dogLifetime <<- dogLifetimeIn
  bugLifetime <<- bugLifetimeIn
  probSuperCuy <<- probSuperCuyIn
  cuyPersistTime <<- cuyPersistTimeIn
  cuyPrefFactor <<- cuyPrefFactorIn
  dogPrefFactor <<- dogPrefFactorIn
  chickenPrefFactor <<- chickenPrefFactorIn
  humanLifetime <<- humanLifetimeIn

  # ==========================================================
  # These parameters related to bug (vector) bottlenecks.
  # bugBottleStartTime is the simulation day when the bug bottleneck
  # starts. If set to -1, there is no bug bottleneck.
  bugBottleStartTime <<- bugBottleStartTimeIn
  # bugBottlePopSize is the size to which the vector population declines
  # during a bottleneck.
  bugBottlePopSize <<- bugBottlePopSizeIn
  # bugBottleEndTime is the simulation day when the bug bottleneck ends.
  # After this day, the bug population gradually increases until it
  # regains its maximum level on day bugBottlePopRegainTime
  bugBottleEndTime <<- bugBottleEndTimeIn
  bugBottlePopRegainTime <<- bugBottlePopRegainTimeIn
  # Make sure the bug bottleneck inputs are reasonable:
  if (bugBottleStartTime > 0) {
    if (bugBottleEndTime <= bugBottleStartTime)
      stop("bugBottleEndTime must be greater than bugBottleStartTime")
    if (bugBottlePopRegainTime <= bugBottleEndTime)
      stop("bugBottlePopRegainTime must be greater than bugBottleEndTime")
  }
  # ==========================================================
  
  # Store the integrated human infectivity over the course of the lifetime:
  integratedHumanInfect <<- integrateHumanInfectivity()

  # Set the probability of feeding on each host type
  setFeedingProbs()
}

# Model that includes vectors, guinea pigs, dogs, humans, and chickens.
# Assumes that no humans or dogs are initially infected.
doSim <- function(nsims = 5, simLength = 365, makeOutput=FALSE, outdir="none",
                  vectorDataPrefix="vector_prev", cuyDataPrefix="cuy_prev",
                  dogDataPrefix="dog_prev", humanDataPrefix="human_prev") {

  # If the global parameters have not been set, initialize with the defaults:
  if (!exists("cuysamplesize")) initialize()
  
  # Set some useful standard deviations for our Gaussian distributions:
  sdFeed <- 2 # sd of days between feeding for vectors
  sdBD <- 30 # sd of bug lifetime (days)
  sdCD <- 25 # sd of cuy lifetime (days)
  sdDD <- 365 # sd of dog lifetime (days)

  # Create matrices to store the simulation outputs that we're tracking:
  # outBug, outCuy, outDog, outHuman are the numbers of bugs, cuyes, dogs
  # and humans infected each day.
  outBug <- outCuy <- outDog <- outHuman <-
    matrix(0, nrow=nsims, ncol=simLength)
  
  # Store the total numbers of bugs and hosts that become infected,
  # as well as the lifetime of each (humans assumed to live forever)
  nCuyInfect <- nBugInfect <- nDogInfect <- nHumanInfect <- (1:nsims) * 0
  cuyLife <- bugLife <- dogLife <- 0

  # Store the total number of bug feeding events
  nFeeding <- (1:nsims) * 0

  # Store the number of feeding events on each host type:
  nCuyFeeding <- nDogFeeding <- nChickenFeeding <- nHumanFeeding <- 0
  
  # Do the requested number of simulations
  for (sim in 1:nsims) {

    cat("Doing simulation",sim,"of",nsims,"\n")

    # Set up the initial infection status and birth/death days
    # of all hosts and vectors:

    # bI = 0/1 whether vector is infected;
    # bID = day number when vector infected
    bI <- bID <- c(rep(1, nbuginfected),
                   rep(0, (vectorsamplesize - nbuginfected)))

    # bFD = day when vector last fed (initialize to less than zero since
    # our start day is zero by definition) -- base this off of the
    # fact that the mean interval between feeds is "feedInterval"
    bFD <- round(rnorm(vectorsamplesize, mean=(-feedInterval/2), sd=sdFeed))
    bFD[which(bFD >= 0)] <- -1

    # Choose the "birth day" of all bugs as a uniform distribution between
    # 0 and "bugLifetime" days before the simulation started:
    bBD <- round(runif(vectorsamplesize, -bugLifetime, -1))

    # Set the death day of each bug based on the expected lifetime:
    bDD <- 1:vectorsamplesize * 0
    for (b in 1:vectorsamplesize)
      bDD[b] <- round(rnorm(1, (bugLifetime + bBD[b]), sdBD))
    bDD[which(bDD <= 0)] <- 1

    if (cuysamplesize > 0) {
      # cI, cID, cBD are similar to above but for cuyes. Make sure any
      # initially infected cuyes have a birth date of 0 so that they will
      # live long enough to infect the expected number of vectors.
      cI <- cID <- c(rep(1, ncuyinfected),
                     rep(0, (cuysamplesize - ncuyinfected)))
      cBD <- c(rep(0, ncuyinfected),
               round(runif((cuysamplesize - ncuyinfected), -cuyLifetime, -1)))

      # Store whether or not each cuy is a supershedder:
      cSuper <- (runif(cuysamplesize) < probSuperCuy)

      # Set the death day of each cuy based on the expected lifetime:
      cDD <- 1:cuysamplesize * 0
      for (cuy in 1:cuysamplesize)
        cDD[cuy] <- round(rnorm(1, (cuyLifetime + cBD[cuy]), sdCD))
      cDD[which(cDD <= 0)] <- 1
    }
    
    if (dogsamplesize > 0) {  
      # dI, dID, dBD are similar to above but for dogs. Initially,
      # no dogs are infected.
      dI <- dID <- rep(0, dogsamplesize)
      dBD <- round(runif(dogsamplesize, -dogLifetime, -1))

      # Set the death day of each dog based on the expected lifetime:
      dDD <- 1:dogsamplesize * 0
      for (dog in 1:dogsamplesize)
        dDD[dog] <- round(rnorm(1, (dogLifetime + dBD[dog]), sdDD))
      dDD[which(dDD <= 0)] <- 1
    }

    if (humansamplesize > 0) {
      # Set up initial infection status and day of infection for humans.
      hI <- hID <- rep(0, humansamplesize)
    }
      
    # NOTE: Chickens do not become infected so we don't need to keep
    # track of their infection status or birth/death days. 
    
    # Now loop over all days in the simulation:
    for (day in 1:simLength) {

      # Loop over bugs and determine if they are feeding today:
      for (bug in 1:vectorsamplesize) {

        daysSinceFeed <- day - bFD[bug]
        if (runif(1) < pnorm(daysSinceFeed, feedInterval, sdFeed)) {
          
          # This bug will feed today. Reset day of last feed.
          bFD[bug] <- day

          nFeeding[sim] <- nFeeding[sim] + 1
          
          # Determine on which host type this bug will feed
          hostType <- getHostType()
          
          if (hostType == 1) {

            # Feeding on a cuy.
            cuyToFeed <- ceiling(runif(1, 0, cuysamplesize))
            nCuyFeeding <- nCuyFeeding + 1
            
            if ((bI[bug] == 1) & (cI[cuyToFeed] == 0)) {
              # Vector is infected and cuy is uninfected.
              # Determine whether to infect the cuy (vector is only
              # infectious if it has been infected longer than the
              # incubation period):
              if (((day - bID[bug]) > incubation) &
                  (runif(1) < BtoCInfectProb)) {

                # Infect this cuy today:
                cI[cuyToFeed] <- 1
                cID[cuyToFeed] <- day

                # Increment the number of cuyes that became infected:
                nCuyInfect[sim] <- nCuyInfect[sim] + 1
              }
            }
            else if ((bI[bug] == 0) & (cI[cuyToFeed] == 1)) {
              # Cuy is infected and vector is uninfected.
              # Determine whether to infect the vector.
              if (runif(1) < CtoBInfectProb(cI[cuyToFeed], cID[cuyToFeed],
                                            cSuper[cuyToFeed], day,
                                            cuyPersistTime)) {
                # Infect this bug today:
                bI[bug] <- 1
                bID[bug] <- day

                # Increment the number of bugs that became infected:
                nBugInfect[sim] <- nBugInfect[sim] + 1
              }
            }
          }
          else if (hostType == 2) {

            # Feeding on a dog.
            dogToFeed <- ceiling(runif(1, 0, dogsamplesize))
            nDogFeeding <- nDogFeeding + 1
            if ((bI[bug] == 1) & (dI[dogToFeed] == 0)) {
              # Vector is infected and dog is uninfected.
              # Determine whether to infect the dog (vector is only
              # infectious if it has been infected longer than the
              # incubation period):
              if (((day - bID[bug]) > incubation) &
                  (runif(1) < BtoDInfectProb)) {

                # Infect this dog today:
                dI[dogToFeed] <- 1
                dID[dogToFeed] <- day

                # Increment the number of dogs that became infected:
                nDogInfect[sim] <- nDogInfect[sim] + 1
              }
            }
            else if ((bI[bug] == 0) & (dI[dogToFeed] == 1)) {
              # Dog is infected and vector is uninfected.
              # Determine whether to infect the vector.
              if (runif(1) < DtoBInfectProb(dI[dogToFeed], dID[dogToFeed],
                                            day)) {
                # Infect this bug today:
                bI[bug] <- 1
                bID[bug] <- day
              }
            }
          }
          else if (hostType == 3) {

            # Feeding on a chicken. Do nothing because bugs cannot
            # infect chickens, and vice versa.
            nChickenFeeding <- nChickenFeeding + 1
          }
          else {

            # Feeding on a human.
            humanToFeed <- ceiling(runif(1, 0, humansamplesize))
            nHumanFeeding <- nHumanFeeding + 1
            if ((bI[bug] == 1) & (hI[humanToFeed] == 0)) {
              # Vector is infected and person is uninfected.
              # Determine whether to infect the person (vector is only
              # infectious if it has been infected longer than the
              # incubation period):
              if (((day - bID[bug]) > incubation) &
                  (runif(1) < BtoHInfectProb)) {

                # Infect this person today:
                hI[humanToFeed] <- 1
                hID[humanToFeed] <- day

                # Increment the number of humans that became infected:
                nHumanInfect[sim] <- nHumanInfect[sim] + 1
              }
            }
            else if ((bI[bug] == 0) & (hI[humanToFeed] == 1)) {
              # Human is infected and vector is uninfected.
              # Determine whether to infect the vector.
              if (runif(1) < HtoBInfectProb(hI[humanToFeed], hID[humanToFeed],
                                            day)) {
                # Infect this bug today:
                bI[bug] <- 1
                bID[bug] <- day
              }
            }
          }
        }
        
        # Now determine if this bug should die today.
        # If so, replace it with a new uninfected bug.
        if (day == bDD[bug]) {
          # Reset this bug
          bugLife <- c(bugLife, (day-bBD[bug]))
          #cat("bug",bug,"dying on day",day,"at age",(day-bBD[bug]),"\n")
          bBD[bug] <- bFD[bug] <- day
          bI[bug] <- bID[bug] <- 0
          
          # Choose the death day:
          bDD[bug] <- round(rnorm(1, (bugLifetime + day), sdBD))
        }
      }

      # Now determine if any cuyes or dogs should die today:
      if (cuysamplesize > 0) {
        for (cuy in 1:cuysamplesize) {
          if (day == cDD[cuy]) {
            # Reset this cuy
            cuyLife <- c(cuyLife, (day-cBD[cuy]))
            cI[cuy] <- cID[cuy] <- 0
            cBD[cuy] <- day
            cSuper[cuy] <- (runif(1) < probSuperCuy)
            cDD[cuy] <- round(rnorm(1, (cuyLifetime + day), sdCD))
          }
        }
      }

      if (dogsamplesize > 0) {
        for (dog in 1:dogsamplesize) {
          if (day == dDD[dog]) {
            # Reset this dog
            dogLife <- c(dogLife, (day-dBD[dog]))
            dI[dog] <- dID[dog] <- 0
            dBD[dog] <- day
            dDD[dog] <- round(rnorm(1, (dogLifetime + day), sdDD))
          }
        }
      }

      # Store the number of infected bugs and hosts:
      outBug[sim, day] <- sum(bI)
      if (cuysamplesize > 0)
        outCuy[sim, day] <- sum(cI)
      else
        outCuy[sim, day] <- 0
      if (dogsamplesize > 0)
        outDog[sim, day] <- sum(dI)
      else
        outDog[sim, day] <- 0
      if (humansamplesize > 0)
        outHuman[sim, day] <- sum(hI)
      else
        outHuman[sim, day] <- 0
    }
  }

  # Now make a plot of the simulation output (infected bugs and hosts
  # as a function of time):
  par(mfrow=c(4,1))
  color <- rainbow(nsims)
  # Plot vectors infected
  plot(outBug[1,], ylim=c(0, vectorsamplesize), main="Number of Vectors Infected Each Day", cex=.2, pch="-", ylab="Number Infected", xlab="Day")
  for(phi in 2:nsims)
    points(outBug[phi,], cex=.2, pch="-", col=color[phi])

  # Plot cuyes infected
  plot(outCuy[1,], ylim=c(0, cuysamplesize), main="Number of Cuyes Infected Each Day", cex=.2, pch="-", ylab="Number Infected", xlab="Day")
  for(phi in 2:nsims)
    points(outCuy[phi,], cex=.2, pch="-", col=color[phi])
  
  # Plot dogs infected
  plot(outDog[1,], ylim=c(0, dogsamplesize), main="Number of Dogs Infected Each Day", cex=.2, pch="-", ylab="Number Infected", xlab="Day")
  for(phi in 2:nsims)
    points(outDog[phi,], cex=.2, pch="-", col=color[phi])

  # Plot humans infected
  plot(outHuman[1,], ylim=c(0, humansamplesize), main="Number of Humans Infected Each Day", cex=.2, pch="-", ylab="Number Infected", xlab="Day")
  for(phi in 2:nsims)
    points(outHuman[phi,], cex=.2, pch="-", col=color[phi])

  par(mfrow=c(1,1))
  
  # Now print some statistics
  cat("Number of cuyes infected mean/median:",mean(nCuyInfect),
      median(nCuyInfect),"\n")
  cat("Number of bugs infected mean/median:",mean(nBugInfect),
      median(nBugInfect),"\n")
  cat("Number of dogs infected mean/median:",mean(nDogInfect),
      median(nDogInfect),"\n")
  cat("Number of humans infected mean/median:",mean(nHumanInfect),
      median(nHumanInfect),"\n")
  cat("Mean number of days between vector feedings:",
      (simLength*vectorsamplesize)/mean(nFeeding),"\n")
  cat("Mean lifespan of",(length(bugLife)-1),"vectors:",mean(bugLife[2:length(bugLife)]),"\n")
  cat("Mean lifespan of",(length(cuyLife)-1),"cuyes:",mean(cuyLife[2:length(cuyLife)]),"\n")
  cat("Mean lifespan of",(length(dogLife)-1),"dog(s):",mean(dogLife[2:length(dogLife)]),"\n")

  # Print probabilities of feeding on each host type:
  totalFeeding <- nCuyFeeding + nDogFeeding + nChickenFeeding + nHumanFeeding
  cat("Cuy feedings:",nCuyFeeding,(nCuyFeeding/totalFeeding*100),"%\n")
  cat("Dog feedings:",nDogFeeding,(nDogFeeding/totalFeeding*100),"%\n")
  cat("Chicken feedings:",nChickenFeeding,(nChickenFeeding/totalFeeding*100),"%\n")
  cat("Human feedings:",nHumanFeeding,(nHumanFeeding/totalFeeding*100),"%\n")

  # Store the simulation output
  if (makeOutput) {
    
    if (outdir == "none") {
      outfolder<-"./output/AT"
      timestamp<- gsub(":","-",gsub(" ","_",Sys.time()))
      outdir<-paste(outfolder,timestamp, sep="_")
    }

    # change to directory where files should be created, and save
    # csv files with daily prevalence data for hosts and vectors.
    dir.create(outdir, recursive=TRUE)
    setwd(outdir)
    write.csv(outBug, file = paste(vectorDataPrefix,".csv",sep=""))
    if (cuysamplesize > 0)
      write.csv(outCuy, file = paste(cuyDataPrefix,".csv",sep=""))
    if (dogsamplesize > 0)
      write.csv(outDog, file = paste(dogDataPrefix,".csv",sep=""))
    if (humansamplesize > 0)
      write.csv(outHuman, file = paste(humanDataPrefix,".csv",sep=""))

    # Save the input parameters to a spreadsheet
    parameter_names<-c("nsims", "vectorsamplesize", "cuysamplesize", "dogsamplesize", "chickensamplesize", "humansamplesize", "ncuyinfected", "nbuginfected", "feedInterval", "cuyLifetime", "bugLifetime", "dogLifetime", "probSuperCuy","incubation", "BtoCInfectProb", "BtoDInfectProb", "BtoHInfectProb", "cuyPrefFactor", "dogPrefFactor", "chickenPrefFactor", "cuyPersistTime")
    parameter_values<-c(nsims, vectorsamplesize, cuysamplesize, dogsamplesize, chickensamplesize, humansamplesize, ncuyinfected, nbuginfected, feedInterval, cuyLifetime, bugLifetime, dogLifetime, probSuperCuy, incubation, BtoCInfectProb, BtoDInfectProb, BtoHInfectProb, cuyPrefFactor, dogPrefFactor, chickenPrefFactor, cuyPersistTime)
    parameter_output<-data.frame(parameter_names, parameter_values)
    write.csv(parameter_output, file = "runparameters.csv")

    # Save a spreadsheet with information about how many hosts and
    # vectors became infected during the simulations
    outData <- data.frame(nBugInfect, nCuyInfect, nDogInfect, nHumanInfect)
    write.csv(outData, file="infection_numbers.csv")
    
    # Save plots of the infectivity of hosts to vectors.
    saveInfectivityPlots(cuyPersistTime)
    
    # Return to previous working directory
    setwd(inFolder)
  }
}

# Function returns the probability that the given cuy will be infectious
# to a vector that feeds on this cuy on the given day.
CtoBInfectProb <- function(cI, cID, isSuper, dayNum, cuyPersistTime) {

  if (cI == 0) return(0) # this cuy is not infected
  else {

    # Infected. Find out how many days have passed since infection:
    daysSinceInfection <- dayNum - cID
    if (isSuper) {
      # This is a supershedder cuy.
      # Assume 100% chance of infecting bug during the "cuyPersistTime", then
      # 20% for rest of life:
      if (daysSinceInfection < cuyPersistTime) return(1)
      else return(0.2)
    }
    else {
      # This is a normal, non-supershedder cuy. Assume 100% chance of
      # infecting bug during the "cuyPersistTime", then 0% for rest of life:
      if (daysSinceInfection < cuyPersistTime) return(1)
      else return(0)
    }
  }
}

# Function returns the probability that the given dog will be infectious
# to a vector that feeds on this dog on the given day.
DtoBInfectProb <- function(dI, dID, dayNum) {
  
  if (dI == 0) return(0) # this dog is not infected

  # Based on Gurtler 1996, the probability that an infected dog
  # transmits infection to an uninfected vector with a single bite is
  # 48.7% overall. This probability doesn't seem to vary much with the
  # age of the dog, so assume it is constant over lifetime.

  return(0.487)
}

# Function returns the probability that the human will be infectious
# to a vector feeding on this day.
HtoBInfectProb <- function(hI, hID, dayNum) {

  if (hI == 0) return (0)

  # Gurtler 1996 gives an overall probability of 2.6% that an infected
  # human will infect a vector after a single bite. This probability
  # decreased from 4.2% in children to 0.48% in adults.
  # My guess is that this decrease is due to the fact that children were
  # infected more recently. So let's model it as 4.2% in the first three
  # years of infection, then 0.48% afterward:
  daysSinceInfection <- dayNum - hID
  if (daysSinceInfection < (3*365)) return (0.042)
  else return (0.0048)
}

# Given the number of different types of hosts, and the factor by which
# vectors prefer other hosts over humans, determine on which host type
# the vector will feed.
# Return values are:
# 1 = cuy
# 2 = dog
# 3 = chicken
# 4 = human
getHostType <- function() {

  # If the total weight and feeding probabilities have not been set,
  # set them here:
  if (!exists("totalWeight")) setFeedingProbs()
  
  # Stochastically determine host type for feeding:
  randNum <- runif(1)
  if (randNum <= probCuy) return (1)
  else if (randNum <= (probCuy + probDog)) return (2)
  else if (randNum <= (probCuy + probDog + probChicken)) return (3)
  else return (4)
}

# Set the probability that the vector will feed on each host type.
setFeedingProbs <- function() {
  
  totalWeight <<- humansamplesize + (cuysamplesize * cuyPrefFactor) +
    (dogsamplesize * dogPrefFactor) + (chickensamplesize * chickenPrefFactor)
  probCuy <<- (cuysamplesize * cuyPrefFactor) / totalWeight
  probDog <<- (dogsamplesize * dogPrefFactor) / totalWeight
  probChicken <<- (chickensamplesize * chickenPrefFactor) / totalWeight
  probHuman <<- humansamplesize / totalWeight
  cat("In setFeedingProbs:probCuy=",probCuy,"probDog=",probDog,"probChicken=",probChicken,"probHuman=",probHuman,"\n")
  cat("humansamplesize=",humansamplesize,"cuysamplesize=",cuysamplesize,"dogsamplesize=",dogsamplesize,"chickensamplesize=",chickensamplesize,"\n")
}

# Save plots of the infectivity of each host to vectors.
saveInfectivityPlots <- function(cuyPersistTime) {
  
  #Do vectors
  color <- grey(uniform(nsims)
  # Plot vectors infected
  plot(outBug[1,], ylim=c(0, vectorsamplesize), main="Number of Vectors Infected Each Day", cex=.2, pch="-", ylab="Number Infected", xlab="Day")
  for(phi in 2:nsims)
    points(outBug[phi,], cex=.2, pch="-", col=color[phi])
  
  
  # Do guinea pigs:
  jpeg(file="cuy_infectivity.jpg")
  cuyInfectivity <- 1:120 * 0
  for (day in 1:120) {
    cuyInfectivity[day] <- CtoBInfectProb(1, 1, 0, day, cuyPersistTime)
  }
  plot(1:120, cuyInfectivity, type="l", xlab="Days since infection", ylab="Probability of transmitting infection to vector", main="Guinea Pigs")
  invisible(dev.off()) # close the jpg plot file

  # Do humans:
  jpeg(file="human_infectivity.jpg")
  humanInfectivity <- 1:2555 * 0
  for (day in 1:2555) {
    humanInfectivity[day] <- HtoBInfectProb(1, 1, day)
  }
  plot(1:2555, humanInfectivity, type="l", xlab="Days since infection", ylab="Probability of transmitting infection to vector", main="Humans")
  invisible(dev.off()) # close the jpg plot file

  # Do dogs
  jpeg(file="dog_infectivity.jpg")
  dogInfectivity <- 1:2555 * 0
  for (day in 1:2555) {
    dogInfectivity[day] <- DtoBInfectProb(1, 1, day)
  }
  plot(1:2555, dogInfectivity, type="l", xlab="Days since infection", ylab="Probability of transmitting infection to vector", main="Dogs")
  invisible(dev.off()) # close the jpg plot file
}


cat("Doing simulation with all hosts. Initial 1 infected cuy.\n")
initialize(cuysamplesizeIn = 10, vectorsamplesizeIn = 1000,
           dogsamplesizeIn = 1, humansamplesizeIn = 0,
           chickensamplesizeIn = 0,
           ncuyinfectedIn = 1, nbuginfectedIn = 0)
doSim(nsims = 2, simLength = 365, makeOutput=TRUE, outdir="~/output/multiple_host2/", vectorDataPrefix="vector_prev", cuyDataPrefix="cuy_prev", dogDataPrefix="dog_prev", humanDataPrefix="human_prev")

###end here for no NGM





# See Roberts (2003) for a description of how to create the next generation
# matrix.
createNextGenMatrix <- function() {

  # Our matrix will be at most 4x4 (cuyes, dogs, humans, vectors) since
  # chickens cannot become infected. But it may be smaller if one of the
  # host types has a zero population size. Test for that here:
  if (!exists("cuysamplesize")) {
    initialize()
    setFeedingProbs()
  }
  if (!exists("probCuy")) setFeedingProbs()

  # Initialize storage for the coefficients
  vectorToHost <- NULL
  hostToVector <- NULL
  
  # Now get the coefficients in the matrix for each host type:
  if (cuysamplesize > 0) {
    vectorToHost <- (1/feedInterval)*probCuy*BtoCInfectProb*bugLifetime*exp(-incubation/bugLifetime)
    hostToVector <- (1/feedInterval)*probCuy*(vectorsamplesize/cuysamplesize)*CtoBInfectProb(1, 1, 0, 1, cuyPersistTime)*cuyPersistTime
  }

  if (humansamplesize > 0) {
    vectorToHost <- c(vectorToHost, (1/feedInterval)*probHuman*BtoHInfectProb*bugLifetime*exp(-incubation/bugLifetime))
    hostToVector <- c(hostToVector, (1/feedInterval)*probHuman*(vectorsamplesize/humansamplesize)*integratedHumanInfect)
  }

  if (dogsamplesize > 0) {
    vectorToHost <- c(vectorToHost, (1/feedInterval)*probDog*BtoDInfectProb*bugLifetime*exp(-incubation/bugLifetime))
    hostToVector <- c(hostToVector, (1/feedInterval)*probDog*(vectorsamplesize/dogsamplesize)*DtoBInfectProb(1, 1, 1)*dogLifetime)
  }

  # Now create the next generation matrix:
  numRows <- length(vectorToHost) + 1
  nextGen <- matrix(rep(0, numRows^2), nrow=numRows, ncol=numRows, byrow=TRUE)
  for (host in 1:(numRows - 1)) {
    nextGen[host, numRows] <- vectorToHost[host]
    nextGen[numRows, host] <- hostToVector[host]
  }

  return(nextGen)
}

# Note: this function returns what is sometimes referred to as
# sqrt(R0). The number of secondary infections resulting from a
# single primary infection is given by R0^2.
# Return R0 of the given next generation matrix.
getR0 <- function(nextGenMatrix) {

  eigens <- eigen(nextGenMatrix)
  R0 <- max(eigens$values)
  return(R0)  
}

# Since humans have two different infectivities (immediately after becoming
# infected, versus later), we need to integrate the infectivity over the
# entire infectious period in order to calculate things like R0:
integrateHumanInfectivity <- function() {

  integral <- 0
  for (day in 1:(humanLifetime))
    integral <- integral + HtoBInfectProb(1, 1, day)

  return(integral)
}

# Commands for running on the cluster:
#cat("Doing simulation with cuyes and vectors. Initial 1 infected cuy.\n")
#initialize(cuysamplesizeIn = 10, vectorsamplesizeIn = 1000,
#           dogsamplesizeIn = 0, humansamplesizeIn = 0,
#           chickensamplesizeIn = 0,
#           ncuyinfectedIn = 1, nbuginfectedIn = 0)
#doSim(nsims = 5, simLength = 365, makeOutput=TRUE, outdir="~/output/multiple_host/cuyes_vectors_1cuy_infected/", vectorDataPrefix="vector_prev", cuyDataPrefix="cuy_prev", dogDataPrefix="dog_prev", humanDataPrefix="human_prev")

cat("Doing simulation with all hosts. Initial 1 infected cuy.\n")
initialize(cuysamplesizeIn = 10, vectorsamplesizeIn = 1000,
           dogsamplesizeIn = 1, humansamplesizeIn = 0,
           chickensamplesizeIn = 0,
           ncuyinfectedIn = 1, nbuginfectedIn = 0)
doSim(nsims = 2, simLength = 365, makeOutput=TRUE, outdir=, vectorDataPrefix="vector_prev", cuyDataPrefix="cuy_prev", dogDataPrefix="dog_prev", humanDataPrefix="human_prev")

#cat("Doing simulation with all hosts. Initial 1 infected vector.\n")
#initialize(cuysamplesizeIn = 10, vectorsamplesizeIn = 1000,
#           dogsamplesizeIn = 2, humansamplesizeIn = 0,
#           chickensamplesizeIn = 0,
#           ncuyinfectedIn = 0, nbuginfectedIn = 1)
#doSim(nsims = 5, simLength = 365, makeOutput=TRUE, outdir="~/output/multiple_host/all_hosts_1vector_infected/", vectorDataPrefix="vector_prev", cuyDataPrefix="cuy_prev", dogDataPrefix="dog_prev", humanDataPrefix="human_prev")

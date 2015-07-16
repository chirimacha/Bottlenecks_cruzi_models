# Updated version of the multiple_host_sim.R with new functionality:
# 1) Time steps are now in weeks, rather than days, to make the code
# run faster.
# 2) Vector population size now depends on whether there are guinea pigs
# present. This is based on discussion with Ricardo that in the field we
# generally observe many more vectors when the domicile contains guinea pigs,
# compared to houses that contain only humans and dogs.
# 3) Code can now model guinea pig bottlenecks (which represent events such
# as mass kill-offs for national feasts) and vector bottlenecks (representing
# a pesticide spraying event). User can vary the bottleneck duration and the
# number of hosts/vectors that return after a bottleneck.

# Remaining TO DO items:
# - Consider adding congenital transmission in hosts

# Store the path to this source code:
inFolder <- getwd()
#inFolder <- "~/Dropbox/Aaron_work/cuy_migration/models/"

# Function to set the important parameters as global variables.
initialize <- function(cuyStartPopulationIn = 10, # initial guinea pig population
                       dogPopulationIn = 2, # dog population size
                       humanPopulationIn = 4, # human population size
                       chickenPopulationIn = 3, # chicken population size
                       ncuyinfectedIn = 1,
                       nbuginfectedIn = 0,
                       feedIntervalIn = 2, # weeks between feedings for bugs
                       incubationIn = 6, # weeks before exposed bug is infectious
                       BtoCInfectProbIn = 0.00068, # prob of infecting a cuy w/ one bite
                       BtoDInfectProbIn = 0.00068, # prob of infecting a dog w/ one bite
                       BtoHInfectProbIn = 0.00058, # prob of infecting a human w/ one bite
                       cuyLifetimeIn = 17, # weeks a gp lives, on average (4 mo)
                       dogLifetimeIn = 260, # 5 years (260 weeks)
                       bugLifetimeIn = 28, # weeks
                       humanLifetimeIn = 2600, # 50 years
                       probSuperCuyIn = 0.33, # prob that each gp is a superspreader
                       cuyPersistTimeIn = 8, # weeks before infection clears in normal (non-superspreader) gp
                       cuyPrefFactorIn = 2.5,
                       dogPrefFactorIn = 2.5,
                       chickenPrefFactorIn = 4.8,
                       cuyBottleStartTimeIn = 20,
                       cuyBottlePopSizeIn = 2,
                       cuyBottleEndTimeIn = 40,
                       cuyToReintroduceIn = 8,
                       infectedCuyToReintroduceIn = 0,
                       # The following inputs control the vector population dynamics
                       vectorPopBaseline = 10, # bug population size when no gp present
                       vectorPopCuyMultiplier = 100, # factor by which bug population increases when gp are present
                       
                       bugBottleStartTimeIn = 30, # set to -1 for no bottleneck
                       bugBottlePopSizeIn = 0,
                       bugBottleEndTimeIn = 40,
                       newBugsPerDayAfterBottle = 5, # rate by which bug population increases after bottleneck
                       infectedBugReintroduceAfterBottle = 1 # number of new bugs that should be infected on the week after the bottleneck.
                       ) {

  # Store all inputs as global variables
  cuyStartPopulation <<- cuyStartPopulationIn 
  dogPopulation <<- dogPopulationIn
  humanPopulation <<- humanPopulationIn
  chickenPopulation <<- chickenPopulationIn
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

  # ========================================================
  # These parameters are related to guinea pig bottlenecks, to model
  # the effect of festivals where the gp population declines precipitously.
  # if cuyBottleStartTime is -1, there is no gp bottleneck.
  cuyBottleStartTime <<- cuyBottleStartTimeIn
  cuyBottlePopSize <<- cuyBottlePopSizeIn
  cuyBottleEndTime <<- cuyBottleEndTimeIn
  if (cuyBottleStartTime > 0 & cuyStartPopulation > 0) {
      doCuyBottleneck <<- TRUE
      if (cuyBottleEndTime <= cuyBottleStartTime)
          stop("cuyBottleEndTime must be greater than cuyBottleStartTime")
  }
  else
      doCuyBottleneck <<- FALSE
  
  # The number of new guinea pigs to add to the population at the end
  # of a guinea pig bottleneck.
  cuyToReintroduce <<- cuyToReintroduceIn
  infectedCuyToReintroduce <<- infectedCuyToReintroduceIn
  
  # ==========================================================
  # These parameters related to bug (vector) bottlenecks. Vector bottlenecks
  # can be used to model pesticide spraying or other causes of sharp declines
  # in vector populations.
  # bugBottleStartTime is the simulation week when the bug bottleneck
  # starts. If set to -1, there is no bug bottleneck.
  bugBottleStartTime <<- bugBottleStartTimeIn
  bugBottlePopSize <<- bugBottlePopSizeIn  # size to which the vector population declines during a bottleneck.
  bugBottleEndTime <<- bugBottleEndTimeIn  # simulation week when bug bottleneck ends.
  vectorPopBaseline <<- vectorPopBaseline # bug population size when no cuyes are present
  vectorPopCuyMultiplier <<- vectorPopCuyMultiplier # factor by which bug population increases when gp are present
  newBugsPerDayAfterBottle <<- newBugsPerDayAfterBottle # rate by which bug population increases after bottleneck
  infectedBugReintroduceAfterBottle <<- infectedBugReintroduceAfterBottle # number of infected bugs that should come back into the population after the bug bottleneck ends.
  
  # Make sure the bug bottleneck inputs are reasonable:
  if (bugBottleStartTime > 0 & bugBottleEndTime <= bugBottleStartTime)
      stop("bugBottleEndTime must be greater than bugBottleStartTime")
  # ==========================================================
  
  # Store the integrated human infectivity over the course of the lifetime:
  integratedHumanInfect <<- integrateHumanInfectivity()
}

# Model that includes vectors, guinea pigs, dogs, humans, and chickens.
# Assumes that no humans or dogs are initially infected.
doSim <- function(nsims = 5,
                  simLength = 52, # weeks
                  makeOutput=FALSE, outdir="none",
                  vectorDataPrefix="vector_prev", cuyDataPrefix="cuy_prev",
                  dogDataPrefix="dog_prev", humanDataPrefix="human_prev") {

  # If the global parameters have not been set, initialize with the defaults:
  if (!exists("cuyStartPopulation")) initialize()
  
  # Set some useful standard deviations for our Gaussian distributions:
  sdFeed <- 1 # sd of weeks between feeding for vectors
  sdBD <- 5 # sd of bug lifetime (weeks)
  sdCD <- 4 # sd of cuy lifetime (weeks)
  sdDD <- 50 # sd of dog lifetime (weeks)

  # Create matrices to store the simulation outputs that we're tracking:
  # infectBug, infectCuy, infectDog, infectHuman are the total number of
  # infected bugs, cuyes, dogs and humans each week.
  # bugPop, cuyPop, dogPop, humanPop are the total number of each host
  # type each week. NOTE: Dog and human pops are assumed to be stable.
  # Cuy and bug populations can vary.
  infectBug <- infectCuy <- infectDog <- infectHuman <-
      bugPop <- cuyPop <- dogPop <- humanPop <- bugPopNew <- matrix(0, nrow=nsims, ncol=simLength)
  
  # Store the total numbers of bugs and hosts that become infected,
  # as well as the lifetime of each (humans assumed to live forever)
  #nCuyInfect <- nBugInfect <- nDogInfect <- nHumanInfect <- (1:nsims) * 0
  #cuyLife <- bugLife <- dogLife <- 0

  # Store the total number of bug feeding events
  #nFeeding <- (1:nsims) * 0

  # Store the number of feeding events on each host type:
  #nCuyFeeding <- nDogFeeding <- nChickenFeeding <- nHumanFeeding <- 0
  
  # Do the requested number of simulations
  for (sim in 1:nsims) {

    cat("Doing simulation",sim,"of",nsims,"\n")

    # Set up the initial infection status and birth/death weeks
    # of all hosts and vectors:

    # First get the starting vector population size -- it depends on whether
    # guinea pigs are present.
    if (cuyStartPopulation > 0)
        vectorPopulation <- vectorPopBaseline*vectorPopCuyMultiplier
    else
        vectorPopulation <- vectorPopBaseline
    bugPopNew[sim, 1] <- vectorPopulation
    
    # bI = 0/1 whether each vector is infected;
    # bID = week number when vector became infected (set to 0 if
    # the bug is not infected)
    bI <- bID <- c(rep(1, nbuginfected),
                   rep(0, (vectorPopulation - nbuginfected)))

    # bFD = week when vector last fed (initialize to less than zero since
    # our start week is zero by definition) -- base this off of the
    # fact that the mean interval between feeds is "feedInterval"
    bFD <- round(rnorm(vectorPopulation, mean=(-feedInterval/2), sd=sdFeed))
    bFD[which(bFD >= 0)] <- -1

    # Choose the "birth week" of all bugs as a uniform distribution between
    # 0 and "bugLifetime" weeks before the simulation started:
    bBD <- round(runif(vectorPopulation, -bugLifetime, -1))

    # Set the death week of each bug based on the expected lifetime:
    bDD <- 1:vectorPopulation * 0
    for (b in 1:vectorPopulation)
      bDD[b] <- round(rnorm(1, (bugLifetime + bBD[b]), sdBD))
    bDD[which(bDD <= 0)] <- 1

    if (cuyStartPopulation > 0) {
      # cI, cID, cBD are similar to above but for cuyes. Make sure any
      # initially infected cuyes have a birth date of 0 so that they will
      # live long enough to infect the expected number of vectors.
      cI <- cID <- c(rep(1, ncuyinfected),
                     rep(0, (cuyStartPopulation - ncuyinfected)))
      cBD <- c(rep(0, ncuyinfected),
               round(runif((cuyStartPopulation - ncuyinfected), -cuyLifetime, -1)))

      # Store whether or not each cuy is a supershedder:
      cSuper <- (runif(cuyStartPopulation) < probSuperCuy)

      # Set the death week of each cuy based on the expected lifetime:
      cDD <- 1:cuyStartPopulation * 0
      for (cuy in 1:cuyStartPopulation)
        cDD[cuy] <- round(rnorm(1, (cuyLifetime + cBD[cuy]), sdCD))
      cDD[which(cDD <= 0)] <- 1

      # Keep track of the current guinea pig population, because it varies...
      cuyCurrentPopulation <<- cuyStartPopulation
    }
    else
        cuyCurrentPopulation <<- 0
    
    if (dogPopulation > 0) {  
      # dI, dID, dBD are similar to above but for dogs. Initially,
      # no dogs are infected.
      dI <- dID <- rep(0, dogPopulation)
      dBD <- round(runif(dogPopulation, -dogLifetime, -1))

      # Set the death week of each dog based on the expected lifetime:
      dDD <- 1:dogPopulation * 0
      for (dog in 1:dogPopulation)
        dDD[dog] <- round(rnorm(1, (dogLifetime + dBD[dog]), sdDD))
      dDD[which(dDD <= 0)] <- 1
    }

    if (humanPopulation > 0) {
      # Set up initial infection status and week of infection for humans.
      hI <- hID <- rep(0, humanPopulation)
    }
      
    # NOTE: Chickens do not become infected so we don't need to keep
    # track of their infection status or birth/death weeks. 
    
    # Now loop over all weeks in the simulation:
    for (week in 1:simLength) {

        # Adjust the guinea pig population size if this is the first
        # week of the guinea pig bottleneck
        if (doCuyBottleneck & week == cuyBottleStartTime) {
            # Reduce the current guinea pig population size:
            cuyCurrentPopulation <<- cuyBottlePopSize
            # Randomly determine which guinea pigs to keep:
            cuyToKeep <- sample(1:length(cI), cuyCurrentPopulation, replace=FALSE)
            cI <- cI[cuyToKeep]
            cID <- cID[cuyToKeep]
            cBD <- cBD[cuyToKeep]
            cSuper <- cSuper[cuyToKeep]
            cDD <- cDD[cuyToKeep]
        }

        # Adjust the guinea pig population if we've reached the
        # end of the guinea pig bottleneck.
        if (doCuyBottleneck & week == cuyBottleEndTime) {
            cuyCurrentPopulation <<- cuyCurrentPopulation + cuyToReintroduce
            # Make all the new cuyes born this week and set death week
            # randomly:
            cBD <- c(cBD, rep(week, cuyToReintroduce))
            cDD <- c(cDD, round(rnorm(cuyToReintroduce, (cuyLifetime + week), sdCD)))
            cDD[which(cDD <= week)] <- (week + 1)
            cI <- c(cI, rep(1, infectedCuyToReintroduce), rep(0, (cuyToReintroduce - infectedCuyToReintroduce)))
            cID <- c(cID, rep(week, infectedCuyToReintroduce), rep(0, (cuyToReintroduce - infectedCuyToReintroduce)))
            cSuper <- c(cSuper, runif(cuyToReintroduce) < probSuperCuy)
        }  

        # Calculate the current vector population size:
        if (bugBottleStartTime > 0 & week >= bugBottleStartTime & week < bugBottleEndTime) {
            # we're within a vector bottleneck.
            currentVectorPop <- bugBottlePopSize
        }
        else {
            # Not within a vector bottleneck
            if (cuyCurrentPopulation > 0)
                currentVectorPop <- vectorPopBaseline*vectorPopCuyMultiplier
            else
                currentVectorPop <- vectorPopBaseline
        }

        # If the current vector population differs from the previous
        # vector population, we must either increase or decrease the
        # number of bugs:
        if (week > 1) {
            if (currentVectorPop > bugPopNew[sim, week-1]) {
                # Increase the number of vectors. Let all new vectors
                # be born this week.
                numNewBugs <- currentVectorPop - bugPopNew[sim, week-1]
                bBD <- c(bBD, rep(week, numNewBugs))
                bDD <- c(bDD, round(rnorm(numNewBugs, (bugLifetime + week), sdBD)))
                bDD[which(bDD <= week)] <- (week + 1)
                bI <- c(bI, rep(1, infectedBugReintroduceAfterBottle), rep(0, (numNewBugs - infectedBugReintroduceAfterBottle)))
                bID <- c(bID, rep(week, infectedBugReintroduceAfterBottle), rep(0, (numNewBugs - infectedBugReintroduceAfterBottle)))
                bFD <- c(bFD, rep((week-1), numNewBugs)) # assume all new bugs have just fed one week ago.
            }
            else if (currentVectorPop < bugPopNew[sim, week-1]) {
                # Decrease the number of vectors.
                if (currentVectorPop == 0) {
                    # Eliminate all vectors during this bug bottleneck.
                    bI <- bID <- bFD <- bBD <- bDD <- 0
                }
                else {
                    # The number of bugs is reduced, but not zero:
                    bugsToKeep <- sample(bugPopNew[sim, week-1], currentVectorPop, replace=FALSE)
                    bI <- bI[bugsToKeep]
                    bID <- bID[bugsToKeep]
                    bFD <- bFD[bugsToKeep]
                    bBD <- bBD[bugsToKeep]
                    bDD <- bDD[bugsToKeep]
            }
            }                
        }
            
        # Set the probability of feeding on each host type
        setFeedingProbs()
        
      # Loop over bugs and determine if they are feeding this week:
        if (currentVectorPop > 0) {
            for (bug in 1:currentVectorPop) {

        weeksSinceFeed <- week - bFD[bug]

        if (runif(1) < pnorm(weeksSinceFeed, feedInterval, sdFeed)) {
          
          # This bug will feed this week. Reset week of last feed.
          bFD[bug] <- week

          #nFeeding[sim] <- nFeeding[sim] + 1
          
          # Determine on which host type this bug will feed
          hostType <- getHostType()
          
          if (hostType == 1) {

            # Feeding on a cuy.
            cuyToFeed <- ceiling(runif(1, 0, cuyCurrentPopulation))
            #nCuyFeeding <- nCuyFeeding + 1
            
            if ((bI[bug] == 1) & (cI[cuyToFeed] == 0)) {
              # Vector is infected and cuy is uninfected.
              # Determine whether to infect the cuy (vector is only
              # infectious if it has been infected longer than the
              # incubation period):
              if (((week - bID[bug]) > incubation) &
                  (runif(1) < BtoCInfectProb)) {

                # Infect this cuy this week:
                cI[cuyToFeed] <- 1
                cID[cuyToFeed] <- week

                # Increment the number of cuyes that became infected:
                #nCuyInfect[sim] <- nCuyInfect[sim] + 1
              }
            }
            else if ((bI[bug] == 0) & (cI[cuyToFeed] == 1)) {
              # Cuy is infected and vector is uninfected.
              # Determine whether to infect the vector.
              if (runif(1) < CtoBInfectProb(cI[cuyToFeed], cID[cuyToFeed],
                                            cSuper[cuyToFeed], week,
                                            cuyPersistTime)) {
                # Infect this bug this week:
                bI[bug] <- 1
                bID[bug] <- week

                # Increment the number of bugs that became infected:
                #nBugInfect[sim] <- nBugInfect[sim] + 1
              }
            }
          }
          else if (hostType == 2) {

            # Feeding on a dog.
            dogToFeed <- ceiling(runif(1, 0, dogPopulation))

            #nDogFeeding <- nDogFeeding + 1
            if ((bI[bug] == 1) & (dI[dogToFeed] == 0)) {
              # Vector is infected and dog is uninfected.
              # Determine whether to infect the dog (vector is only
              # infectious if it has been infected longer than the
              # incubation period):
              if (((week - bID[bug]) > incubation) &
                  (runif(1) < BtoDInfectProb)) {

                # Infect this dog this week:
                dI[dogToFeed] <- 1
                dID[dogToFeed] <- week

                # Increment the number of dogs that became infected:
                #nDogInfect[sim] <- nDogInfect[sim] + 1
              }
            }
            else if ((bI[bug] == 0) & (dI[dogToFeed] == 1)) {
              # Dog is infected and vector is uninfected.
              # Determine whether to infect the vector.
              if (runif(1) < DtoBInfectProb(dI[dogToFeed], dID[dogToFeed],
                                            week)) {
                # Infect this bug this week:
                bI[bug] <- 1
                bID[bug] <- week
              }
            }
          }
          else if (hostType == 3) {

            # Feeding on a chicken. Do nothing because bugs cannot
            # infect chickens, and vice versa.
              #nChickenFeeding <- nChickenFeeding + 1
          }
          else {

            # Feeding on a human.
            humanToFeed <- ceiling(runif(1, 0, humanPopulation))
            #nHumanFeeding <- nHumanFeeding + 1
            if ((bI[bug] == 1) & (hI[humanToFeed] == 0)) {
              # Vector is infected and person is uninfected.
              # Determine whether to infect the person (vector is only
              # infectious if it has been infected longer than the
              # incubation period):
              if (((week - bID[bug]) > incubation) &
                  (runif(1) < BtoHInfectProb)) {

                # Infect this person this week:
                hI[humanToFeed] <- 1
                hID[humanToFeed] <- week

                # Increment the number of humans that became infected:
                #nHumanInfect[sim] <- nHumanInfect[sim] + 1
              }
            }
            else if ((bI[bug] == 0) & (hI[humanToFeed] == 1)) {
              # Human is infected and vector is uninfected.
              # Determine whether to infect the vector.
              if (runif(1) < HtoBInfectProb(hI[humanToFeed], hID[humanToFeed],
                                            week)) {
                # Infect this bug this week:
                bI[bug] <- 1
                bID[bug] <- week
              }
            }
          }
        }

        # Now determine if this bug should die this week.
        # If so, replace it with a new uninfected bug.
        if (week == bDD[bug]) {
          # Reset this bug
            #bugLife <- c(bugLife, (week-bBD[bug]))
          #cat("bug",bug,"dying on week",week,"at age",(week-bBD[bug]),"\n")
          bBD[bug] <- bFD[bug] <- week
          bI[bug] <- bID[bug] <- 0
          
          # Choose the death week:
          bDD[bug] <- round(rnorm(1, (bugLifetime + week), sdBD))
        }
    }
        }

      # Now determine if any cuyes or dogs should die this week:
      if (cuyCurrentPopulation > 0) {
        for (cuy in 1:cuyCurrentPopulation) {
          if (week == cDD[cuy]) {
            # Reset this cuy
              #cuyLife <- c(cuyLife, (week-cBD[cuy]))
              cI[cuy] <- cID[cuy] <- 0
              cBD[cuy] <- week
              cSuper[cuy] <- (runif(1) < probSuperCuy)
              cDD[cuy] <- round(rnorm(1, (cuyLifetime + week), sdCD))
          }
      }
    }

      if (dogPopulation > 0) {
        for (dog in 1:dogPopulation) {
          if (week == dDD[dog]) {
            # Reset this dog
            #dogLife <- c(dogLife, (week-dBD[dog]))
            dI[dog] <- dID[dog] <- 0
            dBD[dog] <- week
            dDD[dog] <- round(rnorm(1, (dogLifetime + week), sdDD))
          }
        }
      }

      # Store the number of infected bugs and hosts, as well as the total
        # number of bugs and hosts:
        infectBug[sim, week] <- sum(bI)
        bugPopNew[sim, week] <- currentVectorPop
        if (cuyCurrentPopulation > 0)
            infectCuy[sim, week] <- sum(cI)
        else
            infectCuy[sim, week] <- 0
        cuyPop[sim, week] <- cuyCurrentPopulation
        if (dogPopulation > 0)
            infectDog[sim, week] <- sum(dI)
        else
            infectDog[sim, week] <- 0
        dogPop[sim, week] <- dogPopulation
        if (humanPopulation > 0)
            infectHuman[sim, week] <- sum(hI)
        else
            infectHuman[sim, week] <- 0
        humanPop[sim, week] <- humanPopulation
    }
}

  # Now make a plot of the simulation output (total and
  # infectious bugs and hosts as a function of time):
  par(mfrow=c(2,4))

  # Plot total vectors
  matplot(t(bugPopNew), main="Total Vectors", type="l", xlab="Week", ylab="N", col="grey")
  points(apply(bugPopNew, 2, quantile, 0.5), type="l", lwd=2, col="red")
  points(apply(bugPopNew, 2, quantile, 0.975), type="l", lwd=2, col="red")
  points(apply(bugPopNew, 2, quantile, 0.025), type="l", lwd=2, col="red")

  # Plot total guinea pigs:
  matplot(t(cuyPop), main="Total Guinea Pigs", type="l", xlab="Week", ylab="N", col="grey")
  points(apply(cuyPop, 2, quantile, 0.5), type="l", lwd=2, col="red")
  points(apply(cuyPop, 2, quantile, 0.975), type="l", lwd=2, col="red")
  points(apply(cuyPop, 2, quantile, 0.025), type="l", lwd=2, col="red")

  # Plot total dogs
  matplot(t(dogPop), main="Total Dogs", type="l", xlab="Week", ylab="N", col="grey")
  points(apply(dogPop, 2, quantile, 0.5), type="l", lwd=2, col="red")
  points(apply(dogPop, 2, quantile, 0.975), type="l", lwd=2, col="red")
  points(apply(dogPop, 2, quantile, 0.025), type="l", lwd=2, col="red")

  # Plot total humans:
  matplot(t(humanPop), main="Total Humans", type="l", xlab="Week", ylab="N", col="grey")
  points(apply(humanPop, 2, quantile, 0.5), type="l", lwd=2, col="red")
  points(apply(humanPop, 2, quantile, 0.975), type="l", lwd=2, col="red")
  points(apply(humanPop, 2, quantile, 0.025), type="l", lwd=2, col="red")
  
  # Plot vectors infected
  matplot(t(infectBug), main="Infected Vectors", type="l", xlab="Week", ylab="N", col="grey")
  points(apply(infectBug, 2, quantile, 0.5), type="l", lwd=2, col="red")
  points(apply(infectBug, 2, quantile, 0.975), type="l", lwd=2, col="red")
  points(apply(infectBug, 2, quantile, 0.025), type="l", lwd=2, col="red")
  
  # Plot cuyes infected
  matplot(t(infectCuy), main="Infected Guinea Pigs", type="l", xlab="Week", ylab="N", col="grey")
  points(apply(infectCuy, 2, quantile, 0.5), type="l", lwd=2, col="red")
  points(apply(infectCuy, 2, quantile, 0.975), type="l", lwd=2, col="red")
  points(apply(infectCuy, 2, quantile, 0.025), type="l", lwd=2, col="red")

  # Plot dogs infected
  matplot(t(infectDog), main="Infected Dogs", type="l", xlab="Week", ylab="N", col="grey")
  points(apply(infectDog, 2, quantile, 0.5), type="l", lwd=2, col="red")
  points(apply(infectDog, 2, quantile, 0.975), type="l", lwd=2, col="red")
  points(apply(infectDog, 2, quantile, 0.025), type="l", lwd=2, col="red")

  # Plot humans infected
  matplot(t(infectHuman), main="Infected Humans", type="l", xlab="Week", ylab="N", col="grey")
  points(apply(infectHuman, 2, quantile, 0.5), type="l", lwd=2, col="red")
  points(apply(infectHuman, 2, quantile, 0.975), type="l", lwd=2, col="red")
  points(apply(infectHuman, 2, quantile, 0.025), type="l", lwd=2, col="red")
  
  par(mfrow=c(1,1))
  
  # Now print some statistics
  #cat("Number of cuyes infected mean/median:",mean(nCuyInfect),
  #    median(nCuyInfect),"\n")
  #cat("Number of bugs infected mean/median:",mean(nBugInfect),
  #    median(nBugInfect),"\n")
  #cat("Number of dogs infected mean/median:",mean(nDogInfect),
  #    median(nDogInfect),"\n")
  #cat("Number of humans infected mean/median:",mean(nHumanInfect),
  #    median(nHumanInfect),"\n")
  #cat("Mean number of weeks between vector feedings:",
   #   (simLength*vectorPopulation)/mean(nFeeding),"\n")
  #cat("Mean lifespan of",(length(bugLife)-1),"vectors:",mean(bugLife[2:length(bugLife)]),"\n")
  #cat("Mean lifespan of",(length(cuyLife)-1),"cuyes:",mean(cuyLife[2:length(cuyLife)]),"\n")
  #cat("Mean lifespan of",(length(dogLife)-1),"dog(s):",mean(dogLife[2:length(dogLife)]),"\n")

  # Print probabilities of feeding on each host type:
  #totalFeeding <- nCuyFeeding + nDogFeeding + nChickenFeeding + nHumanFeeding
  #cat("Cuy feedings:",nCuyFeeding,(nCuyFeeding/totalFeeding*100),"%\n")
  #cat("Dog feedings:",nDogFeeding,(nDogFeeding/totalFeeding*100),"%\n")
  #cat("Chicken feedings:",nChickenFeeding,(nChickenFeeding/totalFeeding*100),"%\n")
  #cat("Human feedings:",nHumanFeeding,(nHumanFeeding/totalFeeding*100),"%\n")

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
    write.csv(infectBug, file = paste(vectorDataPrefix,".csv",sep=""))
    if (cuyStartPopulation > 0)
      write.csv(infectCuy, file = paste(cuyDataPrefix,".csv",sep=""))
    if (dogPopulation > 0)
      write.csv(infectDog, file = paste(dogDataPrefix,".csv",sep=""))
    if (humanPopulation > 0)
      write.csv(infectHuman, file = paste(humanDataPrefix,".csv",sep=""))

    # Save the input parameters to a spreadsheet
    parameter_names<-c("nsims", "cuyStartPopulation", "dogPopulation", "chickenPopulation", "humanPopulation", "ncuyinfected", "nbuginfected", "feedInterval", "cuyLifetime", "bugLifetime", "dogLifetime", "probSuperCuy","incubation", "BtoCInfectProb", "BtoDInfectProb", "BtoHInfectProb", "cuyPrefFactor", "dogPrefFactor", "chickenPrefFactor", "cuyPersistTime")
    parameter_values<-c(nsims, cuyStartPopulation, dogPopulation, chickenPopulation, humanPopulation, ncuyinfected, nbuginfected, feedInterval, cuyLifetime, bugLifetime, dogLifetime, probSuperCuy, incubation, BtoCInfectProb, BtoDInfectProb, BtoHInfectProb, cuyPrefFactor, dogPrefFactor, chickenPrefFactor, cuyPersistTime)
    parameter_output<-data.frame(parameter_names, parameter_values)
    write.csv(parameter_output, file = "runparameters.csv")

    # Save a spreadsheet with information about how many hosts and
    # vectors became infected during the simulations
    #outData <- data.frame(nBugInfect, nCuyInfect, nDogInfect, nHumanInfect)
    #write.csv(outData, file="infection_numbers.csv")
    
    # Save plots of the infectivity of hosts to vectors.
    saveInfectivityPlots(cuyPersistTime)
    
    # Return to previous working directory
    setwd(inFolder)
  }
}

# Function returns the probability that the given cuy will be infectious
# to a vector that feeds on this cuy on the given week.
CtoBInfectProb <- function(cI, cID, isSuper, weekNum, cuyPersistTime) {

  if (cI == 0) return(0) # this cuy is not infected
  else {

    # Infected. Find out how many weeks have passed since infection:
    weeksSinceInfection <- weekNum - cID
    if (isSuper) {
      # This is a supershedder cuy.
      # Assume 100% chance of infecting bug during the "cuyPersistTime", then
      # 20% for rest of life:
      if (weeksSinceInfection < cuyPersistTime) return(1)
      else return(0.2)
    }
    else {
      # This is a normal, non-supershedder cuy. Assume 100% chance of
      # infecting bug during the "cuyPersistTime", then 0% for rest of life:
      if (weeksSinceInfection < cuyPersistTime) return(1)
      else return(0)
    }
  }
}

# Function returns the probability that the given dog will be infectious
# to a vector that feeds on this dog on the given week.
DtoBInfectProb <- function(dI, dID, weekNum) {
  
  if (dI == 0) return(0) # this dog is not infected

  # Based on Gurtler 1996, the probability that an infected dog
  # transmits infection to an uninfected vector with a single bite is
  # 48.7% overall. This probability doesn't seem to vary much with the
  # age of the dog, so assume it is constant over lifetime.

  return(0.487)
}

# Function returns the probability that the human will be infectious
# to a vector feeding on this week.
HtoBInfectProb <- function(hI, hID, weekNum) {

  if (hI == 0) return (0)

  # Gurtler 1996 gives an overall probability of 2.6% that an infected
  # human will infect a vector after a single bite. This probability
  # decreased from 4.2% in children to 0.48% in adults.
  # My guess is that this decrease is due to the fact that children were
  # infected more recently. So let's model it as 4.2% in the first three
  # years of infection, then 0.48% afterward:
  weeksSinceInfection <- weekNum - hID
  if (weeksSinceInfection < (3*365)) return (0.042)
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
  
  totalWeight <<- humanPopulation + (cuyCurrentPopulation * cuyPrefFactor) +
    (dogPopulation * dogPrefFactor) + (chickenPopulation * chickenPrefFactor)
  probCuy <<- (cuyCurrentPopulation * cuyPrefFactor) / totalWeight
  probDog <<- (dogPopulation * dogPrefFactor) / totalWeight
  probChicken <<- (chickenPopulation * chickenPrefFactor) / totalWeight
  probHuman <<- humanPopulation / totalWeight
  #cat("In setFeedingProbs:probCuy=",probCuy,"probDog=",probDog,"probChicken=",probChicken,"probHuman=",probHuman,"\n")
  #cat("humanPopulation=",humanPopulation,"cuyCurrentPopulation=",cuyCurrentPopulation,"dogPopulation=",dogPopulation,"chickenPopulation=",chickenPopulation,"\n")
}

# Save plots of the infectivity of each host to vectors.
saveInfectivityPlots <- function(cuyPersistTime) {
  
  # Do guinea pigs:
  jpeg(file="cuy_infectivity.jpg")
  cuyInfectivity <- 1:120 * 0
  for (week in 1:120) {
    cuyInfectivity[week] <- CtoBInfectProb(1, 1, 0, week, cuyPersistTime)
  }
  plot(1:120, cuyInfectivity, type="l", xlab="Weeks since infection", ylab="Probability of transmitting infection to vector", main="Guinea Pigs")
  invisible(dev.off()) # close the jpg plot file

  # Do humans:
  jpeg(file="human_infectivity.jpg")
  humanInfectivity <- 1:2555 * 0
  for (week in 1:2555) {
    humanInfectivity[week] <- HtoBInfectProb(1, 1, week)
  }
  plot(1:2555, humanInfectivity, type="l", xlab="Weeks since infection", ylab="Probability of transmitting infection to vector", main="Humans")
  invisible(dev.off()) # close the jpg plot file

  # Do dogs
  jpeg(file="dog_infectivity.jpg")
  dogInfectivity <- 1:2555 * 0
  for (week in 1:2555) {
    dogInfectivity[week] <- DtoBInfectProb(1, 1, week)
  }
  plot(1:2555, dogInfectivity, type="l", xlab="Weeks since infection", ylab="Probability of transmitting infection to vector", main="Dogs")
  invisible(dev.off()) # close the jpg plot file
}

# See Roberts (2003) for a description of how to create the next generation
# matrix.
createNextGenMatrix <- function() {

  # Our matrix will be at most 4x4 (cuyes, dogs, humans, vectors) since
  # chickens cannot become infected. But it may be smaller if one of the
  # host types has a zero population size. Test for that here:
  if (!exists("cuyStartPopulation")) {
    initialize()
    setFeedingProbs()
  }
  if (!exists("probCuy")) setFeedingProbs()

  # Initialize storage for the coefficients
  vectorToHost <- NULL
  hostToVector <- NULL
  
  # Now get the coefficients in the matrix for each host type:
  if (cuyStartPopulation > 0) {
    vectorToHost <- (1/feedInterval)*probCuy*BtoCInfectProb*bugLifetime*exp(-incubation/bugLifetime)
    hostToVector <- (1/feedInterval)*probCuy*(vectorPopulation/cuyStartPopulation)*CtoBInfectProb(1, 1, 0, 1, cuyPersistTime)*cuyPersistTime
  }

  if (humanPopulation > 0) {
    vectorToHost <- c(vectorToHost, (1/feedInterval)*probHuman*BtoHInfectProb*bugLifetime*exp(-incubation/bugLifetime))
    hostToVector <- c(hostToVector, (1/feedInterval)*probHuman*(vectorPopulation/humanPopulation)*integratedHumanInfect)
  }

  if (dogPopulation > 0) {
    vectorToHost <- c(vectorToHost, (1/feedInterval)*probDog*BtoDInfectProb*bugLifetime*exp(-incubation/bugLifetime))
    hostToVector <- c(hostToVector, (1/feedInterval)*probDog*(vectorPopulation/dogPopulation)*DtoBInfectProb(1, 1, 1)*dogLifetime)
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
  for (week in 1:(humanLifetime))
    integral <- integral + HtoBInfectProb(1, 1, week)

  return(integral)
}

# Commands for running on the cluster:
#cat("Doing simulation with cuyes and vectors. Initial 1 infected cuy.\n")
#initialize(cuyStartPopulationIn = 10, vectorPopulationIn = 1000,
#           dogPopulationIn = 0, humanPopulationIn = 0,
#           chickenPopulationIn = 0,
#           ncuyinfectedIn = 1, nbuginfectedIn = 0)
#doSim(nsims = 5, simLength = 365, makeOutput=TRUE, outdir="~/output/multiple_host/cuyes_vectors_1cuy_infected/", vectorDataPrefix="vector_prev", cuyDataPrefix="cuy_prev", dogDataPrefix="dog_prev", humanDataPrefix="human_prev")

#cat("Doing simulation with all hosts. Initial 1 infected cuy.\n")
#initialize(cuyStartPopulationIn = 10, vectorPopulationIn = 1000,
#           dogPopulationIn = 2, humanPopulationIn = 4,
#           chickenPopulationIn = 3,
#           ncuyinfectedIn = 1, nbuginfectedIn = 0)
#doSim(nsims = 5, simLength = 365, makeOutput=TRUE, outdir="~/output/multiple_host/all_hosts_1cuy_infected_v2/", vectorDataPrefix="vector_prev", cuyDataPrefix="cuy_prev", dogDataPrefix="dog_prev", humanDataPrefix="human_prev")

cat("Doing simulation...\n")
initialize(cuyStartPopulationIn = 10,
           dogPopulationIn = 2,
           humanPopulationIn = 4,
           chickenPopulationIn = 3,
           ncuyinfectedIn = 1,
           nbuginfectedIn = 0,
           cuyBottleStartTimeIn = 40,
           cuyBottleEndTimeIn = 90,
           cuyBottlePopSizeIn = 0,
           cuyToReintroduceIn = 10,
           infectedCuyToReintroduce = 1,
           vectorPopBaseline = 100, 
           vectorPopCuyMultiplier = 10,
           bugBottleStartTimeIn = -1,
           bugBottleEndTimeIn = 150,
           bugBottlePopSizeIn = 0,
           newBugsPerDayAfterBottle = 5,
           infectedBugReintroduceAfterBottle = 1
           )
doSim(nsims = 100, simLength = 300, makeOutput=TRUE, outdir="../output/all_hosts_1vector_1cuy/", vectorDataPrefix="vector_prev", cuyDataPrefix="cuy_prev", dogDataPrefix="dog_prev", humanDataPrefix="human_prev")

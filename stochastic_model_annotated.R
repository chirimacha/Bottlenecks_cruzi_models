inFolder <- getwd() # store the current working directory

#==============================================
# Main function. Input parameters are assigned default values that
# can be changed by the user.
# Possible values of pairtype: "bottle" -- guinea pigs experience a
# bottleneck; "multiplebottle" -- guinea pigs and vectors experience a
# population bottleneck; "nobottle" -- no bottlenecks of either hosts or
# vectors..
doSim <- function(pairtype = "bottle",
                  nsims = 1, # number of simulations to run
                  cuysamplesize = 10, # guinea pig population size
                  vectorsamplesize = 1000, # triatomine population size
                  incubation = 30, # vector incubation period in days
                  numCuyDuringBottleneck = 2, # size of guinea pig population during bottlenecks
                  bottleneckStart = 120, # day when guinea pig population bottleneck begins
                  numBugDuringBuggleneck = 200, # size of vector population size during vector bottlenecks
                  bottleneckLength = 100, # length of bottleneck in days
                  probSupershedder = 0.5, # probability that a given guinea pig is a "supershedder" with prolonged parasitemia
                  simLength = 3650, # total length of simulation in days
                  makeOutput=FALSE, # boolean determines whether detailed output files should be created.
                  outdir="none", # directory for storage of output files
                  plotPrefix="prevalence_vs_time", # prefix for filename of prevalence plots
                  vectorDataPrefix="vector_data", # prefix for filename for storage of vector infection data
                  hostDataPrefix="host_data" # prefix for filename for storage of host infection data
                  ) {

  # Set desired input parameters
  parameters(cuysamplesize, vectorsamplesize, incubation, simLength)

  # Allocate storage for output data structures
  out<-matrix(NA, nrow=EndTime,ncol=2)
  outVector<-matrix(NA, nrow=nsims,ncol=EndTime)
  outHost<-matrix(NA, nrow=nsims,ncol=EndTime)
 
  for (simulationNumber in 1:nsims) {

    # Reset the global variables and matrices.
    initialize(numCuyDuringBottleneck, bottleneckStart, bottleneckLength, cuysamplesize, probSupershedder, numBugDuringBuggleneck)
        
    Feed<-setFeed()     #Sets feed type and makes feeding matrix.
    
    # Run the model and store outputs of each simulation
    out<-runall(Feed, pairtype)
    outVector[simulationNumber,]<-out[,1]
    outHost[simulationNumber,]<-out[,2]
    print(simulationNumber)
  }

  if (makeOutput) {
    # Plot and store simulation output.
    
    #==============================================
    #This set of plotting commands creates two figures and saves them as a 
    #jpeg.  The figures show the infection prevalence in vectors and
    # hosts on each day of the simulation.
    #==============================================

    # Use desired output directory, or create a new one.
    if (outdir == "none") {
      outfolder<-"./output/SIM"
      timestamp<- gsub(":","-",gsub(" ","_",Sys.time()))
      outdir<-paste(outfolder,timestamp, sep="_")
    }

    # change to directory where files should be created, and save plots
    # and csv files with prevalence data for hosts and vectors.
    dir.create(outdir, recursive=TRUE)
    setwd(outdir)
    plotName <- paste(plotPrefix,".jpg", sep="")
    color<-rainbow(nsims)
    jpeg(file=plotName)
    par(mfrow=c(2,1))
    plot(outVector[1,], ylim=c(0, Nbug), main="Number of Vectors Infected Each 
      Day", cex=.2, pch="-", ylab="Number Infected", xlab="Day")

    if (pairtype != "nobottle") {
      shifts<-seq(0,EndTime,BNrep)
      abline(v=c(startBottle+shifts,startBottle+lengthBottle+shifts),col=4)
    }

    for(phi in 2:nsims) {
      points(outVector[phi,], cex=.2, pch="-", col=color[phi])
    }

    plot(outHost[1,], ylim=c(0, Ncuy), main="Number of Hosts Infected Each Day", cex=.2, pch="-", ylab="Number Infected", xlab="Day")
    for(phi in 2:nsims) {
      points(outHost[phi,], cex=.2, pch="-", col=color[phi])
    }
    par(mfrow=c(1,1))
    dev.off() # close the jpeg plot file

    # Save the prevalence data in spreadsheets for both vectors and hosts.
    write.csv(outVector, file = paste(vectorDataPrefix,".csv",sep=""))
    write.csv(outHost, file = paste(hostDataPrefix,".csv",sep=""))
    setwd(inFolder)

    # Save the input parameters to a spreadsheet
    parameter_names<-c("Ncuy", "NcuyI", "Nbug", "NbugI", "l", "f", "cLife", "bLife", "L", "startBottle", "lengthBottle", "bottleNeck", "BNrep", "startBottlebug", "lengthBottlebug", "bottleNeckbug", "BNrepbug", "probSupershed","incPeriod") 
    parameter_values<-c(Ncuy, NcuyI, Nbug, NbugI, l, f, cLife, bLife, L, startBottle, lengthBottle, bottleNeck, BNrep, startBottlebug, lengthBottlebug, bottleNeckbug, BNrepbug, probSupershed, incPeriod)
    parameter_output<-data.frame(parameter_names, parameter_values)
    setwd(outdir)
    write.csv(parameter_output, file = "runparameters.csv")
    setwd(inFolder)
  }

  return(outVector)
} # end doSim()

# set all the initial parameters as global variables
parameters<-function(cuysamplesize, vectorsamplesize, incubation, simLength) {
  #==============================================
  #Set up guinea pig population (guinea pig is "cuy" in Spanish)
  #==============================================
  Ncuy<<-cuysamplesize    #number of hosts (outside of bottleneck times)
  NcuyI<<-1               #initially infected hosts at start of simulation
                                                          
  #============================================= =
  #Set up triatomine vector ("bug") population
  #==============================================
  Nbug<<-vectorsamplesize  #number of bugs (apart from buggleneck times)
  NbugI<<-0       #initially infected bugs at start of simulation
      
  #==============================================
  # Set up vector feeding parameters
  #==============================================
  l <<- 14               #Number of days between bloodmeals for each vector
  EndTime <<- simLength  # Total number of days in the simulation
  f <<- EndTime/l          #Number of times bugs feed during simulation
      
  #==============================================
  #Parameters for the frequency with which vectors and hosts die
  #==============================================
  cLife<<-120     #frequency with which hosts (guinea pigs) are reset (days)
  bLife<<-180     #frequency with which vectors are reset (days)

  # Incubation period - vectors must be infected for this many days
  # before they can pass infection to hosts.
  incPeriod <<- incubation
      
  #==============================================
  #Set up vector of vector->host infection probabilities with a
  # logistic function, constructed such that maximum infectiousness
  # is reached around 45 days  
  #==============================================
  da<-seq(1, EndTime)
  L<<-0.00068  #L is the maximum infectiousness - taken from literature (Catala 1992 -> .00068)
  k<- 0.15     #k is steepness of curve (increase to reach L faster)
  c<- 20       #c is number of times initial output must grow to reach L
  InfectProbBugs<<- L/(1+c*exp(-k*da))  #Logistic function --> Vector of bug infectiousness to hosts

  #==============================================
  #Create the matrices that will hold information regarding 
  #which vector infected which guinea pig (and vice versa) and on what day
  #==============================================
  bugToCuy<<-matrix(0, nrow=0, ncol=3)  #col 1= vector number, col 2= guinea pig number, col 3= day
  cuyToBug<<-matrix(0, nrow=0, ncol=3)  #col 1= guinea pig number, col 2= vector number, col 3= day
}   #end parameters
 
 
#==============================================
#this function is called at the start of every run to set the 
#run-specific parameters and global variables.
#==============================================
initialize<-function(numCuyDuringBottleneck, bottleneckStart,
                     bottleneckLength, cuysamplesize, probSupershedder,
                     numBugDuringBuggleneck) {

  #==============================================
  # Sets initial parameters for the host infection status (cI), host infection 
  # day (cID), which hosts are active/alive (activeHost), and the
  # total number of active hosts that are infected.
  #==============================================
  cI<<-rep(0, Ncuy)     #create a 0/1 vector to record host infection status
  cID<<-rep(0, Ncuy)    #create a vector to record day when host became infected
  activeHost<<- seq(1, Ncuy)  #create a vector of available/living hosts
  cuyTime<<-rep(0, EndTime)   #A matrix to hold the daily sum of hosts infected
  phi<-sample(1:Ncuy, NcuyI, replace=FALSE)  # host(s) to be infected at start
  cID[phi]<<-1   # Assign the infection day to day 1 for infected hosts
  cI[phi]<<-1 # assign the infection status to 1 for infected hosts
                                  
  #==============================================
  #Sets the initial parameters for the vector infection status (bI), the vector
  #infection day (bID), which vectors are active/alive and the total number of 
  #active vectors that are infected on a given day
  #==============================================
  bI<<-rep(0, Nbug)
  bID<<-rep(0, Nbug)
  activeVector<<-seq(1, Nbug)
  bugTime<<-rep(0, EndTime) #A matrix to hold the daily sum of vectors infected
  phi<-sample(1:Nbug, NbugI , replace=FALSE)  # vectors to be infected at start
  bID[phi]<<-1   #Assign the infection day to day 1 for all infected vectors
  bI[phi]<<-1

  # Set up the host-to-vector (C->B) infection probabilities.
  # The following section probabilistically determines whether each host
  # is a supershedder guinea pig, or not. It then assigns the correct
  # transmission probability (from guinea pig -> vector) to the host
  # for each of 52 weeks (assumption is that guinea pigs will never live
  # longer than one year).
  cuyType <- runif(cuysamplesize, min=0, max=1)
  CBInfectMatrix<<-matrix(NA, nrow=cuysamplesize, ncol=52)
  for (i in 1:52) {
    for(j in 1:cuysamplesize) {
      # Set probability that this host is infectious for each of 52 weeks:
      if (cuyType[j] < probSupershedder)
        # this is a supershedder host
        CBInfectMatrix[j,i] <<- probInfectionSupershedder(i)
      else
        # this is a "normal" (non-supershedder) host
        CBInfectMatrix[j,i] <<- probInfectionNormalCuy(i)
    }
  }

  # Store the probSupershedder so that it can be used when new hosts
  # enter the simulation:
  probSupershed <<- probSupershedder

  #==============================================
  # Store parameters for host bottlenecks. 
  #==============================================
  startBottle<<-bottleneckStart        #Day to start the bottleneck
  lengthBottle<<-bottleneckLength      #Number of days in in the bottleneck
  endBottle<<-startBottle+lengthBottle #Day to end the bottleneck 
  bottleNeck<<-numCuyDuringBottleneck  #Number of hosts that stay active during a bottleneck
  BNrep<<-365      #Frequency with which bottleneck should be repeated
  
  # Store days when bottlenecks start and end
  startDays <<- seq(startBottle, EndTime, BNrep)
  endDays <<- seq(endBottle, EndTime, BNrep)
  # If there are more startDays than endDays, bottleneck is still
  # active when the simulation ends. In this case, add an
  # extra element onto endDays so that the getBottleneckStatus() method
  # will work properly:
  if (length(startDays) > length(endDays))
    endDays <<- c(endDays, EndTime + 1)
  
  #==============================================
  # Store parameters for vector bottlenecks - they
  # start and end at the same times that host bottlenecks occur, so
  # that the same functions can be used to determine whether
  # a vector bottleneck ("buggleneck") is currently active.
  #==============================================      
  startBottlebug <<- bottleneckStart
  lengthBottlebug <<-bottleneckLength
  endBottlebug<<- startBottlebug+lengthBottlebug
  bottleNeckbug<<-numBugDuringBuggleneck
  BNrepbug<<-365
      
  #==============================================
  # create matrices of assigned death days for each host and vector
  #==============================================
  bDeath<<-bugDeath()
  cDeath<<-cuyDeath()

} #end initialize()

# This function returns the probability that a supershedder guinea pig
# will infect a vector, on a given week. Supershedders infect 100% of
# vectors during the first 8 weeks, then 20% of vectors for the remainder
# of the host's lifetime.
probInfectionSupershedder <- function(weekNumber) {
  if (weekNumber <= 8) return(1)
  else return(0.2)
}

# Function to return the probability that a normal non-supershedder guinea
# pig will infect a vector on a given week -- these hosts infect 100%
# of biting vectors during the first 8 weeks of infection, then control their
# parasitemia and do not infect vectors thereafter.
probInfectionNormalCuy <- function(weekNumber) {
  if (weekNumber <= 8) return(1)
  else return(0)
}

#==============================================
#Creates a matrix that pre-determines the death day of each vector.  Death 
#frequency is determined by the parameter bLife. This function returns
# a matrix bDeath with a 1 assigned to the vector's death day.  
#==============================================
bugDeath<-function() {
  #Determines the number of times each vector could die considering 
  #the length of a vector's life (bLife) and number of days in the simulation
  #It then produces a vector of last possible death day for each death cycle 
  deathFreq<-1:ceiling((EndTime)/bLife)*bLife
  bDeath<-matrix(0, nrow=Nbug, ncol=max(deathFreq)) #Matrix with columns = days and vectors = rows
  for (k in 1:Nbug) {
    eL<-deathFreq-sample(0:(bLife-1),1) #creates a vector of days vector k will die
    bDeath[k, eL]<-1                  #assigns vector of days to bDeath matrix
  } #end outer for loop
  return(bDeath)
} #end bugDeath

#==============================================
#function creates a matrix that pre-determines the death day 
#of each host.  The start day is random but hosts are reset every cLife 
#number of days, therefore cLife controls the frequency of death.  
#==============================================
cuyDeath<-function() {
  deathFreq<-1:ceiling((EndTime)/cLife)*cLife  #vector of maximum death days.
  cDeath<-matrix(0, nrow=Ncuy, ncol=max(deathFreq))
  for (k in 1:Ncuy) {
    eL<-deathFreq-sample(0:(cLife-1),1) #creates a vector of when host will die, staggering death days
    cDeath[k, eL]<-1   #Assigns vector of death days to cDeath matrix
  }
  return(cDeath)
}
  
 #==============================================
 #resetBugs is a function that looks at column day of the bDeath matrix to
# see which vectors are due to die and be reset on that day. The function
# sets bI and bID to 0 for each reset vector.
 #==============================================
resetBugs<-function(day) {
  if (sum(bDeath[,day])==0) 
    {
      return          # no vectors died on this day
    } #end if
  else
    {
        # reset vectors that have died
      reset<-which(bDeath[,day]==1) #Rows in bDeath of vectors to die that day
      bID[reset]<<-0
      bI[reset]<<-0
    } #end else
} #end resetBugs()
        
#==============================================
# function that resets cI and cID for hosts that died on the given day.
# creates new replacement hosts (either supershedders or not)
#==============================================
resetCuy<-function(day) {
  if (sum(cDeath[,day])==0) 
    return
  else {
    # Reset the infection day and infection status of the new hosts
    # (which are all assumed to be uninfected).
    reset<-which(cDeath[,day]==1) #Rows in cDeath of hosts to die that day
    cID[reset] <<- 0
    cI[reset] <<- 0

    # Update the C->B infection probability matrix as hosts regenerate.
    for (cuy in 1:length(reset)) {
      randNum <- runif(1, min=0, max=1)
      if (randNum < probSupershed) {
        # This new host is a supershedder
        for (j in 1:52) {
          # Set infection probability for each week
          CBInfectMatrix[reset[cuy], j] <<- probInfectionSupershedder(j)
        }
      }
      else {
        # This new host is a normal (non-supershedder) host
        for (j in 1:52) {
          CBInfectMatrix[reset[cuy], j] <<- probInfectionNormalCuy(j)
        }
      }
    }
  }
}

 #==============================================
 #Function to handle the "buggleNeck" = vector bottleneck.
# resize the vector population at the edges of bottlenecks, when
# the population may either increase or decrease. Returns
 # a list of all vectors active on the given day.
 #==============================================
buggleNeck<-function(day, pairtype, activeVector) {
  
  if (pairtype=="multipleBottle") {

    # Get the status of today's day compared to when bugglenecks occur:
    # Return values: 1 = between bottlenecks; 2 = start day of bottleneck;
    # 3 = within a bottleneck; 4 = end day of bottleneck
    bottleneckStatus <- getBottleneckStatus(day)
    
    if (bottleneckStatus == 4) {
      # Return to normal vector population size
      
      cleanVect<-resumePop(bI,bID,activeVector)   # On the last day of the bottleneck, bI and bID are reset
      activeVector <- seq(1:Nbug)
      bI<<-cleanVect$I
      bID <<- cleanVect$ID
          
      return(activeVector)
    }
    else if (bottleneckStatus == 2) {
      # This is the start day of the bottleneck. Reduce vector population size.
      # And select the subset of vectors that remains "active".
      activeVector<-sample(activeVector, bottleNeckbug, replace=FALSE) 
      return(activeVector)
    }
    else if ((bottleneckStatus == 3) | (bottleneckStatus == 1))
      # within or between bottlenecks. Do nothing.
      return(activeVector)
  } #end if pairtype "multiplebottle"
  else {
    # There is no buggleneck. No need to modify activeVector:
    return(activeVector)
  }
} #end buggleneck()

 #==============================================
 #Function whichFeed determines which vectors are feeding on the given day.
 #whichFeed returns FT, or returns NA if no vectors are feeding
 #==============================================
whichFeed<-function(Feed, day, pairtype, activeVector) {
  total<-sum(Feed[,day])
  if (total==0) return(NA)       
  activeVector <<- buggleNeck(day, pairtype, activeVector)
  FT<-intersect(activeVector,which(Feed[,day]==1)) #determines the rows in column day that have a 1
  return(FT)  #return the vector of vectors feeding that day
} #end whichFeed
        
 #==============================================
 #Function pairfunction matches vectorss feeding on a given day
 #(in input object FT) with potential hosts. pairfunction returns a matrix
 #with column 1= vectors and column 2= the host that vector will feed on.
 # The input parameter pairtype should be given as "nobottle",
 # "bottle" (host bottleneck only),
 # or "multipleBottle" (both host and vector bottlenecks)
 #==============================================
pairfunction<-function(FT, day, pairtype) {

  BugHost<-matrix(NA, nrow=length(FT), ncol=2) #Create matrix of vectors (col 1) and hosts (col 2)
  BugHost[,1]<-FT   #Assign vectors feeding on day to column 1
  
  #==============================================
  #No bottleneck or special pairing option
  #all pairing is random and even
  #==============================================
  if (pairtype=="nobottle") {
    BugHost[,2]<-sample(activeHost, length(FT), replace=TRUE)  #Assigns host to the second column
    return(BugHost)
  }
  
  #==============================================
  #Bottleneck option (host or vector)
  #==============================================
  else if (pairtype=="bottle" | pairtype=="multipleBottle") {

    # Get the status of today's day compared to when bottlenecks occur:
    # Return values: 1 = between bottlenecks; 2 = start day of bottleneck;
    # 3 = within a bottleneck; 4 = end day of bottleneck
    bottleneckStatus <- getBottleneckStatus(day)
    
    if (bottleneckStatus == 1) { # between bottlenecks
      BugHost[,2]<-sample(1:Ncuy,length(FT), replace=TRUE)  #Assigns hosts to the second column randomly
      return(BugHost)
    }
    else if (bottleneckStatus ==2) { # start day of bottleneck
      hostBottle<<-sample(1:Ncuy, bottleNeck, replace=FALSE)  # determines which hosts will remain active during the bottleneck. All others have died.  
      BugHost[,2]<-sample(hostBottle,length(FT), replace=TRUE)
      return(BugHost)
    }
    else if (bottleneckStatus == 3) { # within bottleneck
      BugHost[,2]<-sample(hostBottle,length(FT), replace=TRUE)
      return(BugHost)
    }
    else if (bottleneckStatus == 4) { # end day of bottleneck
      BugHost[,2]<-sample(1:Ncuy,length(FT), replace=TRUE)  #Assigns host to the second column randomly

      # Reset cI and cID to 0 for all hosts that have just regenerated,
      # and also determine the host's infection probability by deciding
      # if each new host is a supershedder or not.
      for (i in 1:Ncuy) {
        if (!(i %in% hostBottle)) {
          # This is a newly regenerated host that needs to be reset.
          cI[i] <<- 0
          cID[i] <<- 0
          randNum <- runif(1, min=0, max=1)
          if (randNum < probSupershed) {
            # This new host is a supershedder
            for (j in 1:52) {
              # Set infection probability for each week
              CBInfectMatrix[i, j] <<- probInfectionSupershedder(j)
            }
          }
          else {
            # This new host is a normal (non-supershedder) guinea pig
            for (j in 1:52) {
              CBInfectMatrix[i, j] <<- probInfectionNormalCuy(j)
            }
          }
        }
      }
      
      return(BugHost)
    } #end if
    
  } #end if pairtype "bottle"
  else {
    cat("invalid pairtype", pairtype,"\n")
    stop()
  }
} #pairfuntion()	 

# Function to determine where this day falls with respect to host bottlenecks.
# Return values: 1 = between bottlenecks; 2 = start day of bottleneck;
# 3 = within a bottleneck; 4 = end day of bottleneck
getBottleneckStatus <- function(day) {
  
  if (day %in% startDays)
    return(2) # start day of bottleneck
  else if (day %in% endDays)
    return (4) # end day of bottleneck
  else {
    # check each bottleneck to see whether current day is within bottleneck:
    for (i in 1:length(startDays)) {
      if ((day > startDays[i]) && (day < endDays[i]))
        return(3) # within a bottleneck
    }
  }

  # between bottlenecks - return 1
  return(1)
}
        
#==============================================
#resumePop is a function that is called at the end of a vector bottleneck.
#It puts the vector infection day (bID) and infection status (bI) back to 0
# for all vectors that were removed during the "buggleneck"
#==============================================
resumePop<-function(InfStatus,InfDay,activeOnes) {
  cleanInfStatus<-rep(0, length(InfStatus))  #place holder vector for bI
  cleanInfDay<-rep(NA, length(InfDay))       #place holder vector for bID
  for (w in 1:length(activeOnes))            #Loops through the active vectors
    {
      cleanInfStatus[activeOnes[w]]<- InfStatus[activeOnes[w]] # retain bI status for active vectors
      if (is.na(InfDay[activeOnes[w]])==FALSE)    # retain bID for active vectors
        {
          cleanInfDay[activeOnes[w]]<- InfDay[activeOnes[w]]
        } #end if
    } #end for loop

  return(list(I=cleanInfStatus,ID=cleanInfDay))
} #end resumePop		
        
 #==============================================
 #transmissionSort has inputs BugHost and day.  It 
 #determines if the vector/host pair should be passed to the InfectBtoC, 
 #InfectCtoB or neither
 #==============================================
transmissionSort<-function(BugHost, day) {
  DAY<<-day
  for(alpha in 1:length(BugHost[,1])) {
    
    vector<- BugHost[alpha, 1]
    host<- BugHost[alpha, 2]
              
    if ((bI[vector] == 1) & (cI[host] == 0) & (DAY>bID[vector])) {
      #Vector that has been infected for at least one day and uninfected host
      # Pass the vector, host, and day to InfectBtoC()
      InfectBtoC(vector, host, day)
    }
    
    if((bI[vector]==0) & (cI[host]==1) & (DAY>cID[host])) {
     #Host that has been infected for at least one day and an uninfected vector
      # Pass the vector, host, and day to InfectCtoB()
      InfectCtoB(vector, host, day, CBInfectMatrix)
    }
  }
}

#==============================================
#InfectBtoC is a function that determines if infection is passed from an 
#infected vector to a non-infected host. Inputs are the day, vector number and 
#host number.  The function determines if the infection was passed between 
#vector and host.  If so, then cI, cID, and bugToCuy are altered.
#==============================================
InfectBtoC<-function(vector, host, day) {
      
  daysSinceInfection<-(day-bID[vector])  #Days since vector was initially infected

  # If still in the incubation period, vector cannot infect host:
  if (daysSinceInfection < incPeriod) return
  
  InfectProb<-InfectProbBugs[daysSinceInfection]  #InfectProbBugs is a vector of solutions to the vector infectiousness equations

  randomNum<-runif(1, min=0, max=1)
  if (randomNum <= InfectProb) {
    #If the random number is less than the vector's infectiousness, then the host becomes infected
    cID[host]<<- day   # store day (cID) when host became infected
    cI[host]<<- 1      # specify that this host is now infected
    bugToCuy<<-rbind(bugToCuy, c(vector, host, day)) #Columns are: vector number, host number, day when vector infected host
  } #end if
} #end InfectBtoC

#==============================================
#InfectCtoB is a function that determines if infection is passed from an 
#infected host to a non-infected vector. Inputs are the day, vector number and 
#host number.  The function determines if the infection was passed between 
#vector and host.  If so, then bI, bID cuyToBug are altered. 
#==============================================
InfectCtoB<-function(vector, host, day, CBInfectMatrix) {
  
  daysSinceInfection<-(day-cID[host])  #Days since host was initially infected
  randomNum<-runif(1, min=0, max=1)
  probaInfect<-CBInfectMatrix[host, ceiling(daysSinceInfection/7)]
  
  if(randomNum <= probaInfect) {
    bID[vector]<<-day  #If the vector becomes infected, that vector's infection day is altered
    bI[vector]<<-1     #The vector's infection status is chaged from 0 to 1
    cuyToBug<<-rbind(cuyToBug, c(host, vector, day)) #Columns are: host number, vector number, day on which host infected vector
  }
}

#==============================================
# function creates a matrix "Feed" that stores the pre-determined days
# on which each vector will feed.
#==============================================
setFeed<-function() {
  Feed <- matrix(0,nrow=Nbug, ncol=EndTime)  #allocate storage
  everyl <- 1:f * l

  for (k in 1:Nbug) {
    eL <- everyl-sample(0:(l-1),length(everyl), replace=TRUE)
    for (j in 1:f) {
      Feed[k,eL[j]] <- 1
    }
  }

  return(Feed)
}
        
#==============================================
# Main loop for the simulation. 
#==============================================
runall<-function(Feed, pairtype) {

  for (day in 1:EndTime) {
      # step through time and call feeding functions on each day

    day<<-day
    FT<<- whichFeed(Feed, day, pairtype, activeVector)  #Returns vectors feeding on day i as a vector of vector numbers
    
    if(is.na(FT[1]) == TRUE)
      return     # no bugs are feeding on day i
    else {
      BugHost<<-pairfunction(FT, day, pairtype)    #Pairs vectors and hosts in a matrix, column 1= bugs, column 2= host

      transmissionSort(BugHost, day)  #Determines if host-vector pair 
                                        #should be sent to InfectBtoC, 
                                        #InfectCtoB or neither. 
    }

    # Track infection prevalence in guinea pigs and vectors
    stat <- getBottleneckStatus(day)
    if ((pairtype=="nobottle") || (stat != 3)) {
      # Either there are no bottlenecks, or current day does not fall
        # within a bottleneck
      cuyTime[day]<-sum(cI[activeHost])
      bugTime[day]<-sum(bI[activeVector])
    }
    else {
      # The current day falls within a bottleneck:
      cuyTime[day]<-sum(cI[hostBottle])
      bugTime[day]<-sum(bI[activeVector])
    }
    
    resetBugs(day)
    resetCuy(day)
  } #end for(day in 1:EndTime)
  
  return(cbind(bugTime,cuyTime))
} #end runall()


# ============================================================
# What follows are two example commands that demonstrate how to run the code.
# Either of these lines may be uncommented to run the examples.
# ============================================================

# Example 1: Simulation with no bottleneck and no supershedder guinea pigs.
doSim(pairtype = "nobottle", nsims = 100, cuysamplesize = 10, vectorsamplesize = 1000, incubation = 45, numCuyDuringBottleneck = 2, bottleneckStart = 0, bottleneckLength = 0, probSupershedder = 0, simLength = 3650, makeOutput=TRUE, outdir="./output/test1", plotPrefix="prevalence_vs_time", vectorDataPrefix="vector_prev_vs_time", hostDataPrefix="host_prev_vs_time")

# ============================================================

# Example 2: Simulation with a guinea pig bottleneck that lasts 90 days
# (guinea pig population drops from 10 -> 2 hosts during bottleneck)
#doSim(pairtype = "bottle", nsims = 10, cuysamplesize = 10, vectorsamplesize = 1000, incubation = 30, numCuyDuringBottleneck = 2, bottleneckStart = 120, bottleneckLength = 90, probSupershedder = 0, simLength = 1000, makeOutput=TRUE, outdir="./output/test2", plotPrefix="prevalence_vs_time", vectorDataPrefix="vector_prev_vs_time", hostDataPrefix="host_prev_vs_time")


#experimentfunctions.r - script containing functions to run experiments

# establish burn-in period - i.e. number of generations for within community diversity to stabilise.
runBurnIn <- function(populationSize = 200, stepSize = 50, iterations = 100, endpoint = 2000){
  #initialisations
  generations <- seq(stepSize, endpoint, by = stepSize) #generation vector
  populations <- as.factor(rep('a', each=populationSize)) #for pi calculations
  #load output file to restart where last experiment left off if interrupted
  outputfilepath <- paste0('./burninresults/pop', populationSize, 'step', stepSize, '.rds')
  if (file.exists(outputfilepath)){
    out <- readRDS(outputfilepath)
  } else {
    out <- data.frame("Iteration"=1, "Generation"=stepSize, "Pi"=NA, "seed"=NA)
  }
  
  for (iteration in as.numeric(tail(out$Iteration),1):iterations){ #loop through iterations
    #RunSimulation
    if (tail(out$Generation, 1) != endpoint){
      seed <- tail(out$seed, 1)
      genstart = tail(out$Generation, 1) + stepSize
    } else {
      seed <- sample(1:1e12,1) #set seed for repeatability
      system(paste0("slim -s ",seed," -d steps=",stepSize," -d endpoint=",endpoint," -d popsize=",populationSize," ./SLiMscripts/BurnIn.slim"))
      genstart = stepSize
    }
    for (g in seq(genstart, endpoint, by = stepSize)){ #loop through generations
      #read in vcf file
      vcfFilepath = paste0("vcffiles/burnin",seed,"-",g,".vcf")
      vcf <- read.vcfR(vcfFilepath, verbose = FALSE)
      #calculate heterozygosities for pi calculation
      gendif <- genetic_diff(vcf, populations, methoed = 'nei')
      #add zeros so means can be compared
      zeros <- as.data.frame(matrix(0,nrow = 8000-nrow(gendif), ncol = ncol(gendif)))
      names(zeros) <- names(gendif)
      gendif <- suppressWarnings(rbind(gendif, zeros))
      
      pi = mean(gendif$Ht, na.rm = TRUE)
      
      out <- rbind(out, c(iteration, g, pi, seed))
      #save output in case of interruptions
      saveRDS(out, outputfilepath)
    }
  }
}
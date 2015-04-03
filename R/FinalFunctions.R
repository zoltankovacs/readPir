PirInfo <- function(dataFile) {
  fCon <- file(dataFile, open = "rb")
  nByts <- readBin(fCon, integer(), n = 5, size = 4)[5]
  nColSpect <- readBin(fCon, integer(), n = 1, size = 4)
  nColClass <- readBin(fCon, integer(), n = 1, size = 4)
  nColDepen <- readBin(fCon, integer(), n = 1, size = 4)
  nColAll <- readBin(fCon, integer(), n = 30, size = 4)[30]
  nRow <- readBin(fCon, integer(), n = 1, size = 4)
  close(fCon)
  out <- list(NoColAll = nColAll, NoRow = nRow, NoColSpect = nColSpect, 
              NoColClass = nColClass , NoColDepen = nColDepen, NoByts = nByts) 
} #Eof functio

PirIndexes <- function(dataFile, nByts) { #could be quicker with for like in header function
  fCon <- file(dataFile, open = "rb")
  Raw <- readBin(fCon, raw(), n = nByts * 128)
  close(fCon)
  nrow<-floor(length(Raw)/4)
  AAA <- matrix(Raw, ncol=nrow)
  AAA<-t(AAA)
  
  posi <- NA
  while(TRUE) {
  for (i in 1:nrow){
    match <- which(AAA[i,]=="ff")
    if (length(match)==4){
      posi <- c(posi,i)
    }
  }
  if(length(posi)>6){
  	break
  }
  }
  out <- posi <- 4*(posi[-c(1,2)])+136
}#Eof function

getPirTable <- function(dataFile, fromPosi, toPosi, nRow, nColAll) {
  fCon <- file(dataFile, open = "rb")
  a <- readBin(fCon, raw(), fromPosi + 4) # jump to the start point
  
  DatatoRead <- toPosi - fromPosi
  a <- readBin(fCon, raw(), DatatoRead) # read all the data as raw
  
  IndNoData <- seq(116, length(a), 128)
  IndNoData <- sort(c(IndNoData, IndNoData-1, IndNoData-2, IndNoData-3))
  
  a <- a[-IndNoData] # eliminate the byte indexes (one block contains 128 units but always starts with FF FF FF FF, we do not need)
  close(fCon)
  a <- readBin(a, double(), size=4, DatatoRead) # translet it into numbers
  a <- a[1:(nRow*nColAll)]
  a[which(a==9.99999968028569246551e+37)] <- NA
  
  out <- NIR <- matrix(a, nRow, nColAll)
}#Eof function

getPirHead <- function(dataFile, fromPosi, toPosi, lengthHead) {
  fCon <- file(dataFile, open = "rb")
  a<-readBin(fCon, raw(), fromPosi)
    
  b<-readBin(fCon, raw(), toPosi-fromPosi)#135
  close(fCon)
  if (length(b) > 128){
    IndNoData <- seq(120, length(b), 128)
    IndNoData <- sort(c(IndNoData, IndNoData-1, IndNoData-2, IndNoData-3))
    b <- b[-IndNoData]
  }
  colNames <- readBin(b, character(), length(b))
  out <- colNames[1:(lengthHead)]  
}#Eof function

creatPirMatrix <- function(PirTable, colNames, rowNames, nColSpect, nColClass, noColDepen) {
  if(nColSpect+nColClass+noColDepen==1){
    PirTable <- matrix(PirTable, ncol=1)
  }
  colnames(PirTable) <- colNames
  rownames(PirTable) <- rowNames
  if (nColSpect > 0){
    spectCol <- 1:nColSpect
    IndepenVar = matrix(PirTable[,spectCol], ncol = nColSpect)
    colnames(IndepenVar) <- paste0("X", colNames[spectCol])
    rownames(IndepenVar) <- rowNames
    colnames(PirTable)[spectCol] <- colnames(IndepenVar)
  } else {
    IndepenVar = NA
  }
  if (nColClass > 0){
    classCol <- (nColSpect+1):(nColSpect+nColClass)
    ClassVar = PirTable[,classCol]
  } else {
    ClassVar = NA
  }
  if (noColDepen > 0){
    depenCol <- c((ncol(PirTable) - noColDepen) : ncol(PirTable))
    DepenVar = PirTable[,depenCol]
  } else {
    DepenVar = NA
  }
  
  out <- list(FullData = PirTable, IndepenVar = IndepenVar, ClassVar = ClassVar, DepenVar = DepenVar)
}

#'	@title Import Pirouette File.
#'	@description Imports directly from a Pirouette file and gives back all the 
#'		data in a list consisting of 4 elements.
#'	@details 'importPir' function imports directly from a Pirouette file and gives back all the 
#'  	data in a list consisting of the following 4 elements:
#'    FullData: matrix containing all the 3 types of variables 
#'    in the same order as in Pirouette software (IndepenVar, ClassVar, DepenVar)
#'    IndepenVar: matrix containing only the independent variables (i.e. spectra)
#'    with the wavelengths as the header starting with 'X'
#'    ClassVar: matrix containing only the variables which were set as class variable in Pirouette
#'    DepenVar: matrix containing only the variables which were set as dependent variable in Pirouette  
#'    The original rownames from the Pirouette file are used in all the 4 matrixes unless they are not unique.
#'    This case the function force to generate uniqe names adding dot and numbering.
#'	@param dataFile Character, a valid path to a Pirouette data file.
#'	@return A list consisting  of 4 elements (FullData, IndepenVar, ClassVar, DepenVar).
#'	@examples \dontrun{
#'	pirData <- importPir("~/Documents/pathToData.pir")
#'	}
#'	@export
importPir <- function(dataFile) {
  Info <- PirInfo(dataFile)
  Indexes <- PirIndexes(dataFile, nByts = Info$NoByts)
  
  PirData <- getPirTable(dataFile, fromPosi=Indexes[1], toPosi=Indexes[2], #toPosi is still not in use
                         nRow=Info$NoRow, nColAll=Info$NoColAll)
  headAll <- getPirHead(dataFile, fromPosi = Indexes[2], toPosi = Indexes[3], lengthHead = Info$NoColAll)
  
  rowNames <- getPirHead(dataFile, fromPosi = Indexes[3], toPosi = Indexes[4], lengthHead = Info$NoRow)
  
  out <- PirMatrix <- creatPirMatrix(PirTable=PirData, colNames=headAll, rowNames=rowNames, nColSpect=Info$NoColSpect, nColClass=Info$NoColClass, noColDepen=Info$NoColDepen)
} #Eof function
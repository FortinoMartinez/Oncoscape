###
#
#       This Script Executes Basic Processing On TCGA Files
#       Specifically It Types, Uppercases and In Cases Enforces Enumeration Types
#       
###

# Configuration -----------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)
os.data.batch.inputFile <- "tcga.filename.manifest.txt"
os.data.batch.outputDir <- "../tcga.clean/"

os.data.batch.inputFile.fileCols <- c("pt", "drug", "rad","f1","f2", "f3","nte","omf","nte_f1","nte_f2")
os.data.batch.inputFile.studyCol <- "study"
os.data.batch.inputFile.dirCol   <- "directory"

dir.create(file.path(os.data.batch.outputDir), showWarnings = FALSE)


# Library Imports ---------------------------------------------------------
library(RUnit)
library(R.utils)
library(stringr)
library(plyr)
library(jsonlite)


#os.tcga.batch.inputFile    <- fromJSON("os.tcga.file.manifest.json")
os.tcga.field.enumerations  <- fromJSON("os.tcga.field.enumerations.json")
os.tcga.column.enumerations <- fromJSON("os.tcga.column.enumerations.json")

# Class Definitions :: Enumerations -------------------------------------------------------
os.enum.na <- c("", "NA", "[NOTAVAILABLE]","[UNKNOWN]","UNKOWN","[NOT AVAILABLE]","[NOT EVALUATED]","UKNOWN","[DISCREPANCY]","[DISCREPANCY]|[DISCREPANCY]",
"NOT LISTED IN MEDICAL RECORD","[NOT APPLICABLE]","[PENDING]","PENDING", "[NOT AVAILABLE]","[PENDING]","[NOTAVAILABLE]","NOT SPECIFIED","N/A",
"[NOT AVAILABLE]|[NOT AVAILABLE]", "[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]", "[NOT APPLICABLE]|[NOT APPLICABLE]","[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]",
"[UNKNOWN]|[UNKNOWN]|[UNKNOWN]|[UNKNOWN];[UNKNOWN]|[UNKNOWN]|[UNKNOWN]", "[UNKNOWN]|[UNKNOWN]|[UNKNOWN]", "[UNKNOWN]|[UNKNOWN]", "[UNKNOWN]|[UNKNOWN]|[UNKNOWN]|[UNKNOWN]",
"[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]","[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]",
"[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]")
#os.enum.other <- c( "OTHER","OTHER: SPECIFY IN NOTES","OTHER (SPECIFY BELOW)","SPECIFY","OTHER REPORTING SCALE")
os.enum.logical.true  <- c("TRUE","YES","1","Y")
os.enum.logical.false <- c("FALSE","NO","0","N")
os.tcga.ignore.columns <- c("bcr_patient_uuid","bcr_drug_uuid","bcr_drug_barcode",
                            "bcr_followup_uuid", "bcr_followup_barcode",
                            "bcr_radiation_uuid","bcr_radiation_barcode","bcr_omf_uuid", "bcr_omf_barcode",
                            "informed_consent_verified", "form_completion_date","project_code", "patient_id") 

# aggregate list of unmapped data & cde id mapping
unmapped.List <- list()
cde.List <- data.frame()
                            
                                         
Map( function(key, value, env=parent.frame()){
  setClass(key)
  setAs("character", key, function(from){ 
    # Convert To Upper + Set NAs  
    from<-toupper(from) 
    from.na<-which(from %in% os.enum.na)
    from[from.na]<-NA    
    
    # Return Enum or NA
    standardVals <- names(os.tcga.field.enumerations[[key]])
    for(fieldName in standardVals){
      values <-os.tcga.field.enumerations[[key]][[fieldName]]
      from[ which(from %in% values)] <- fieldName
    }
    # !!!! need to convert so all standardized names stored for each raw value
    
    if(all(from %in% c(standardVals, NA)))
      return(from)
    
    # Kill If Not In Enum or Na
    stop(paste(key, " not set due to: ", paste(setdiff(from,c(standardVals, NA)), collapse=";"), " not belonging to ", paste(standardVals, collapse=";")))
  })
}, names(os.tcga.field.enumerations), os.tcga.field.enumerations);

# Class Definitions :: TCGA [ID | DATE | CHAR | NUM | BOOL] -------------------------------------------------------

### TCGA ID
setClass("os.class.tcgaId")
setAs("character","os.class.tcgaId", function(from) {
  as.character(str_replace_all(from,"-","." )) 
})

### TCGA Date
setClass("os.class.tcgaDate");
setAs("character","os.class.tcgaDate", function(from){
  
  # Convert Input Character Vector To Uppercase
  from<-toupper(from) 
  
  # Validate Format + Convert Day-Month to 1-1
  if ((str_length(from)==4) && !is.na(as.integer(from) ) ){
    return(format(as.Date(paste(from, "-1-1", sep=""), "%Y-%m-%d"), "%m/%d/%Y"))
  }
  
  # Return NA If Validation Fails
  return(NA)
})

### TCGA Character
setClass("os.class.tcgaCharacter");
setAs("character","os.class.tcgaCharacter", function(from){
  
  # Convert Input Character Vector To Uppercase
  from<-toupper(from) 
  
  # Get Indexes Of Fram Where Value Is In NA
  from.na<-which(from %in% os.enum.na)
  
  # Set From Indexes Values To NA
  from[from.na]<-NA 
  
  return(from)
})

### TCGA Numeric Radiation
setClass("os.class.tcgaNumeric.radiation");
setAs("character","os.class.tcgaNumeric.radiation", function(from){
  
  # Convert Input Character Vector To Uppercase
  from<-toupper(from) 
  
  # Get Indexes Of Fram Where Value Is In NA
  from.na<-which(from %in% os.enum.na)
  
  # Set From Indexes Values To NA
  from[from.na]<-NA 
  
  from<- gsub("MCI|MILLICURIES|-MILLICURIE|MCI (3730 MBQ)|MILLICURIES 131-IODINE", "", from)
  trim(from)
  
  from <- as.numeric(from)
  
  if(all(is.numeric(from))) return (from)
  
  # Kill If Not In Enum or Na
  stop(paste("os.class.tcgaNumeric.radiation not properly set: ", from[!is.numeric(from)], collapse=";"))
  
})


### TCGA Numeric
setClass("os.class.tcgaNumeric");
setAs("character","os.class.tcgaNumeric", function(from){
  
  # Convert Input Character Vector To Uppercase
  from<-toupper(from) 
  
  # Get Indexes Of Fram Where Value Is In NA
  from.na<-which(from %in% os.enum.na)
  
  # Set From Indexes Values To NA
  from[from.na]<-NA 
  
  from <- as.numeric(from)
  
  if(all(is.numeric(from))) return (from)
  
  # Kill If Not In Enum or Na
  stop(paste("os.class.tcgaNumeric not properly set: ", from[!is.numeric(from)], collapse=";"))
  
})

### TCGA Boolean
setClass("os.class.tcgaBoolean");
setAs("character","os.class.tcgaBoolean", function(from){
  
  from<-toupper(from) 
  
  from.na<-which(from %in% os.enum.na)
  from[from.na]<-NA  
  
  from.true <- which( from %in% os.enum.logical.true )
  from[from.true] <- "TRUE"
  
  from.false <- which(from %in% os.enum.logical.false )
  from[from.false] <- "FALSE"
  
  from <- as.logical(from)
  
  # Return Enum or NA        
  if( all(from %in% c( TRUE, FALSE, NA))) return( from )
  
  # Kill If Not In Enum or Na
  stop(paste("os.class.tcgaBoolean not properly set: ", setdiff(from,c( TRUE, FALSE, NA )), collapse=";"))
})

# IO Utility Functions :: [Batch, Load, Save]  -------------------------------------------------------

### Save Function Takes A matrix/data.frame + Base File Path (w/o extension) & Writes to Disk In Multiple (optionally specified) Formats
os.data.save <- function(df, file, format = c("tsv", "csv", "RData","json")){
  
  # Write Tab Delimited
   if("tsv" %in% format)
     write.table(df, file=paste(file,".tsv", sep = ""), quote=F, sep="\t")
  
  # # Write CSV Delimited
   if("csv" %in% format)
     write.csv(df, file=paste(file,".csv",sep = ""), quote = F)
  
  # # Write RData File
   if("RData" %in% format)
     save(df, file=paste(file,".RData", sep = "") )
  
  if("json" %in% format)
   write( toJSON(df, pretty = TRUE), file=paste(file,".json", sep = "") )
   
  # Return DataFrame For Chaining
  return(df)
}

### Load Function Takes An Import File + Column List & Returns A DataFrame
os.data.load <- function(inputFile, checkEnumerations=FALSE, checkClassType = "character"){
  
  # Columns :: Create List From Url
  header <- readLines(inputFile, n=3)
  columns <- unlist(strsplit(header[1],'\t'));
  cde_ids <- unlist(strsplit(header[3],'\t'));

  unMappedData <- list();
  
  if(grepl("clinical_patient_skcm.txt",inputFile)){
  	columns[match("submitted_tumor_site", columns)] = "skcm_tissue_site"
  	columns[match("submitted_tumor_site", columns)] = "skcm_tumor_type"
  }
  if(grepl("follow_up_v2.0_skcm.txt",inputFile)){
  	columns[match("new_tumor_event_type", columns)] = "skcm_tumor_event_type"
  }
  if(grepl("clinical_patient_thca.txt",inputFile)){
    columns[columns=="metastatic_dx_confirmed_by_other"] = "thca_metastatic_dx_confirmed_by_other"
  }
  if(grepl("clinical_patient_kirp.txt",inputFile)){
    columns[columns=="tumor_type"] = "disease_subtype"
  }
  
  # if checkEnumerations - all columns will be read in and assigned 'character' class by default
  # otherwise only classes with defined enumerations will be stored in the mapped table
  if(checkEnumerations) { column.type <- rep("character", length(columns))}
  else                  { column.type <- rep("NULL", length(columns)) }
  

  # assign class types for recognized columns
  #   for each enumerated class type, 
  #     rename matching column to mapped name and assign appropriate type
  os.tcga.classes <- names(os.tcga.column.enumerations)
  for(class.type in os.tcga.classes){
    for(colName in names(os.tcga.column.enumerations[[class.type]])){
      values <-os.tcga.column.enumerations[[class.type]][[colName]]
      matching.values <- which(columns %in% values)
      columns[matching.values ] <- colName
      column.type[ matching.values] <- class.type
    }
  }
  
  # Table :: Read Table From URL
  mappedTable<-read.delim(inputFile,
                          header = FALSE, 
                          skip = 3,
                          dec = ".", 
                          sep = "\t",
                          strip.white = TRUE,
                          check.names=FALSE,
                          numerals = "warn.loss",
                          col.names = columns,
                          colClasses = column.type
  );
  
  if(checkEnumerations) {
    
    # Grab columns matching class type and remove those within the ignore list
    headerWithData <- columns[column.type == checkClassType]
    ignoreCols <- which(headerWithData %in% os.tcga.ignore.columns)
    if(length(ignoreCols > 0))       headerWithData <- headerWithData[- ignoreCols ]
    if(length(headerWithData) == 0)  return(list(mapped=mappedTable, unmapped=unMappedData));
    
    # Discard columns where all values are NA
    DataIndicator <- sapply(headerWithData, function(colName){!all(toupper(mappedTable[,colName]) %in% os.enum.na)})
    headerWithData <- headerWithData[DataIndicator]
    if(length(headerWithData) == 0) return(list(mapped=mappedTable, unmapped=unMappedData));
    
    # Print list of unique values for each column
    unMappedData <- lapply(headerWithData, function(colName){ unique(toupper(mappedTable[,colName]))})
    names(unMappedData) <- headerWithData
    print("---Unused columns")
    print(unMappedData)

  }
  return(list(mapped=mappedTable, unmapped = unMappedData))
  #return(list("mapped"=mappedTable, "unmapped" = unMappedData, "cde"=cbind(columns,cde_ids)))
}

### Batch Is Used To Process Multiple TCGA Files Defined 
os.data.batch <- function(inputFile, outputDirectory, tables,checkEnumerations, ...){
  
  # Load Input File 
  inputFiles <- read.delim(inputFile, sep="\t", header=TRUE)
  
  # Loop Column Wise: for each file type
  for (currentTable in tables)
  {
    # Loop Row Wise: for each disease type
    for (rowIndex in 1:nrow(inputFiles))
    {
      currentDisease   <- inputFiles[ rowIndex, os.data.batch.inputFile.studyCol ];
      currentDirectory <- inputFiles[ rowIndex, os.data.batch.inputFile.dirCol ]
      currentDataFile  <- inputFiles[ rowIndex, currentTable]
      if (is.na(currentDataFile)) next()
      cat(currentDisease, currentTable,"\n")
      inputFile <- paste(currentDirectory, currentDataFile, sep = "")
      outputFile <- paste(outputDirectory, currentDisease, "_", currentTable, sep="")
      
      # Load Data Frame - map and filter by named columns
      MapData <- os.data.load( inputFile = inputFile, ...)
      df <- MapData$mapped
	  unmapped.List <- appendList(unmapped.List, MapData$unmapped)
	  
	  #cde.List <- rbind(cde.List, MapData$cde)

      # Save Data Frame
      os.data.save(
        df = df,
        file = outputFile,
        format = "json")
      
      # Remove Df From Memory
      rm(df)
    }
  }
  if(checkEnumerations){
      os.data.save(
          df = unmapped.List,
          file = "unmapped.List",
          format = "json")
  }
   
}

# Aggregate unmapped column names and classes into a single list  ------------------
appendList <- function (x, val) 
{
    if(!is.list(x) && !is.list(val)) return(x)
    xnames <- names(x)
    for (v in names(val)) {
        x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
            appendList(x[[v]], val[[v]])
        else unique(c(x[[v]], val[[v]]))
    }
    x
}


# Run Block  -------------------------------------------------------
os.data.batch(
  inputFile = os.data.batch.inputFile,
  outputDirectory = os.data.batch.outputDir,
  tables = os.data.batch.inputFile.fileCols, 
  checkEnumerations = FALSE,
  checkClassType = "character");

#tables = os.data.batch.inputFile.fileCols, 

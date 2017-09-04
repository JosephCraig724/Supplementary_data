#######################################################################
##############     Supplementary Data - 1   ###########################
#######################################################################

#######################################################################
##############     R code          ####################################
#######################################################################

#######################################################################
###############       MATRIX CREATION      ############################
#######################################################################

source("../Code/ReadMorphNexus.txt")

if(!require(devtools)) install.packages("devtools")
library(devtools)

install_github("TGuillerme/Claddis")
library(Claddis)

if(!require(gdata)) install.packages("gdata")
library(gdata)

## Read in the data ##

Beck_Nexus<- ReadMorphNexus("../In/2014-Beck-ProcB-matrix-morpho.nex", gap.as.missing=FALSE)
Halliday_Nexus<- ReadMorphNexus("../In/2015-Halliday-LinSoc-matrix-morpho.nex", gap.as.missing=FALSE)

## Extract the matrices ##

Beck_matrix<- Beck_Nexus$matrix
Halliday_matrix<- Halliday_Nexus$matrix

## Extract the row names ##

Beck_species<-row.names(Beck_matrix)
Halliday_species<-row.names(Halliday_matrix)

## Find the matching and unique species ##

Match_Beck<-match(Beck_species, Halliday_species)

Unique_Beck_species<-Beck_species[is.na(Match_Beck)]
Common_Beck_species<-Beck_species[!is.na(Match_Beck)]

Match_Halliday<-match(Halliday_species, Beck_species)

Unique_Halliday_species<-Halliday_species[is.na(Match_Halliday)]
Common_Halliday_species<-Halliday_species[!is.na(Match_Halliday)]

## Becuase it alphebetized them so use function sort ##

sort(Common_Halliday_species)==sort(Common_Beck_species)

## Upload the combination list ##

Combo_Excel<-read.csv("../In/Combining Sheet.csv")

## Create a matrix of becks unique species and Halidays characters ## 

Draft_Beck_matrix<-matrix(NA,nrow=length(Unique_Beck_species),ncol = (ncol(Halliday_matrix)))

## Name the new matrix rows ##

row.names(Draft_Beck_matrix)<-Unique_Beck_species

## Create a loop to place all of the correct characters in the correct columns ##

for(character in 1:nrow(Combo_Excel)){
  cat(paste("replace character" ,Combo_Excel[character, 3], "into", Combo_Excel[character, 2], "\n"))
  Draft_Beck_matrix[,Combo_Excel[character,3]]<-Beck_matrix[match(Unique_Beck_species,row.names(Beck_matrix)),Combo_Excel[character,2]]
}

## Combine the matrix of unique Beck taxa with all of Hallidays taxa and characters ##

Beck_Halliday_matrix_draft_first<-rbind(Halliday_matrix,Draft_Beck_matrix)

## Create a matrix of Becks unique characters ##

New_matrix_beck_unique_characters<-matrix(NA,ncol = 13,nrow= nrow(Halliday_matrix))

## Name the new matrix with Hallidays species ##

row.names(New_matrix_beck_unique_characters)<-row.names(Halliday_matrix)

## Create a Matrix of Becks unique characters with Becks taxa ##

Unique_Beck_mini_matrix<-Beck_matrix[match(Unique_Beck_species,row.names(Beck_matrix)),409:421]

## Combine the two new smaller matrices ##

Beck_characters_matrix<-rbind(New_matrix_beck_unique_characters,Unique_Beck_mini_matrix)

## Combine all into the BIG matrix ##

Big_matrix<-cbind(Beck_Halliday_matrix_draft_first,Beck_characters_matrix)


############################################################################
###############      LOOPING SECTION      ##################################
############################################################################

source("flip.characters.R")

## Read in Flip excel ##

Flip_Loop<-read.csv("../In/FlippingLOOPdaLOOP.csv")

# Flip_Loop[,2] <- as.character(Flip_Loop[,2])
# Flip_Loop[,3] <- as.character(Flip_Loop[,3])

## replacing the question marks (?) in Flip_Loop to be NAs
# Flip_Loop[which(Flip_Loop[,3] == "?"), 3] <- NA


## Extract which species need flipping (all Beck species) ##

Species_list<- which(rownames(Big_matrix) %in% Beck_species)
rownames(Big_matrix)
which(rownames(Big_matrix) %in% Beck_species)
rownames(Big_matrix)[Species_list]

## Create a new Big matrix to compare and protect ##

Big_matrix_flipped<-Big_matrix

## Where are character numbers = Every row, colum one ##

Character_numbers<-unique(Flip_Loop[,1])

## Loop for each character (unique) read left to right flip characters using function ##
for(one_character_number in 1:length(Character_numbers)){ ## For as long as Flip sheet (ie. character numbers) ##
  Character_number_one<-Character_numbers[one_character_number] ## Indiviudual (unique) character ##
  Character <- Big_matrix[,Character_number_one]  ## Extract from Big Matrix ##
  Character_rows<-which(Flip_Loop[,1] %in% Character_number_one) ## Which characters to flip ##
  tmp_list <- t(Flip_Loop[Character_rows,2:3]) ## Flip states in row 2 to be states in row 3 ##
  conv_list <- list() ## Create an empty conversion list ##
  for(tmpcolumn in 1:ncol(tmp_list)){ ## Sub-loop using all columns to fill conv_list and print as character ##
    conv_list[[tmpcolumn]]<-as.character(tmp_list[,tmpcolumn])
  }
  Big_matrix_flipped[,Character_number_one]<-flip.characters(Character, Species_list,conv_list) ## Insert the flipped characters into a New flipped Matrix ##
}

#TG: replacing the question marks (?) introduced from the Flip_Loop into proper NAs
Big_matrix_flipped[which(Big_matrix_flipped == "?")] <- NA

###################################################################################
##########################       WRITING TO NEXUS       ###########################
###################################################################################

## Modified write.nexus.data function to cope with standard data
# write.nexus.std <- ape::write.nexus.data
# body(write.nexus.std)[[2]] <- substitute(format <- match.arg(toupper(format), c("DNA", "PROTEIN", "STANDARD")))
# body(write.nexus.std)[[26]][[3]][[4]] <- substitute(fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, " ", GAP, 
#                                                          " ", INTERLEAVE, " symbols=\"0123456789\";\n"))

####################################################################################
##########################       CREATING CONSTRAINTS       #######################
####################################################################################

install.packages("devtools")
library(devtools)
install_github("dwbapst/paleotree")
library(paleotree)

Hall_Tree_Constraint<-read.tree("../In/Halliday_topology.con.tre")
Beck_Tree_Constraint<-read.tree("../In/Beck_topology.con.tre")
Super_Tree_Constraint<-read.tree("../In/SuperTree.tre")

Hall_Constaint<-createMrBayesConstraints(Hall_Tree_Constraint, partial = TRUE, file = NULL)
Beck_Constraint<-createMrBayesConstraints(Beck_Tree_Constraint, partial = TRUE, file= NULL)
Super_Constraint<-createMrBayesConstraints(Super_Tree_Constraint, partial = TRUE, file= NULL)

# To get the full list of node constraints (ie. without "<Truncated>"), this code must be executed in R as apposed to Rstudio.
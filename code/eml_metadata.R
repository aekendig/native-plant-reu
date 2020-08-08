## Goal: create EML metadata for project using the Environmental Data Initiative's EML Assembly Line (https://github.com/EDIorg/EMLassemblyline)

# Put data and code in a folder together to be grabbed by make_eml
# Generate metadata files by editing current ones or generating them (see Github page for tutorial)
# Edit and run this script


#### set up ####

# clear environment
rm(list=ls())

# load libraries
library(EMLassemblyline)
library(tidyverse)
library(knitr)


#### import templates ####

# list of data files
dlist <- list.files(path = "./data",
                    pattern = ".csv")

# create high level text files
template_core_metadata(path = "./metadata",
                       license = "CCBY")

# create an attribute file for each data table
template_table_attributes(path = "./metadata",
                          data.path = "./data",
                          data.table = dlist)

# create a categorical value file for each data table
template_categorical_variables(path = "./metadata",
                               data.path = "./data")

# look at units
view_unit_dictionary()


#### data file descriptors ####

dlist

# description list
ddlist <- NA

ddlist[1] <- "biomass data collected from Microstegium and native plants"
ddlist[2] <- "infection data collected from Microstegium and native plants"

# name list
dnlist <- c("Biomass data", 
            "Infection data")

# print table
# dtable <- data.frame(data = dlist, description = ddlist)
# kable(dtable)


#### code descriptors ####

# list of code files
clist <- list.files(path = "./code",
                    pattern = ".R")

# remove the eml code
clist <- clist[-2]

# code descripions
cdlist <- c(
  "code to analyze biomass data",
  "code to analyze infection data"
)

# name list
cnlist <- c("Biomass analysis", 
            "Infection analysis")

# print table
# ctable <- data.frame(code = clist, desription = cdlist)
# kable(ctable)


#### make eml ####

# copied data and code from the respective folders and put into edi folder

make_eml(path = "./metadata",
         data.path = "./edi",
         dataset.title = "Emerging fungal pathogen on an invasive grass differentially affects native species",
         data.table = dlist,
         data.table.name = dnlist,
         data.table.description = ddlist,
         data.table.quote.character = rep("\"", length(dlist)),
         other.entity = clist,
         other.entity.name = cnlist,
         other.entity.description = cdlist,
         temporal.coverage = c("2019-06-26", "2019-09-12"),
         geographic.description = "Gainesville, FL, USA",
         geographic.coordinates = c(29.64, -82.36, 29.64, -82.36),
         maintenance.description = "completed", 
         user.id = "aekendig",
         user.domain = "EDI",
         package.id = "edi.63.1")


#### check warnings ####

eml <- EML::read_eml("./metadata/edi.63.1.xml")
EML::eml_validate(eml)



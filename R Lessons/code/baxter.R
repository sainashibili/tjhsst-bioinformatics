library(tidyverse)
library(readxl)

get_metadata <- function(){
  metadata <- read_excel(path="raw_data/baxter.metadata.xlsx",
                         col_types=c(sample = "text", fit_result = "numeric", Site = "text",
                                     Dx_Bin = "text", dx = "text", Hx_Prev = "logical",
                                     Hx_of_Polyps = "logical", Age = "numeric", Gender = "text",
                                     Smoke = "logical", Diabetic = "logical", Hx_Fam_CRC = "logical",
                                     Height = "numeric", Weight = "numeric", NSAID = "logical",
                                     Diabetes_Med = "logical", stage = "text")
  )
  metadata <- mutate(metadata, Height = na_if(Height, 0))
  metadata <- mutate(metadata, Weight = na_if(Weight, 0))
  metadata <- mutate(metadata, Site = recode(.x=Site, "U of Michigan"="U Michigan"))
  metadata <- mutate(metadata, Dx_Bin = recode(.x=Dx_Bin, "Cancer."="Cancer"))
  metadata <- mutate(metadata, Gender = recode(.x=Gender, "f"="female", "m"="male"))
  
  metadata <- rename_all(.tbl=metadata, .funs=tolower)
  metadata <- rename(.data=metadata,
                     previous_history=hx_prev,
                     history_of_polyps=hx_of_polyps,
                     family_history_of_crc=hx_fam_crc,
                     diagnosis_bin=dx_bin,
                     diagnosis=dx,
                     sex=gender)
  
  metadata <- mutate(metadata, diagnosis = factor(diagnosis, levels=c("normal", "adenoma", "cancer")))
  
  #metadata <- mutate(metadata, bmi = get_bmi(weight_kg=weight, height_cm = height)) 
 metadata <- metadata %>%
             mutate(bmi = get_bmi(weight_kg=weight, height_cm=height),
                    bmi_category = get_bmi_category(weight_kg=weight, height_cm=height),
                    obese = is_obese(weight_kg=weight, height_cm=height)
                    )
  
  return(metadata)
}

get_bmi <- function(weight_kg, height_cm){
  return(weight_kg / (height_cm/100) ^ 2)
}

get_bmi_category <- function(weight_kg, height_cm){
  bmi <- get_bmi(weight_kg, height_cm)
  
  bmi_cat <- case_when(bmi >= 30 ~ "obese",
                       bmi >= 25 ~ "overweight",
                       bmi >= 18.5 ~ "normal",
                       bmi > 0 ~ "underweight",
                       TRUE ~ NA_character_
  )
  
  return(bmi_cat)
}

is_obese <- function(weight_kg, height_cm){
  bmi_category <- get_bmi_category(weight_kg, height_cm)
  return(bmi_category == "obese")
}


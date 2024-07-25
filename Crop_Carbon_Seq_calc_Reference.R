library(dplyr)
library(rgcam)
library(tidyr)
library(ggplot2)

##################
#df <- read.csv("L2252.LN5_MgdCarbon_crop.csv")


#India_df <- df %>%
#  filter(region == "India")
### to check if the columns are matching or not..

#same <- India_df$hist.veg.carbon.density == India_df$veg.carbon.density
#counts <- table(same)
###########################

setwd("D:/GCAM_Output_Analysis/GCAM_Crop_CarbonSeq")

df<- read.csv("Crop_C_density_Information.csv")

df <- df %>%
  select(-c(X,year,yield,veg.carbon.density))



################## Establishing Connection ########

conn <- localDBConn('C:/GCAM_V6/GCAM_New/output', 'database_basexdb')

############# Fetching detail land allocation data ###

Reference_crop_area <- addScenario(conn, 'Detail_Land_Allocation.dat','Reference-Default_2100','Land_query.xml')
Reference_crop_area <- loadProject('Detail_Land_Allocation.dat')
scenario_crop_area <- listScenarios(Reference_crop_area)
queries1 <- listQueries(Reference_crop_area, 'Reference-Default_2100')

Detail_LA <- getQuery(Reference_crop_area, 'detailed land allocation') # land allocation data

#### Fetching Agricultural yield data ###

Reference_crop_yield <- addScenario(conn, 'Crop_Yield.dat','Reference-Default_2100','Agri_Tech_Yield.xml')
Reference_crop_yield <- loadProject('Crop_Yield.dat')
scenario_yield <- listScenarios(Reference_crop_yield)
queries2 <- listQueries(Reference_crop_yield, 'Reference-Default_2100')

Crop_Yield <- getQuery (Reference_crop_yield,'ag tech yield')


########### Calculation of Vegetative Carbon Sequestration Following GCAM Methods  ##########


# Function to process data for a given year
process_year <- function(yr) {
  Crop_Yield %>%
    filter(region=="India") %>%
    filter(year==yr) %>%
    rename(LandLeaf=technology) %>%
    select(LandLeaf,year,value)->
    Crop_Yield_year

  India_df <- df %>%
    left_join(Crop_Yield_year, by = 'LandLeaf')

  Detail_LA %>%
    filter(year==yr) %>%
    rename(LandLeaf=landleaf, Area = value) %>%
    mutate(Area = Area*10^9) %>%  # This is to convert thousand km2 to m2
    mutate(Units='m2') %>%
    left_join(India_df, by='LandLeaf') %>%
    rename(Region=region.x) %>%
    mutate(C_density = round(value/ (HarvestIndex) * (1 + Root_Shoot) * (1 - WaterContent) *
                                0.45 * 0.126,  3), # ideally the multiplication factor should be 0.5, however, we have used 0.126 to match with BUR reported value 
           # For tree crops, replace the values calculated above with tree-specific carbon contents elsewhere calculated
           C_density = if_else(is.na(Tree_Cdensity_kgm2), C_density,
                               round(Tree_Cdensity_kgm2, 3))) %>%
    mutate(Veg_carbon_total_Mtonnes=((Area*C_density)/10^9)*3.67) %>%  # Carbon (kg to Million tonnes CO2 eq)
    filter(!is.na(Veg_carbon_total_Mtonnes)) ->
    Detail_LA_year

  Cropwise_Seq <- Detail_LA_year %>%
    separate(LandLeaf, into = c("Crops"), sep = "_", remove = TRUE) %>%
    group_by(Crops) %>%
    summarise(Carbon_Seq = sum(Veg_carbon_total_Mtonnes)) %>%
    mutate(Year = yr)

  return(Cropwise_Seq)
}

# Initialize an empty dataframe to store all results

all_results <- data.frame()

# Loop over every 5 years from 2010 to 2100

for (yr in seq(2010, 2100, by = 5)) {
  yearly_result <- process_year(yr)
  all_results <- bind_rows(all_results, yearly_result)
}

Yearly_Cropland_Seq  <- all_results %>%
  group_by(Year) %>%
  summarise(Cropland_VegC_Seq=sum(Carbon_Seq))


################### Final Calculation of Cropland Sequestration #########

crop_soil_seq <- read.csv("D:/GCAM_Output_Analysis/GCAM_Crop_CarbonSeq/Cropland_Soil_Seq/Output/Cropland_Soil_Emissions_Sequestration_Reference.csv")
#Other_Arable_crop_soil_seq <- read.csv("D:/GCAM_Output_Analysis/GCAM_Crop_CarbonSeq/Cropland_Soil_Seq/Output/Other_Arable_Soil_Emissions_Sequestration_Reference.csv")

years <- seq(2010, 2100, by = 5)

crop_soil_seq_5_years <- crop_soil_seq[crop_soil_seq$Year %in% years, ]
#Other_Arable_soil_seq_5_years <- Other_Arable_crop_soil_seq[Other_Arable_crop_soil_seq$Year %in% years, ]


Cropland_Seq <- Yearly_Cropland_Seq %>%
  left_join(crop_soil_seq_5_years, by='Year') %>%
  select(-c('Carbon_Seq','Carbon_Stock','Carbon_Stock_MtCO2')) %>%
  mutate(Cropland_Seq = ((Cropland_VegC_Seq)*-1) + Carbon_Seq_MtCO2 ) 

write.csv(Cropland_Seq,"Final_Cropland_Sequestration_Reference.csv", row.names = FALSE)
###########################

ggplot(Cropland_Seq, aes(x = Year, y = Cropland_Seq)) +
  geom_bar(stat = "identity", fill = "forestgreen", alpha=0.8) +
  scale_x_continuous(breaks = seq(min(Cropland_Seq$Year), max(Cropland_Seq$Year), by = 10)) +
  labs(title = "Yearly Cropland Vegetation Carbon Sequestration",
       x = "Year",
       y = "Cropland VegC Sequestration (MtCO2)") +
  theme_minimal()

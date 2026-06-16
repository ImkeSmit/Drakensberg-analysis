####Identify growing season from TOMSt data####

#import raw tomst data 
tms <- read.csv("C:\\Users\\imke6\\Documents\\PhD 2025\\Ch2 niche modelling\\All_data\\Clean_data\\microclimate\\Pekka_cleaned_data\\TOMST data\\transfer_402058_files_56a24364\\tomst_data_imputed.csv") |> 
  mutate(area = case_when(area == "BOK" ~ "BK", 
                          area == "GOL" ~ "GG", 
                          area == "WIT" ~ "WH"), 
         grid = str_split_i(grid, "_", 2), 
         cell = toupper(str_split_i(site, "_", 3)), 
         Cell_ID = paste0(area, "_G", grid, "_", cell)) |> 
  select(!c(site, impprop)) |> 
  rename(site = "area") |> 
  #pick variables we are interested in
  filter(stat %in% c("avg", "fdh", "gdh5"), sensor %in% c("T1", "moist"), period %in% c("annual", "summer")) |>  #T1 is the soil temperature, this is the most relaible temperature as it isn't excesively heated by solar radiation
  mutate(variable = paste(sensor, stat, period, sep = "_")) |> 
  select(!c(sensor, stat, period)) |> 
  pivot_wider(names_from = variable, values_from = value) |> 
  arrange(site, grid)
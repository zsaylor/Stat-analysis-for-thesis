library(terra)
library(stats)

library(trend)   # pentru Sen's slope
library(Kendall) # pentru Mann–Kendal



#load rasters
#current_wy <- rast("wyield_current.tif")
#fu_wy_2081_2100 <- rast("wyield_2081_2100_ensemble.tif")

#compareGeom(current_wy, fu_wy_2081_2100, stopOnError = TRUE)

#ext(current_wy)
#ext(fu_wy_2081_2100)

#res(current_wy)
#res(fu_wy_2081_2100)

#crs(current_wy)
#crs(fu_wy_2081_2100)

#fu_aligned <- resample(
#  fu_wy_2081_2100,
#  current_wy,
#  method = "bilinear"
#)

#compareGeom(current_wy, fu_aligned, stopOnError = TRUE)



#v_current <- values(current_wy)
#v_future <- values(fu_aligned)

#idx <- complete.cases(v_current, v_future)
#v_current <- v_current[idx]
#v_future <- v_future[idx]

# Test principal
#t_test <- t.test(v_current, v_future, paired = TRUE)
#boxplot(v_current, v_future)

# Alternativ non-parametric
#W_test <- wilcox.test(v_current, v_future, paired = TRUE)


#diff_r <- fu_aligned - current_wy

#eps <- global(diff_r, "sd", na.rm = TRUE)[1,1]





#change_r <- classify(
#  diff_r,
#  rbind(
#    c(-Inf,    -eps,  -1),  # scade
#    c(-eps,     eps,   0),  # neschimbat
#    c( eps,   Inf,    1)    # creste
#  )
#)

#cols <- c("blue", "grey80", "red")

#plot(
#  change_r,
#  col = cols,
#  breaks = c(-1.5, -0.5, 0.5, 1.5),
#  legend = FALSE,
#  main = "Schimbare future vs current"
#)

#legend(
#  "topright",
#  legend = c("Scădere", "Fără schimbare", "Creștere"),
#  fill = cols,
#  bty = "n"
#)



###########################################

setwd("F:/01_Zane/WYM_zane_new/stat_analysis")

#load rasters
wy_current  <- rast("wyield_current.tif")
wy_21_40   <- rast("wyield_2021_2040_ensemble.tif")
wy_41_60   <- rast("wyield_2041_2060_ensemble.tif")
wy_61_80   <- rast("wyield_2061_2080_ensemble.tif")
wy_81_100   <- rast("wyield_2081_2100_ensemble.tif")

#allign rasters

r_list <- list(wy_current, wy_21_40, wy_41_60, wy_61_80, wy_81_100)

r_aligned <- lapply(
  r_list,
  function(r) resample(r, wy_current, method = "bilinear")
)

lapply(r_aligned, function(r) compareGeom(wy_current, r))

# extragerea valorilor

vals <- lapply(r_aligned, values)

# combinăm într-un singur data.frame
df <- data.frame(
  wy_current = vals[[1]],
  wy_21_40 = vals[[2]],
  wy_41_60 = vals[[3]],
  wy_61_80 = vals[[4]],
  wy_81_100 = vals[[5]]
)

# eliminăm NA-uri simultan
df <- df[complete.cases(df), ]


colnames(df) <- c(
  "1990-2000",
  "2021–2040",
  "2041–2060",
  "2061–2080",
  "2081–2100"
)


# box plot

boxplot(
  df,
  outline = FALSE,
  col = c("grey", "lightblue", "lightgreen", "orange", "red"),
  ylab = "Water Yield (mm)",
  main = "Water yield distibution"
)

# takes too long and doesn't bring to much info
# boxplot(
#  log1p(df),  # log(1 + x)
#  col = c("grey", "lightblue", "lightgreen", "orange", "red"),
#  ylab = "log(Water Yield + 1)",
#  main = "Water yield (log)distibution"
# )


#esantionam setul de date
set.seed(123)

n <- 5000   # 3000–10000 ideal
idx <- sample(seq_len(nrow(df)), n)

df_s <- df[idx, ]


# non_param ANOVA, Friedman test
friedman <- friedman.test(as.matrix(df_s))
friedman # significant difference between periods

#post hoc test to see which groups are different

wilcox <- pairwise.wilcox.test(
  x = as.vector(as.matrix(df_s)),
  g = rep(colnames(df_s), each = nrow(df_s)),
  paired = TRUE,
  p.adjust.method = "holm"
)
 wilcox # all the groups are different 
 

 summary_stats <- apply(df, 2, function(x) {
   c(
     mean = mean(x),
     median = median(x),
     q25 = quantile(x, 0.25),
     q75 = quantile(x, 0.75),
     range = range(x, na.rm = TRUE)
   )
 })
 
 t(summary_stats)

# let s do the stack
 r_stack <- rast(r_aligned) 
  
# put names to bands in the stack
 names(r_stack) <- c(
   "1990_2000",
   "2021_2040",
   "2041_2060",
   "2061_2080",
   "2081_2100"
 )
 
 # allocate times
 time <- c(2000, 2030, 2050, 2070, 2090)
 
 mk_sen_fun <- function(x) {
   
   # elimină NA
   x <- x[!is.na(x)]
   
   # prea puține valori
   if (length(x) < 4) {
     return(c(NA_real_, NA_real_))
   }
   
   # seria este constantă (fără trend posibil)
   if (sd(x) == 0) {
     return(c(0, 1))  # pantă 0, p-value = 1
   }
   
   mk <- Kendall::MannKendall(x)
   sen <- trend::sens.slope(x)$estimates
   
   c(
     sen_slope = sen,
     p_value   = mk$sl
   )
 }

 
 trend_rasters <- app(
   r_stack,
   fun = mk_sen_fun
 )
 
 names(trend_rasters)
 # [1] "sen_slope" "p_value"
 
 writeRaster(
   trend_rasters$`sen_slope.Sen's slope`,
   "sen_slope.tif",
   overwrite = TRUE
 )
 
 
 plot(
   trend_rasters$`sen_slope.Sen's slope`,
   col = hcl.colors(30, "YIOrBr", rev = TRUE),
   main = "Sen's slope (trend temporal WY)"
 )
 
 sig_trend <- trend_rasters$`sen_slope.Sen's slope`
 sig_trend[trend_rasters$p_value >= 0.05] <- NA
 
 writeRaster(
   sig_trend,
   "sig_trend.tif",
   overwrite = TRUE
 )
 
 plot(
   sig_trend,
   col = hcl.colors(30, "RdBu", rev = TRUE),
   main = "Trend temporal Sen – semnificativ statistic (p < 0.05)"
 )
 
 
 #####
 library(sf)
 library(dplyr)
 library(exactextractr)
 library(tidyr)
 
 clc <- st_read("F:/01_Zane/input_data/U2018_CLC2018_V2020_20u1.gpkg", layer = "U2018_CLC2018_V2020_20u1")
 clc_3844 <- st_transform(clc, 3844)
 
 romania <- st_read("F:/01_Zane/ro_admin_country_polygon/ro_admin_country_polygon.shp")
 
 clc_ro <- st_intersection(clc_3844, romania)
 
 st_agr(clc_ro) <- "constant"
 st_crs(clc_ro)
 
 clc_ro$Area_Ha <- as.numeric(st_area(clc_ro)) / 10000
 
 str(clc_ro)
 
 ######
 
 ############################################################
 ## ANALIZA WATER YIELD PE CORINE LEVEL 1
 ## Script complet, adaptat la structura:
 ## Code_18  | Shape (geometry)
 ############################################################
 
 ############################
 ## 0. PACHETE NECESARE
 ############################
 
 library(sf)
 library(terra)
 library(exactextractr)
 library(dplyr)
 library(tidyr)
 
 ############################
 ## 1. VERIFICARE OBIECTE
 ############################
 
 stopifnot(
   inherits(clc_ro, "sf"),
   inherits(wy_current, "SpatRaster"),
   inherits(wy_21_40,  "SpatRaster"),
   inherits(wy_41_60,  "SpatRaster"),
   inherits(wy_61_80,  "SpatRaster"),
   inherits(wy_81_100, "SpatRaster")
 )
 
 ############################
 ## 2. DEFINIRE CORINE LEVEL 1 (CLC_L1)
 ##
 ## Structura reală:
 ## - coloană cod CORINE: Code_18
 ## - geometrie: Shape
 ############################
 
 # Conversie sigură la character
 clc_ro$Code_18 <- as.character(clc_ro$Code_18)
 
 # Definire CORINE Level 1 (regula oficială)
 clc_ro$CLC_L1 <- factor(
   substr(clc_ro$Code_18, 1, 1),
   levels = c("1", "2", "3", "4", "5"),
   labels = c(
     "Artificial",
     "Agricultural",
     "Forest_semi-natural",
     "Wetlands",
     "Water"
   )
 )
 
 # Verificare obligatorie
 print(table(clc_ro$CLC_L1, useNA = "ifany"))
 
 ############################
 ## 3. ZONAL STATISTICS
 ##    MEDIANA Water Yield pe fiecare poligon CORINE
 ##
 ## IMPORTANT:
 ## exact_extract se aplică direct pe obiectul sf
 ############################
 
 clc_ro$WY_current <- exact_extract(wy_current, clc_ro, "median")
 clc_ro$WY_21_40   <- exact_extract(wy_21_40,   clc_ro, "median")
 clc_ro$WY_41_60   <- exact_extract(wy_41_60,   clc_ro, "median")
 clc_ro$WY_61_80   <- exact_extract(wy_61_80,   clc_ro, "median")
 clc_ro$WY_81_100  <- exact_extract(wy_81_100,  clc_ro, "median")
 
 ############################
 ## 4. PREGĂTIRE DATE PENTRU BOXPLOT (FORMAT LONG)
 ############################
 
 df_box <- clc_ro %>%
   st_drop_geometry() %>%
   select(
     CLC_L1,
     WY_current,
     WY_21_40,
     WY_41_60,
     WY_61_80,
     WY_81_100
   ) %>%
   pivot_longer(
     cols = starts_with("WY_"),
     names_to = "Period",
     values_to = "WY"
   )
 
 # Etichete clare ale perioadelor
 df_box$Period <- factor(
   df_box$Period,
   levels = c(
     "WY_current",
     "WY_21_40",
     "WY_41_60",
     "WY_61_80",
     "WY_81_100"
   ),
   labels = c(
     "Current",
     "2021–2040",
     "2041–2060",
     "2061–2080",
     "2081–2100"
   )
 )
 
 ############################
 ## 5. BOXPLOT‑URI
 ############################
 
 lc_cols <- c(
   "Artificial"          = "grey60",
   "Agricultural"        = "gold",
   "Forest_semi-natural" = "forestgreen",
   "Wetlands"            = "turquoise",
   "Water"               = "royalblue"
 )
 
 par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
 
 for (lc in levels(df_box$CLC_L1)) {
   
   sub_df <- df_box[df_box$CLC_L1 == lc, ]
   
   boxplot(
     WY ~ Period,
     data = sub_df,
     outline = FALSE,
     col = lc_cols[lc],
     main = lc,
     ylab = "Annual Water Yield",
     xlab = "Period"
   )
 }
 
 par(mfrow = c(1, 1))
 
 ############################
 ## 6. EXEMPLU TEST STATISTIC
 ##    Friedman – Agricultural
 ############################
 
 ############################################################
 ## ANALIZA WATER YIELD PE CORINE LEVEL 1
 ## Script complet, adaptat la structura:
 ## Code_18  | Shape (geometry)
 ############################################################
 
 ############################
 ## 0. PACHETE NECESARE
 ############################
 
 library(sf)
 library(terra)
 library(exactextractr)
 library(dplyr)
 library(tidyr)
 
 ############################
 ## 1. VERIFICARE OBIECTE
 ############################
 
 stopifnot(
   inherits(clc_ro, "sf"),
   inherits(wy_current, "SpatRaster"),
   inherits(wy_21_40,  "SpatRaster"),
   inherits(wy_41_60,  "SpatRaster"),
   inherits(wy_61_80,  "SpatRaster"),
   inherits(wy_81_100, "SpatRaster")
 )
 
 ############################
 ## 2. DEFINIRE CORINE LEVEL 1 (CLC_L1)
 ##
 ## Structura reală:
 ## - coloană cod CORINE: Code_18
 ## - geometrie: Shape
 ############################
 
 # Conversie sigură la character
 clc_ro$Code_18 <- as.character(clc_ro$Code_18)
 
 # Definire CORINE Level 1 (regula oficială)
 clc_ro$CLC_L1 <- factor(
   substr(clc_ro$Code_18, 1, 1),
   levels = c("1", "2", "3", "4", "5"),
   labels = c(
     "Artificial",
     "Agricultural",
     "Forest_semi-natural",
     "Wetlands",
     "Water"
   )
 )
 
 # Verificare obligatorie
 print(table(clc_ro$CLC_L1, useNA = "ifany"))
 
 ############################
 ## 3. ZONAL STATISTICS
 ##    MEDIANA Water Yield pe fiecare poligon CORINE
 ##
 ## IMPORTANT:
 ## exact_extract se aplică direct pe obiectul sf
 ############################
 
 clc_ro$WY_current <- exact_extract(wy_current, clc_ro, "median")
 clc_ro$WY_21_40   <- exact_extract(wy_21_40,   clc_ro, "median")
 clc_ro$WY_41_60   <- exact_extract(wy_41_60,   clc_ro, "median")
 clc_ro$WY_61_80   <- exact_extract(wy_61_80,   clc_ro, "median")
 clc_ro$WY_81_100  <- exact_extract(wy_81_100,  clc_ro, "median")
 
 ############################
 ## 4. PREGĂTIRE DATE PENTRU BOXPLOT (FORMAT LONG)
 ############################
 
 df_box <- clc_ro %>%
   st_drop_geometry() %>%
   select(
     CLC_L1,
     WY_current,
     WY_21_40,
     WY_41_60,
     WY_61_80,
     WY_81_100
   ) %>%
   pivot_longer(
     cols = starts_with("WY_"),
     names_to = "Period",
     values_to = "WY"
   )
 
 # Etichete clare ale perioadelor
 df_box$Period <- factor(
   df_box$Period,
   levels = c(
     "WY_current",
     "WY_21_40",
     "WY_41_60",
     "WY_61_80",
     "WY_81_100"
   ),
   labels = c(
     "Current",
     "2021–2040",
     "2041–2060",
     "2061–2080",
     "2081–2100"
   )
 )
 
 ############################
 ## 5. BOXPLOT‑URI
 ############################
 
 lc_cols <- c(
   "Artificial"          = "grey60",
   "Agricultural"        = "gold",
   "Forest_semi-natural" = "forestgreen",
   "Wetlands"            = "turquoise",
   "Water"               = "royalblue"
 )
 
 par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
 
 for (lc in levels(df_box$CLC_L1)) {
   
   sub_df <- df_box[df_box$CLC_L1 == lc, ]
   
   boxplot(
     WY ~ Period,
     data = sub_df,
     outline = FALSE,
     col = lc_cols[lc],
     main = lc,
     ylab = "Annual Water Yield",
     xlab = "Period"
   )
 }
 
 par(mfrow = c(1, 1))
 
 ############################
 ## 6. EXEMPLU TEST STATISTIC
 ##    Friedman – Agricultural
 ############################
 
 df_agri <- clc_ro %>%
   st_drop_geometry() %>%
   filter(CLC_L1 == "Agricultural") %>%
   select(
     WY_current,
     WY_21_40,
     WY_41_60,
     WY_61_80,
     WY_81_100
   )
 
 friedman.test(as.matrix(df_agri))
 
 ############################
 ## SFÂRȘIT SCRIPT
 ############################

 
 ############################################################
 ## BOXPLOT: CORINE Level 1 – perioada CURRENT
 ############################################################
 
 ############################################################
 ## BOXPLOT: CORINE Level 1 + toate perioadele
 ############################################################
 
 ############################################################
 ## BOXPLOT GRUPAT:
 ## - X = perioade (o singură etichetă per perioadă)
 ## - 5 boxplot-uri per perioadă (CORINE Level 1)
 ## - culoare = CLC_L1
 ############################################################
 
 # Ordini explicite (evită haosul)
 periods <- levels(df_box$Period)
 lc_levels <- levels(df_box$CLC_L1)
 
 # Culori pentru CORINE Level 1
 lc_cols <- c(
   "Artificial"          = "lightgrey",
   "Agricultural"        = "gold",
   "Forest_semi-natural" = "forestgreen",
   "Wetlands"            = "turquoise",
   "Water"               = "royalblue"
 )
 
 # Lățime și offset pentru box-uri
 n_lc <- length(lc_levels)
 box_width <- 0.12
 offsets <- seq(-2, 2) * box_width
 
 # Inițializare plot "gol"
 plot(
   NA,
   xlim = c(0.5, length(periods) + 0.5),
   ylim = range(df_box$WY, na.rm = TRUE),
   xaxt = "n",
   xlab = "Period",
   ylab = "Annual Water Yield",
   main = "Water Yield by Land Cover and Period"
 )
 
 axis(1, at = seq_along(periods), labels = periods)
 
 # Desenarea boxplot-urilor
 for (j in seq_along(lc_levels)) {
   lc <- lc_levels[j]
   
   for (i in seq_along(periods)) {
     p <- periods[i]
     
     sub <- df_box$WY[
       df_box$Period == p &
         df_box$CLC_L1 == lc
     ]
     
     if (length(sub) > 0) {
       boxplot(
         sub,
         at = i + offsets[j],
         add = TRUE,
         boxwex = box_width,
         outline = FALSE,
         col = lc_cols[lc],
         border = "black"
       )
     }
   }
 }
 
 # Legendă corectă (acum chiar corespunde desenului)
 legend(
   "topright",
   legend = names(lc_cols),
   fill = lc_cols,
   title = "CORINE Level 1",
   bty = "n"
 )

 
 ############################################################
 ## BOXPLOT GRUPAT – axa Y fixată la 700
 ############################################################
 
 # Ordini
 periods <- levels(df_box$Period)
 lc_levels <- levels(df_box$CLC_L1)
 
 # Culori CORINE
 lc_cols <- c(
   "Artificial"          = "lightgrey",
   "Agricultural"        = "gold",
   "Forest_semi-natural" = "forestgreen",
   "Wetlands"            = "turquoise",
   "Water"               = "royalblue"
 )
 
 # Offset pentru box-uri
 box_width <- 0.12
 offsets <- seq(-2, 2) * box_width
 
 # Plot de bază cu Y limitat la 700
 plot(
   NA,
   xlim = c(0.5, length(periods) + 0.5),
   ylim = c(0, 700),          # 👈 AICI este ajustarea
   xaxt = "n",
   xlab = "Period",
   ylab = "Annual Water Yield",
   main = "Water Yield by Land Cover and Period"
 )
 
 axis(1, at = seq_along(periods), labels = periods)
 
 # Boxplot-uri
 for (j in seq_along(lc_levels)) {
   lc <- lc_levels[j]
   for (i in seq_along(periods)) {
     p <- periods[i]
     sub <- df_box$WY[
       df_box$Period == p &
         df_box$CLC_L1 == lc
     ]
     if (length(sub) > 0) {
       boxplot(
         sub,
         at = i + offsets[j],
         add = TRUE,
         boxwex = box_width,
         outline = FALSE,
         col = lc_cols[lc],
         border = "black"
       )
     }
   }
 }
 
 # Legendă coborâtă
 legend(
   "topright",
   inset = c(0.02, 0.08),     # control fin pe poziție
   legend = names(lc_cols),
   fill = lc_cols,
   title = "CORINE Level 1",
   bty = "n",
   cex = 0.9
 )
 
 ############################################################
 ## TEST KRUSKAL–WALLIS
 ## Diferențe între CORINE Level 1
 ## pentru fiecare perioadă
 ############################################################
 
 library(dplyr)
 
 kw_results <- df_box %>%
   group_by(Period) %>%
   summarise(
     KW_p_value = kruskal.test(WY ~ CLC_L1)$p.value
   )
 
 print(kw_results)
 
 ############################################################
 ## POST-HOC WILCOXON
 ## pentru fiecare perioadă
 ############################################################
 
 posthoc_results <- lapply(
   levels(df_box$Period),
   function(p) {
     cat("\n============================\n")
     cat("Perioada:", p, "\n")
     cat("============================\n")
     
     sub <- df_box[df_box$Period == p, ]
     
     print(
       pairwise.wilcox.test(
         sub$WY,
         sub$CLC_L1,
         p.adjust.method = "holm"
       )
     )
   }
 )
 
 
 ################## analiza pe clase de altitudine
 
 ############################################################
 ## EXTRAGERE DATE – DEM class + WY
 ############################################################
 
 ############################################################
 ## 1. SETARE REFERINȚĂ
 ############################################################
 
 ref <- wy_current   # raster de referință
 
 ############################################################
 ## 2. ALINIERE DEM LA WY
 ############################################################
 
 # DEM continuu → bilinear
 dem_r <- resample(dem_clip, ref, method = "bilinear")
 
 ############################################################
 ## 3. CLASIFICARE ALTIMETRICĂ
 ############################################################
 
 rcl <- matrix(c(
   -Inf,  200, 1,
   200,  500, 2,
   500,  800, 3,
   800, 1500, 4,
   1500,  Inf, 5
 ), ncol = 3, byrow = TRUE)
 
 dem_class <- classify(dem_r, rcl)
 dem_class <- as.factor(dem_class)
 
 ############################################################
 ## 4. ALINIERE WY (siguranță)
 ############################################################
 
 wy_current_r <- resample(wy_current, ref, method = "bilinear")
 wy_21_40_r   <- resample(wy_21_40,  ref, method = "bilinear")
 wy_41_60_r   <- resample(wy_41_60,  ref, method = "bilinear")
 wy_61_80_r   <- resample(wy_61_80,  ref, method = "bilinear")
 wy_81_100_r  <- resample(wy_81_100, ref, method = "bilinear")
 
 ############################################################
 ## 5. STACK RASTERE
 ############################################################
 
 r_stack <- c(
   dem_class,
   wy_current_r,
   wy_21_40_r,
   wy_41_60_r,
   wy_61_80_r,
   wy_81_100_r
 )
 
 names(r_stack) <- c(
   "elev_class",
   "WY_current",
   "WY_2021_2040",
   "WY_2041_2060",
   "WY_2061_2080",
   "WY_2081_2100"
 )
 
 
 
 ############################################################
 ## 6. EXTRAGERE VALORI
 ############################################################
 
 df <- as.data.frame(values(r_stack))
 df <- df[complete.cases(df), ]
 
 ############################################################
 ## 7. LONG FORMAT
 ############################################################
 
 df_long <- pivot_longer(
   df,
   cols = starts_with("WY"),
   names_to = "Period",
   values_to = "WY"
 )
 
 df_long$Period <- factor(
   df_long$Period,
   levels = c(
     "WY_current",
     "WY_2021_2040",
     "WY_2041_2060",
     "WY_2061_2080",
     "WY_2081_2100"
   ),
   labels = c(
     "Current",
     "2021–2040",
     "2041–2060",
     "2061–2080",
     "2081–2100"
   )
 )
 
 ############################################################
 ## 8. ETICHETE CLASE ALTIMETRICE
 ############################################################
 
 df_long$elev_class <- factor(
   df_long$elev_class,
   levels = 1:5,
   labels = c(
     "Plains (<200 m)",
     "Hills (200–500 m)",
     "High hills (500–800 m)",
     "Mountains (800–1500 m)",
     "Alpine (>1500 m)"
   )
 )
 
 ############################################################
 ## 9. BOXPLOT FINAL (GRUPAT CORECT)
 ############################################################
 
 periods <- levels(df_long$Period)
 classes <- levels(df_long$elev_class)
 
 # culori (gradient logic)
 cols <- c(
   "#4daf4a",  # plains
   "#a1d76a",  # hills
   "#ffffbf",  # high hills
   "#fdae61",  # mountains
   "#d73027"   # alpine
 )
 
 # layout
 box_width <- 0.12
 offsets <- seq(-2, 2) * box_width
 
 # plot gol
 plot(
   NA,
   xlim = c(0.5, length(periods) + 0.5),
   ylim = c(0, 700),
   xaxt = "n",
   xlab = "Period",
   ylab = "Annual Water Yield",
   main = "Water Yield by Elevation Class and Period"
 )
 
 axis(1, at = seq_along(periods), labels = periods)
 
 # desen boxploturi
 for (j in seq_along(classes)) {
   for (i in seq_along(periods)) {
     
     sub <- df_long$WY[
       df_long$Period == periods[i] &
         df_long$elev_class == classes[j]
     ]
     
     if (length(sub) > 0) {
       boxplot(
         sub,
         at = i + offsets[j],
         add = TRUE,
         boxwex = box_width,
         outline = FALSE,
         col = cols[j],
         border = "black"
       )
     }
   }
 }
 
 # legendă
 legend(
   "topright",
   inset = c(0.02, 0.08),
   legend = classes,
   fill = cols,
   title = "Elevation classes",
   bty = "n"
 )
 
################################## 
# analysis per hydrobasins
 
# library(sf)
# library(terra)
# library(exactextractr)
# library(dplyr)
# library(tidyr)
 
# load shapefile
basins <- st_read("F:/01_Zane/WYM_zane_new/stat_analysis/bazinehidrografice/BazineHidrografice.shp")

basins <- st_transform(basins, crs(wy_current))

basins$WY_current <- exact_extract(wy_current, basins, "median")
basins$WY_21_40   <- exact_extract(wy_21_40,   basins, "median")
basins$WY_41_60   <- exact_extract(wy_41_60,   basins, "median")
basins$WY_61_80   <- exact_extract(wy_61_80,   basins, "median")
basins$WY_81_100  <- exact_extract(wy_81_100,  basins, "median")

df_basins <- basins %>%
  st_drop_geometry()%>%
select(
  denumire,
  WY_current,
  WY_21_40,
  WY_41_60,
  WY_61_80,
  WY_81_100
)

df_long <- df_basins %>%
  pivot_longer(
    cols = starts_with("WY"),
    names_to = "Period",
    values_to = "WY"
  )


df_long$basin <- df_long$denumire


############################################################
## BOXPLOTURI SEPARATE PE BAZINE
############################################################

# listă bazine
basins_names <- unique(df_long$denumire)

# layout: ajustezi în funcție de numărul de bazine
par(mfrow = c(3, 5), mar = c(4, 4, 3, 1))

for (b in basins_names) {
  
  sub <- df_long[df_long$denumire == b, ]
  
  boxplot(
    WY ~ Period,
    data = sub,
    outline = FALSE,
    col = "lightblue",
    main = b,
    ylab = "Water Yield",
    xlab = "Period"
  )
}

par(mfrow = c(1,1))
``

 
 
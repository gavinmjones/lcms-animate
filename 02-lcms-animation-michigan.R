# ==========================================================
#  LCMS Cumulative Disturbances 1985–2024 (Upper Peninsula, MI — Land Only)
# ==========================================================

# --- Libraries ---
library(terra)
library(sf)
library(tigris)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(magick)

options(tigris_use_cache = TRUE)

# ==========================================================
# --- Step 1: Define Upper Peninsula land area (no Great Lakes) ---
# ==========================================================

mi_counties <- counties(state = "MI", year = 2023)
up_counties <- mi_counties %>%
  filter(NAME %in% c(
    "Gogebic", "Ontonagon", "Houghton", "Keweenaw", "Iron", "Baraga",
    "Marquette", "Dickinson", "Menominee", "Alger", "Delta",
    "Schoolcraft", "Luce", "Mackinac", "Chippewa"
  ))

up_shape <- st_union(up_counties)

# Remove large lakes
lakes <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass = "sf")
up_shape <- st_make_valid(up_shape)
lakes <- st_make_valid(lakes)
lakes_union <- st_union(lakes)
lakes_union <- st_transform(lakes_union, st_crs(up_shape))
up_land <- st_difference(up_shape, lakes_union)
up_land <- st_make_valid(up_land)

# ==========================================================
# --- Step 2: Setup file paths and template raster ---
# ==========================================================

years <- 1985:2024
template_path <- "C:/Users/gmjon/Downloads/LCMS_CONUS_v2024-10_Change_Annual_1985/LCMS_CONUS_v2024-10_Change_1985.tif"

r_first <- rast(gsub("1985", years[1], template_path))
up_land_proj <- st_transform(up_land, crs(r_first))
up_vect <- vect(up_land_proj)

# Crop and mask base raster to UP land
cum_r <- crop(r_first, up_vect)
cum_r <- mask(cum_r, up_vect)
values(cum_r) <- NA  # start empty

# ==========================================================
# --- Step 3: Get disturbance category table and colors ---
# ==========================================================

cat_table <- cats(r_first)[[1]]
coltab_first <- coltab(r_first)[[1]]

# Identify columns
name_col <- grep("class", names(cat_table), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(name_col)) {
  legend_labels <- cat_table[[name_col]]
} else {
  legend_labels <- paste("Class", cat_table$ID)
}

plot_codes <- cat_table$ID

# Codes to exclude (non-disturbance)
exclude_codes <- c(14, 15, 16)  # adjust if needed
keep_idx <- !(plot_codes %in% exclude_codes)

plot_codes <- plot_codes[keep_idx]
legend_labels <- legend_labels[keep_idx]

# Match LCMS colors to codes
if (!is.null(coltab_first)) {
  plot_colors <- coltab_first[match(plot_codes, coltab_first$ID), "value"]
} else {
  plot_colors <- terrain.colors(length(plot_codes))
}

# Sanity check
data.frame(ID = plot_codes, Label = legend_labels, Color = plot_colors)

# ==========================================================
# --- Step 4: Output folder ---
# ==========================================================
out_dir <- "C:/Users/gmjon/LCMS_UP_Cumulative"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ==========================================================
# --- Step 5: Loop through years and build cumulative rasters ---
# ==========================================================
pb <- txtProgressBar(min = 0, max = length(years), style = 3)

for (i in seq_along(years)) {
  yr <- years[i]
  r_path <- gsub("1985", yr, template_path)
  
  r <- rast(r_path)
  r <- crop(r, up_vect)
  r <- mask(r, up_vect)
  
  # Remove excluded codes
  r[r %in% exclude_codes] <- NA
  
  # Update cumulative raster: overwrite where new disturbance exists
  r_vals <- values(r)
  cum_vals <- values(cum_r)
  cum_vals[!is.na(r_vals)] <- r_vals[!is.na(r_vals)]
  values(cum_r) <- cum_vals
  
  # Save cumulative raster
  writeRaster(cum_r,
              filename = file.path(out_dir, paste0("Cumulative_UP_", yr, ".tif")),
              datatype = "INT1U",
              overwrite = TRUE)
  
  setTxtProgressBar(pb, i)
}
close(pb)

library(terra)
library(magick)

# ----------------------------
# Step 6: Compute cumulative rasters (safe overwrite)
# ----------------------------

# Remove old cumulative files (start fresh)
existing_files <- list.files(out_dir, pattern = "Cumulative_UP_.*\\.tif$", full.names = TRUE)
if(length(existing_files) > 0) file.remove(existing_files)

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = length(years), style = 3)

for (i in seq_along(years)) {
  yr <- years[i]
  
  # Load annual raster and crop/mask to UP land
  r_path <- gsub("1985", yr, template_path)
  r <- rast(r_path)
  r <- crop(r, up_vect)
  r <- mask(r, up_vect)
  
  # Remove excluded codes (stable, veg growth, etc.)
  r[r %in% exclude_codes] <- NA
  
  # ----------------------------
  # Cumulative calculation:
  # overwrite cum_r only where new disturbance exists
  # ----------------------------
  r_vals <- values(r)
  cum_vals <- values(cum_r)
  cum_vals[!is.na(r_vals)] <- r_vals[!is.na(r_vals)]
  values(cum_r) <- cum_vals
  
  # Save annual cumulative raster
  out_file <- file.path(out_dir, paste0("Cumulative_UP_", yr, ".tif"))
  writeRaster(
    cum_r,
    filename = out_file,
    datatype = "INT1U",
    overwrite = TRUE  # ensure safe overwrite
  )
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Clean memory
  rm(r, r_vals, cum_vals)
  gc()
}

close(pb)


# ----------------------------
# Step 7: Generate GIF (optimized)
# ----------------------------
library(magick)
frames <- list()

# Pre-downsampled for faster plotting
agg_factor <- 4

pb <- txtProgressBar(min = 0, max = length(years), style = 3)

for (i in seq_along(years)) {
  yr <- years[i]
  tif_file <- file.path(out_dir, paste0("Cumulative_UP_", yr, ".tif"))
  r <- rast(tif_file)
  
  # Downsample for plotting only
  r_plot <- aggregate(r, fact = agg_factor, fun = "first")
  
  pngfile <- tempfile(fileext = ".png")
  png(pngfile, width = 1600, height = 900, res = 150, bg = "black")
  par(mar=c(2,2,3,2), bg="black")
  
  plot(r_plot,
       col = plot_colors,
       levels = plot_codes,
       axes = FALSE,
       box = FALSE,
       legend = FALSE,
       main = ""
  )
  
  # Title (upper-right)
  text(
    x = ext(r_plot)[2]-15000,
    y = ext(r_plot)[4]-40000,
    labels = paste0("Cumulative Disturbances: ", yr),
    col="white", cex=1.5, pos=2, font=2
  )
  
  # Legend (all classes, fixed)
  legend(
    "bottomright",
    inset = c(0.15, 0.15),
    legend = legend_labels,
    fill = plot_colors,
    border = NA,
    text.col = "white",
    bg = rgb(0,0,0,0.55),
    cex = 0.9,
    pt.cex = 1.2,
    y.intersp = 0.8,
    xpd = TRUE
  )
  
  dev.off()
  frames[[as.character(yr)]] <- image_read(pngfile)
  
  # Clean memory
  rm(r, r_plot)
  gc()
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Animate and save
animation <- image_animate(image_join(frames), fps = 4, dispose="previous")
image_write(animation, file.path(out_dir, "Cumulative_Disturbances_UP_1985_2024.gif"))

cat("\n✅ GIF complete! Saved to:", file.path(out_dir, "Cumulative_Disturbances_UP_1985_2024.gif"), "\n")
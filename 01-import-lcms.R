# ==========================================================
#  LCMS Raster Cropping for Michigan's Upper Peninsula (Land Only)
# ==========================================================

# --- Install and load necessary packages ---


library(terra)
library(sf)
library(tigris)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

options(tigris_use_cache = TRUE)  # speeds up repeated downloads


# --- Step 1: Get Michigan counties from TIGER/Line ---
mi_counties <- counties(state = "MI", year = 2023)

# --- Step 2: Subset to the Upper Peninsula ---
up_counties <- mi_counties %>%
  filter(NAME %in% c(
    "Gogebic", "Ontonagon", "Houghton", "Keweenaw", "Iron", "Baraga",
    "Marquette", "Dickinson", "Menominee", "Alger", "Delta",
    "Schoolcraft", "Luce", "Mackinac", "Chippewa"
  ))

# Merge all UP counties into a single polygon
up_shape <- st_union(up_counties)

# --- Step 3: Remove large water bodies (e.g., Great Lakes) ---
# Download high-resolution lakes
lakes <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass = "sf")

# Fix invalid geometries for both polygons
up_shape <- st_make_valid(up_shape)
lakes <- st_make_valid(lakes)

# Union all lake polygons into one geometry
lakes_union <- st_union(lakes)

# Ensure both layers share the same CRS
lakes_union <- st_transform(lakes_union, st_crs(up_shape))

# Perform safe difference (land minus lakes)
up_land <- st_difference(up_shape, lakes_union)

# Make geometry valid again after difference
up_land <- st_make_valid(up_land)

# --- Step 4: Read the LCMS raster ---
r_path <- "C:/Users/gmjon/Downloads/LCMS_CONUS_v2024-10_Change_Annual_2024/LCMS_CONUS_v2024-10_Change_2024.tif"
r <- rast(r_path)

# --- Step 5: Reproject AOI and convert to terra vector ---
up_land_proj <- st_transform(up_land, crs(r))
up_vect <- vect(up_land_proj)

# --- Step 6: Crop & mask the raster efficiently ---
r_up <- crop(r, up_vect)
r_up <- mask(r_up, up_vect)

# --- Step 7: Save cropped raster (land only) ---
writeRaster(r_up, "LCMS_UP_Change_2024_LandOnly.tif", datatype = "INT1U", overwrite = TRUE)

# --- Step 8: Check memory use & summary ---
summary(r_up)

# --- Step 9: Minimalist plot (clean, embedded legend) ---
# Plot raster with no outer box or title
tiff("lcms_2024.tif", width=13.3, height=7.5, units="in", res=300, compression="lzw")

par(mar = c(0, 0, 0, 0), bg = "black")
# Plot the raster without axes, box, or title
plot(
  r_up,
  axes = FALSE,
  box = FALSE,
  main = "",
  legend = FALSE,
)
dev.off()


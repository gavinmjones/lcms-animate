# ==============================================================================
# FULL WORKFLOW — LCMS Keweenaw 1985–1995
# ==============================================================================

# -------------------------------
# Load/install packages
# -------------------------------
required <- c("terra", "sf", "fs")
missing_pkg <- setdiff(required, installed.packages()[, "Package"])
if(length(missing_pkg)) install.packages(missing_pkg)

library(terra)
library(sf)
library(fs)

# -------------------------------
# Directories
# -------------------------------
base_dir <- "C:/Users/gmjon/Downloads"   # parent folder with LCMS folders
out_dir <- file.path(base_dir, "LCMS_clipped")
dir_create(out_dir)

# -------------------------------
# Years to process
# -------------------------------
years <- 1987:1995

# Build list of input TIFs
tifs <- file.path(
  base_dir,
  paste0("LCMS_CONUS_v2024-10_Change_Annual_", years),
  paste0("LCMS_CONUS_v2024-10_Change_", years, ".tif")
)
tifs <- normalizePath(tifs, winslash = "/", mustWork = FALSE)

# Keep only existing files
exists <- file_exists(tifs)
if(!all(exists)){
  warning("Some expected TIFs are missing:\n", paste(years[!exists], collapse = ", "))
  tifs <- tifs[exists]
  years <- years[exists]
}
message("Found ", length(tifs), " annual LCMS change rasters.")

# -------------------------------
# Keweenaw polygon
# -------------------------------

keweenaw_bbox <- st_as_sfc(st_bbox(c(
  xmin = -89.7, ymin = 46.2, xmax = -87.2, ymax = 47.6), crs = st_crs(4326)))
keweenaw_vect <- st_sf(geometry = keweenaw_bbox, crs = 4326)


# -------------------------------
# Clipping helpers
# -------------------------------
clip_with_gdal <- function(infile, outfile, cutline, overwrite = FALSE){
  if(file.exists(outfile) && !overwrite) return(TRUE)
  cl_tmp <- tempfile(fileext = ".geojson")
  st_write(cutline, cl_tmp, driver = "GeoJSON", quiet = TRUE)
  args <- c(
    "-cutline", shQuote(cl_tmp),
    "-crop_to_cutline",
    "-of", "GTiff",
    "-dstnodata", "0",
    shQuote(infile), shQuote(outfile)
  )
  system2("gdalwarp", args = args, stdout = TRUE, stderr = TRUE)
  invisible(TRUE)
}

clip_with_terra <- function(infile, outfile, cutline, overwrite = FALSE){
  if(file.exists(outfile) && !overwrite) return(TRUE)
  r <- rast(infile)
  cutline_prj <- st_transform(cutline, crs = st_crs(r))
  v <- vect(cutline_prj)
  r_clipped <- crop(r, v)
  r_clipped <- mask(r_clipped, v)
  writeRaster(r_clipped, outfile, overwrite = TRUE, gdal=c("COMPRESS=LZW"), datatype="INT1U")
  invisible(TRUE)
}

use_gdalwarp <- nchar(Sys.which("gdalwarp")) > 0

# -------------------------------
# Clip rasters
# -------------------------------
clipped_files <- list()
for(i in seq_along(tifs)){
  infile <- tifs[i]
  year <- years[i]
  outname <- file.path(out_dir, paste0("LCMS_Change_", year, "_Keweenaw.tif"))
  
  message("\nProcessing year ", year)
  if(use_gdalwarp){
    tryCatch({
      clip_with_gdal(infile, outname, keweenaw_vect)
    }, error=function(e){
      message("gdalwarp failed, using terra::crop")
      clip_with_terra(infile, outname, keweenaw_vect)
    })
  } else {
    clip_with_terra(infile, outname, keweenaw_vect)
  }
  
  clipped_files[[as.character(year)]] <- outname
}

# -------------------------------
# Water mask using rnaturalearthdata
# -------------------------------
if(!requireNamespace("rnaturalearth", quietly=TRUE)) install.packages("rnaturalearth")
if(!requireNamespace("rnaturalearthdata", quietly=TRUE)) install.packages("rnaturalearthdata")
library(rnaturalearth)
library(rnaturalearthdata)

# Download lakes as sf
lakes <- ne_download(
  scale = "large",
  type = "lakes",
  category = "physical",
  returnclass = "sf"
)

# Transform lakes to match raster CRS
r_ref <- rast(clipped_files[[1]])           # reference raster
lakes_prj <- st_transform(lakes, crs=st_crs(r_ref))

# Transform Keweenaw polygon to same CRS as lakes
keweenaw_prj <- st_transform(keweenaw_vect, st_crs(lakes_prj))

# Intersection while both are sf
lakes_crop_sf <- st_intersection(lakes_prj, keweenaw_prj)

# Convert to terra SpatVector
lakes_vect <- vect(lakes_crop_sf)

# Loop through clipped rasters and mask water
masked_files <- list()
for(yr in names(clipped_files)){
  r <- rast(clipped_files[[yr]])
  
  # Rasterize lakes: cells overlapping water get NA, other cells stay as-is
  r_mask <- rasterize(lakes_vect, r, field=1)   # water cells = 1, others NA
  r[r_mask == 1] <- NA                          # set water to NA
  
  # Save masked raster
  outname <- file.path(out_dir, paste0("LCMS_Change_", yr, "_Keweenaw_masked.tif"))
  writeRaster(r, outname, overwrite=TRUE, gdal=c("COMPRESS=LZW"), datatype="INT1U")
  
  masked_files[[yr]] <- outname
}



# -------------------------------
# Fixed categories and colors
# -------------------------------
all_categories <- c(
  "Stable",
  "Mechanical Land Transformation",
  "Insect, Disease, or Drought Stress",
  "Vegetation Successional Growth",
  "Other Loss",
  "Non-Processing Area Mask",
  "Wildfire",
  "Tree Removal",
  "Desiccation",
  "Inundation",
  "Other/Unknown"
)

colors <- c(
  "Stable" = "#b2b2b2",
  "Mechanical Land Transformation" = "#e41a1c",
  "Insect, Disease, or Drought Stress" = "#377eb8",
  "Vegetation Successional Growth"= "#b2b2b2", #"#4daf4a",
  "Other Loss" = "#ff7f00",
  "Non-Processing Area Mask" = "#984ea3",
  "Wildfire" = "#f781bf",
  "Tree Removal" = "#a65628",
  "Desiccation" = "#ffff33",
  "Inundation" = "#00ffff",
  "Other/Unknown" = "grey"
)

# -------------------------------
# Plot function
# -------------------------------
# # -------------------------------
# # Plot function
# # -------------------------------
plot_lcms_consistent <- function(raster_file, year){
  r <- rast(raster_file)
  if(!is.factor(r)) r <- as.factor(r)
  
  lev_table <- levels(r)[[1]]
  
  # Ensure there is a 'category' column; if missing, create it
  if(!"category" %in% names(lev_table)){
    lev_table$category <- paste0("Class_", lev_table$value)
  }
  
  lev_table$category <- trimws(as.character(lev_table$category))
  
  # Map unknown categories to 'Other/Unknown'
  lev_table$category[!lev_table$category %in% all_categories] <- "Other/Unknown"
  lev_table$col <- colors[lev_table$category]
  
  # Ensure all colors are assigned
  lev_table$col[is.na(lev_table$col)] <- "grey"
  
  levels(r)[[1]] <- lev_table
  
  # Remove margins
  par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  
  # Plot without padding, NA (water) as black
  plot(r,
       col=lev_table$col,
       colNA="black",        # <--- add this line
       legend=FALSE,
       axes=FALSE,
       box=FALSE,
       asp=NA,
       plg=list(shrink=FALSE))
}


# -------------------------------
# Plot all years 1985–1995
# -------------------------------
for(yr in names(masked_files)){
  plot_lcms_consistent(masked_files[[yr]], yr)
}






# -------------------------------
# Directory to store PNG frames
# -------------------------------
frames_dir <- file.path(tempdir(), "lcms_frames_cumulative")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# Helper function: write PNG for a raster with progress bar
# -------------------------------
plot_lcms_to_png <- function(raster_file, year, year_index, total_years, out_file, w = 1300, h = 800, res = 150){
  png(out_file, width = 13.3, height = 7.5, res = 300, units = "in", bg = "black")
  par(mar = c(0,0,4,0), oma = c(0,0,0,0))  # extra top margin for year
  
  # -------------------------------
  # Plot the raster
  # -------------------------------
  r <- rast(raster_file)
  if(!is.factor(r)) r <- as.factor(r)
  
  lev_table <- levels(r)[[1]]
  if(!"category" %in% names(lev_table)) lev_table$category <- paste0("Class_", lev_table$value)
  lev_table$category <- trimws(as.character(lev_table$category))
  lev_table$category[!lev_table$category %in% all_categories] <- "Other/Unknown"
  lev_table$col <- colors[lev_table$category]
  lev_table$col[is.na(lev_table$col)] <- "grey"
  levels(r)[[1]] <- lev_table
  
  plot(r,
       col = lev_table$col,
       colNA = "black",
       legend = FALSE, axes = FALSE, box = FALSE,
       asp = NA, plg = list(shrink = FALSE))
  
  # -------------------------------
  # Year label
  # -------------------------------
  mtext(paste0("Year: ", year), side = 3, line = -1, adj = 0.02, cex = 1.7, font = 2, col = "white")
  
  # -------------------------------
  # Draw progress bar
  # -------------------------------
  bar_height <- 0.03  # fraction of device height
  bar_bottom <- 0.02
  bar_top <- bar_bottom + bar_height
  
  # Full bar (grey background)
  rect(grconvertX(0, "npc"), grconvertY(bar_bottom, "npc"),
       grconvertX(1, "npc"), grconvertY(bar_top, "npc"),
       col = "grey20", border = NA)
  
  # White tick for current year
  tick_fraction <- (year_index - 1) / (total_years - 1)  # 0..1
  tick_width <- 0.02
  rect(grconvertX(tick_fraction - tick_width/2, "npc"), grconvertY(bar_bottom, "npc"),
       grconvertX(tick_fraction + tick_width/2, "npc"), grconvertY(bar_top, "npc"),
       col = "white", border = NA)
  
  dev.off()
}

# -------------------------------
# 1) Generate PNG frames for each year
# -------------------------------
years <- names(masked_files)
total_years <- length(years)
frames_dir <- file.path(tempdir(), "lcms_frames_progress")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)
png_files <- c()

for(i in seq_along(years)){
  yr <- years[i]
  out_png <- file.path(frames_dir, paste0("lcms_", yr, ".png"))
  plot_lcms_to_png(masked_files[[yr]], yr, year_index = i, total_years = total_years, out_file = out_png)
  png_files <- c(png_files, out_png)
  message("Wrote frame: ", out_png)
}


# -------------------------------
# 3) Animate and save GIF
# -------------------------------
library(magick)

img <- image_read(png_files)        # read all frames
img_animated <- image_animate(img, fps = 4)   # adjust fps as desired
image_write(img_animated, "lcms_animation.gif")




library(stars)
library(terra)

# focus area with buffer
fa = st_read("AIRREG/focus_area.gpkg")
fa_buff = st_buffer(fa, 30000)


# land cover
clc_proxy = read_stars("/palma/scratch/tmp/jheisig/aq/supplementary/static/lcv_landcover.clc_corine_c_100m_0..0cm_2018_eumap_epsg3035_v2020.tif")
clc = st_as_stars(clc_proxy[fa_buff]) |> 
  setNames("Land_Cover") |> 
  mutate(Land_Cover = as.factor(Land_Cover),
         nat_areas_bin = ifelse(Land_Cover %in% c(23:34), 1, 0))
clc$CLC_Nat = clc["nat_areas_bin"] |> 
           rast() |> 
           focal(w = 5, fun = mean) |> 
           st_as_stars() |> pull()
clc$nat_areas_bin = NULL
plot(merge(clc))

# elevation
cop_dem_proxy = read_stars("supplementary/static/COP-DEM/COP_DEM_Europe_100m_epsg3035.tif")
cop_dem = st_as_stars(cop_dem_proxy[fa_buff]) |> setNames("Elevation")
plot(cop_dem)

# population
pop_proxy = read_stars("supplementary/static/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R4_C19.tif") 
pop_proxy = pop_proxy[st_transform(st_buffer(fa_buff, 5000), st_crs(pop_proxy)), drop=T] 
pop_proxy = st_as_stars(pop_proxy)
pop_warp = st_warp(pop_proxy, clc)
pop = pop_warp[fa_buff] |> setNames("Population")
pop
plot(pop)

# traffic
roi = st_as_text(st_as_sfc(st_bbox(fa_buff)))

grip_1 = st_read("supplementary/static/GRIP/grip_type_1_buff_chunks.gpkg",
                  geometry_column="geom", wkt_filter=roi)
grip_2 = st_read("supplementary/static/GRIP/grip_type_2_buff_chunks.gpkg",
                 geometry_column="geom", wkt_filter=roi)
grip_3 = st_read("supplementary/static/GRIP/grip_type_3_buff_chunks.gpkg",
                 geometry_column="geom", wkt_filter=roi)

mapview::mapview(list(grip_1, grip_2, grip_3), color = list("blue", "red", "purple"))


tmplt = clc["Land_Cover"]
tmplt$Land_Cover = NULL
tmplt$g1 = NA

r1 = rasterizeGeom(vect(grip_1), rast(tmplt), fun = "area", unit = "m") |> st_as_stars()
r2 = rasterizeGeom(vect(grip_2), rast(tmplt), fun = "area", unit = "m") |> st_as_stars()
r3 = rasterizeGeom(vect(grip_3), rast(tmplt), fun = "area", unit = "m") |> st_as_stars()

# merge
roads = c(r1,r2,r3) |> setNames(c("Highway", "Major_Roads", "Minor_Roads"))
other = c(clc, cop_dem, pop)
st_dimensions(other) = st_dimensions(roads)

cov = c(other, roads)
saveRDS(cov, "AIRREG/supplementary/static_covariates_buffered.rds")
saveRDS(cov[fa[1,]], "AIRREG/supplementary/static_covariates_focus_area_1.rds")
saveRDS(cov[fa[2,]], "AIRREG/supplementary/static_covariates_focus_area_2.rds")

cov_3d = merge(cov, name="bands")
write_stars(cov_3d, "AIRREG/supplementary/static_covariates_buffered.tif")
write_stars(cov_3d[fa[1,]], "AIRREG/supplementary/static_covariates_focus_area_1.tif")
write_stars(cov_3d[fa[2,]], "AIRREG/supplementary/static_covariates_focus_area_2.tif")




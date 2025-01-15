"""
Implements the three hierarchical decision trees from:
Miller, B. A., & Schaetzl, R. J. (2015).
"Digital Classification of Hillslope Position."
Soil Science Society of America Journal, 79(1), 132–145.

- Each tree (SlgP, PrcP, ReeP) follows the logical flow from Figure 3.
- Thresholds are user-defined (no field-based calibration).
- Relative Elevation can be computed using two approaches:
    1) 'mean' => DEM - local_mean
    2) 'minmax' => DEM - ((focal_min + focal_max) - DEM)

After classification, we run a final step that masks out any pixels where the
DEM is NoData—ensuring that "outside" remains NoData instead of being labeled.

Requires: ArcPy + Spatial Analyst + 3D Analyst extension
"""

import os
import arcpy
from arcpy.sa import (
    Fill,
    FocalStatistics,
    NbrCircle,
    NbrRectangle,
    Raster,
    Con,
    Minus,
    IsNull,
    SetNull
)
from arcpy.ddd import SurfaceParameters
import numpy as np

#-----------------------------------------------------------------------
# FUNCTION: compute_terrain_derivatives
#-----------------------------------------------------------------------
def compute_terrain_derivatives(
    dem_path,
    out_folder,
    slope_radius=9,
    curvature_radius=63,
    elev_radius=135,
    rel_elev_method="mean",
    nodata_val=-9999,
    overwrite=True
):
    """
    Computes slope, curvature, and relative elevation using the author's approach:

      1) Fill the DEM to remove depressions (optional step from the original code).
      2) Slope => use SurfaceParameters(..., parameter_type="SLOPE", ...)
         - By default: 9m, QUADRATIC, PERCENT_RISE
      3) Profile Curvature => from SurfaceParameters(..., 63m, QUADRATIC)
      4) Relative Elevation => either 'mean' (DEM - focal_mean)
         or 'minmax' (DEM - ((focal_min + focal_max) - DEM))

    Returns the paths to the *raw* rasters for slope, curvature, and relative elevation.

    Parameters
    ----------
    dem_path : str
        Path to the input DEM raster.
    out_folder : str
        Folder to store intermediate rasters.
    slope_radius : int
        Neighborhood distance (in meters) for slope calculation.
    curvature_radius : int
        Neighborhood distance (in meters) for curvature calculation.
    elev_radius : int
        Neighborhood distance (in cells for 'MAP') for relative elevation.
    rel_elev_method : str
        Method to compute relative elevation: "mean" or "minmax".
    nodata_val : int or float
        Value used as NoData in the output NumPy arrays.
    overwrite : bool
        Whether to overwrite existing output rasters.

    Returns
    -------
    slope_raw_path : str
        Path to the slope raster.
    profc_raw_path : str
        Path to the profile curvature raster.
    rel_raw_path : str
        Path to the relative elevation raster.
    """
    #---
    # 1) Check out required ArcGIS extensions
    #---
    print("[DEBUG] Checking out Spatial and 3D Analyst extensions.")
    arcpy.CheckOutExtension("Spatial")
    arcpy.CheckOutExtension("3D")

    #---
    # 2) Set environment overwrite
    #---
    print(f"[DEBUG] Setting arcpy.env.overwriteOutput = {overwrite}")
    arcpy.env.overwriteOutput = overwrite

    #---
    # 3) Create output folder if needed
    #---
    if not os.path.isdir(out_folder):
        print(f"[DEBUG] Creating directory: {out_folder}")
        os.makedirs(out_folder)
    else:
        print(f"[DEBUG] Output folder already exists: {out_folder}")

    #---
    # 4) Fill DEM to remove depressions
    #   (Some workflows skip fill if the DEM is already processed,
    #    but we'll follow the original approach.)
    #---
    print("[DEBUG] Filling DEM to remove any depressions.")
    dem_filled_path = os.path.join(out_folder, "DEM_Filled.tif")
    fill_dem = Fill(dem_path)
    fill_dem.save(dem_filled_path)

    #---
    # 5) Compute Slope (PERCENT_RISE, QUADRATIC)
    #---
    print("[DEBUG] Computing slope (percent rise).")
    slope_raw_path = os.path.join(out_folder, "Slope_Raw.tif")
    SurfaceParameters(
        in_raster=dem_filled_path,
        out_raster=slope_raw_path,
        parameter_type="SLOPE",
        local_surface_type="QUADRATIC",
        neighborhood_distance=f"{slope_radius} Meters",
        use_adaptive_neighborhood="FIXED_NEIGHBORHOOD",
        z_unit="Meter",
        output_slope_measurement="PERCENT_RISE"
    )

    #---
    # 6) Compute Profile Curvature (63 m, QUADRATIC)
    #---
    print("[DEBUG] Computing profile curvature.")
    profc_raw_path = os.path.join(out_folder, "ProfileCurv_Raw.tif")
    SurfaceParameters(
        in_raster=dem_filled_path,
        out_raster=profc_raw_path,
        parameter_type="PROFILE_CURVATURE",
        local_surface_type="QUADRATIC",
        neighborhood_distance=f"{curvature_radius} Meters",
        use_adaptive_neighborhood="FIXED_NEIGHBORHOOD",
        z_unit="Meter",
        output_slope_measurement="PERCENT_RISE"
    )

    #---
    # 7) Compute Relative Elevation
    #       "mean"   => dem_filled - focal_mean
    #       "minmax" => dem_filled - ((focal_min + focal_max) - dem_filled)
    #---
    print(f"[DEBUG] Computing relative elevation: method='{rel_elev_method}'.")
    rel_raw_path = os.path.join(out_folder, "RelElev_Raw.tif")
    dem_filled = Raster(dem_filled_path)

    if rel_elev_method.lower() == "minmax":
        #--- minmax approach
        print("[DEBUG] Using 'minmax' => DEM - ((min + max) - DEM).")
        rect_nbr = NbrRectangle(elev_radius, elev_radius, "MAP")
        max_focal = FocalStatistics(dem_filled, rect_nbr, "MAXIMUM", "")
        min_focal = FocalStatistics(dem_filled, rect_nbr, "MINIMUM", "")
        rel_calc = dem_filled - ((min_focal + max_focal) - dem_filled)
    else:
        #--- 'mean' approach
        print("[DEBUG] Using 'mean' => DEM - focal_mean.")
        circle_nbr = NbrCircle(elev_radius, "MAP")
        focal_mean = FocalStatistics(dem_filled, circle_nbr, "MEAN", "")
        rel_calc = dem_filled - focal_mean

    rel_calc.save(rel_raw_path)

    #---
    # 8) Return paths to the derived rasters
    #---
    return slope_raw_path, profc_raw_path, rel_raw_path


#-----------------------------------------------------------------------
# FUNCTION: classify_hillslopes
#-----------------------------------------------------------------------
def classify_hillslopes(
    dem_path,
    out_folder,
    method="SlgP",
    slope_break_low=2.4,
    slope_break_high=5.1,
    slope_break=10.0,
    curvature_neg_tol=0.0,
    curvature_pos_tol=0.0,
    rel_elev_break=0.0,
    rel_elev_high=10.0,
    rel_elev_low=-10.0,
    slope_radius=9,
    curvature_radius=63,
    elev_radius=135,
    rel_elev_method="mean",
    no_data_value=-9999,
    overwrite=True
):
    """
    Classifies each pixel into one of five hillslope positions (Summit=1, 
    Shoulder=2, Backslope=3, Footslope=4, Toeslope=5) using one of three 
    hierarchical methods from Miller & Schaetzl (2015).

    The final classification is written to a raster, and then a masked version
    (based on the original DEM's NoData) is also created. The masked version is
    typically the final output used by the workflow.

    Parameters
    ----------
    dem_path : str
        Path to the input DEM raster.
    out_folder : str
        Folder where outputs will be stored.
    method : str
        One of {"SlgP", "PrcP", "ReeP"} specifying which hierarchy to use.
    slope_break_low : float
        Lower slope threshold for the slope-based reclass.
    slope_break_high : float
        Upper slope threshold for the slope-based reclass.
    slope_break : float
        A single threshold used by other classification hierarchies.
    curvature_neg_tol : float
        Tolerance for negative curvature (currently not used explicitly, 
        but can be integrated if needed).
    curvature_pos_tol : float
        Tolerance for positive curvature (currently not used explicitly).
    rel_elev_break : float
        A threshold for relative elevation (not used in this example, can be extended).
    rel_elev_high : float
        Possible upper threshold for relative elevation (not used in this example).
    rel_elev_low : float
        Possible lower threshold for relative elevation (not used in this example).
    slope_radius : int
        Neighborhood size (meters) for slope.
    curvature_radius : int
        Neighborhood size (meters) for curvature.
    elev_radius : int
        Neighborhood size (in 'MAP' units, typically cells or meters) for relative elevation.
    rel_elev_method : str
        Either "mean" or "minmax" for relative elevation approach.
    no_data_value : int or float
        Value to treat as NoData in the NumPy arrays. 
    overwrite : bool
        If True, allows overwriting existing rasters.

    Returns
    -------
    masked_out_path : str
        Path to the final classification raster, masked to the original DEM's NoData.
    """

    #---
    # 1) Check out extensions & set env
    #---
    print("[DEBUG] Checking out Spatial & 3D extensions, setting overwrite...")
    arcpy.CheckOutExtension("Spatial")
    arcpy.CheckOutExtension("3D")
    arcpy.env.overwriteOutput = overwrite

    #---
    # 2) Compute terrain derivatives (slope, curvature, relative elev)
    #---
    print("[DEBUG] Computing terrain derivatives for classification.")
    slope_raw, curv_raw, rel_raw = compute_terrain_derivatives(
        dem_path=dem_path,
        out_folder=out_folder,
        slope_radius=slope_radius,
        curvature_radius=curvature_radius,
        elev_radius=elev_radius,
        rel_elev_method=rel_elev_method,
        nodata_val=no_data_value,
        overwrite=overwrite
    )

    #---
    # 3) Load the derived rasters + DEM as NumPy arrays
    #   NOTE: Using no_data_value=-9999 by default to avoid conflict 
    #         with real slope=0 or curvature=0.
    #---
    print("[DEBUG] Converting slope, curvature, relative elev to NumPy arrays.")
    slope_arr = arcpy.RasterToNumPyArray(slope_raw, nodata_to_value=no_data_value)
    curv_arr  = arcpy.RasterToNumPyArray(curv_raw,  nodata_to_value=no_data_value)
    rel_arr   = arcpy.RasterToNumPyArray(rel_raw,   nodata_to_value=no_data_value)

    #---
    # 4) Also load the filled DEM for final references (extent, cell size).
    #---
    print("[DEBUG] Reading filled DEM for final references.")
    filled_dem_path = os.path.join(out_folder, "DEM_Filled.tif")
    filled_dem      = Raster(filled_dem_path)
    dem_arr         = arcpy.RasterToNumPyArray(filled_dem, nodata_to_value=no_data_value)

    #--- Collect geo-referencing info
    lower_left   = arcpy.Point(filled_dem.extent.XMin, filled_dem.extent.YMin)
    cell_width   = filled_dem.meanCellWidth
    cell_height  = filled_dem.meanCellHeight
    spatial_ref  = filled_dem.spatialReference

    #---
    # 5) Prepare final classification array + valid mask
    #   We ONLY classify where slope_arr, curv_arr, rel_arr != no_data_value.
    #---
    print("[DEBUG] Creating classification array (5 classes) and valid_mask.")
    final_class = np.zeros_like(slope_arr, dtype=np.int32)  # default 0 => unclassified

    valid_mask = (
        (slope_arr != no_data_value) &
        (curv_arr  != no_data_value) &
        (rel_arr   != no_data_value)
    )

    #-------------------------------------------------------------------
    # 6) Reclassify slope, curvature, relative elev => intermediate codes
    #    - slope_code: -1 (low), 0 (middle), +1 (high)
    #    - curvature_code: -1 (negative), +1 (positive)
    #    - rel_elev_code: -1 (below threshold), +1 (above threshold)
    #
    #    Note: We ignore "tolerance" values for curvature in this example,
    #    but you can incorporate them by checking, e.g.:
    #      curv_arr < curvature_neg_tol => -1
    #      curv_arr > curvature_pos_tol => +1
    #      else => 0 (flat/near 0 curvature)
    #-------------------------------------------------------------------
    slope_code = np.where(
        slope_arr < slope_break_low, 
        -1, 
        np.where(slope_arr > slope_break_high, 1, 0)
    )

    #--- Example (binary curvature code):
    #    you could do: curv_arr < -0.001 => -1, curv_arr > +0.001 => +1, else => 0
    #    For now, we'll do a strict negative/positive split at 0.0:
    curvature_code = np.where(curv_arr < 0.0, -1, 1)

    #--- For relative elev, we'll just do <0 => -1, >=0 => +1
    rel_elev_code = np.where(rel_arr < 0.0, -1, 1)

    #-------------------------------------------------------------------
    # 7) Define the hierarchical classification methods
    #
    #    SlgP => Slope Gradient Priority:
    #       slope==1 => Backslope(3)
    #       slope==0 => if curvature==1 => Shoulder(2) else => Footslope(4)
    #       slope==-1 => if rel_elev==1 => Summit(1) else => Toeslope(5)
    #
    #    PrcP => Profile Curvature Priority:
    #       [Based on curvature first, then slope, then rel_elev]
    #
    #    ReeP => Relative Elevation Priority:
    #       [Based on relative elevation first, then slope, then curvature]
    #-------------------------------------------------------------------
    def slope_gradient_priority():
        #---
        # Only classify where valid_mask is True
        #---
        sp = (slope_code == 1) & valid_mask
        sz = (slope_code == 0) & valid_mask
        sn = (slope_code == -1) & valid_mask

        #--- slope==1 => Backslope=3
        final_class[sp] = 3

        #--- slope==0 => Shoulder if curvature==1 else Footslope
        zero_shoulder = sz & (curvature_code == 1)
        final_class[zero_shoulder] = 2
        zero_foot = sz & (final_class == 0)
        final_class[zero_foot] = 4

        #--- slope==-1 => Summit if rel_elev==1 else Toeslope
        neg_sum = sn & (rel_elev_code == 1)
        final_class[neg_sum] = 1
        neg_toe = sn & (final_class == 0)
        final_class[neg_toe] = 5

    def profile_curvature_priority():
        #---
        # Priority: curvature => slope => relative elev
        #---
        mask_all = valid_mask & (final_class == 0)

        #--- curvature==+ => Shoulder=2, curvature==- => Footslope=4
        curv_plus = mask_all & (curvature_code == 1)
        final_class[curv_plus] = 2
        curv_neg = mask_all & (curvature_code == -1)
        final_class[curv_neg] = 4

        #--- next step: slope
        pass_mask = (final_class == 0) & valid_mask
        high_slope = pass_mask & (slope_arr > slope_break)
        final_class[high_slope] = 3  # Backslope=3

        #--- next step: rel elev => Summit=1, else Toeslope=5
        pass_mask2 = (final_class == 0) & valid_mask
        sum_mask = pass_mask2 & (rel_elev_code == 1)
        final_class[sum_mask] = 1
        toe_mask = pass_mask2 & (final_class == 0)
        final_class[toe_mask] = 5

    def relative_elevation_priority():
        #---
        # Priority: relative elevation => slope => curvature
        #---
        mask_all = valid_mask & (final_class == 0)

        #--- rel_elev==+ => Summit=1
        re_high = mask_all & (rel_elev_code == 1)
        final_class[re_high] = 1

        #--- rel_elev==- => Toeslope=5
        re_low = mask_all & (rel_elev_code == -1)
        final_class[re_low] = 5

        #--- next step: slope => Backslope if slope > slope_break
        pass_mask = (final_class == 0) & valid_mask
        high_slope = pass_mask & (slope_arr > slope_break)
        final_class[high_slope] = 3

        #--- next step: curvature => Shoulder=2 if +, else Footslope=4
        pass_mask2 = (final_class == 0) & valid_mask
        conv_mask = pass_mask2 & (curvature_code == 1)
        final_class[conv_mask] = 2

        foot_mask = pass_mask2 & (curvature_code == -1) & (final_class == 0)
        final_class[foot_mask] = 4

        #--- anything leftover => Shoulder=2 by default
        leftover = pass_mask2 & (final_class == 0)
        final_class[leftover] = 2

    #---
    # 8) Run the chosen hierarchy
    #---
    print(f"[DEBUG] Selected hierarchy method: {method}")
    if method.lower() == "slgp":
        slope_gradient_priority()
    elif method.lower() == "prcp":
        profile_curvature_priority()
    elif method.lower() == "reep":
        relative_elevation_priority()
    else:
        raise ValueError(f"Invalid method '{method}'. Use 'SlgP', 'PrcP', or 'ReeP'.")

    #---
    # 9) (Optional) Debug: Count how many cells ended up in each class
    #---
    for c_id in [1, 2, 3, 4, 5]:
        count_c = np.sum(final_class == c_id)
        print(f"[DEBUG] Class {c_id} => {count_c} pixels")

    #---
    # 10) Convert final_class => Raster and define projection
    #---
    final_ras = arcpy.NumPyArrayToRaster(
        final_class,
        lower_left,
        cell_width,
        cell_height,
        value_to_nodata=no_data_value
    )
    arcpy.management.DefineProjection(final_ras, spatial_ref)

    out_ras = os.path.join(out_folder, f"Hillslope_{method}.tif")
    final_ras.save(out_ras)
    print(f"[INFO] Unmasked classification saved to: {out_ras}")

    #---
    # 11) Mask out NoData from the original DEM 
    #   => so that outside remains NoData instead of receiving a class.
    #---
    print("[DEBUG] Masking out NoData from the *original* DEM.")
    original_dem = Raster(dem_path)
    classification_masked = SetNull(IsNull(original_dem), final_ras)

    #---
    # 12) Save the masked classification
    #---
    masked_out_path = os.path.join(out_folder, f"Hillslope_{method}_masked.tif")
    classification_masked.save(masked_out_path)

    #---
    # 13) Build RAT on masked version & label classes
    #---
    arcpy.management.BuildRasterAttributeTable(masked_out_path, "Overwrite")
    arcpy.management.AddField(masked_out_path, "ClassName", "TEXT", field_length=25)

    code_labels = {
        1: "Summit",
        2: "Shoulder",
        3: "Backslope",
        4: "Footslope",
        5: "Toeslope"
    }
    with arcpy.da.UpdateCursor(masked_out_path, ["Value", "ClassName"]) as cursor:
        for row in cursor:
            row[1] = code_labels.get(row[0], "NoData")
            cursor.updateRow(row)

    print(f"[INFO] Final classification masked to DEM NoData => {masked_out_path}")
    return masked_out_path


#-----------------------------------------------------------------------
# USAGE: if called as a script
#-----------------------------------------------------------------------
if __name__ == "__main__":
    #--- Input parameters
    dem_input = r"U:/t_StandardProcessGuides/b_DEM_Processes/Test_DEMS/DEM_3m.tif"
    output_dir = r"U:/t_StandardProcessGuides/b_DEM_Processes/Test_DEMS/Output"

    #--- Call: SlgP method, using 'minmax' relative elev approach
    classify_hillslopes(
        dem_path=dem_input,
        out_folder=output_dir,
        method="SlgP",
        slope_break_low=2.4,
        slope_break_high=5.1,
        slope_break=10.0,
        curvature_neg_tol=0.0,
        curvature_pos_tol=0.0,
        rel_elev_break=0.0,
        rel_elev_high=10.0,
        rel_elev_low=-10.0,
        slope_radius=9,
        curvature_radius=63,
        elev_radius=135,
        rel_elev_method="minmax",  # or "mean"
        no_data_value=-9999,       # more common than 0 to avoid confusion with real data=0
        overwrite=True
    )

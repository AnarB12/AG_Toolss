"""
This code was modified from https://github.com/asub-sandwich/relief-analysis-2
and is based on the algorithm described in the following paper:

Miller, B. A., & Schaetzl, R. J. (2015). Digital classification of hillslope position.
Soil Science Society of America Journal, 79(1), 132-145.
"""

# -*- coding: utf-8 -*-
import arcpy
from arcpy import Raster, Parameter, EnvManager
from arcpy.sa import *
from arcpy.ddd import SurfaceParameters
import os
from glob import glob

# ----------------------------------------------------------------------
# HELPER FUNCTION: Determine the correct name/extension for intermediate files.
# If the workspace is a folder, we append ".tif". If it's a file GDB (LocalDatabase
# or RemoteDatabase), we do not add the .tif extension. This way, intermediate
# files can be saved in both folder or gdb without issues.
# ----------------------------------------------------------------------
def get_intermediate_name(base_name, workspace):
    """
    Returns a full path for 'base_name' with .tif extension if workspace
    is a folder, or no extension if it's a GDB.

    :param base_name:   (str) The core file name, e.g. "slope" or "profc".
    :param workspace:   (str) The workspace path; can be a folder or GDB.
    :return:            (str) The appropriate path.
    """
    # If no workspace is provided or if it's not recognized as a GDB, assume a folder -> .tif
    if not workspace:
        return os.path.join(r"C:\tmp", f"{base_name}.tif")

    desc = arcpy.Describe(workspace)
    if hasattr(desc, "workspaceType") and desc.workspaceType == "FileSystem":
        # It's a folder -> use .tif
        return os.path.join(workspace, f"{base_name}.tif")
    elif hasattr(desc, "workspaceType") and desc.workspaceType in ["LocalDatabase", "RemoteDatabase"]:
        # It's a file GDB -> omit .tif
        return os.path.join(workspace, base_name)
    else:
        # Fallback: treat as folder with .tif
        return os.path.join(workspace, f"{base_name}.tif")


class Toolbox:
    def __init__(self) -> None:
        self.label = "Relief Analysis 2.0"
        self.alias = "Relief_Analysis_2"

        # Add all tools here ...
        self.tools = [
            # The order in which tools are listed here is the order
            # in which they appear in ArcGIS.
            # Include your entire list of tools:
            HillslopeAutomatic,
            HillslopeManual,
            Reclassify_2,
            Reclassify_3,
            RelativeElevation,
            ZoneCleanup
        ]


class HillslopeAutomatic:
    def __init__(self) -> None:
        self.label = "Hillslope Position (Automatic)"
        self.description = "Calculates hillslope position from DEM automatically."
    
    def getParameterInfo(self) -> list[Parameter]:
        """
        Define the tool parameters for HillslopeAutomatic.
        """
        # ----------------------------------------------------------------------
        # 1. Input DEM
        # ----------------------------------------------------------------------
        param0 = Parameter(
            displayName="Input DEM",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )

        # ----------------------------------------------------------------------
        # 2. Output Raster
        # ----------------------------------------------------------------------
        param1 = Parameter(
            displayName="Output Raster",
            name="output",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output"
        )

        # ----------------------------------------------------------------------
        # 3. Boolean parameter for saving or deleting intermediate files
        # ----------------------------------------------------------------------
        param2 = Parameter(
            displayName="Save Intermediate Rasters?",
            name="save_intermediates",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        )
        param2.value = False  # Default is False

        # ----------------------------------------------------------------------
        # 4. Workspace for storing intermediate files (folder or GDB)
        #    if save_intermediates == True.
        # ----------------------------------------------------------------------
        param3 = Parameter(
            displayName="Intermediate Output Workspace (folder or GDB)",
            name="intermediate_ws",
            datatype="DEWorkspace",  
            parameterType="Optional", 
            direction="Input"
        )
        # No default—user can leave blank if not saving.

        # ----------------------------------------------------------------------
        # 5. Slope Neighborhood Distance (meters)
        # ----------------------------------------------------------------------
        param4 = Parameter(
            displayName="Slope Neighborhood Distance (meters)",
            name="slope_nhd",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        param4.value = 9  # Default was "9 Meters"

        # ----------------------------------------------------------------------
        # 6. Profile Curvature Neighborhood Distance (meters)
        # ----------------------------------------------------------------------
        param5 = Parameter(
            displayName="Profile Curvature Neighborhood Distance (meters)",
            name="profc_nhd",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        param5.value = 63  # Default was "63 Meters"

        # ----------------------------------------------------------------------
        # 7. Lower Break (Slope Gradient)
        # ----------------------------------------------------------------------
        param6 = Parameter(
            displayName="Lower Break (Slope Gradient)",
            name="lbr",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        param6.value = 2.4  # Default if slope type is PERCENT_RISE

        # ----------------------------------------------------------------------
        # 8. Upper Break (Slope Gradient)
        # ----------------------------------------------------------------------
        param7 = Parameter(
            displayName="Upper Break (Slope Gradient)",
            name="ubr",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        param7.value = 5.1  # Default if slope type is PERCENT_RISE

        # ----------------------------------------------------------------------
        # 9. Relative Elevation Scale
        # ----------------------------------------------------------------------
        param8 = Parameter(
            displayName="Relative Elevation Scale",
            name="rel_elev_scale",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input"
        )
        param8.value = 135  # Default was 135

        # ----------------------------------------------------------------------
        # 10. NEW Parameter: Slope Analysis Type (Percent Rise or Degrees)
        # ----------------------------------------------------------------------
        param9 = Parameter(
            displayName="Slope Analysis Type",
            name="slope_analysis_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        )
        # Allowed options: "PERCENT_RISE" or "DEGREES".
        param9.filter.list = ["PERCENT_RISE", "DEGREE"]
        param9.value = "PERCENT_RISE"  # Default

        return [
            param0, param1, param2, param3,
            param4, param5, param6, param7,
            param8, param9
        ]
    
    def isLicensed(self) -> bool:
        return True
    
    def execute(self, parameters, messages) -> None:
        """
        Main execution of the Hillslope Automatic classification.
        """
        # ----------------------------------------------------------------------
        # 1. Extract parameters from user input (or defaults)
        # ----------------------------------------------------------------------
        dem = Raster(parameters[0].valueAsText)       # DEM
        final_output = parameters[1].valueAsText      # Final DHP output
        save_intermediates = parameters[2].value      # Boolean (True/False)
        user_ws = parameters[3].valueAsText           # Folder or GDB (optional)
        slope_nhd_val = float(parameters[4].value)    # Slope neighborhood distance
        profc_nhd_val = float(parameters[5].value)    # Profile curvature neighborhood
        lbr_user = float(parameters[6].value)         # Lower break (user-supplied)
        ubr_user = float(parameters[7].value)         # Upper break (user-supplied)
        rel_elev_scale_val = int(parameters[8].value) # Relative elevation scale
        slope_type = parameters[9].value              # "PERCENT_RISE" or "DEGREES"

        # ----------------------------------------------------------------------
        # 2. If slope type is DEGREES, adjust the default slope breaks from
        #    2.4 and 5.1 to 1.4 and 2.9 if the user hasn't manually changed them.
        # ----------------------------------------------------------------------
        if slope_type == "DEGREE":
            # If user is still using the default 2.4, switch to 1.4
            if lbr_user == 2.4:
                lbr_user = 1.4
            # If user is still using the default 5.1, switch to 2.9
            if ubr_user == 5.1:
                ubr_user = 2.9

        # ----------------------------------------------------------------------
        # 3. Determine where to store intermediates
        #    If user wants to save intermediates and gave a workspace, use it,
        #    else default to C:\tmp.
        # ----------------------------------------------------------------------
        if save_intermediates and user_ws:
            tmpdir = user_ws  # Store in user-chosen location
        else:
            tmpdir = r"C:\tmp"
            if not os.path.isdir(tmpdir):
                os.mkdir(tmpdir, 0o0777)
            else:
                os.chmod(tmpdir, 0o0777)

        arcpy.env.workspace = tmpdir

        # ----------------------------------------------------------------------
        # 4. Calculate Slope (SurfaceParameters -> SLOPE)
        #    We set output_slope_measurement based on the chosen slope_type.
        #    -> If slope_type == "DEGREES" => "DEGREE"
        #    -> Else => "PERCENT_RISE"
        # ----------------------------------------------------------------------
#        if slope_type == "DEGREES":
#            slope_measurement = "DEGREE"
#        else:
#            slope_measurement = "PERCENT_RISE"

        # Slope neighborhood distance as string
        slope_nhd_str = f"{slope_nhd_val} Meters"

        # Build path for slope intermediate (with or without .tif)
        slope_tif = get_intermediate_name("slope", tmpdir)

        with EnvManager(parallelProcessingFactor="75%"):
            SurfaceParameters(
                in_raster=dem,
                out_raster=slope_tif,
                parameter_type="SLOPE",
                local_surface_type="QUADRATIC",
                neighborhood_distance=slope_nhd_str,
                use_adaptive_neighborhood="FIXED_NEIGHBORHOOD",
                z_unit="Meter",
                output_slope_measurement=slope_type,  # <-- key change
                project_geodesic_azimuths="GEODESIC_AZIMUTHS",
                use_equatorial_aspect="NORTH_POLE_ASPECT",
                in_analysis_mask=None
            )
        
        # ----------------------------------------------------------------------
        # 5. Reclassify Slope Based on lbr_user / ubr_user
        # ----------------------------------------------------------------------
        slope_class_tif = get_intermediate_name("slope_class", tmpdir)
        slope_ras = Raster(slope_tif)
        # Create a new raster: -1 where slope < lower break, 1 where slope > upper break, else 0
        out_slope_class = Con(slope_ras < lbr_user, -1, Con(slope_ras > ubr_user, 1, 0))
        out_slope_class.save(slope_class_tif)

        # ----------------------------------------------------------------------
        # 6. Calculate Profile Curvature & reclassify
        #    We'll keep "PERCENT_RISE" for profile curvature measurement, as in the original code.
        # ----------------------------------------------------------------------
        profc_tif = get_intermediate_name("profc", tmpdir)
        profc_nhd_str = f"{profc_nhd_val} Meters"

        with EnvManager(parallelProcessingFactor="75%"):
            SurfaceParameters(
                in_raster=dem,
                out_raster=profc_tif,
                parameter_type="PROFILE_CURVATURE",
                local_surface_type="QUADRATIC",
                neighborhood_distance=profc_nhd_str,
                use_adaptive_neighborhood="FIXED_NEIGHBORHOOD",
                z_unit="Meter",
#                output_slope_measurement="PERCENT_RISE",
#                project_geodesic_azimuths="GEODESIC_AZIMUTHS",
#                use_equatorial_aspect="NORTH_POLE_ASPECT",
                in_analysis_mask=None
            )

        profc_class_tif = get_intermediate_name("profc_class", tmpdir)
        profc_ras = Raster(profc_tif)
        br = 0  # Hard-coded threshold for curvature classification
        out_profc_class = Con(profc_ras >= br, 1, -1)
        out_profc_class.save(profc_class_tif)

        # ----------------------------------------------------------------------
        # 7. Calculate Relative Elevation & reclassify
        #    (Calls the RelativeElevation tool)
        # ----------------------------------------------------------------------
        # Build path for relative elevation intermediate
        rel = RelativeElevation()
        rel_out = Parameter()
        rel_out.value = get_intermediate_name("rel_elev", tmpdir)

        # Provide the user-defined relative elevation scale
        scale_param = Parameter()
        scale_param.value = rel_elev_scale_val

        # Execute RelativeElevation
        rel.execute([parameters[0], rel_out, scale_param], None)

        relel_class_tif = get_intermediate_name("rel_elev_class", tmpdir)
        relel_ras = Raster(rel_out.value)
        out_relel_class = Con(relel_ras >= 0, 1, -1)
        out_relel_class.save(relel_class_tif)

        # ----------------------------------------------------------------------
        # 8. Run the "HillslopeManual" tool to get final classification
        # ----------------------------------------------------------------------
        dhp = HillslopeManual()
        slope_param = Parameter()
        slope_param.value = slope_class_tif
        profc_param = Parameter()
        profc_param.value = profc_class_tif
        relel_param = Parameter()
        relel_param.value = relel_class_tif
        out_param = Parameter()
        out_param.value = final_output
        dhp.execute([slope_param, profc_param, relel_param, out_param], None)

        return None
    
    def postExecute(self, parameters) -> None:
        """
        Cleanup or keep the intermediate outputs depending on user choice.
        If NOT saving, we remove intermediate rasters from c:\tmp (or the user
        workspace if that was used).
        We do NOT remove the workspace folder or GDB itself.
        """
        save_intermediates = parameters[2].value
        user_ws = parameters[3].valueAsText if parameters[3].value else None
        tmpdir = user_ws if (save_intermediates and user_ws) else r"C:\tmp"

        # If NOT saving intermediates, remove them from the chosen tmpdir.
        if not save_intermediates:
            # We'll look for the known patterns (with or without .tif).
            files_to_delete = [
                "slope*", "slope_class*",
                "profc*", "profc_class*",
                "rel_elev*", "rel_elev_class*"
            ]
            for pattern in files_to_delete:
                for f in glob(os.path.join(tmpdir, pattern)):
                    try:
                        os.remove(f)
                    except:
                        pass

        return


class Reclassify_2:
    def __init__(self) -> None:
        self.label = "Reclassify 2"
        self.description = ("Reclassify a raster into 2 classes by "
                            "specifying a break (Defaults to 0).")

    def getParameterInfo(self) -> list[Parameter]:
        """Define the tool parameters."""
        param0 = Parameter(
            displayName="Input Raster",
            name="in_ras",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        param1 = Parameter(
            displayName="Output Raster",
            name="out_ras",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output"
        )
        param2 = Parameter(
            displayName="Break",
            name="br",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        param2.value = 0

        return [param0, param1, param2]
            
    def isLicensed(self) -> bool:
        return True

    def execute(self, parameters, messages) -> None:
        in_ras = Raster(parameters[0].valueAsText)
        out_ras = parameters[1].valueAsText
        br = float(parameters[2].value)

        out = Con(in_ras >= br, 1, -1)
        out.save(out_ras)

        return None


class Reclassify_3:
    def __init__(self) -> None:
        self.label = "Reclassify 3"
        self.description = ("Reclassify a raster into 3 classes by "
                            "specifying 2 breaks.")
    
    def getParameterInfo(self) -> list[Parameter]:
        param0 = Parameter(
            displayName="Input Raster",
            name="in_ras",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        param1 = Parameter(
            displayName="Output Raster",
            name="out_ras",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output"
        )
        param2 = Parameter(
            displayName="Lower Break",
            name="lower_break",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input"
        )
        param3 = Parameter(
            displayName="Upper Break",
            name="upper_break",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input"
        )

        return [param0, param1, param2, param3]
    
    def isLicensed(self) -> bool:
        return True
    
    def execute(self, parameters, messages) -> None:
        in_ras = Raster(parameters[0].valueAsText)
        out_ras = parameters[1].valueAsText
        lbr = float(parameters[2].value)
        ubr = float(parameters[3].value)

        out = Con(in_ras < lbr, -1, Con(in_ras > ubr, 1, 0))
        out.save(out_ras)

        return None


class HillslopeManual:
    def __init__(self) -> None:
        self.label = "Hillslope Position (Manual)"
        self.description = (
            "Calculate hillslope positions given classified slope, profile "
            "curvature, and relative elevation. Based on the algorithm "
            "by Miller & Schaetzl (2015)."
        )
    
    def getParameterInfo(self) -> list[Parameter]:
        param0 = Parameter(
            displayName="Slope Gradient (3 Classes)",
            name="slope",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        param1 = Parameter(
            displayName="Profile Curvature (2 Classes)",
            name="profc",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        param2 = Parameter(
            displayName="Relative Elevation (2 Classes)",
            name="rel_elev",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        param3 = Parameter(
            displayName="Output Raster",
            name="output",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output"
        )
        return [param0, param1, param2, param3]
    
    def isLicensed(self) -> bool:
        return True
    
    def execute(self, parameters, messages) -> None:
        """
        Execute the manual hillslope classification:
          slope == 1 => Backslope (3)
          slope == 0 => if profc==1 => Shoulder (2), else => Footslope (4)
          slope == -1 => if rel_elev==1 => Summit (1), else => Toeslope (5)
        """
        slope = Raster(parameters[0].valueAsText)
        profc = Raster(parameters[1].valueAsText)
        relel = Raster(parameters[2].valueAsText)
        out_ras = parameters[3].valueAsText

        out = Con(
            slope == 1, 3, 
            Con(
                slope == 0,
                Con(profc == 1, 2, 4),
                Con(relel == 1, 1, 5)
            )
        )
        out.save(out_ras)

        # Build RAT, add "ClassName" field, and populate classification names
        arcpy.management.BuildRasterAttributeTable(out_ras, "Overwrite")
        arcpy.management.AddField(out_ras, "ClassName", "TEXT", field_length=25)

        code_labels = {
            1: "Summit",
            2: "Shoulder",
            3: "Backslope",
            4: "Footslope",
            5: "Toeslope"
        }

        with arcpy.da.UpdateCursor(out_ras, ["Value", "ClassName"]) as cursor:
            for row in cursor:
                row[1] = code_labels.get(row[0], "NoData")
                cursor.updateRow(row)

        return None
        

class RelativeElevation:
    def __init__(self) -> None:
        self.label = "Relative Elevation"
        self.description = "Calculates Relative Elevation of an input DEM."
    
    def isLicensed(self) -> bool:
        return True
    
    def getParameterInfo(self) -> list[Parameter]:
        param0 = Parameter(
            displayName="Input DEM",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        param1 = Parameter(
            displayName="Output Raster",
            name="out_ras",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output"
        )
        param2 = Parameter(
            displayName="Analysis Scale",
            name="scale",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input"
        )
        param2.value = 135

        return [param0, param1, param2]
    
    def execute(self, parameters, messages) -> None:
        """
        Calculate a simple relative elevation model by subtracting the
        difference between min and max from the DEM in a neighborhood window.
        """
        dem = Raster(parameters[0].valueAsText)
        out_ras = parameters[1].valueAsText
        scale = int(parameters[2].value)

        # Create rectangular neighborhood
        # window = NbrRectangle(scale, scale, "MAP")
        window = NbrCircle(scale, "MAP")
        # Compute local max and min
        max_focal = FocalStatistics(dem, window, "MAXIMUM", "")
        min_focal = FocalStatistics(dem, window, "MINIMUM", "")

        # Subtract difference from DEM
        out = dem - ((min_focal + max_focal) - dem)
        out.save(out_ras)

        return None


class ZoneCleanup:
    def __init__(self) -> None:
        self.label = "Zone Cleanup"
        self.description = "Uses a moving window to clean up delineations."
    
    def isLicensed(self) -> bool:
        return True
    
    def getParameterInfo(self) -> list[Parameter]:
        param0 = Parameter(
            displayName="Input DHP Raster",
            name="in_dhp",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        param1 = Parameter(
            displayName="Output DHP Raster",
            name="out_dhp",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output"
        )
        
        return [param0, param1]

    def execute(self, parameters, messages) -> None:
        """
        Expands, majority-filters (twice), re-expands, and boundary cleans
        the input raster. Useful for removing 'salt-and-pepper' noise
        from classification outputs.
        """
        in_dhp = Int(parameters[0].valueAsText)
        out_dhp = parameters[1].valueAsText

        # Expand certain zones
        expanded_dhp = Expand(
            in_raster=in_dhp, 
            number_cells=1, 
            zone_values=[5, 4, 1, 2],
            expand_method="MORPHOLOGICAL"
        )

        # First majority filter
        first_filter = FocalStatistics(
            in_raster=expanded_dhp,
            neighborhood=NbrCircle(5, "CELL"),
            statistics_type="MAJORITY",
            ignore_nodata="DATA"
        )

        # Second majority filter
        second_filter = FocalStatistics(
            in_raster=Int(first_filter),
            neighborhood=NbrCircle(5, "CELL"),
            statistics_type="MAJORITY",
            ignore_nodata="DATA"
        )

        # Expand again
        expanded_filtered = Expand(
            in_raster=second_filter,
            number_cells=2,
            zone_values=[5, 4, 1, 2],
            expand_method="MORPHOLOGICAL"
        )

        # Boundary clean
        out = BoundaryClean(
            in_raster=expanded_filtered,
            sort_type="ASCEND",
            number_of_runs="TWO_WAY"
        )
        out.save(out_dhp)

        return None



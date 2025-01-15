import arcpy
import os
import requests

class Toolbox(object):
    """Python Toolbox for Downloading and Clipping a DEM based on an AOI"""
    def __init__(self):
        self.label = "Download and Clip DEM"
        self.alias = "DEM_Tool"
        self.tools = [DownloadAndClipDEM]


class DownloadAndClipDEM(object):
    def __init__(self):
        self.label = "Download and Clip DEM"
        self.description = (
            "Takes a field boundary (AOI), optionally buffers it, projects it to WGS84, "
            "uses its extent to download a DEM, and clips the DEM to that AOI. "
            "The user provides the full output raster path (in a folder or GDB). "
            "Intermediate data is stored in a scratch folder and scratch GDB, which can be cleaned up at the end."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        # AOI feature class
        param_aoi = arcpy.Parameter(
            displayName="AOI Feature Class",
            name="aoi",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input"
        )
        # The input feature class representing the field boundary or AOI.

        # Apply Buffer checkbox
        param_do_buffer = arcpy.Parameter(
            displayName="Apply Buffer",
            name="do_buffer",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input"
        )
        param_do_buffer.value = False
        # Check to apply a buffer to the AOI before processing.

        # Buffer distance
        param_buffer_distance = arcpy.Parameter(
            displayName="Buffer Distance (meters)",
            name="buffer_distance",
            datatype="Double",
            parameterType="Optional",
            direction="Input",
            enabled=False
        )
        # Buffer distance in meters to apply to the AOI.

        # Output Raster full path
        param_output_raster = arcpy.Parameter(
            displayName="Output Raster (Full Path)",
            name="output_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output"
        )
        # Full path to the output raster. If pointing to a folder, use a .tif extension.
        # If pointing to a GDB, just provide a raster name without extension, e.g. C:\\path\\my.gdb\\rastername.

        # Scratch Location (Optional)
        param_scratch_location = arcpy.Parameter(
            displayName="Scratch Location (Optional)",
            name="scratch_location",
            datatype="Folder",
            parameterType="Optional",
            direction="Input"
        )
        # A folder where 'scratch' folder and 'scratch.gdb' will be created.
        # If not provided, they will be created in the same directory as the output raster or
        # in the parent directory if output is inside a GDB.

        # Clean Intermediate
        param_clean_intermediate = arcpy.Parameter(
            displayName="Clean Intermediate Layers",
            name="clean_intermediate",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input"
        )
        param_clean_intermediate.value = False
        # Check to delete intermediate layers after processing.

        param_buffer_distance.parameterDependencies = [param_do_buffer.name]

        return [
            param_aoi,
            param_do_buffer,
            param_buffer_distance,
            param_output_raster,
            param_scratch_location,
            param_clean_intermediate
        ]

    def updateParameters(self, parameters):
        # Enable buffer distance if do_buffer is True
        do_buffer = parameters[1].value
        parameters[2].enabled = bool(do_buffer)
        return

    def execute(self, parameters, messages):
        # Retrieve parameters
        aoi = parameters[0].valueAsText
        do_buffer = parameters[1].value
        buffer_distance = parameters[2].value
        output_raster = parameters[3].valueAsText
        scratch_location = parameters[4].valueAsText
        clean_intermediate = parameters[5].value

        # Validate the output path and determine base directory
        # If output is inside a .gdb (e.g. C:\path\my.gdb\rastername),
        # base_dir = C:\path\my.gdb. For scratch, we should use C:\path as scratch location if not provided.
        # If output is a .tif file (e.g. C:\path\myraster.tif),
        # base_dir = C:\path.
        
        # Check if output_raster points to a GDB dataset:
        # A GDB path typically ends in .gdb at some point:
        # Example: C:\path\my.gdb\rastername
        # We can detect if it's a GDB by searching for ".gdb" in the path.
        
        lower_out = output_raster.lower()
        if ".gdb" in lower_out:
            # Split at ".gdb"
            gdb_index = lower_out.rfind(".gdb")
            # gdb_index gives the position of ".gdb" start
            # base_dir includes the ".gdb"
            base_dir = output_raster[:gdb_index+4]
            # Ensure that the part after ".gdb" is actually a raster name (i.e. not empty)
            if len(output_raster) <= gdb_index + 4:
                arcpy.AddError("Output raster inside a GDB requires a raster name after the GDB path.")
                return
        else:
            # Assume it's a file-based raster (like a .tif)
            # or a folder path if user made a mistake
            base_dir = os.path.dirname(output_raster)
            # If no extension, and not inside gdb, must provide a .tif or a known raster type
        
        if not scratch_location:
            # If output is inside a GDB, use the parent folder of the GDB for scratch
            if ".gdb" in lower_out:
                parent_folder = os.path.dirname(base_dir)
                scratch_location = parent_folder
            else:
                # If it's a .tif or similar, just use the same folder as output
                scratch_location = base_dir

        # Validate that scratch_location is a folder (not a GDB)
        if scratch_location.lower().endswith(".gdb"):
            arcpy.AddError("Scratch location cannot be a GDB. Please provide a folder.")
            return

 
        scratch_folder = os.path.join(scratch_location, "scratch")
        scratch_gdb = os.path.join(scratch_location, "scratch.gdb")


        if not os.path.exists(scratch_folder):
            os.makedirs(scratch_folder)


        if not arcpy.Exists(scratch_gdb):
            arcpy.management.CreateFileGDB(scratch_location, "scratch.gdb")

        arcpy.AddMessage(f"Scratch folder: {scratch_folder}")
        arcpy.AddMessage(f"Scratch GDB: {scratch_gdb}")

        # Intermediate outputs
        projected_aoi = os.path.join(scratch_gdb, "AOI_Projected")
        buffered_aoi = os.path.join(scratch_gdb, "AOI_Buffered")
        downloaded_dem = os.path.join(scratch_folder, "dem.tif")

        # Project AOI to WGS84
        arcpy.AddMessage("Projecting AOI to WGS84...")
        arcpy.management.Project(
            in_dataset=aoi,
            out_dataset=projected_aoi,
            out_coor_system=arcpy.SpatialReference(4326)
        )

        final_aoi_for_extent = projected_aoi

        # If buffering is requested:
        if do_buffer and buffer_distance and buffer_distance > 0:
            arcpy.AddMessage(f"Buffering AOI by {buffer_distance} meters...")
            arcpy.analysis.Buffer(
                in_features=projected_aoi,
                out_feature_class=buffered_aoi,
                buffer_distance_or_field=f"{buffer_distance} Meters",
                dissolve_option="ALL"
            )
            final_aoi_for_extent = buffered_aoi

        # Describe final AOI to get bounding box
        desc = arcpy.Describe(final_aoi_for_extent)
        extent = desc.extent
        xmin, ymin, xmax, ymax = extent.XMin, extent.YMin, extent.XMax, extent.YMax
        arcpy.AddMessage(f"Bounding Box (lat/lon): xmin={xmin}, ymin={ymin}, xmax={xmax}, ymax={ymax}")

        # Download DEM
        arcpy.AddMessage("Downloading DEM from USGS...")
        api_key = "demoapikeyot2022"  # API key
        dem_download_url = (
            f"https://portal.opentopography.org/API/usgsdem?"
            f"datasetName=USGS1m&south={ymin}&north={ymax}&west={xmin}&east={xmax}"
            f"&outputFormat=GTiff&API_Key={api_key}"
        )

        response = requests.get(dem_download_url, stream=True)
        if response.status_code == 200:
            with open(downloaded_dem, 'wb') as dem_file:
                for chunk in response.iter_content(chunk_size=8192):
                    dem_file.write(chunk)
            arcpy.AddMessage(f"DEM downloaded successfully: {downloaded_dem}")
        else:
            arcpy.AddError(f"Error downloading DEM: {response.status_code}")
            return

        # Clip the DEM to AOI
        arcpy.AddMessage("Clipping DEM to AOI...")
        arcpy.management.Clip(
            in_raster=downloaded_dem,
            rectangle=f"{xmin} {ymin} {xmax} {ymax}",
            out_raster=output_raster,
            in_template_dataset=final_aoi_for_extent,
            clipping_geometry="ClippingGeometry",
            maintain_clipping_extent="NO_MAINTAIN_EXTENT"
        )
        arcpy.AddMessage(f"Clipped DEM saved at: {output_raster}")

        # Clean intermediate layers if requested
        if clean_intermediate:
            arcpy.AddMessage("Cleaning intermediate layers...")
            # Delete intermediate feature classes
            to_delete = [projected_aoi, buffered_aoi]
            for fc in to_delete:
                if arcpy.Exists(fc):
                    arcpy.management.Delete(fc)

            # Delete downloaded DEM
            if os.path.exists(downloaded_dem):
                try:
                    os.remove(downloaded_dem)
                except Exception as e:
                    arcpy.AddWarning(f"Could not remove DEM file: {e}")

        arcpy.AddMessage("Processing complete.")



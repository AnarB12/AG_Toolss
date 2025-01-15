import arcpy
import csv
import os
import shutil

# -------------------------------------------------------------------------
# 1) Function for the first tool: Join CSV to Feature Class
# -------------------------------------------------------------------------

def clean_temp_dir(temp_dir):
    """Deletes the specified temporary directory and its contents."""
    try:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            arcpy.AddMessage(f"Deleted temporary directory: {temp_dir}")
    except Exception as e:
        arcpy.AddWarning(f"Failed to delete temporary directory {temp_dir}: {str(e)}")


def join_csv_to_feature_class(
    shapefile,
    csv_file,
    output_fc,
    csv_time_field,
    longitude_field,
    latitude_field,
    order_field=None,
    coord_system=None,
    search_tolerance=0.0
):
    """
    Joins a CSV field (e.g., "GPS Time") to an existing point shapefile by:
      1) Creating a temp directory (C:\\temp).
      2) Copying the user CSV into that directory.
      3) Using Python's csv module to remove unwanted columns from the CSV,
         keeping only longitude, latitude, and time fields.
      4) Creating a scratch GDB in C:\\temp if needed, then XYTableToPoint
         on the pruned CSV.
      5) Spatially joining the shapefile with the new point FC (keeping shapefile
         fields + time).
      6) Optionally sorting by that time field and assigning an order in 'order_field'.

    Parameters
    ----------
    shapefile : str
        Path to the existing point shapefile/feature class (target).
    csv_file : str
        Path to the original CSV file.
    output_fc : str
        Path to the output feature class (.shp or in a GDB).
    csv_time_field : str
        Name of the time field in the CSV, e.g., "GPS Time".
    longitude_field : str
        Name of the longitude field in the CSV.
    latitude_field : str
        Name of the latitude field in the CSV.
    order_field : str, optional
        If provided, creates/updates a field to store a 1..N order after sorting by the time field.
    coord_system : arcpy.SpatialReference, optional
        Spatial reference for the CSV points. Defaults to WGS84 if None.
    search_tolerance : float
        If > 0, uses 'INTERSECT' for Spatial Join with that distance (meters).
        If 0, uses 'INTERSECT' (exact overlap).
    """

    arcpy.AddMessage("Starting CSV-to-feature join by copying CSV to C:\\temp and removing unwanted columns...")

    # ---------------------------------------------------------------------
    # 1) Ensure C:\temp exists
    # ---------------------------------------------------------------------
    temp_dir = r"C:\tempDir_Yield"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
        arcpy.AddMessage(f"Created temp directory: {temp_dir}")

    # Paths for the CSV copy and the reduced CSV
    csv_copy = os.path.join(temp_dir, "temp_copy.csv")
    csv_reduced = os.path.join(temp_dir, "temp_reduced.csv")

    # ---------------------------------------------------------------------
    # 2) Copy the user CSV into C:\temp
    # ---------------------------------------------------------------------
    shutil.copy2(csv_file, csv_copy)
    arcpy.AddMessage(f"Copied CSV to: {csv_copy}")

    # ---------------------------------------------------------------------
    # 3) Use Python's csv module to remove unwanted columns, saving a reduced CSV
    #    that keeps only [longitude_field, latitude_field, csv_time_field].
    # ---------------------------------------------------------------------
    keep_fields_lower = {longitude_field.lower(), latitude_field.lower(), csv_time_field.lower()}

    with open(csv_copy, "r", newline="", encoding="utf-8") as fin, \
         open(csv_reduced, "w", newline="", encoding="utf-8") as fout:
        
        reader = csv.DictReader(fin)
        # Build a list of fields we'll keep, preserving the same case as in the CSV's header
        keep_fields_actual = [f for f in reader.fieldnames if f.lower() in keep_fields_lower]

        writer = csv.DictWriter(fout, fieldnames=keep_fields_actual)
        writer.writeheader()

        for row in reader:
            # Build a reduced dict with only the keep_fields
            reduced_row = {k: row[k] for k in keep_fields_actual}
            writer.writerow(reduced_row)

    arcpy.AddMessage(f"Reduced CSV saved as: {csv_reduced}")
    arcpy.AddMessage(f"Kept only fields: {keep_fields_actual}")

    # ---------------------------------------------------------------------
    # 4) Create or reuse a scratch GDB in C:\temp, then XYTableToPoint
    # ---------------------------------------------------------------------
    # scratch_gdb = os.path.join(temp_dir, "scratch.gdb")
    # if not arcpy.Exists(scratch_gdb):
    #     arcpy.management.CreateFileGDB(temp_dir, "scratch.gdb")
    #     arcpy.AddMessage(f"Created scratch GDB: {scratch_gdb}")
    # else:
    #     arcpy.AddMessage(f"Using existing scratch GDB: {scratch_gdb}")

    if coord_system is None:
        coord_system = arcpy.SpatialReference(4326)  # WGS84 as default

    csv_points_fc = arcpy.CreateUniqueName("CSVPoints", arcpy.env.scratchGDB)
    # csv_points_fc = os.path.join(scratch_gdb, "ReducedCSVPoints")
    arcpy.management.XYTableToPoint(
        in_table=csv_reduced,
        out_feature_class=csv_points_fc,
        x_field=longitude_field,
        y_field=latitude_field,
        coordinate_system=coord_system
    )
    arcpy.AddMessage(f"Created point feature class from reduced CSV: {csv_points_fc}")

    # ---------------------------------------------------------------------
    # 5) Identify the final time field name in the CSV point FC
    #    (ArcGIS may have changed "GPS Time" -> "GPS_Time", etc.)
    # ---------------------------------------------------------------------
    fc_fields = arcpy.ListFields(csv_points_fc)
    final_time_field = None

    # We'll try to find a field that matches csv_time_field ignoring underscores/spaces
    # or if there's only one leftover user field besides the OID/Shape, we take that.
    for fld in fc_fields:
        # Check aliasName or a special property to see if it was derived from csv_time_field
        # Often, XYTableToPoint might keep the alias or do a safe name, e.g. "GPS_Time".
        # We'll do a simple check: if either the name or the alias (case-insensitive)
        # is the same as csv_time_field, or if ArcGIS replaced space with underscore, etc.
        if fld.type in ("OID", "Geometry"):
            continue
        # We'll do a simple "normalize" of the name:
        normalized_name = fld.name.replace("_", "").replace(" ", "").lower()
        normalized_target = csv_time_field.replace("_", "").replace(" ", "").lower()
        if normalized_name == normalized_target:
            final_time_field = fld.name
            break
        # Or if the alias is the same as the original csv_time_field:
        elif fld.aliasName.lower() == csv_time_field.lower():
            final_time_field = fld.name
            break


    if not final_time_field:
        user_fields = [f.name for f in fc_fields if f.type not in ("OID", "Geometry")]
        if len(user_fields) == 1:
            final_time_field = user_fields[0]
            arcpy.AddWarning(
                f"Couldn't match '{csv_time_field}' exactly, but found a single field: '{final_time_field}'. "
                "Using that as the time field."
            )
        else:
            err_msg = (
                f"Could not find a matching time field for '{csv_time_field}' in {csv_points_fc}. "
                "Please check your CSV field names."
            )
            arcpy.AddError(err_msg)
            raise ValueError(err_msg)

    arcpy.AddMessage(f"Final time field in point FC is: {final_time_field}")

    # ---------------------------------------------------------------------
    # 6) Spatial Join with shapefile (target), keeping shapefile fields + final_time_field
    # ---------------------------------------------------------------------
    # field_mappings = arcpy.FieldMappings()
    # field_mappings.addTable(shapefile)       # Target
    # field_mappings.addTable(csv_points_fc)   # Join

    # # Remove all CSV fields except final_time_field
    # for i in reversed(range(field_mappings.fieldCount)):
    #     f_map = field_mappings.getFieldMap(i)
    #     src_name = f_map.getInputFieldName(0)
    #     if f_map.findInputFieldIndex(csv_points_fc) != -1:
    #         # It's from CSV
    #         if src_name.lower() != final_time_field.lower():
    #             field_mappings.removeFieldMap(i)

    if search_tolerance > 0:
        match_option = "INTERSECT"
        distance_str = f"{search_tolerance} Meters"
    else:
        match_option = "INTERSECT"
        distance_str = None

    arcpy.analysis.SpatialJoin(
        target_features=shapefile,
        join_features=csv_points_fc,
        out_feature_class=output_fc,
        join_operation="JOIN_ONE_TO_ONE",
        join_type="KEEP_ALL",
        match_option=match_option,
        search_radius=distance_str,
        # field_mapping=field_mappings
    )
    arcpy.AddMessage(f"Spatial Join complete. Output: {output_fc}")

    fields_to_remove = ["Join_Count", 
                        "TARGET_FID", 
                        "Field"]
    
    # List all fields in the feature class
    fields = [field.name for field in arcpy.ListFields(output_fc)]

    # Filter fields_to_remove to include only those that exist in the feature class
    existing_fields_to_remove = [field for field in fields_to_remove if field in fields]

    # Delete only the fields that exist
    if existing_fields_to_remove:
        arcpy.management.DeleteField(output_fc, existing_fields_to_remove)
        print(f"Deleted fields: {existing_fields_to_remove}")
    else:
        print("No matching fields to delete.")
    
    # ---------------------------------------------------------------------
    # 7) (Optional) Sort by final_time_field and assign order
    # ---------------------------------------------------------------------
    if order_field:
        existing_fields = [f.name for f in arcpy.ListFields(output_fc)]
        if order_field not in existing_fields:
            arcpy.management.AddField(output_fc, order_field, "LONG")

        features = []
        with arcpy.da.SearchCursor(output_fc, ["OID@", final_time_field]) as s_cursor:
            for row in s_cursor:
                oid_val, time_val = row
                if time_val not in (None, ""):
                    features.append((oid_val, time_val))

        # Sort by time_val
        features.sort(key=lambda x: x[1])
        oid_to_rank = {oid: i + 1 for i, (oid, _) in enumerate(features)}

        with arcpy.da.UpdateCursor(output_fc, ["OID@", order_field]) as u_cursor:
            for row in u_cursor:
                oid = row[0]
                if oid in oid_to_rank:
                    row[1] = oid_to_rank[oid]
                    u_cursor.updateRow(row)

        arcpy.AddMessage(f"Assigned order (1..N) in field: {order_field}")

    arcpy.AddMessage("Join completed successfully (CSV reduced in C:\\temp).")

    clean_temp_dir(temp_dir)
    arcpy.management.Delete(csv_points_fc)

    
# -------------------------------------------------------------------------
# 2) Function for the second tool: Create Buffer Pass
# -------------------------------------------------------------------------
def create_buffer_pass(in_shapefile, 
                        gdb_path, 
                        working_fd, 
                        yield_fd, 
                        outFC_name, 
                        sortFC=False, 
                        sort_field=None, 
                        sortedField_name=None, 
                        remove_duplicates=False, 
                        pass_threshold=10, 
                        min_pass_length=50, 
                        buffer_distance=4.5):
    """
    Processes a strip trial workflow by copying an input shapefile into a geodatabase,
    optionally sorting features, removing duplicates, buffering lines, and more.

    :param in_shapefile: Path to the input shapefile.
    :param gdb_path: Path to the file geodatabase.
    :param working_fd: Name of the working feature dataset.
    :param yield_fd: Name of the yield feature dataset.
    :param outFC_name: Name of the output feature class to create in the yield FD.
    :param sortFC: Boolean indicating whether to sort features.
    :param sort_field: Field name used for sorting if sortFC is True.
    :param sortedField_name: New field name to store sorted order if sortFC is True.
    :param remove_duplicates: Boolean indicating whether to remove duplicate geometry.
    :param swath_width: Swath width used in point-to-line processes.
    :param min_pass_length: Minimum length for a pass line feature.
    :param buffer_distance: Buffer distance for line buffering.
    """
    try:
        arcpy.AddMessage("Starting the strip trial processing...")

        # Set the coordinate system for the new feature datasets
        coordinate_system = arcpy.SpatialReference("NAD 1983 UTM Zone 15N")

        # Create the geodatabase if it doesn't exist
        if not arcpy.Exists(gdb_path):
            arcpy.management.CreateFileGDB(
                out_folder_path=os.path.dirname(gdb_path),
                out_name=os.path.basename(gdb_path)
            )
            arcpy.AddMessage(f"Created Geodatabase: {gdb_path}")

        # Create the working feature dataset if it doesn't exist
        working_fd_path = os.path.join(gdb_path, working_fd)
        if not arcpy.Exists(working_fd_path):
            arcpy.management.CreateFeatureDataset(
                out_dataset_path=gdb_path,
                out_name=working_fd,
                spatial_reference=coordinate_system
            )
            arcpy.AddMessage(f"Created Working Feature Dataset: {working_fd_path}")

        # Create the yield feature dataset if it doesn't exist
        yield_fd_path = os.path.join(gdb_path, yield_fd)
        if not arcpy.Exists(yield_fd_path):
            arcpy.management.CreateFeatureDataset(
                out_dataset_path=gdb_path,
                out_name=yield_fd,
                spatial_reference=coordinate_system
            )
            arcpy.AddMessage(f"Created Yield Feature Dataset: {yield_fd_path}")

        # Set environment settings
        arcpy.env.workspace = gdb_path
        arcpy.env.overwriteOutput = True

        # Copy input shapefile to the yield feature dataset
        output_feature_class = os.path.join(yield_fd_path, outFC_name)
        arcpy.management.CopyFeatures(in_features=in_shapefile, out_feature_class=output_feature_class)
        arcpy.AddMessage(f"Copied input shapefile to: {output_feature_class}")

        # Remove duplicates if specified
        if remove_duplicates:
            arcpy.management.DeleteIdentical(output_feature_class, ["Shape"])
            arcpy.AddMessage("Removed duplicate features based on geometry.")

        # Sort and create IDs if specified
        if sortFC and sort_field:
            # If a new sortedField_name is provided and doesn't exist, add it
            existing_fields = [f.name for f in arcpy.ListFields(output_feature_class)]
            if sortedField_name and sortedField_name not in existing_fields:
                arcpy.management.AddField(output_feature_class, sortedField_name, "LONG")

            # Sort features by sort_field
            sorted_object_ids = [
                row[0] for row in arcpy.da.SearchCursor(
                    output_feature_class,
                    ["OID@", sort_field],
                    sql_clause=(None, f"ORDER BY {sort_field} ASC"))
            ]

            # Create a dictionary: OID -> new ID
            oid_to_id = {oid: i + 1 for i, oid in enumerate(sorted_object_ids)}

            # Update the sortedField_name with the new ID
            if sortedField_name:
                with arcpy.da.UpdateCursor(output_feature_class, ["OID@", sortedField_name]) as cursor:
                    for row in cursor:
                        row[1] = oid_to_id[row[0]]
                        cursor.updateRow(row)
                arcpy.AddMessage(f"Sorted features and assigned IDs in field: {sortedField_name}")

        # ---------------------------------------------------------------------
        # The following lines create intermediate outputs and transform the data
        # ---------------------------------------------------------------------

        # Define paths for intermediate layers
        base_name = outFC_name
        sms_min = output_feature_class
        sms_line_base = os.path.join(gdb_path, working_fd, f"{base_name}_line_base")
        sms_split_line = os.path.join(gdb_path, working_fd, f"{base_name}_split_line")
        sms_select_line = os.path.join(gdb_path, working_fd, f"{base_name}_select_line")
        sms_pass_line = os.path.join(gdb_path, working_fd, f"{base_name}_pass_line")
        sms_pass_line_F = os.path.join(gdb_path, working_fd, f"{base_name}_pass_line_F")
        sms_min_pass = os.path.join(gdb_path, yield_fd, f"{base_name}_pass")
        sms_pass_line_buffer = os.path.join(gdb_path, working_fd, f"{base_name}_line_buffer")
        sms_pass_line_union = os.path.join(gdb_path, working_fd, f"{base_name}_line_union")
        sms_line_identical = os.path.join(gdb_path, f"{base_name}_Line_Identical")
        line_delete = os.path.join(gdb_path, working_fd, f"{base_name}_Line_Delete")
        line_dissolve = os.path.join(gdb_path, working_fd, f"{base_name}_Line_Dissolve")

        # Convert points to line
        arcpy.management.PointsToLine(
            Input_Features=sms_min,
            Output_Feature_Class=sms_line_base,
            Line_Field=None,
            Sort_Field=sort_field,
            Close_Line="NO_CLOSE",
            Line_Construction_Method="CONTINUOUS",
            Attribute_Source="NONE",
            Transfer_Fields=None
        )
        arcpy.AddMessage("Converted points to line.")

        # Split lines at intersections
        arcpy.management.SplitLine(
            in_features=sms_line_base,
            out_feature_class=sms_split_line
        )
        arcpy.AddMessage("Split lines at intersections.")

        # Select lines shorter than half the swath width
        # half_width = swath_width / 2
        arcpy.management.MakeFeatureLayer(
            in_features=sms_split_line,
            out_layer="split_line_layer",
            where_clause=f"Shape_Length < {pass_threshold}"
        )
        arcpy.management.CopyFeatures(
            in_features="split_line_layer",
            out_feature_class=sms_select_line
        )
        arcpy.AddMessage(f"Selected lines shorter than {pass_threshold} meters.")

        # Unsplit selected lines
        arcpy.management.UnsplitLine(
            in_features=sms_select_line,
            out_feature_class=sms_pass_line,
            dissolve_field=None,
            statistics_fields=None,
            concatenation_separator=""
        )
        arcpy.AddMessage("Unsplit selected lines.")

        # Select lines longer than min_pass_length
        arcpy.management.MakeFeatureLayer(
            in_features=sms_pass_line,
            out_layer="pass_line_layer",
            where_clause=f"Shape_Length > {min_pass_length}"
        )
        arcpy.management.CopyFeatures(
            in_features="pass_line_layer",
            out_feature_class=sms_pass_line_F
        )
        arcpy.AddMessage(f"Selected lines longer than {min_pass_length} meters.")

        # Add Pass_ID field and populate with OID
        arcpy.management.AddField(
            in_table=sms_pass_line_F,
            field_name="Pass_ID",
            field_type="SHORT",
            field_is_nullable="NULLABLE"
        )
        arcpy.management.CalculateField(
            in_table=sms_pass_line_F,
            field="Pass_ID",
            expression="!OBJECTID!",
            expression_type="PYTHON3"
        )
        arcpy.AddMessage("Added Pass_ID field and calculated values.")

        # Perform spatial join to attach pass info to points
        arcpy.analysis.SpatialJoin(
            target_features=sms_min,
            join_features=sms_pass_line_F,
            out_feature_class=sms_min_pass,
            join_operation="JOIN_ONE_TO_ONE",
            join_type="KEEP_ALL",
            match_option="WITHIN_A_DISTANCE",
            search_radius="1 Meters"
        )
        arcpy.AddMessage("Performed spatial join.")

        # Use an update cursor to iterate over the rows and delete those with Pass_ID as null
        with arcpy.da.UpdateCursor(sms_min_pass, ["Pass_ID"]) as cursor:
            for row in cursor:
                if row[0] is None:  # Check if Pass_ID is null
                    cursor.deleteRow()

        arcpy.AddMessage("Rows with null Pass_ID have been deleted.")

        # Buffer lines
        # arcpy.analysis.PairwiseBuffer
        arcpy.analysis.Buffer(
            in_features=sms_pass_line_F,
            out_feature_class=sms_pass_line_buffer,
            buffer_distance_or_field=f"{buffer_distance} Meters",
            line_side="FULL",
            line_end_type="FLAT",
            dissolve_option="NONE",
            dissolve_field=None,
            method="GEODESIC"
        )
        arcpy.AddMessage("Buffered lines.")

        # Union buffered lines
        arcpy.analysis.Union(
            in_features=[sms_pass_line_buffer],
            out_feature_class=sms_pass_line_union,
            join_attributes="ALL",
            gaps="GAPS"
        )
        arcpy.AddMessage("Performed union on buffered lines.")

        # Find identical features
        arcpy.management.FindIdentical(
            in_dataset=sms_pass_line_union,
            out_dataset=sms_line_identical,
            fields=["Shape_Length", "Shape_Area"],
            output_record_option="ONLY_DUPLICATES"
        )
        arcpy.AddMessage("Found identical features.")

        # Join identical features back
        addJoin_result = arcpy.management.AddJoin(
            in_layer_or_view=sms_pass_line_union,
            in_field="OBJECTID",
            join_table=sms_line_identical,
            join_field="IN_FID",
            join_type="KEEP_ALL"
        )
        arcpy.management.CopyFeatures(
            in_features=addJoin_result,
            out_feature_class=line_delete
        )
        arcpy.AddMessage("Joined identical features back.")

        # Delete identical features
        arcpy.management.DeleteIdentical(
            in_dataset=line_delete,
            fields=["FEAT_SEQ", "Shape_Area"]
        )
        arcpy.AddMessage("Deleted identical features.")

        # Dissolve line features
        # Construct the field name dynamically from the joined output
        # This uses the naming convention from AddJoin (tableName_fieldName).
        dissolveField = f"{os.path.basename(sms_pass_line_union)}_FID_{os.path.basename(sms_pass_line_buffer)}"
        arcpy.management.Dissolve(
            in_features=line_delete,
            out_feature_class=line_dissolve,
            dissolve_field=[dissolveField],
            multi_part="SINGLE_PART",
            unsplit_lines="DISSOLVE_LINES"
        )
        arcpy.AddMessage("Dissolved line features.")

        # Clean up temporary layers
        arcpy.management.Delete(["split_line_layer", "pass_line_layer", addJoin_result])
        arcpy.AddMessage("Cleaned up temporary layers.")

        arcpy.AddMessage("Process finished successfully!")

    except arcpy.ExecuteError as e:
        arcpy.AddError(f"ExecuteError: {e}")
    except Exception as e:
        arcpy.AddError(f"Unexpected error: {str(e)}")


# -------------------------------------------------------------------------
# 3) Function for the third tool: Clean Points by Pass ID
# -------------------------------------------------------------------------
def clean_points_by_passid(
    input_points,
    input_polygons,
    output_points,
    search_distance,
    pass_id_field
):
    """
    Copies the input points (so the original remains unmodified), then:
      - Iterates through polygons
      - Selects points by location (search distance)
      - Finds min PASS_ID
      - Deletes points with PASS_ID above that min
      - Appends survivors to a new output FC

    input_points (str): Path to original points (not modified).
    input_polygons (str): Path to polygons.
    output_points (str): Path to new output feature class (created).
    search_distance (str): e.g. "5 Meters".
    pass_id_field (str): Name of PASS_ID field in points.
    """
    # Add a message to show tool started
    arcpy.AddMessage("Starting 'clean_points_by_passid'...")

    # We'll wrap the entire process in try/except to handle any errors gracefully
    try:
        # ---------------------------------------------------------------------
        # 1) Copy the input points to an in-memory dataset, so we don't modify original
        # ---------------------------------------------------------------------
        # temp_points = r"in_memory\temp_points"
        temp_points = arcpy.CreateUniqueName("temp_points", arcpy.env.scratchGDB)
        if arcpy.Exists(temp_points):
            arcpy.management.Delete(temp_points)

        arcpy.AddMessage("Copying input points to temporary in-memory dataset...")
        arcpy.management.CopyFeatures(input_points, temp_points)

        # Make a layer from this temporary copy
        arcpy.management.MakeFeatureLayer(temp_points, "points_lyr")
        arcpy.AddMessage("Created feature layer from the temporary points copy.")


        # ---------------------------------------------------------------------
        # 2) Remove any points that 'BOUNDARY_TOUCH' polygons
        # ---------------------------------------------------------------------
        # We'll select by location with 'BOUNDARY_TOUCHES' and delete them
        arcpy.management.SelectLayerByLocation(
            in_layer="points_lyr",
            overlap_type="BOUNDARY_TOUCHES",
            select_features=input_polygons,
            search_distance=None,
            selection_type="NEW_SELECTION"
        )
        boundary_count = int(arcpy.management.GetCount("points_lyr").getOutput(0))
        if boundary_count > 0:
            arcpy.AddMessage(f"Found {boundary_count} point(s) on polygon boundaries. Deleting them...")
            with arcpy.da.UpdateCursor("points_lyr", ["OID@"]) as ucur:
                for row in ucur:
                    ucur.deleteRow()
        else:
            arcpy.AddMessage("No points found on polygon boundaries.")

        # Clear selection
        arcpy.management.SelectLayerByAttribute("points_lyr", "CLEAR_SELECTION")


        # ---------------------------------------------------------------------
        # 3) Create an empty output feature class (schema only)
        # ---------------------------------------------------------------------
        if arcpy.Exists(output_points):
            arcpy.management.Delete(output_points)
            arcpy.AddMessage(f"Deleted existing {output_points} to recreate it.")

        # Copy schema from temp_points, then delete rows to make it empty
        arcpy.management.CopyFeatures(temp_points, output_points)
        arcpy.management.DeleteRows(output_points)
        arcpy.AddMessage(f"Created empty output feature class: {output_points}")


        # ---------------------------------------------------------------------
        # 4) Iterate through each polygon, select points, find min PASS_ID, delete > min
        # ---------------------------------------------------------------------
        with arcpy.da.SearchCursor(input_polygons, ["OID@"]) as poly_cursor:
            for poly_row in poly_cursor:
                polygon_oid = poly_row[0]

                # 4a) Make a single-polygon layer for the current polygon
                where_clause = f"OBJECTID = {polygon_oid}"
                arcpy.management.MakeFeatureLayer(
                    input_polygons,
                    "single_poly_lyr",
                    where_clause=where_clause
                )

                # 4b) Select points that intersect (within search distance if > 0)
                arcpy.management.SelectLayerByLocation(
                    in_layer="points_lyr",
                    overlap_type="INTERSECT",
                    select_features="single_poly_lyr",
                    search_distance=search_distance,  # e.g. "5 Meters"
                    selection_type="NEW_SELECTION"
                )

                # Count selected
                count_sel = int(arcpy.management.GetCount("points_lyr").getOutput(0))
                arcpy.AddMessage(
                    f"Polygon OID={polygon_oid}, selected {count_sel} point(s)."
                )

                if count_sel > 0:
                    # 4c) Find the minimum PASS_ID among these selected points
                    pass_ids = []
                    with arcpy.da.SearchCursor("points_lyr", [pass_id_field]) as s_cursor:
                        for s_row in s_cursor:
                            pass_ids.append(s_row[0])

                    if pass_ids:
                        min_passid = min(pass_ids)
                        # 4d) Delete any selected points with PASS_ID > min_passid
                        with arcpy.da.UpdateCursor("points_lyr", [pass_id_field]) as u_cursor:
                            for u_row in u_cursor:
                                if u_row[0] > min_passid:
                                    u_cursor.deleteRow()

                        arcpy.AddMessage(
                            f"Deleted points in selection with {pass_id_field} > {min_passid}."
                        )

                    # 4e) Append survivors to the output FC
                    arcpy.management.Append("points_lyr", output_points, "NO_TEST")
                    arcpy.AddMessage("Appended cleaned selection to output.")

                # Clear selection, delete single_poly_lyr
                arcpy.management.SelectLayerByAttribute("points_lyr", "CLEAR_SELECTION")
                arcpy.management.Delete("single_poly_lyr")

        arcpy.AddMessage(
            f"Finished iterating polygons. Cleaned points saved to: {output_points}"
        )

        arcpy.management.Delete(temp_points)

    except arcpy.ExecuteError:
        # ArcPy-specific error
        arcpy.AddError(f"ArcPy ExecuteError: {arcpy.GetMessages(2)}")

    except Exception as ex:
        # General Python error
        arcpy.AddError(f"Unexpected error: {str(ex)}")


# -------------------------------------------------------------------------
# 3) Function for the fourth tool: Create Thiessen Polygons
# -------------------------------------------------------------------------
def create_combined_clipped_thiessen(
    in_polygons,
    in_points,
    out_thiessens
):
    """
    Iterates over each polygon in 'in_polygons'. For each polygon:
      1. Selects points in 'in_points' that fall within the polygon.
      2. Creates Thiessen polygons for those selected points.
      3. Clips the Thiessen polygons to the polygon boundary.
      4. Appends the clipped polygons to a single output feature class.
         - If 'out_thiessens' doesn't exist yet, it is created using CopyFeatures.
         - Otherwise, new polygons are appended.

    Parameters
    ----------
    in_polygons : str
        Path to the input polygon feature class.
    in_points : str
        Path to the input point feature class.
    out_thiessens : str
        Path to the final polygon feature class (combined output). This is created
        automatically on the first iteration that produces polygons, and appended
        on subsequent iterations.
    """

    arcpy.AddMessage("Starting creation of combined clipped Thiessen polygons...")


    try:
        # ---------------------------------------------------------------------
        # 1) Create layers from the polygon and point feature classes
        # ---------------------------------------------------------------------
        arcpy.management.MakeFeatureLayer(in_polygons, "polygons_lyr")
        arcpy.management.MakeFeatureLayer(in_points,   "points_lyr")
        arcpy.AddMessage("Created in-memory layers for polygons and points.")

        # A flag to track if we've created the 'out_thiessens' dataset yet
        output_created = False

        # ---------------------------------------------------------------------
        # 2) Iterate over each polygon with a SearchCursor
        # ---------------------------------------------------------------------
        with arcpy.da.SearchCursor(in_polygons, ["OID@"]) as poly_cursor:
            for row in poly_cursor:
                poly_oid = row[0]

                # Make a layer for just this one polygon (filter by OID)
                where_clause = f"OBJECTID = {poly_oid}"
                arcpy.management.MakeFeatureLayer(
                    in_features="polygons_lyr",
                    out_layer="single_poly_lyr",
                    where_clause=where_clause
                )

                # Select points that are WITHIN the single polygon
                arcpy.management.SelectLayerByLocation(
                    in_layer="points_lyr",
                    overlap_type="WITHIN",
                    select_features="single_poly_lyr",
                    selection_type="NEW_SELECTION"
                )

                # Count the selected points
                sel_points_count = int(arcpy.management.GetCount("points_lyr").getOutput(0))
                arcpy.AddMessage(f"Polygon OID={poly_oid} => selected {sel_points_count} point(s).")

                if sel_points_count > 0:
                    # -----------------------------------------------------------------
                    # 2a) Create Thiessen polygons for the selected points
                    #     We'll set the environment extent to the bounding box
                    #     of the single polygon layer to focus the Thiessen creation.
                    # -----------------------------------------------------------------
                    thiessen_temp = arcpy.CreateUniqueName("ThiessenTemp", "in_memory")
                    thiessen_clipped = arcpy.CreateUniqueName("ThiessenClip", "in_memory")

                    # Get the extent of the single polygon (for env)
                    poly_extent = arcpy.Describe("single_poly_lyr").extent

                    # Use EnvManager to limit the processing extent
                    with arcpy.EnvManager(extent=poly_extent):
                        arcpy.AddMessage(f"Creating Thiessen polygons with extent={poly_extent}.")
                        arcpy.analysis.CreateThiessenPolygons(
                            in_features="points_lyr",
                            out_feature_class=thiessen_temp,
                            fields_to_copy="ALL"
                        )

                    # -----------------------------------------------------------------
                    # 2b) Clip those Thiessen polygons to the polygon boundary
                    # -----------------------------------------------------------------
                    arcpy.analysis.Clip(
                        in_features=thiessen_temp,
                        clip_features="single_poly_lyr",
                        out_feature_class=thiessen_clipped
                    )

                    thiessen_count = int(arcpy.management.GetCount(thiessen_clipped).getOutput(0))
                    arcpy.AddMessage(f"Created and clipped {thiessen_count} Thiessen polygon(s).")

                    # Only proceed if we actually got polygons
                    if thiessen_count > 0:
                        # -----------------------------------------------------------------
                        # 2c) Combine (create or append) to out_thiessens
                        # -----------------------------------------------------------------
                        if not arcpy.Exists(out_thiessens):
                            # If out_thiessens doesn't exist yet, we copy the first
                            # set of polygons to create it
                            arcpy.management.CopyFeatures(thiessen_clipped, out_thiessens)
                            output_created = True
                            arcpy.AddMessage(
                                f"Created {out_thiessens} from the first polygon's Thiessen polygons."
                            )
                        else:
                            # Otherwise, append to the existing dataset
                            arcpy.management.Append(
                                inputs=thiessen_clipped,
                                target=out_thiessens,
                                schema_type="NO_TEST"
                            )
                            arcpy.AddMessage(
                                f"Appended clipped Thiessen polygons to existing {out_thiessens}."
                            )

                    # Clean up the in_memory temp FCs
                    arcpy.management.Delete(thiessen_temp)
                    arcpy.management.Delete(thiessen_clipped)

                # Clear the point selection for next iteration
                arcpy.management.SelectLayerByAttribute("points_lyr", "CLEAR_SELECTION")

                # Delete single polygon layer
                arcpy.management.Delete("single_poly_lyr")

        # ---------------------------------------------------------------------
        # If no polygons produced any Thiessen polygons, we won't have created
        # out_thiessens. But that's okay; it just means we had no data to save.
        # ---------------------------------------------------------------------
        if arcpy.Exists(out_thiessens):
            final_count = int(arcpy.management.GetCount(out_thiessens).getOutput(0))
            arcpy.AddMessage(f"Final result has {final_count} polygons in {out_thiessens}.")
        else:
            arcpy.AddWarning("No Thiessen polygons were generated, so no output dataset was created.")

    except arcpy.ExecuteError:
        # Handle ArcPy-specific errors
        arcpy.AddError(f"ArcPy ExecuteError: {arcpy.GetMessages(2)}")

    except Exception as ex:
        # Handle all other Python errors
        arcpy.AddError(f"Unexpected error: {str(ex)}")

    finally:
        # Clean up any leftover layers
        if arcpy.Exists("polygons_lyr"):
            arcpy.management.Delete("polygons_lyr")
        if arcpy.Exists("points_lyr"):
            arcpy.management.Delete("points_lyr")
        if arcpy.Exists("single_poly_lyr"):
            arcpy.management.Delete("single_poly_lyr")

        arcpy.AddMessage("Done with 'create_combined_clipped_thiessen' function.")


# -------------------------------------------------------------------------
# Toolbox Class
# -------------------------------------------------------------------------
class Toolbox(object):
    """
    The main Python Toolbox class. This defines the toolbox and lists the tools within.
    """
    def __init__(self):
        """Initialize the toolbox (the name is the .pyt filename)."""
        self.label = "Yield Process Toolbox"
        self.alias = "YieldProcessTool"
        # The toolbox contains four tools:
        self.tools = [JoinCSVToFeatureClass, 
                      ProcessBufferPass, 
                      CleanPointsByPassIDTool,
                      ThiessenByPolygonTool]


# -------------------------------------------------------------------------
# Tool 1) JoinCSVToFeatureClass
# -------------------------------------------------------------------------
class JoinCSVToFeatureClass(object):
    """
    Converts CSV -> in-memory table -> XY point FC (with only user-supplied fields),
    then does a Spatial Join with a shapefile, optionally sorts by the joined time field.
    """
    def __init__(self):
        self.label = "1. Join CSV to Feature Class"
        self.description = (
            "Creates a new feature class from a shapefile by joining CSV data. "
            "Removes unwanted CSV fields before XYTableToPoint, keeps only coordinate/time fields, "
            "then optionally sorts by the joined time field."
        )
        self.canRunInBackground = True

    def getParameterInfo(self):
        # params = []

        param_shapefile = arcpy.Parameter(
            displayName="Input Shapefile",
            name="in_shapefile",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )
        # params.append(param_shapefile)

        param_csv = arcpy.Parameter(
            displayName="CSV File",
            name="in_csv_file",
            datatype="DEFile",
            parameterType="Required",
            direction="Input"
        )
        # params.append(param_csv)

        param_output_fc = arcpy.Parameter(
            displayName="Output Feature Class",
            name="out_feature_class",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output"
        )
        # params.append(param_output_fc)

        param_csv_time_field = arcpy.Parameter(
            displayName="CSV Time Field (Exact name in CSV)",
            name="csv_time_field",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        # params.append(param_csv_time_field)

        param_long_field = arcpy.Parameter(
            displayName="CSV Longitude Field",
            name="longitude_field",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        # params.append(param_long_field)

        param_lat_field = arcpy.Parameter(
            displayName="CSV Latitude Field",
            name="latitude_field",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        # params.append(param_lat_field)

        param_order_field = arcpy.Parameter(
            displayName="Order Field (Optional)",
            name="order_field",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        )
        # params.append(param_order_field)

        param_coord_system = arcpy.Parameter(
            displayName="Coordinate System for CSV Points (Optional)",
            name="coord_system",
            datatype="GPSpatialReference",
            parameterType="Optional",
            direction="Input"
        )
        # params.append(param_coord_system)

        param_search_tolerance = arcpy.Parameter(
            displayName="Search Tolerance (Meters) (Optional)",
            name="search_tolerance",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        param_search_tolerance.value = 0.0
        # params.append(param_search_tolerance)

        # return params
        return [
            param_shapefile,
            param_csv,
            param_output_fc,
            param_csv_time_field,
            param_long_field,
            param_lat_field,
            param_order_field,
            param_coord_system,
            param_search_tolerance
        ]
    
    def execute(self, parameters, messages):
        # Retrieve parameter values
        shapefile = parameters[0].valueAsText
        csv_file = parameters[1].valueAsText
        output_fc = parameters[2].valueAsText
        csv_time_field = parameters[3].valueAsText
        longitude_field = parameters[4].valueAsText
        latitude_field = parameters[5].valueAsText
        order_field = parameters[6].valueAsText
        coord_system = parameters[7].value
        search_tolerance = parameters[8].value

        if search_tolerance is None:
            search_tolerance = 0.0

        # Call the main function
        join_csv_to_feature_class(
            shapefile=shapefile,
            csv_file=csv_file,
            output_fc=output_fc,
            csv_time_field=csv_time_field,
            longitude_field=longitude_field,
            latitude_field=latitude_field,
            order_field=order_field,
            coord_system=coord_system,
            search_tolerance=search_tolerance
        )

    def isLicensed(self):
        return True

    def updateMessages(self, parameters):
        return

# -------------------------------------------------------------------------
# Tool 2) ProcessBufferPass
# -------------------------------------------------------------------------
class ProcessBufferPass(object):
    """
    Tool that implements the strip trial workflow:
      - Copies an input shapefile into a geodatabase feature dataset
      - Optionally sorts features and removes duplicates
      - Converts points to lines, splits/unsplits them
      - Buffers and dissolves lines
      - Performs a spatial join to associate pass lines with point features
    """
    def __init__(self):
        self.label = "2. Create Buffer Passes"
        self.description = ("Copies an input shapefile to a file GDB, optionally sorts features, "
                            "removes duplicates, performs line operations, and buffers/dissolves lines.")
        self.canRunInBackground = True

    def getParameterInfo(self):
        """
        Define parameter definitions for the Process Strip Trial tool.
        """
        # 0) Input shapefile
        param_in_shp = arcpy.Parameter(
            displayName="Input Shapefile",
            name="in_shapefile",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )

        # 1) File GDB Path
        param_gdb_path = arcpy.Parameter(
            displayName="File GDB Path",
            name="gdb_path",
            datatype="DEWorkspace",  # File GDB is a workspace
            parameterType="Required",
            direction="Input"
        )

        # 2) Working Feature Dataset Name
        param_working_fd = arcpy.Parameter(
            displayName="Working Feature Dataset Name",
            name="working_fd",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )

        # 3) Yield Feature Dataset Name
        param_yield_fd = arcpy.Parameter(
            displayName="Yield Feature Dataset Name",
            name="yield_fd",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )

        # 4) Output Feature Class Name
        param_outFC_name = arcpy.Parameter(
            displayName="Output Feature Class Name",
            name="outFC_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )

        # 5) Sort Feature Class (Boolean)
        param_sortFC = arcpy.Parameter(
            displayName="Sort Feature Class? (Optional!! If already sorted, skip this step)",
            name="sortFC",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        )
        param_sortFC.value = False

        # 6) Sort Field (Used if sortFC = True)
        param_sort_field = arcpy.Parameter(
            displayName="Field to Sort Features by",
            name="sort_field",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        )

        # 7) Sorted Field Name (Store new IDs)
        param_sortedField_name = arcpy.Parameter(
            displayName="Create a new field name for 'Sorted Values' if necessary",
            name="sortedField_name",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        )

        # 8) Remove Duplicates? (Boolean)
        param_remove_duplicates = arcpy.Parameter(
            displayName="Remove Duplicate Geometries?",
            name="remove_duplicates",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        )
        param_remove_duplicates.value = False

        # 9) Swath Width
        param_pass_threshold = arcpy.Parameter(
            displayName="Pass Threshold",
            name="pass_threshold",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        param_pass_threshold.value = 5

        # 10) Min Pass Length
        param_min_pass_length = arcpy.Parameter(
            displayName="Minimum Pass Length",
            name="min_pass_length",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        param_min_pass_length.value = 30.0

        # 11) Buffer Distance
        param_buffer_distance = arcpy.Parameter(
            displayName="Buffer Distance",
            name="buffer_distance",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        param_buffer_distance.value = 6.096

        return [
            param_in_shp,
            param_gdb_path,
            param_working_fd,
            param_yield_fd,
            param_outFC_name,
            param_sortFC,
            param_sort_field,
            param_sortedField_name,
            param_remove_duplicates,
            param_pass_threshold,
            param_min_pass_length,
            param_buffer_distance
        ]

    def execute(self, parameters, messages):
        """
        Execute runs when the tool is launched from ArcGIS.
        """
        # Retrieve parameter values (matching the index in getParameterInfo)
        in_shapefile = parameters[0].valueAsText
        gdb_path = parameters[1].valueAsText
        working_fd = parameters[2].valueAsText
        yield_fd = parameters[3].valueAsText
        outFC_name = parameters[4].valueAsText

        sortFC = parameters[5].value  # Boolean
        sort_field = parameters[6].valueAsText
        sortedField_name = parameters[7].valueAsText
        remove_duplicates = parameters[8].value  # Boolean
        pass_threshold = parameters[9].value
        min_pass_length = parameters[10].value
        buffer_distance = parameters[11].value

        # Call the function
        create_buffer_pass(
            in_shapefile=in_shapefile,
            gdb_path=gdb_path,
            working_fd=working_fd,
            yield_fd=yield_fd,
            outFC_name=outFC_name,
            sortFC=sortFC,
            sort_field=sort_field,
            sortedField_name=sortedField_name,
            remove_duplicates=remove_duplicates,
            pass_threshold=pass_threshold,
            min_pass_length=min_pass_length,
            buffer_distance=buffer_distance
        )
        return

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateMessages(self, parameters):
        """
        Modify messages created by internal validation if needed.
        For example, you could require sort_field if sortFC == True.
        """
        # If 'sortFC' is True but no 'sort_field' is provided, show a warning
        if parameters[5].value and not parameters[6].value:
            parameters[6].setWarningMessage(
                "Sort Field is required if 'Sort Feature Class?' is True."
            )
        return


# -------------------------------------------------------------------------
# Tool 3) CleanPointsByPassIDTool
# -------------------------------------------------------------------------
class CleanPointsByPassIDTool(object):
    """
    This tool copies the input points to a temporary dataset (in_memory),
    iterates through polygons, selects points within a search distance,
    finds the minimum PASS_ID, deletes points with PASS_ID > min, then
    appends the survivors into a new output FC. The original is untouched.
    """
    def __init__(self):
        """Set up the tool (label, description, etc.)."""
        self.label = "3. Clean Points By Pass ID"
        self.description = (
            "Iterates through polygons, selects points by location, finds minimum PASS_ID, "
            "deletes points above that PASS_ID, and appends the cleaned points to a new "
            "output feature class. The original points are not modified."
        )
        self.canRunInBackground = True

    def getParameterInfo(self):
        """
        Define parameters that the user provides in the tool interface.
        """
        # 1) Input point feature class
        param_in_points = arcpy.Parameter(
            displayName="Input Points",
            name="in_points",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )

        # 2) Input polygon feature class
        param_in_polygons = arcpy.Parameter(
            displayName="Input Polygons",
            name="in_polygons",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )

        # 3) Output point feature class (cleaned result)
        param_out_points = arcpy.Parameter(
            displayName="Output Cleaned Points",
            name="out_points",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output"
        )

        # 4) Search distance (e.g. "5 Meters")
        param_search_distance = arcpy.Parameter(
            displayName="Search Distance",
            name="search_distance",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        )
        param_search_distance.value = "0 Meters"  # default is no distance

        # 5) PASS ID Field
        param_pass_id_field = arcpy.Parameter(
            displayName="PASS ID Field",
            name="pass_id_field",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )

        return [
            param_in_points,
            param_in_polygons,
            param_out_points,
            param_search_distance,
            param_pass_id_field
        ]

    def execute(self, parameters, messages):
        """
        The source code of the tool, run when the user clicks 'OK'.
        """
        # Extract parameter values
        in_points = parameters[0].valueAsText
        in_polygons = parameters[1].valueAsText
        out_points = parameters[2].valueAsText
        search_distance = parameters[3].valueAsText
        pass_id_field = parameters[4].valueAsText


        clean_points_by_passid(
            in_points, in_polygons, out_points, search_distance, pass_id_field
        )

        return
    
    def isLicensed(self):
        return True

    def updateMessages(self, parameters):
        return
    

# -------------------------------------------------------------------------
# Tool 4) ThiessenByPolygonTool
# -------------------------------------------------------------------------
class ThiessenByPolygonTool(object):
    """
    Tool that:
      - Iterates polygons,
      - Selects points in each polygon,
      - Creates Thiessen polygons (with env extent),
      - Clips to polygon,
      - Appends to output FC.
    """
    def __init__(self):
        self.label = "4. Create Thiessen per Polygon"
        self.description = (
            "For each polygon, create Thiessen polygons from the points within it, "
            "clip them to the polygon boundary, and append them to an existing empty "
            "polygon feature class."
        )
        self.canRunInBackground = True

    def getParameterInfo(self):
        """
        Define tool parameters:
          1) Input Polygons
          2) Input Points
          3) Output Polygon FC
        """
        # params = []

        # 1) Input polygon feature class
        param_polygons = arcpy.Parameter(
            displayName="Input Polygon Feature Class",
            name="in_polygons",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )


        # 2) Input point feature class
        param_points = arcpy.Parameter(
            displayName="Input Point Feature Class",
            name="in_points",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )


        # 3) Output polygon feature class (existing empty FC to store results)
        param_out_thiessens = arcpy.Parameter(
            displayName="Output Thiessen Feature Class",
            name="out_thiessens",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output"
        )


        return [
            param_polygons,
            param_points,
            param_out_thiessens,
        ]
    
    def execute(self, parameters, messages):
        """
        The main execution code, run when the tool is launched.
        """
        # Extract parameter values
        in_polygons = parameters[0].valueAsText
        in_points   = parameters[1].valueAsText
        out_thiessens = parameters[2].valueAsText

        # Call our function
        create_combined_clipped_thiessen(
            in_polygons=in_polygons,
            in_points=in_points,
            out_thiessens=out_thiessens
        )
        return

    def isLicensed(self):
        return True

    def updateMessages(self, parameters):
        return
        
        
        

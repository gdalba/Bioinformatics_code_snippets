def filter_excel(identifiers, input_file):
    """
    Filter Excel data based on a list of identifiers.
    Searches for identifiers in columns B and C of the Excel file with flexible matching.
    
    Args:
        identifiers: List of identifiers to filter by
        input_file: Path to the Excel file
        
    Returns:
        Filtered pandas DataFrame
    """
    import pandas as pd
    import re
    
    print(f"Reading Excel file: {input_file}")
    print(f"Excel file columns: {df.columns.tolist()}")
    
    clean_identifiers = []
    for identifier in identifiers:
        # Trailing tabs or spaces need to be removed and not empty
        clean_id = identifier.strip()
        if clean_id:
            clean_identifiers.append(clean_id)
    
    identifiers_set = set(clean_identifiers)
    print(f"Looking for {len(identifiers_set)} unique identifiers")
    
    # Debug
    sample_ids = list(identifiers_set)[:5]
    print(f"Sample identifiers: {sample_ids}")
    
    # Create a mask that will be True for rows that contain any identifier in columns B and C
    mask = pd.Series(False, index=df.index)
    
    # Get the column names for columns B and C (index 1 and 2)
    search_columns = []
    if len(df.columns) > 1:
        search_columns.append(df.columns[1])  # Column B
    if len(df.columns) > 2:
        search_columns.append(df.columns[2])  # Column C
    
    print(f"Searching only in columns: {search_columns}")
    
    # Check columns B and C for matches
    for column in search_columns:
        print(f"Checking column: {column}")
        column_values = df[column].astype(str)
        
        # Check for matches using a more efficient approach
        for identifier in identifiers_set:
            # Get the base identifier (part after pipe)
            base_id = identifier.split('|')[1] if '|' in identifier else identifier
            
            # Simpler approach to avoid the warning
            # Look for base_id surrounded by non-alphanumeric chars or at start/end
            # Instead of using regex groups, we'll just search for the base_id
            # This is less precise but avoids the warning
            id_mask = column_values.str.contains(re.escape(base_id), na=False)
            
            if id_mask.any():
                print(f"Found matches for {base_id} in column {column}")
                mask = mask | id_mask
    
    filtered_df = df[mask]
    
    print(f"Original data: {len(df)} rows")
    print(f"Filtered data: {len(filtered_df)} rows")
    
    if len(filtered_df) == 0:
        print("No matches found! Here are the first few rows of the Excel file:")
        print(df.head().to_string())
        print("\nAnd here are a few sample identifiers we're looking for:")
        for i in range(min(5, len(identifiers))):
            print(f"  - {identifiers[i]}")

    return filtered_df

def save_filtered_data(filtered_data, output_file):
    """
    Save the filtered DataFrame to an Excel file.
    
    Args:
        filtered_data: Pandas DataFrame with filtered data
        output_file: Path where to save the Excel file
    """
    filtered_data.to_excel(output_file, index=False)
    
    if len(filtered_data) == 0:
        print(f"Warning: No matching data found. Output file {output_file} contains only headers.")
    else:
        print(f"Successfully saved {len(filtered_data)} rows to {output_file}")
import os
import json
import argparse

def copy_attributes_to_building_parts(cityjson_data):
    """Copy attributes from parent Building to its BuildingParts."""
    # Get the main objects
    city_objects = cityjson_data.get('CityObjects', {})
    
    # Keep track of processed buildings
    for obj_id, obj_data in city_objects.items():
        if obj_data.get('type') == 'Building':
            # Get parent building attributes
            parent_attributes = obj_data.get('attributes', {})
            # Get children (BuildingParts)
            children_ids = obj_data.get('children', [])
            
            # Copy attributes to each child
            for child_id in children_ids:
                if child_id in city_objects:
                    child_obj = city_objects[child_id]
                    if child_obj.get('type') == 'BuildingPart':
                        # Create attributes dict if it doesn't exist
                        if 'attributes' not in child_obj:
                            child_obj['attributes'] = {}
                        # Update child attributes with parent attributes
                        child_obj['attributes'].update(parent_attributes)
    
    return cityjson_data

def process_cityjson_files(input_folder, output_folder, verbose):
    """Recursively find and process CityJSON files."""
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Count processed files
    processed_count = 0
    
    # Walk through directory tree
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            if filename.lower().endswith('.city.jsonl'):
                input_path = os.path.join(root, filename)
                
                # Create corresponding output path
                relative_path = os.path.relpath(root, input_folder)
                output_dir = os.path.join(output_folder, relative_path)
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, filename)
                
                try:
                    # Read the input file
                    with open(input_path, 'r', encoding='utf-8') as f:
                        cityjson_data = json.load(f)
                    
                    # Process the data
                    modified_data = copy_attributes_to_building_parts(cityjson_data)
                    
                    # Write the modified data
                    with open(output_path, 'w', encoding='utf-8') as f:
                        json.dump(modified_data, f)
                    
                    if verbose:
                      processed_count += 1
                      print(f"Processed: {input_path} -> {output_path}")
                    
                except Exception as e:
                    print(f"Error processing {input_path}: {str(e)}")
    
    if verbose:
      print(f"\nProcessing complete. Total files processed: {processed_count}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Copy attributes from Buildings to descendant BuildingParts in CityJSON files"
    )
    parser.add_argument(
        "input_folder",
        help="Path to the input folder containing CityJSON files"
    )
    parser.add_argument(
        "output_folder",
        help="Path to the output folder for processed files"
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose output"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate input folder
    if not os.path.isdir(args.input_folder):
        print(f"Error: Input folder '{args.input_folder}' does not exist!")
        return
    
    if args.verbose:
        print(f"Starting processing...")
        print(f"Input folder: {args.input_folder}")
        print(f"Output folder: {args.output_folder}")
    
    # Process the files
    process_cityjson_files(args.input_folder, args.output_folder, args.verbose)

if __name__ == "__main__":
    main()
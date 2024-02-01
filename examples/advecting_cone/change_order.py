import argparse
import json

def modify_polynomial_order(json_file_path, new_polynomial_order):
    # Open the JSON file for reading
    with open(json_file_path, 'r') as file:
        data = json.load(file)

    # Update the polynomial_order value under the specified path
    path = ["case", "numerics", "polynomial_order"]
    current_node = data
    for key in path[:-1]:
        current_node = current_node[key]
    current_node[path[-1]] = new_polynomial_order

    # Save the modified data back to the JSON file
    with open(json_file_path, 'w') as file:
        json.dump(data, file, indent=2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modify polynomial order in the case file")
    parser.add_argument("new_polynomial_order", type=int, help="New polynomial order value")

    args = parser.parse_args()

    new_polynomial_order = args.new_polynomial_order

    modify_polynomial_order("advecting_cone.case", new_polynomial_order)

import sys
import os
"""
Author: Michael Hathaway
Description: script writes the absolute paths for all files in a directory to
an output file
"""


def get_directory_contents(dir_name):
    paths_list = []
    root_path = os.getcwd()
    for root, _, files in os.walk(dir_name):
        for file in files:
            paths_list.append(os.path.abspath(os.path.join(root, file)))

    return paths_list


def write_paths_to_file(paths_list, output_name):
    with open(output_name, "w") as f:
        for path in paths_list:
            f.write(f"{path}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python3 {sys.argv[0]} <dir> <output_file>")
        sys.exit()

    dir_name = sys.argv[1]
    output_name = sys.argv[2]

    paths_list = get_directory_contents(dir_name)
    write_paths_to_file(paths_list, output_name)

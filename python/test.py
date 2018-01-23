import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

import parameters

a = Parameters("parameters.xml")

def main():
    print(a.config_path)


if __name__ == "__main__":
    main()

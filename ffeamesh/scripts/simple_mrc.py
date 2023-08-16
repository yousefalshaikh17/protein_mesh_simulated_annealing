import numpy as np
import pathlib
from ffeamesh.mrc_utility import (CellSize, CellAngles, CellProps, write_mrcfile)

def fill_test(image):
    for i in range(0, 4):
        for j in range(0, 4):
            for k in range(0, 4):
                image[i][j][k] = 0.8

    image[2][2][2] = 1.2

def main():
    cell_size = 0.5
    cell_angle = 90.0
    label="Simple test"
    test_image = np.full((5, 5, 5), dtype=np.float32, fill_value=0.2)

    props = CellProps(CellSize(cell_size, cell_size, cell_size),
                      CellAngles(cell_angle, cell_angle, cell_angle))

    write_mrcfile(test_image, props, pathlib.Path("simple_test.mrc"), label, True)

if __name__ == "__main__":
    main()

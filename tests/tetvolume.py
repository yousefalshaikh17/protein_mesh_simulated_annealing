from ffeamesh.fivetets import (Coordinate, CoordTransform, tet_volume)

def main():
    """
    test using (0,0,0) (18,0,0) (9,18,0) (9,9,9)
    """
    target = 486.0

    coords = []
    coords.append(CoordTransform(Coordinate(0,0,0), Coordinate(0.0, 0.0, 0.0)))
    coords.append(CoordTransform(Coordinate(0,0,0), Coordinate(18.0, 0.0, 0.0)))
    coords.append(CoordTransform(Coordinate(0,0,0), Coordinate(9.0, 18.0, 0.0)))
    coords.append(CoordTransform(Coordinate(0,0,0), Coordinate(9.0, 9.0, 9.0)))

    volume = tet_volume(coords)
    message = "{}: volume calculated was {} should be {}."

    if abs(target-volume) > 0.0:
        print(message.format("Fail", volume, target))
    else:
        print(message.format("Pass", volume, target))

if __name__ == '__main__':
    main()

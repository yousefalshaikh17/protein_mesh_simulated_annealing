from ffeamesh.fivetets import VolumeBins

def main():
    vb = VolumeBins()

    for i in [1.0, 1.23, 1.7, 2.3, 1.1, 2.7]:
        vb.add(i)

    print("Should be:")
    print("1.0 => 3")
    print("2.0 => 2")
    print("3.0 => 1")

    print("Was:")
    for a, b in vb.bin_counts.items():
        print(f"{a} => {b}")

if __name__ == '__main__':
    main()
import numpde

def main():
    c = int(input("c: "))
    dt = float(input("dt: "))
    dx = float(input("dx: "))
    T = int(input("T: "))
    L = int(input("L: "))

    numpde.ex1(c, dt, T, dx, L)

if __name__ == '__main__':
    main()

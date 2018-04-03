import sys

def print_stats(f):
    line = f.readline()
    s_co = 0
    p_co = 0
    while line:
        if line[0] == 'S':
            if line.find('*') == -1:
                s_co += 1
        elif line[0] == 'P':
            p_co += 1
        line = f.readline()
    print("S: " + str(s_co))
    print("P: " + str(p_co))

def main():
    f1_name = str(sys.argv[1])   #GFA1 file
    f2_name = str(sys.argv[2])   #Binary Decoded file

    f1 = open(f1_name, 'r')
    f2 = open(f2_name, 'r')

    print_stats(f1)
    print_stats(f2)

    f1.close()
    f2.close()

if __name__ == '__main__':
    main()

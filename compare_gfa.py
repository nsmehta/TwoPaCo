import sys
import mmap
import multiprocessing
from functools import partial
from textwrap import dedent
from itertools import zip_longest

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def process_chunk(file_two, d):
    if d == None:
        return
    f_2 = open(file_two, 'rb', 0)
    file_map = mmap.mmap(f_2.fileno(), 0, access=mmap.ACCESS_READ)
    trans = d.strip()
    if file_map.find(trans.encode()) != -1:
        print("Match Found")
    else:
        reverse = reverse_complement(trans)
        if file_map.find(reverse.encode()) != -1:
            print("Reverse Complement Match found")
        else:
            print("No match found for : " + trans)

def grouper(n, iterable, padvalue=None):
    return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)


if __name__ == "__main__":
    print("First file : " + sys.argv[1])
    print("Second file : " + sys.argv[2])
    print()
    
    file_one = sys.argv[1]
    file_two = sys.argv[2]

    #f_2 = open(file_two, 'rb', 0)
    #file_map = mmap.mmap(f_2.fileno(), 0, access=mmap.ACCESS_READ)


    p = multiprocessing.Pool(10)
    input_data = open(file_one, 'r')

    for chunk in grouper(500, input_data):
        func = partial(process_chunk, file_two)
        p.map(func, chunk)



    print("Completed")

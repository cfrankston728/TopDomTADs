import argparse
import os
import csv
import subprocess
import multiprocessing

import numpy as np
import inspect

import math
import hicstraw
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.pyplot as plt
#import seaborn as sns

from bisect import insort

from pytadbit import tadbit
from pytadbit import Chromosome
from pytadbit.parsers.hic_parser import load_hic_data_from_bam

# Get the relevant putative TAD data from the TAD .bed file
#   Get: Start index, end index
#   Retrieve .hic sub-matrix [distance normalized?] 
#       --> investigate boundary and interior characteristics 
#               -normalized signal maxima along each boundary and normalized signal maximum from interior
#   Decide whether putative TAD satisfies the strong boundary property

#   Can we also determine any sub-TADs based on the strong boundary property (SBP)?
#       What happens to enrichments if we restrict to TADs with SBP?
#       Can we distinguish between gaps and TADs according to the putative TAD characteristics?
#           I cannot determine this since I have screened out gaps--I will have to remove the gap filter

#################################################################################
# INITIAL PARAMETERS

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Process Hi-C data and detect TAD boundaries.')
parser.add_argument('--hic_type', type=str, default='IN_SITU',
                    help='Hi-C data type: "IN_SITU" or "IN_TACT" (default: IN_SITU)')
parser.add_argument('--window', type=int, default=5,
                    help='Window size for TopDom TAD boundary detection (default: 5)')
parser.add_argument('--reso', type=int, default=10000,
                    help='Base-pair resolution for the .hic data (default: 10000 = 10 kb)')

args = parser.parse_args()

hic_type = args.hic_type
window = args.window
reso = args.reso

if hic_type == 'IN_SITU':
    hic_path = "https://www.encodeproject.org/files/ENCFF555ISR/@@download/ENCFF555ISR.hic" # in-situ hiC
elif hic_type == 'IN_TACT':
    hic_path = "https://www.encodeproject.org/files/ENCFF318GOM/@@download/ENCFF318GOM.hic" # in-tact hiC (microC)

# Define a useful auxiliary function that will enhance the detection of strong boundary domains.
# In particular, we will want to produce lower bounds on the interior maximum of many domains. Therefore it is sensible to check the locations in the matrix that are shared by as many domains as possible. The non-exhaustive algorithm to prioritize these matrix locations requires determining the divisors of a certain integer within a fixed range (or more generally within a given set).
# There are many different algorithms to accomplish this depending on the size of the integers in question. For now I will do a brute force search over a sieved set based on prime factors of the integer in question (those in p_list).
def div_power(args):
    K, factor = args
    if factor == 1:
        return [1,float('inf')]
    power = 0
    prod = factor
    while K % prod == 0:
        power += 1
        prod = prod * factor
    return [factor,power]

def find_powers(p_list, K_values):
    with multiprocessing.Pool() as pool:
        p_powers = pool.map(div_power, [(K, factor) for K in K_values for factor in p_list])
    return np.array(p_powers).reshape(-1, len(p_list), 2)

def div(args):
    K, factor = args
    is_div = (K % factor == 0)
    return [factor, int(K/factor)*is_div, K*is_div]

def find_divs(p_list, K_values):
    with multiprocessing.Pool() as pool:
        divs = pool.map(div, [(K, factor) for K in K_values for factor in p_list])
    
    return np.array(divs).reshape(-1, len(p_list), 3)

# Define a function isSPD that takes a .hic map and domain coordinates as an input, and returns TRUE if the domain is a strong boundary domain, and FALSE otherwise. Additionally, provide the three features determining the output of the function:
#   (1) The maximum signal and its coordinate within the left boundary.
#   (2) The maximum signal and its coordinate within the right boundary.
#   (3) The maximum signal within the interior.
#       If the maximum signal coordinates of both boundaries are the same, then the domain perhaps can also be called an apex domain.
def isSBD(numpy_matrix):#, domain_start, domain_end):
    #numpy_matrix = mzd.getRecordsAsMatrix(domain_start, domain_end, domain_start, domain_end)

    interior = numpy_matrix[1:-1,1:-1]

    if interior.size > 0:
        interior_max = max(max(interior))
    else:
        return True
    
    A=numpy_matrix[-1,-1]

    return (min(max(numpy_matrix[0,0:]), 
                max(numpy_matrix[0:,-1])) >= interior_max)

# If a domain is an SBD, stop and return its coordinates.
# Otherwise, put its reduced sub-domains on the end of the queue and move to the next domain in the queue.
# check_queue will be a list containing the matrix mzd [matrix zoom data].
# depth_range will indicate the depths between which to search for an SBD
def findSBD(matrix, depth_range=(0,50)):#domain_start, domain_end, mzd, depth_range=(0,100)):
    
    #numpy_matrix = mzd.getRecordsAsMatrix(domain_start, domain_end, domain_start, domain_end)

    for depth in range(depth_range[0], depth_range[1]+1):
        if (isSBD(matrix[depth:,depth:])):
            return depth
        else:
            for h in range(0,depth):
                if (isSBD(matrix[h:h-depth,h:h-depth])):
                    return depth

    return None

# Create a function to take in a "greedy index" and return the next greedy index
# If passed, d_list should be a complete list of divisor pairs d < co(d) of d1*d2 where d2 < d, of the form [[d, co(d), d*co(d)],[d', co(d'), d'*co(d')],...], and where d is increasing in the list.
def next_greedy_index(n, d1=None, d2=None, d_list=None):
    this_d1 = d1
    this_d2 = d2
    this_d_list = d_list

    if this_d1 is None:
        this_d1 = math.floor(n/2)
    if this_d2 is None:
        this_d2 = math.ceil(n/2)
        
    if this_d1 < this_d2:
        return [this_d2, this_d1, this_d_list]
    
    next_K = this_d1*this_d2 - 1
    if next_K == 0:
        return [1, 1, None]
    
    # If we have been passed d_list, and it is not empty, we do not need to compute the d_list and can simply return the next item in the sequence.
    if not (this_d_list is None) and len(this_d_list) > 0:
        index = 0
        for d in this_d_list:
            if d[0] > this_d2:
                return [d[0], d[1], this_d_list[index:]]
            index += 1
        this_d_list = []

    # If no d_list is passed, construct the appropriate one.
    if this_d_list is None:

        # The value of d2 and d1 indicate where in the greedy indexing search we are.
        # We will not need to consider divisors that are less than d2.
        if n % 2 == 0:
            k = n/2
            this_x = k*k - (next_K + 1)
            check_list = range(math.ceil(k-math.sqrt(this_x)),math.floor(math.sqrt(k**2-this_x))+1)
        else:
            k = (n-1)/2
            this_x =  k*(k+1) - (next_K + 1)
            check_list = range(math.ceil(k+1/2-math.sqrt(this_x+1/4)),math.floor(math.sqrt((k+1/2)**2-(this_x+1/4)))+1)

        this_d_list = [ds for ds in find_divs(check_list, [next_K + 1])[0] if ds[2] != 0 and ds[0] > this_d2]
        if len(this_d_list) > 0:
            return [this_d_list[0][0], this_d_list[0][1], this_d_list]
    
    # If d_list is passed or constructed as empty, then we need to move on to the next K.
    if len(this_d_list) == 0:
        while len(this_d_list) == 0:
            if n % 2 == 0:
                k = n/2
                this_x = k*k - next_K
                check_list = range(math.ceil(k-math.sqrt(this_x)),math.floor(math.sqrt(k**2-this_x))+1)
            else:
                k = (n-1)/2
                this_x = k*(k+1) - next_K
                check_list = range(math.ceil(k+1/2-math.sqrt(this_x+1/4)),math.floor(math.sqrt((k+1/2)**2-(this_x+1/4)))+1)
            this_d_list = [ds for ds in find_divs(check_list, [next_K])[0] if ds[2] != 0]
            if len(this_d_list) > 0:
                return [this_d_list[0][0], this_d_list[0][1], this_d_list]
            next_K = next_K - 1
            if next_K == 0:
                return [1, 1, None]

# Create a symmetric matrix object, which can be constructed by passing a binary function F over integers to the constructor, the ranges for horizontal and vertical indices, and set of arguments for evaluation.
# If the function is specified as being symmetric, it will be evaluated only over indices where the vertical index is greater than or equal to the horizontal index. Otherwise it will take the arithmetic mean of both transposed inputs.
#  For fast checking, one may wish to pass h_vals in descending order, and v_vals in ascending order. This makes checking to form the dictionary faster.
# Furthermore, it will not store zero values in its dictionary.
# Finally, one can pass F either as a binary function, as a numpy matrix, or as an already-existing SymmetricMatrix object.
class SymmetricMatrix(object):
    def __init__(self, F, h_vals, v_vals, symm=True, good_order=True, set_order=True):
        super().__init__()
        self.F = F
        self.length = -float("inf")
        self.depth = -float("inf")
        
        self.dict = dict()
        
        self.callable = callable(F)
        self.SymMat = isinstance(self.F, SymmetricMatrix)

        # If F is a SymmetricMatrix, then it is surely symmetric...
        if self.SymMat:
            self.symm = True
        else:
            self.symm = symm

        self.good_order = good_order
        
        if self.good_order:
            self.h_vals = h_vals
            self.v_vals = v_vals
        elif set_order:
            self.h_vals = sorted(h_vals, reverse=True)
            self.v_vals = sorted(v_vals)

        self.make_dict()

    def make_dict(self):
        self.dict = {}
        for h_val in self.h_vals:
            in_dict = {}
            for v_val in self.v_vals:
                # If the values are in a good order, this means h_vals are descending and v_vals are ascending.
                # We do not want to store values where h_val > v_val
                if (self.good_order and (int(v_val) > int(h_val))) or (v_val < 0) or (h_val < 0):
                    break
                if self.symm:
                    if self.callable:
                        value = self.F(h_val, v_val)
                    elif self.SymMat:
                        value = self.F.get(h_val, v_val)
                    else:
                        value = self.F[h_val, v_val]
                else:
                    if self.callable:
                        value = (self.F(h_val, v_val) + self.F(v_val, h_val)) / 2
                    else:
                        value = (self.F[h_val, v_val] + self.F[v_val, h_val]) / 2
                if int(h_val) != 0:
                    if value != 0:
                        in_dict[int(v_val)] = value
                    self.dict[int(h_val)] = in_dict
                else:
                    self.dict[int(h_val)] = value
                if h_val+1 >= self.length:
                        self.length = h_val+1
                if v_val+1 >= self.depth:
                    self.depth = v_val+1


    def get(self, a:int, b:int):
        if (a < b):
            return self.get(b,a)
        if (a < 0) or (b < 0):
            return 0
        if (a not in self.dict or b not in self.dict[a]):
            return 0
        if a == 0:
            return self.dict[a]
        return self.dict[a][b]

def symmetric_compress(matrix):
    h_vals = np.arange(len(matrix[0])-1, -1, -1)
    v_vals = np.arange(0, len(matrix[0]), 1)
    return SymmetricMatrix(matrix, h_vals, v_vals)

if __name__ == "__main__":
    random_matrix = np.random.rand(12, 12)
    sym_mat = symmetric_compress(random_matrix)
    
    this_index = [None, None, None]
    descending_intensity_list = []
    while this_index[0] is None or this_index[0]*this_index[1] != 1:
        next_index = next_greedy_index(12, this_index[0], this_index[1], this_index[2])
        this_index = next_index
        this_weight = sym_mat.get(int(12-this_index[0]), int(this_index[1]))
        # insort(descending_intensity_list, (this_index[0],this_index[1],this_weight), key=lambda x: x[2]) # insort incompatible with pyTAD_analysis env version of python :/
    print(descending_intensity_list)
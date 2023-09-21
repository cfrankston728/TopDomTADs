import argparse

import numpy as np
import math

import hicstraw

from pytadbit import tadbit
from pytadbit import Chromosome

# If a control is desired, one may like to employ the "random" package as well:
# import random

# INPUT: 
#   .hic file parameters for Database extraction (ENCODE);
#   filter parameters for Data quality;
#   TAD caller parameters (TopDom).

# OUTPUT: .bed file of Putative TAD boundaries 

#################################################################################
# OVERVIEW:

# 0. Set initial parameters.
# 1. Get HiC object.
# 2. For each Chromosome in the HiC object: Get Chromosome's Zoom Matrix.
# 3. For each Chunk in the Zoom Matrix: Get Chunk's TAD Boundary List.
# 4. For each Chunk's TAD Boundary List: Write TAD Boundaries to .Bed File

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

print(f'This resolution is {reso}.')
# "hic_path" is a path to a HiC object file containing chromatin conformation data.

if hic_type == 'IN_SITU':
    hic_path = "https://www.encodeproject.org/files/ENCFF555ISR/@@download/ENCFF555ISR.hic" # in-situ hiC
elif hic_type == 'IN_TACT':
    hic_path = "https://www.encodeproject.org/files/ENCFF318GOM/@@download/ENCFF318GOM.hic" # in-tact hiC (microC)

# "window" is the window size for the TopDom TAD boundary detection algorithm employed by tadbit.

# tadbit's TopDom algorithm has a bug with small matrices--it may return a division by zero error due to too low of signal. Therefore I will screen out matrix chunks that are too small to reliably pass into TopDom. If the matrix chunk has a size less than the TopDomScreen parameter, it will be ignored and no TAD boudaries will be reported for the corresponding small region.
# The bug can also possibly be addressed with a throw and catch, I may address that later.
# Alternatively, perhaps if the signal is zero it can be smoothed.
TopDomScreen = 22

# "chunk_size" determines how to break a large matrix into smaller chunks.
# If chunk_size is large, the computation may become too demanding.
# It is not necessarily better to use larger chunk sizes, even if resources allow.
chunk_size = 100*reso

# "gap_radius" determines how much space to give around a detected gap in the HiC data.
gap_radius = 4

# "gap_threshold" is a threshold above which a locus must score in order not to be defined as a gap in the HiC data.
# It's value is based on observations of HiC data at detected TAD boundaries at the given resolution and chunk size.
# It is modified in order to avoid spurious detection of TAD boundaries due to such gaps.
gap_threshold = 70 # 70 is good for 10kb resolution on in-tact hiC (microC)

# "signal_threshold" is a threshold below which an expansive gap will be identified.
signal_threshold = 200 # 200 is good for 10kb resolution and in-tact hiC (microC)

# "bed_name" is the file name in which to write the detected TAD boundary information in .bed format.
# "w" is for window size, "r" is for resolution (in kilo base pairs)
# If writing to a CONTROL file, the file name should be changed accordingly.
bed_name = f"GM12878_TopDom_window_{window}_r{int(reso/1000)}kbp_TAD_boundaries_{hic_type}.bed"

#################################################################################
# AUXILIARY FUNCTIONS
                           
def getFirstGap(L, chunk, threshold):
    print("called getFirstGap.")
    # L = len(chunk[0])
    # G is an index to find the first gap in the chunk along its diagonal.
    G = -1
    
    # "this_cross_score" is the heuristic score used to detect a gap at a given locus along the chunk diagonal.
    this_cross_score = threshold
    
    # If this_cross_score is below the "threshold", we detect a gap at the diagonal locus.
    print(f"L is equal to {L}.")
    while this_cross_score >= threshold and G < L-1:
        G += 1
        # Iterate through the diagonal coordinates in the chunk        
        # The cross_score is a localized measure of the gap near the diagonal. 
        # The following diagram illustrates the loci in the matrix of interest for gap detection:
        
        # [ ][ ][ ][ ][ ][ ][ ][X]
        # [ ][ ][ ][ ][ ][ ][X][ ]
        # [ ][ ][ ][ ][O][X][ ][ ]
        # [ ][ ][ ][O][G][ ][ ][ ]
        # [ ][ ][ ][X][ ][ ][ ][ ]
        # [ ][ ][X][ ][ ][ ][ ][ ]
        # [ ][X][ ][ ][ ][ ][ ][ ]
        # [X][ ][ ][ ][ ][ ][ ][ ]
        
        neighborhood_radius = min(G, 10, L-1-G)
        this_cross_score = 0
        for y in range(neighborhood_radius+1):
            this_cross_score += (chunk[G,G-y] + chunk[G,G+y])
        this_cross_score = this_cross_score/(2*neighborhood_radius + 1)
        print(f'This cross score is {this_cross_score}.')
    return G

def getChunkTADBounds(chunk, chunk_start, chunk_end, gap_threshold, signal_threshold, gap_radius):
    print("called getChunkTADBounds.")
    # First, we need to screen out any gaps that may appear in the chunk using getFirstGap(...)
    # Rather than compute the chunk length twice, we compute once and pass the argument:
    L = chunk.shape[0]
    G = getFirstGap(L, chunk, gap_threshold)
    print(f'This chunk start in the Chromosome is {chunk_start}, and the end is {chunk_end}.')
    print(f'As a matrix its length is {L}, and the gap index is {G}')
    # If a gap is detected, then G will be strictly less than L-1.
    # In this case, we need to split the chunk into two separate chunks on either side of the gap.
    
    # Next, we check to see how extensive the gap is. We do this by summing the diagonal entries over an increasing expanse, and checking the running average signal over the diagonal in the expanse until it passes above a signal threshold. This will determine the end of the gap, and this divides the chunk into left and right sections. The run_length determines the maximum size of the expanse.
    chunkTADBounds = {}
    chunkTADBounds['start']=[]
    chunkTADBounds['tag']=[]
    chunkTADBounds['score']=[]
    
    if G < L-1:
        print(f'The gap index {G} is less than L-1 = {L-1}.')
        gap_left = max(0, G-gap_radius)
        print(f'Therefore we compute the gap_left index as {gap_left}.')
        left_chunk = chunk[0:gap_left, 0:gap_left]
        left_chunk_start = chunk_start
        left_chunk_end = chunk_start + reso*gap_left
        print(f'The left chunk_end in the Chromosome is {left_chunk_end}')
        
        average_signal = 0
        signal_sum = 0
        expanse = 0
        run_length = 50
        while average_signal < signal_threshold and G + expanse < L: # L is the length of the current chunk.
            signal_sum += chunk[G+expanse,G+expanse]
            expanse += 1
            if expanse >= run_length:
                signal_sum -= chunk[G+expanse-run_length,G+expanse-run_length]
                average_signal = signal_sum/run_length
            else:
                average_signal = signal_sum/expanse
            print(f'The average signal is {average_signal}')
        print(f'The gap extends {expanse} indices to {G+expanse} until it reaches an average signal of {average_signal}.')
        gap_right = min(G+expanse+gap_radius, L-1)
        print(f'Now the gap_right is computed as {gap_right}')
        right_chunk = chunk[gap_right:L, gap_right:L]
        right_chunk_start = chunk_start + reso*gap_right
        print(f'The right chunk now starts at {right_chunk_start} in the Chromosome')
        right_chunk_end = chunk_end
        
        # [ ][ ][ ][ ][ ][ ][R][R]
        # [ ][ ][ ][ ][ ][ ][R][ ]
        # [ ][ ][ ][X][O][X][ ][ ]
        # [ ][ ][ ][O][G][ ][ ][ ]
        # [ ][ ][ ][X][ ][ ][ ][ ]
        # [L][L][L][ ][ ][ ][ ][ ]
        # [L][L][ ][ ][ ][ ][ ][ ]
        # [L][ ][ ][ ][ ][ ][ ][ ]
        
        print(f'Now we examine the left chunk:')
        leftChunkTADBounds = getChunkTADBounds(left_chunk,left_chunk_start,left_chunk_end,gap_threshold,signal_threshold, gap_radius)
        
        print(f'Now we examine the right chunk:')
        rightChunkTADBounds = getChunkTADBounds(right_chunk,right_chunk_start,right_chunk_end,gap_threshold,signal_threshold, gap_radius)
        
        print(f'Now we combine the left and right chunks together:')

        chunkTADBounds['start'] = leftChunkTADBounds['start'] + rightChunkTADBounds['start']
        
        chunkTADBounds['tag'] = leftChunkTADBounds['tag'] + rightChunkTADBounds['tag']
        
        chunkTADBounds['score'] = leftChunkTADBounds['score'] + rightChunkTADBounds['score']
        
        print(f'The original chunk from {chunk_start} to {chunk_end} was broken into a left chunk from {left_chunk_start} to {left_chunk_end}, a gap, and a right chunk from {right_chunk_start} to {right_chunk_end}.')
        
        return chunkTADBounds
    
    # TopDom has issues operating on matrices that are too small. It will yield "division by zero" errors. To guard against that, I will not pass chunks with a size less than the TopDomScreen into TopDom.
    if L >= TopDomScreen:
        print("passed TopDomScreen.")
        chunkTADBounds = tadbit([chunk.tolist()], use_topdom=True, topdom_window=window)
        print(f'TopDom was used on this chunk from {chunk_start} to {chunk_end}.')
        chunkTADBounds['start'] = [chunk_start + reso*bound for bound in chunkTADBounds['start']]
        return chunkTADBounds
    return chunkTADBounds
    
    
#################################################################################
# PROCESS
with open(bed_name, "w") as bed_file:
    
    # Initialize an index for each TAD boundary entry in the .bed file
    index = 1 
    
    # Get HiC object.
    hic = hicstraw.HiCFile(hic_path)
    print("used hicstraw.")
    print(f"resolutions available: {hic.getResolutions()}")

    # Get the list of Chromosomes in the HiC object:
    chrms = hic.getChromosomes()
    print("used getChromosomes.")
    
    # For each Chromosome in the HiC object (starting at index 1, since index 0 is "All"):
    for this_chr in chrms[1:3]:#-1]: #MAKE SURE FIRST CHROMOSOME 1 IS COMPLETE by moving from CHR1 to CHR2. 
        
        # Get the Chromosome's Matrix Zoom Data:
        chrname = this_chr.name

        mzd = hic.getMatrixZoomData(chrname, chrname,"observed","NONE","BP",reso)
        print("used getMatrixZoomData.")
        # We will require the Chromosome length for processing:
        chrlength = this_chr.length
        print(f'Chromosome: {chrname}, length: {chrlength}.')
        
        # Now we move on to breaking the Matrix Zoom Data into manageable Chunks. 
        # Each chunk is bounded in size by the chunk_size parameter.
        # TAD boundaries will be computed for each chunk in sequence
        # Each successive chunk will begin at the TAD boundary with the highest index in the previous chunk.
        # Computing the TAD boundaries in each chunk will be facilitated by auxiliary functions:
            # getChunkTADBounds(...), which will be facilitated by
                # getFirstGap(...) in order to recursivley avoid gaps in the HiC data.
        
        # The chunk_start and chunk_end coordinates are initialized to 0 at the beginning of each Chromosome
        chunk_start = 0
        chunk_end = 0
        
        # The chunk will be restricted so as not to extend beyond the confines of the Chromosome matrix zoom data;
        # A resolution-sized buffer is also employed to avoid edge effects.
        while max(chunk_start, chunk_end) < chrlength-reso:
        
            # The chunk_end will either extend to the full chunk_size, or to the edge of the Chromosome matrix:
            chunk_end = min(chunk_start + chunk_size + reso, chrlength-1)
            
            # The chunk itself is then extracted from the zoom data and processed into an ordinary matrix
            this_chunk = mzd.getRecordsAsMatrix(chunk_start, chunk_end, chunk_start, chunk_end) 
            print("used getRecordsAsMatrix.")
            print(this_chunk)
            
            #### this_chunk = [this_chunk.tolist()]

            # Next, we get this chunk's TAD Boundary List.
            # Use function getChunkTADBounds(...)
                # Returns TopDom TAD boundaries detected within the chunk according to chosen parameters, 
                # Recursively avoid gaps in the HiC data to reduce spurious detection
            these_TAD_bounds = getChunkTADBounds(this_chunk, chunk_start, chunk_end, gap_threshold, signal_threshold, gap_radius)
            # print(f'Chunk_TAD_bounds: {these_TAD_bounds}')
            # The next chunk will begin at the highest recorded TAD boundary index "maxTAD"
            # It is initialized as chunk_start, but will be updated as larger TAD boundary indices are found.
            maxTAD = chunk_start            
            
            # Write these_TAD_bounds to the opened .bed file:       
            for tad_num in range(len(these_TAD_bounds['start'])):
                
                # Record only the 'domain' TADs(?)
                if these_TAD_bounds['tag'][tad_num]=='domain':
                    
                    # this_TAD_bound will be the boundary position as an index in the Chromosome matrix
                    this_TAD_bound = these_TAD_bounds['start'][tad_num]
                    
                    # Update the maxTAD as larger TAD boundary indices are found
                    if this_TAD_bound > maxTAD:
                        maxTAD = this_TAD_bound
                    
                    # We also want to record the TopDom score of each boundary in the .bed file
                    # This will enable analysis of algorithmic confidence in reference to enrichment data.
                    Score = these_TAD_bounds['score'][tad_num]

                    bed_file.write(f'{chrname}\t{int(max(1,this_TAD_bound-reso/2))}\t{int(min(chrlength-1,this_TAD_bound+reso/2))}\tregion_{index}\t{Score}\n')
                    # print(f'Region {index} written, on Chromosome {chrname} at locus {int(max(1,this_TAD_bound-reso/2))}.')
                    index += 1
                    
            # After the TAD boundaries within this chunk have been recorded, we move on to the next chunk.
            # The next chunk will start at the TAD boundary in the previous chunk that has the largest index:
            # However, it may be that only a single TAD boundary is detected in the chunk, and that this detected TAD boundary just happens to be at the chunk_start. This will cause an infinite loop, since the movement through the chromosome matrix will not progress. Therefore we need to check for this case, and move the process along ``manually."
            # In this case, there are no other TADs detected in the chunk, so we can move to (pretty much) the end of the chunk: 
            if maxTAD == chunk_start:
                chunk_start = chunk_end - reso
            else:
                chunk_start = maxTAD
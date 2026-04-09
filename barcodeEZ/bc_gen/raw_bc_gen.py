import random
from rapidfuzz import distance
import os

# generate n barcodes
n = 20000  # number of barcodes to generate
bc_length = 60

def genBC(count: int, seq_len: int):#, out_file: str):
    barcodes = {'bc': [],
                'edit_distance': [bc_length]} # set first distance to bc_length 
    while len(barcodes['bc']) < count:
        seq = ''.join(random.choices(['A', 'T', 'C', 'G'], k=seq_len))
        GC = (seq.count('G') + seq.count('C')) / len(seq)

        # GC length check
        if GC > 0.59 or GC < 0.41:
            continue
        # homopolymer check
        hp_flag = False
        for hp in ['AAA', 'CCC', 'GGG', 'TTT']:
            if hp in seq:
                hp_flag = True
                break
        if hp_flag: continue
        
        local_dists = []
        for bc in barcodes['bc']:
            dist = distance.Levenshtein.distance(seq, bc)
            local_dists.append(dist)
        
        if barcodes['bc']: # if not empty
            min_dist = min(local_dists)
            if min_dist < 0.3 * seq_len:
                continue
            else:
                barcodes['edit_distance'].append(min(local_dists))
        
        barcodes['bc'].append(seq)

        if len(barcodes['bc']) % 1000 == 0:
            print(f'{len(barcodes["bc"])} barcodes generated...')
    return barcodes

BCs = genBC(n, bc_length)

fh_out = f'{str(int(n/1000))}k_barcodes_{bc_length}mers.fa'
with open(fh_out, 'w') as out:
    for i,bc in enumerate(BCs['bc']):
        out.write(f'>BC_{i+1}\n{bc}\n')

os.system(f'gzip {fh_out}')
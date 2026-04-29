import math
import Bio.Restriction as RS
import random
import gzip
from importlib.resources import files
import pandas as pd

with gzip.open(str(files('barcodeEZ.bc_gen').joinpath('20k_barcodes_60mers.fa.gz')), 'rt') as f:
    bc_pool = []
    for i,line in enumerate(f):
        line = line.strip()
        if i % 2 == 1:
            bc_pool.append(line)

class Site:

    def __init__(self, site_id, left_enzyme, right_enzyme):
        self.id = site_id
        self.left_enzyme = left_enzyme
        self.right_enzyme = right_enzyme
        self.bc_only = []
        self.bc = []
        self.fixed_left = ''
        self.fixed_right = ''
        self.overhangs = []

    def __repr__(self):
        n_bc = len(self.bc)
        bc_len = len(self.bc[0]) if self.bc else 0
        return (f"Site {self.id}: {self.left_enzyme} -- SITE{self.id} -- {self.right_enzyme} "
                f"| {n_bc} barcodes, {bc_len}bp")

    def add_fixed_sequence(self, seq, side):
        if side == 'left':
            self.fixed_left = seq
        else:
            self.fixed_right = seq
        self._rebuild_bc()

    def _rebuild_bc(self):
        self.bc = [self.fixed_left + bc + self.fixed_right for bc in self.bc_only]


class Barcodes:
    
    def __init__(self, n_sites=None, custom_enzymes=None):
        self.n_sites = n_sites
        self.sites = {}
        # default site enzymes (Assembly Plasmid):
        self.site_enzymes = ['EcoRI', 'BamHI', 'NheI', 'XhoI', 'PlutI', 'AgeI', 'MluI']
        self._validate_enzymes(custom_enzymes)
        self._build_sites(n_sites)
        self.positions = 1
        self.overhangs = ['TGCC', 'GCAA', 'AGGA', 'TGTG',
                          'GAGC', 'ATTC', 'ATAG'] # allows <= 8 internal positions
        self.avoid_seqs = [] # rs and sequences to avoid (sequences, not RE names)
        self._bc_pool = bc_pool.copy()

    def _validate_enzymes(self, enzymes):
        self._custom_enz = False
        if enzymes is not None and not isinstance(enzymes, list): 
            raise TypeError("""
                            Custom enzymes must be a list of enzyme names 
                            (e.g. ['EcoRI', 'BamHI', ...])
                            """)
        elif enzymes:
            check_rs = [x in RS.AllEnzymes for x in enzymes]
            if all(check_rs):
                self.site_enzymes = enzymes # overwrite default site enzymes
                self._custom_enz = True
            else:
                raise ValueError("""
                                 One or more restriction enzymes were not found.
                                 Names are case-sensitive.
                                 """) 
    
    def _build_sites(self, n):
        if not isinstance(self.n_sites, int):
            raise TypeError('n_sites must be provided and a valid integer.')
        if self._custom_enz: # if user provides custom enzyme list
            if len(self.site_enzymes) != (n+1):
                raise ValueError("""
                                 Incorrect combination of n_sites with provided 
                                 number of restriction sites. Number of restriction
                                 sites must be n_sites + 1.
                                 """)
        else:
            if n > 6:
                raise ValueError('For the default restriction sites, n_sites must <= 6')
            self.site_enzymes = self.site_enzymes[0:(n+1)]
        
        for i in range(self.n_sites):
            site_id = i + 1
            self.sites[site_id] = Site(site_id, self.site_enzymes[i], self.site_enzymes[i+1])
    
    def add_positions(self, *, n_per_site):
        print(f"Method to add {n_per_site} positions within sites")
        # method to add n_per_site number of positions within each site
        # use default sticky overhangs for middle positions
        # allow users to customize this eventually

    def print_structure(self):
        print('----- ',end='')
        for i in range(self.n_sites):
            current_site = self.sites[i+1]
            if (i+1) == 1:
                print(
                    f'{current_site.left_enzyme} - \033[31mSITE{i+1}\033[0m - '
                    f'{current_site.right_enzyme} - ', end='')
            elif i == max(range(self.n_sites)):
                print(f'\033[31mSITE{i+1}\033[0m - {current_site.right_enzyme} -----')
            else:
                print(f'\033[31mSITE{i+1}\033[0m - {current_site.right_enzyme} - ',
                      end='')
    
    def _draw_bc(self, bc_len):
        n_segs = math.ceil(bc_len / 60)
        segments = []
        for _ in range(n_segs):
            index = random.randint(0, len(self._bc_pool)-1)
            segments.append(self._bc_pool.pop(index))
        return ''.join(segments)[:bc_len]

    def generate_barcodes(self, bc_len, n_barcodes):
        print('Generating barcodes...')
        for i in range(self.n_sites):
            site_id = i + 1
            site = self.sites[site_id]
            site.bc_only = []
            for _ in range(n_barcodes):
                site.bc_only.append(self._draw_bc(bc_len))
            site.bc = site.bc_only.copy()

    def add_fixed_sequence(self, seq, site, side):
        if side not in ['left', 'right']:
            raise ValueError('Side must be either "left" or "right".')
        if site not in self.sites:
            raise ValueError(f'Site {site} does not exist.')
        self.sites[site].add_fixed_sequence(seq, side)

    def view(self):
        cols = ['site', 'position', 'barcode', 'barcode_assembled']
        rows = []
        for site in self.sites.values():
            for raw, assembled in zip(site.bc_only, site.bc):
                rows.append({
                    'site': site.id,
                    'position': 'A',
                    'barcode': raw,
                    'barcode_assembled': assembled,
                })
        if not rows:
            return pd.DataFrame(columns=cols)
        return pd.DataFrame(rows).copy()
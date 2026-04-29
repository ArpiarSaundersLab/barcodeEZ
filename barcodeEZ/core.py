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

_POSITION_LABELS = list('ABCDEFGH')


class Site:

    def __init__(self, site_id, left_enzyme, right_enzyme):
        self.id = site_id
        self.left_enzyme = left_enzyme
        self.right_enzyme = right_enzyme
        self.fixed_left = ''
        self.fixed_right = ''
        # positions dict: label -> {'bc_only': [], 'bc': [], 'left_oh': str|None, 'right_oh': str|None}
        # left_oh/right_oh are None at enzyme boundaries, 4bp strings for internal overhangs
        self.positions = {'A': {'bc_only': [], 'bc': [], 'left_oh': None, 'right_oh': None}}

    def __repr__(self):
        n_pos = len(self.positions)
        n_bc = sum(len(p['bc']) for p in self.positions.values())
        bc_len = next((len(p['bc'][0]) for p in self.positions.values() if p['bc']), 0)
        return (f"Site {self.id}: {self.left_enzyme} -- SITE{self.id} -- {self.right_enzyme} "
                f"| {n_pos} position(s), {n_bc} barcodes, {bc_len}bp")

    def add_fixed_sequence(self, seq, side):
        if side == 'left':
            self.fixed_left = seq
        else:
            self.fixed_right = seq
        self._rebuild_bc()

    def _rebuild_bc(self):
        for pos_data in self.positions.values():
            pos_data['bc'] = [
                self.fixed_left + bc + self.fixed_right for bc in pos_data['bc_only']
            ]


class Barcodes:
    
    def __init__(self, n_sites=None, custom_enzymes=None):
        self.n_sites = n_sites
        self.sites = {}
        # default site enzymes (Assembly Plasmid):
        self.site_enzymes = ['EcoRI', 'BamHI', 'NheI', 'XhoI', 'PlutI', 'AgeI', 'MluI']
        self._validate_enzymes(custom_enzymes)
        self._build_sites(n_sites)
        self.overhangs = ['TGCC', 'GCAA', 'AGGA', 'TGTG',
                          'GAGC', 'ATTC', 'ATAG'] # 7 overhangs support up to 8 positions
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
        if not isinstance(n_per_site, int) or n_per_site < 1:
            raise TypeError('n_per_site must be a positive integer.')
        if n_per_site > 8:
            raise ValueError('n_per_site must be <= 8.')
        labels = _POSITION_LABELS[:n_per_site]
        for site in self.sites.values():
            site.positions = {}
            for i, label in enumerate(labels):
                site.positions[label] = {
                    'bc_only': [],
                    'bc': [],
                    'left_oh': self.overhangs[i - 1] if i > 0 else None,
                    'right_oh': self.overhangs[i] if i < n_per_site - 1 else None,
                }

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
        for site in self.sites.values():
            for pos_data in site.positions.values():
                pos_data['bc_only'] = [self._draw_bc(bc_len) for _ in range(n_barcodes)]
                pos_data['bc'] = pos_data['bc_only'].copy()

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
            for pos_label, pos_data in site.positions.items():
                for raw, assembled in zip(pos_data['bc_only'], pos_data['bc']):
                    rows.append({
                        'site': site.id,
                        'position': pos_label,
                        'barcode': raw,
                        'barcode_assembled': assembled,
                    })
        if not rows:
            return pd.DataFrame(columns=cols)
        return pd.DataFrame(rows).copy()
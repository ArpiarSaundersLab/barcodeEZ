import Bio.Restriction as RS
import random
import gzip
from importlib.resources import files

with gzip.open(str(files('barcodeEZ.bc_gen').joinpath('20k_barcodes_60mers.fa.gz')), 'rt') as f:
    bc_pool = []
    for i,line in enumerate(f):
        line = line.strip()
        if i % 2 == 1:
            bc_pool.append(line)

class Barcodes:
    
    def __init__(self, n_sites=None, custom_enzymes=None):
        self.n_sites = n_sites
        self.sites = {}
        # default site enzymes (Assembly Plasmid):
        self.site_enzymes = ['EcoRI', 'BamHI', 'NheI', 'XhoI', 'PlutI', 'AgeI', 'MluI']
        self._validate_enzymes(custom_enzymes)
        self._build_sites(n_sites)
        self.positions = 1
        self.overhangs = ['TGCC', 'GCAA', 'AGGA'] # add 6 total for each site
        self.avoid_seqs = [] # rs and sequences to avoid
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
            site_id = i+1
            self.sites[site_id] = {'left': self.site_enzymes[i],
                                   'right': self.site_enzymes[i+1],
                                   'bc': None}
    
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
                    f'{current_site["left"]} - \033[31mSITE{i+1}\033[0m - '
                    f'{current_site["right"]} - ', end='')
            elif i == max(range(self.n_sites)):
                print(f'\033[31mSITE{i+1}\033[0m - {current_site["right"]} -----')
            else:
                print(f'\033[31mSITE{i+1}\033[0m - {current_site["right"]} - ',
                      end='')
    
    def generate_barcodes(self, bc_len, n_barcodes):
        print('Generating barcodes...')
        # generate barcodes for each site
        if bc_len > 60:
            # implement concatenation of bc_pool sequences for longer barcodes
            raise ValueError('Barcode length must be 60 or less.') # avoid for now
        for i in range(self.n_sites):
            site_id = i + 1
            self.sites[site_id]['bc_only'] = []
            for j in range(n_barcodes):
                index = random.randint(0, len(self._bc_pool)-1)
                bc_insert = self._bc_pool.pop(index)[0:bc_len] # remove from pool once used
                self.sites[site_id]['bc_only'].append(bc_insert)
            self.sites[site_id]['bc'] = self.sites[site_id]['bc_only'].copy()

    def add_fixed_sequence(self, seq, site, side):
        if side not in ['left', 'right']:
            raise ValueError('Side must be either "left" or "right".')
        if site not in self.sites:
            raise ValueError(f'Site {site} does not exist.')
        self.sites[site].setdefault('fixed_left', '')
        self.sites[site].setdefault('fixed_right', '')
        if side == 'left':
            self.sites[site]['fixed_left'] = seq
        else:
            self.sites[site]['fixed_right'] = seq
        self._rebuild_bc(site)

    def _rebuild_bc(self, site):
        left = self.sites[site].get('fixed_left', '')
        right = self.sites[site].get('fixed_right', '')
        self.sites[site]['bc'] = [
            left + bc + right for bc in self.sites[site]['bc_only']
        ]
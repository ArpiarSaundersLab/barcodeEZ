import Bio.Restriction as RS
import random

with open('./100k_barcodes.txt', 'r') as f:
    bc_pool = []
    for i,line in enumerate(f):
        line = line.strip()
        if i % 2 == 1:
            bc_pool.append(line)

class barcodeEZ:
    
    def __init__(self, n_sites=None, custom_enzymes=None):
        self.n_sites = n_sites
        self.sites = {}
        # default site enzymes (a.31):
        self.site_enzymes = ['EcoRI', 'BamHI', 'NheI', 'XhoI', 'PlutI', 'AgeI', 'MluI']
        self._validate_enzymes(custom_enzymes)
        self._build_sites(n_sites)
        self.positions = 1
        self.overhangs = ['TGCC', 'GCAA', 'AGGA'] # add 6 total for each site
        self.avoid_seqs = [] # rs and sequences to avoid
        self.bc_pool = bc_pool
    
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
            # implement concatenation for longer barcodes
            raise ValueError('Barcode length must be 60 or less.')
        for i in range(self.n_sites):
            site_id = i + 1
            self.sites[site_id]['bc'] = []
            for j in range(n_barcodes):
                index = random.randint(0, len(self.bc_pool)-1)
                bc_insert = self.bc_pool.pop(index)[0:bc_len]
                self.sites[site_id]['bc'].append(bc_insert) 

    def add_fixed_sequence(self, seq, site, side):
        ## add logic so that users don't add fixed_sequence multiple times
        ## if this function has already been called for left:
        ## use original bc to append to
        ## Same logic for right (needs to be separate in case user calls this twice for left and right)
        if side not in ['left', 'right']:
            raise ValueError('Side must be either "left" or "right".')
        if site not in self.sites:
            raise ValueError(f'Site {site} does not exist.')
        update_BCs = []
        for bc in self.sites[site]['bc']:
            if side == 'left':
                update_BCs.append(seq + bc)
            elif side == 'right':
                update_BCs.append(bc + seq)
        self.sites[site]['bc_only'] = self.sites[site]['bc']
        self.sites[site]['bc'] = update_BCs
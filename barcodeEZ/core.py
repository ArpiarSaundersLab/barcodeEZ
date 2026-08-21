import math
import Bio.Restriction as RS
import random
import gzip
from importlib.resources import files
import pandas as pd

with gzip.open(str(files('barcodeEZ').joinpath('corpus/20k_barcodes_60mers.fa.gz')), 'rt') as f:
    bc_pool = []
    for i,line in enumerate(f):
        line = line.strip()
        if i % 2 == 1:
            bc_pool.append(line)

_POSITION_LABELS = list('ABCDEFGH')

_DEFAULT_AVOID_ENZYMES = [
    'BsiWI', 'MreI', 'FseI', 'EcoRI', 'AvrII', 'BamHI', 'KpnI', 'NheI',
    'PciI', 'XhoI', 'SpeI', 'PluTI', 'NotI', 'AgeI', 'AsiSI', 'MluI',
    'SbfI', 'MauBI',
]


_FASTA_EXTENSIONS = {'.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz'}


def _expand_motif_list(avoid):
    """Expand any FASTA file paths in avoid into their constituent sequences."""
    expanded = []
    for item in avoid:
        lower = item.lower()
        is_fasta = any(lower.endswith(ext) for ext in _FASTA_EXTENSIONS)
        if is_fasta:
            from pathlib import Path
            p = Path(item)
            if not p.exists():
                print(f"Error: FASTA file not found: {item}")
                continue
            opener = gzip.open if lower.endswith('.gz') else open
            with opener(item, 'rt') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('>'):
                        expanded.append(line)
        else:
            expanded.append(item)
    return expanded


def _resolve_motifs(names_or_seqs):
    """Return {recognition_sequence: label} for a list of enzyme names or raw sequences."""
    result = {}
    for item in names_or_seqs:
        if item in RS.AllEnzymes:
            result[getattr(RS, item).site] = item
        elif all(c in 'ACGT' for c in item.upper()):
            result[item.upper()] = item
        else:
            raise ValueError(f"'{item}' is not a recognized enzyme name or valid DNA sequence.")
    return result


def _rc(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def _enzyme_parts(enzyme_name):
    """Return (cut_right, cut_left) for a restriction enzyme.
    cut_right = sequence to the right of the cut (insert side when enzyme is on the left).
    cut_left  = base(s) to the left of the cut (vector side).
    """
    enz = getattr(RS, enzyme_name)
    return enz.site[enz.fst5:], enz.site[:enz.fst5]


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
    """Design a combinatorial barcode library for molecular assembly.

    The only object users need. Methods that configure the library return
    ``self``, so steps can be chained or called sequentially on the same
    variable -- both styles are equivalent.

    Parameters
    ----------
    n_sites : int
        Number of barcoded sites. With the default enzyme panel, must be <= 6.
    custom_enzymes : list of str, optional
        Restriction-enzyme panel defining site boundaries. Must contain exactly
        ``n_sites + 1`` enzyme names (case-sensitive, Biopython-recognized).
        Overrides the default panel (EcoRI, BamHI, NheI, XhoI, PluTI, AgeI, MluI).

    Raises
    ------
    TypeError
        If ``n_sites`` is not an integer.
    ValueError
        If ``n_sites > 6`` with default enzymes, the custom panel length is not
        ``n_sites + 1``, or an enzyme name is not recognized by Biopython.

    Examples
    --------
    Chained:

    >>> b = (Barcodes(n_sites=2)
    ...         .add_positions(n_per_site=4)
    ...         .generate_barcodes(bc_len=25, n_barcodes=96)
    ...         .validate())

    Sequential (equivalent):

    >>> b = Barcodes(n_sites=2)
    >>> b.add_positions(n_per_site=4)
    >>> b.generate_barcodes(bc_len=25, n_barcodes=96)
    >>> b.validate()

    Export:

    >>> b.view()
    >>> b.write_order_form('order.csv')
    """

    def __init__(self, n_sites: int = None, custom_enzymes: list[str] | None = None):
        """Initialize."""
        self.n_sites = n_sites
        self.sites = {}
        # default site enzymes (Assembly Plasmid):
        self.site_enzymes = ['EcoRI', 'BamHI', 'NheI', 'XhoI', 'PluTI', 'AgeI', 'MluI']
        self._validate_enzymes(custom_enzymes)
        self._build_sites(n_sites)
        self.overhangs = ['TGCC', 'GCAA', 'AGGA', 'TGTG',
                          'GAGC', 'ATTC', 'ATAG'] # 7 overhangs support up to 8 positions
        self.avoid_seqs = [] # rs and sequences to avoid (sequences, not RE names)
        self._validated = False
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
    
    def add_positions(self, *, n_per_site: int) -> 'Barcodes':
        """Add multiple barcode positions per site for combinatorial designs.

        Call before ``generate_barcodes()`` -- resets any existing barcodes.
        Positions are labelled A-H; adjacent positions are joined by
        automatically assigned 4 bp overhangs.

        Parameters
        ----------
        n_per_site : int
            Number of positions per site (1-8). Keyword-only.

        Raises
        ------
        TypeError
            If ``n_per_site`` is not a positive integer.
        ValueError
            If ``n_per_site > 8``.
        """
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
        return self

    def print_structure(self) -> None:
        """Print the site-and-enzyme layout to stdout."""
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

    def generate_barcodes(self, bc_len: int, n_barcodes: int) -> 'Barcodes':
        """Draw barcodes from the corpus and assemble forward/reverse oligos.

        Barcodes are drawn without replacement from the ~20,000-member corpus,
        so all barcodes in the library are unique. Values of ``bc_len > 60``
        are built by concatenating multiple 60-mer corpus sequences. Resets
        the validated status of the library.

        Parameters
        ----------
        bc_len : int
            Length of each barcode in bp.
        n_barcodes : int
            Number of barcodes to generate per position.
        """
        print('Generating barcodes...')
        self._bc_len = bc_len
        self._validated = False
        for site in self.sites.values():
            for pos_data in site.positions.values():
                pos_data['bc_only'] = [self._draw_bc(bc_len) for _ in range(n_barcodes)]
                pos_data['bc'] = pos_data['bc_only'].copy()
        self._generate_oligos()
        return self

    def _generate_oligos(self):
        for site in self.sites.values():
            for pos_data in site.positions.values():
                if pos_data['left_oh'] is None:
                    f_pre, r_suf = _enzyme_parts(site.left_enzyme)
                else:
                    f_pre, r_suf = pos_data['left_oh'], ''

                if pos_data['right_oh'] is None:
                    r_pre, f_suf = _enzyme_parts(site.right_enzyme)
                else:
                    f_suf, r_pre = '', _rc(pos_data['right_oh'])

                pos_data['forward_oligos'] = [f_pre + bc + f_suf for bc in pos_data['bc']]
                pos_data['reverse_oligos'] = [r_pre + _rc(bc) + r_suf for bc in pos_data['bc']]

    def add_fixed_sequence(self, seq: str, site: int, side: str) -> 'Barcodes':
        """Attach a constant flanking sequence to one side of a site.

        Applies to all positions in the site, immediately adjacent to the
        barcode. Calling again for the same site and side overwrites the
        previous sequence. Resets the validated status.

        Parameters
        ----------
        seq : str
            Fixed DNA sequence to attach.
        site : int
            Site number to modify (1-indexed).
        side : str
            ``'left'`` or ``'right'``.

        Raises
        ------
        ValueError
            If ``side`` is not ``'left'`` or ``'right'``, or ``site`` does not exist.
        """
        if side not in ['left', 'right']:
            raise ValueError('Side must be either "left" or "right".')
        if site not in self.sites:
            raise ValueError(f'Site {site} does not exist.')
        self.sites[site].add_fixed_sequence(seq, side)
        self._validated = False
        self._generate_oligos()
        return self

    def validate(self, ignore_defaults: bool = False, motifs: list | None = None) -> 'Barcodes':
        """Screen barcodes for unwanted motifs and replace contaminated ones.

        Each assembled barcode is checked for the presence of unwanted
        sequences. Contaminated barcodes are replaced with clean draws from
        the corpus. Prints the number of replacements made, or a confirmation
        that no motifs were found. Marks the library as validated.

        Each item in ``motifs`` may be a restriction-enzyme name (e.g.
        ``'BsaI'``), a raw DNA sequence (e.g. ``'TATAAA'``), or a path to a
        FASTA file (``.fa``, ``.fasta``, ``.fna``, optionally ``.gz``).

        The default motif panel (used unless ``ignore_defaults=True``):
        BsiWI, MreI, FseI, EcoRI, AvrII, BamHI, KpnI, NheI, PciI, XhoI,
        SpeI, PluTI, NotI, AgeI, AsiSI, MluI, SbfI, MauBI.

        Parameters
        ----------
        ignore_defaults : bool, optional
            If True, skip the default restriction-enzyme panel. Default is False.
        motifs : list, optional
            Additional motifs to screen. Items may be enzyme names, raw DNA
            sequences, or FASTA file paths (mixed freely). Default is None.

        Raises
        ------
        ValueError
            If a ``motifs`` item is neither a recognized enzyme name nor a
            valid DNA sequence.
        RuntimeError
            If the corpus is exhausted before a clean replacement is found.
            Usually means the motif comes from a fixed sequence or overhang.
        """
        expanded = _expand_motif_list(motifs or [])
        names = ([] if ignore_defaults else _DEFAULT_AVOID_ENZYMES) + expanded
        self.avoid_seqs = _resolve_motifs(names)
        n_replaced = self._validate_and_replace()
        self._generate_oligos()
        if n_replaced:
            print(f'Validation complete: {n_replaced} barcode(s) replaced.')
        else:
            print('Validation complete: no unwanted motifs found.')
        self._validated = True
        return self

    def _validate_and_replace(self):
        n_replaced = 0
        for site in self.sites.values():
            fixed_left = site.fixed_left
            fixed_right = site.fixed_right
            for pos_data in site.positions.values():
                for i, bc in enumerate(pos_data['bc']):
                    current_bc = bc
                    hits = [label for motif, label in self.avoid_seqs.items()
                            if motif in current_bc]
                    while hits:
                        print(f"Contamination ({', '.join(hits)}): {current_bc!r} — replacing barcode.")
                        if not self._bc_pool:
                            raise RuntimeError(
                                'bc_pool exhausted; could not find clean replacement barcode. '
                                'The unwanted motif may be coming from a fixed sequence or overhang.'
                            )
                        new_bc_only = self._draw_bc(self._bc_len)
                        new_bc = fixed_left + new_bc_only + fixed_right
                        hits = [label for motif, label in self.avoid_seqs.items()
                                if motif in new_bc]
                        current_bc = new_bc
                        if not hits:
                            pos_data['bc_only'][i] = new_bc_only
                            pos_data['bc'][i] = new_bc
                            n_replaced += 1
        return n_replaced

    def view(self) -> 'pd.DataFrame':
        """Return the library as a pandas DataFrame.

        Columns appear conditionally based on what has been configured:

        - ``site``, ``position``, ``barcode``: always present.
        - ``barcode_with_fixed_seq``: after ``add_fixed_sequence()``.
        - ``forward_oligo``, ``reverse_oligo``: after ``generate_barcodes()``.

        Returns
        -------
        pandas.DataFrame
            A copy of the library. Modifying it does not affect the library.
        """
        has_fixed = any(
            site.fixed_left or site.fixed_right
            for site in self.sites.values()
        )
        has_oligos = any(
            'forward_oligos' in pos_data
            for site in self.sites.values()
            for pos_data in site.positions.values()
        )
        cols = ['site', 'position', 'barcode']
        if has_fixed:
            cols += ['barcode_with_fixed_seq']
        if has_oligos:
            cols += ['forward_oligo', 'reverse_oligo']
        rows = []
        for site in self.sites.values():
            for pos_label, pos_data in site.positions.items():
                fwd = pos_data.get('forward_oligos', [])
                rev = pos_data.get('reverse_oligos', [])
                for j, (raw, assembled) in enumerate(zip(pos_data['bc_only'], pos_data['bc'])):
                    row = {
                        'site': site.id,
                        'position': pos_label,
                        'barcode': raw,
                    }
                    if has_fixed:
                        row['barcode_with_fixed_seq'] = assembled
                    if has_oligos:
                        row['forward_oligo'] = fwd[j] if j < len(fwd) else None
                        row['reverse_oligo'] = rev[j] if j < len(rev) else None
                    rows.append(row)
        if not rows:
            return pd.DataFrame(columns=cols)
        return pd.DataFrame(rows).copy()

    def write_order_form(self, file: str, metadata: bool = False) -> None:
        """Export the oligo pool as a CSV ready for synthesis.

        Writes ``opool_name`` and ``oligo_sequence`` columns. Each barcode
        produces two rows -- a forward oligo (``site{N}_f``) and a reverse
        oligo (``site{N}_r``). Prints a warning if ``validate()`` has not
        been run since the last design change.

        Parameters
        ----------
        file : str
            Output CSV path.
        metadata : bool, optional
            If True, include ``site``, ``position``, and ``barcode`` columns
            alongside the oligo sequences. Default is False.

        Raises
        ------
        RuntimeError
            If no oligos exist (call ``generate_barcodes()`` first).
        """
        df = self.view()
        if 'forward_oligo' not in df.columns:
            raise RuntimeError('No oligo sequences found. Run generate_barcodes() first.')
        if not self._validated:
            print('Warning: validate() has not been run. Unwanted restriction sites or motifs '
                  'may be present in the library.')
        order_rows = []
        for _, row in df.iterrows():
            name = f"site{row['site']}"
            fwd = {'opool_name': f'{name}_f', 'oligo_sequence': row['forward_oligo']}
            rev = {'opool_name': f'{name}_r', 'oligo_sequence': row['reverse_oligo']}
            if metadata:
                meta = {'site': row['site'], 'position': row['position'], 'barcode': row['barcode']}
                fwd = {**fwd, **meta}
                rev = {**rev, **meta}
            order_rows.append(fwd)
            order_rows.append(rev)
        pd.DataFrame(order_rows).to_csv(file, index=False)
import pytest
from barcodeEZ import Barcodes

# SDF2: 3-site library, position A only, 256 barcodes/site, fixed left sequence per site
# SDF3: 4-site library, positions A-D, 3 barcodes/site/position — requires add_positions()


class TestSDF2:

    def test_site_count(self):
        b = Barcodes(n_sites=3)
        assert len(b.sites) == 3

    def test_default_enzyme_structure(self):
        b = Barcodes(n_sites=3)
        assert b.sites[1].left_enzyme == 'EcoRI'
        assert b.sites[1].right_enzyme == 'BamHI'
        assert b.sites[2].left_enzyme == 'BamHI'
        assert b.sites[2].right_enzyme == 'NheI'
        assert b.sites[3].left_enzyme == 'NheI'
        assert b.sites[3].right_enzyme == 'XhoI'

    def test_fixed_sequence_assembly(self, sdf2_ref, monkeypatch):
        import random
        fwd = sdf2_ref[sdf2_ref['opool_name'].str.endswith('_f')]

        # Sort into draw order: site 1->3, then by variable_sequence_number
        ordered = fwd.sort_values(['site_number', 'variable_sequence_number'])
        pool = [bc.ljust(60, 'A') for bc in ordered['barcode_only']]

        b = Barcodes(n_sites=3)
        b._bc_pool = pool
        monkeypatch.setattr(random, 'randint', lambda a, b: 0)
        b.generate_barcodes(25, 256)

        for site_num in [1, 2, 3]:
            site_rows = fwd[fwd['site_number'] == site_num].sort_values('variable_sequence_number')
            b.add_fixed_sequence(site_rows['fixed_sequence_left'].iloc[0], site=site_num, side='left')

        df = b.view()
        for site_num in [1, 2, 3]:
            site_rows = fwd[fwd['site_number'] == site_num].sort_values('variable_sequence_number')
            generated = list(df[df['site'] == site_num]['barcode_assembled'])
            expected = list(site_rows['full_barcode_sequence'])
            assert generated == expected, f"Mismatch at site {site_num}"

    def test_barcode_and_fixed_seq_lengths(self, sdf2_ref):
        fwd = sdf2_ref[sdf2_ref['opool_name'].str.endswith('_f')]
        assert (fwd['barcode_only'].str.len() == 25).all()
        assert (fwd['fixed_sequence_left'].str.len() == 25).all()
        assert (fwd['full_barcode_sequence'].str.len() == 50).all()

    def test_barcode_count_per_site(self, sdf2_ref):
        fwd = sdf2_ref[sdf2_ref['opool_name'].str.endswith('_f')]
        counts = fwd.groupby('site_number').size()
        assert (counts == 256).all()

    def test_positions(self, sdf2_ref):
        assert sdf2_ref['site_position'].unique().tolist() == ['A']

    @pytest.mark.skip(reason="requires oligo_sequence generation implementation")
    def test_forward_oligo_sequence(self, sdf2_ref):
        pass

    @pytest.mark.skip(reason="requires oligo_sequence generation implementation")
    def test_reverse_oligo_sequence(self, sdf2_ref):
        pass


class TestSDF3:

    def test_site_count(self):
        b = Barcodes(n_sites=4)
        assert len(b.sites) == 4

    def test_default_enzyme_structure(self):
        b = Barcodes(n_sites=4)
        assert b.sites[1].left_enzyme == 'EcoRI'
        assert b.sites[4].right_enzyme == 'PlutI'

    def test_barcode_lengths(self, sdf3_ref):
        fwd = sdf3_ref[sdf3_ref['opool_name'].str.endswith('_f')]
        assert (fwd['barcode_only'].str.len() == 30).all()

    def test_barcode_count_per_site_per_position(self, sdf3_ref):
        fwd = sdf3_ref[sdf3_ref['opool_name'].str.endswith('_f')]
        counts = fwd.groupby(['site_number', 'site_position']).size()
        assert (counts == 3).all()

    def test_positions_A_through_D(self, sdf3_ref, monkeypatch):
        import random
        fwd = sdf3_ref[sdf3_ref['opool_name'].str.endswith('_f')]

        # Sort into the draw order: site -> position -> variant, matching
        # the iteration order of generate_barcodes() (sites 1-4, positions A-D)
        ordered = fwd.sort_values(['site_number', 'site_position', 'variable_sequence_number'])
        # Pad each 30bp barcode to 60bp so it fits the pool format; _draw_bc slices to 30
        pool = [bc.ljust(60, 'A') for bc in ordered['barcode_only']]

        b = Barcodes(n_sites=4)
        b.add_positions(n_per_site=4)
        b._bc_pool = pool
        monkeypatch.setattr(random, 'randint', lambda a, b: 0)
        b.generate_barcodes(30, 3)

        df = b.view()
        for site_num in [1, 2, 3, 4]:
            for pos in ['A', 'B', 'C', 'D']:
                ref = list(fwd[(fwd['site_number'] == site_num) &
                               (fwd['site_position'] == pos)]
                           .sort_values('variable_sequence_number')['barcode_only'])
                generated = list(df[(df['site'] == site_num) &
                                    (df['position'] == pos)]['barcode'])
                assert generated == ref, f"Mismatch at site {site_num} position {pos}"

    @pytest.mark.skip(reason="requires oligo_sequence generation implementation")
    def test_forward_oligo_sequence(self, sdf3_ref):
        pass

    @pytest.mark.skip(reason="requires oligo_sequence generation implementation")
    def test_reverse_oligo_sequence(self, sdf3_ref):
        pass

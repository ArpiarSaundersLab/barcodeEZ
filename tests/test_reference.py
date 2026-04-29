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

    def test_fixed_sequence_assembly(self, sdf2_ref):
        """
        Inject known barcodes from reference data, apply fixed sequences,
        and verify assembled barcodes match full_barcode_sequence.
        """
        fwd = sdf2_ref[sdf2_ref['opool_name'].str.endswith('_f')]
        b = Barcodes(n_sites=3)
        for site_num in [1, 2, 3]:
            site_rows = fwd[fwd['site_number'] == site_num]
            fixed_seq = site_rows['fixed_sequence_left'].iloc[0]
            barcodes = list(site_rows['barcode_only'])
            b.sites[site_num].bc_only = barcodes
            b.sites[site_num].bc = barcodes.copy()
            b.add_fixed_sequence(fixed_seq, site=site_num, side='left')
            expected = list(site_rows['full_barcode_sequence'])
            assert b.sites[site_num].bc == expected

    def test_barcode_and_fixed_seq_lengths(self, sdf2_ref):
        fwd = sdf2_ref[sdf2_ref['opool_name'].str.endswith('_f')]
        assert (fwd['barcode_only'].str.len() == 25).all()
        assert (fwd['fixed_sequence_left'].str.len() == 25).all()
        assert (fwd['full_barcode_sequence'].str.len() == 50).all()

    def test_barcode_count_per_site(self, sdf2_ref):
        fwd = sdf2_ref[sdf2_ref['opool_name'].str.endswith('_f')]
        counts = fwd.groupby('site_number').size()
        assert (counts == 256).all()

    @pytest.mark.skip(reason="requires add_positions() implementation")
    def test_positions(self, sdf2_ref):
        pass

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

    @pytest.mark.skip(reason="requires add_positions() implementation")
    def test_positions_A_through_D(self, sdf3_ref):
        pass

    @pytest.mark.skip(reason="requires oligo_sequence generation implementation")
    def test_forward_oligo_sequence(self, sdf3_ref):
        pass

    @pytest.mark.skip(reason="requires oligo_sequence generation implementation")
    def test_reverse_oligo_sequence(self, sdf3_ref):
        pass

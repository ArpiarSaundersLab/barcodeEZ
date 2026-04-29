import pytest
from barcodeEZ import Barcodes


class TestInit:

    def test_correct_site_count(self):
        b = Barcodes(n_sites=3)
        assert len(b.sites) == 3

    def test_default_enzymes_assigned(self):
        b = Barcodes(n_sites=3)
        assert b.sites[1].left_enzyme == 'EcoRI'
        assert b.sites[1].right_enzyme == 'BamHI'
        assert b.sites[2].left_enzyme == 'BamHI'
        assert b.sites[2].right_enzyme == 'NheI'
        assert b.sites[3].right_enzyme == 'XhoI'

    def test_custom_enzymes(self):
        b = Barcodes(n_sites=2, custom_enzymes=['EcoRI', 'BamHI', 'NheI'])
        assert b.sites[1].left_enzyme == 'EcoRI'
        assert b.sites[2].right_enzyme == 'NheI'

    def test_invalid_enzyme_raises(self):
        with pytest.raises(ValueError):
            Barcodes(n_sites=2, custom_enzymes=['EcoRI', 'FakeEnzyme', 'NheI'])

    def test_n_sites_not_int_raises(self):
        with pytest.raises(TypeError):
            Barcodes(n_sites='two')

    def test_n_sites_exceeds_default_raises(self):
        with pytest.raises(ValueError):
            Barcodes(n_sites=7)

    def test_custom_enzyme_count_mismatch_raises(self):
        with pytest.raises(ValueError):
            Barcodes(n_sites=2, custom_enzymes=['EcoRI', 'BamHI'])  # needs 3 enzymes for 2 sites


class TestGenerateBarcodes:

    def test_correct_barcode_count_per_site(self):
        b = Barcodes(n_sites=2)
        b.generate_barcodes(30, 10)
        assert len(b.sites[1].bc_only) == 10
        assert len(b.sites[2].bc_only) == 10

    def test_correct_barcode_length(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        assert all(len(bc) == 30 for bc in b.sites[1].bc_only)

    def test_bc_matches_bc_only_before_fixed_seq(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        assert b.sites[1].bc == b.sites[1].bc_only

    def test_long_barcode_over_60(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(90, 5)
        assert all(len(bc) == 90 for bc in b.sites[1].bc_only)

    def test_pool_depleted_after_generate(self):
        b = Barcodes(n_sites=1)
        pool_before = len(b._bc_pool)
        b.generate_barcodes(30, 10)
        assert len(b._bc_pool) == pool_before - 10

    def test_pool_depleted_more_for_long_barcodes(self):
        b = Barcodes(n_sites=1)
        pool_before = len(b._bc_pool)
        b.generate_barcodes(90, 10)  # each barcode needs ceil(90/60)=2 pool draws
        assert len(b._bc_pool) == pool_before - 20

    def test_barcodes_unique_across_sites(self):
        b = Barcodes(n_sites=2)
        b.generate_barcodes(30, 5)
        all_barcodes = b.sites[1].bc_only + b.sites[2].bc_only
        assert len(all_barcodes) == len(set(all_barcodes))


class TestAddFixedSequence:

    def test_fixed_left_prepended(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        raw = b.sites[1].bc_only.copy()
        b.add_fixed_sequence('AAAA', site=1, side='left')
        for original, assembled in zip(raw, b.sites[1].bc):
            assert assembled == 'AAAA' + original

    def test_fixed_right_appended(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        raw = b.sites[1].bc_only.copy()
        b.add_fixed_sequence('TTTT', site=1, side='right')
        for original, assembled in zip(raw, b.sites[1].bc):
            assert assembled == original + 'TTTT'

    def test_fixed_both_sides(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        raw = b.sites[1].bc_only.copy()
        b.add_fixed_sequence('AAAA', site=1, side='left')
        b.add_fixed_sequence('TTTT', site=1, side='right')
        for original, assembled in zip(raw, b.sites[1].bc):
            assert assembled == 'AAAA' + original + 'TTTT'

    def test_fixed_sequence_overwrite(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        b.add_fixed_sequence('AAAA', site=1, side='left')
        b.add_fixed_sequence('CCCC', site=1, side='left')
        assert all(bc.startswith('CCCC') for bc in b.sites[1].bc)
        assert not any(bc.startswith('AAAA') for bc in b.sites[1].bc)

    def test_bc_only_unchanged_after_fixed_seq(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        raw = b.sites[1].bc_only.copy()
        b.add_fixed_sequence('AAAA', site=1, side='left')
        assert b.sites[1].bc_only == raw

    def test_invalid_side_raises(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 3)
        with pytest.raises(ValueError):
            b.add_fixed_sequence('AAAA', site=1, side='middle')

    def test_invalid_site_raises(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 3)
        with pytest.raises(ValueError):
            b.add_fixed_sequence('AAAA', site=99, side='left')


class TestView:

    def test_columns(self):
        b = Barcodes(n_sites=2)
        b.generate_barcodes(30, 5)
        df = b.view()
        assert list(df.columns) == ['site', 'position', 'barcode', 'barcode_assembled']

    def test_row_count(self):
        b = Barcodes(n_sites=2)
        b.generate_barcodes(30, 5)
        df = b.view()
        assert len(df) == 10

    def test_site_values(self):
        b = Barcodes(n_sites=3)
        b.generate_barcodes(30, 4)
        df = b.view()
        assert set(df['site']) == {1, 2, 3}

    def test_position_defaults_to_A(self):
        b = Barcodes(n_sites=2)
        b.generate_barcodes(30, 5)
        df = b.view()
        assert (df['position'] == 'A').all()

    def test_barcode_matches_bc_only(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        df = b.view()
        assert list(df['barcode']) == b.sites[1].bc_only

    def test_barcode_assembled_reflects_fixed_seq(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        b.add_fixed_sequence('GGGG', site=1, side='left')
        df = b.view()
        assert df['barcode_assembled'].str.startswith('GGGG').all()
        assert not df['barcode'].str.startswith('GGGG').any()

    def test_returns_copy(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        df = b.view()
        df['barcode'] = 'MUTATED'
        assert b.sites[1].bc_only[0] != 'MUTATED'

    def test_empty_before_generate(self):
        b = Barcodes(n_sites=2)
        df = b.view()
        assert list(df.columns) == ['site', 'position', 'barcode', 'barcode_assembled']
        assert len(df) == 0


class TestSiteRepr:

    def test_repr_before_barcodes(self):
        b = Barcodes(n_sites=1)
        r = repr(b.sites[1])
        assert 'Site 1' in r
        assert 'EcoRI' in r
        assert '0 barcodes' in r

    def test_repr_after_barcodes(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        r = repr(b.sites[1])
        assert '5 barcodes' in r
        assert '30bp' in r

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

    def test_default_single_position(self):
        b = Barcodes(n_sites=1)
        assert list(b.sites[1].positions.keys()) == ['A']


class TestAddPositions:

    def test_correct_position_labels(self):
        b = Barcodes(n_sites=2)
        b.add_positions(n_per_site=4)
        assert list(b.sites[1].positions.keys()) == ['A', 'B', 'C', 'D']
        assert list(b.sites[2].positions.keys()) == ['A', 'B', 'C', 'D']

    def test_single_position(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=1)
        assert list(b.sites[1].positions.keys()) == ['A']

    def test_max_eight_positions(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=8)
        assert list(b.sites[1].positions.keys()) == list('ABCDEFGH')

    def test_exceeds_max_raises(self):
        b = Barcodes(n_sites=1)
        with pytest.raises(ValueError):
            b.add_positions(n_per_site=9)

    def test_invalid_type_raises(self):
        b = Barcodes(n_sites=1)
        with pytest.raises(TypeError):
            b.add_positions(n_per_site='four')

    def test_outer_positions_have_no_internal_overhangs(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=4)
        assert b.sites[1].positions['A']['left_oh'] is None
        assert b.sites[1].positions['D']['right_oh'] is None

    def test_internal_overhangs_assigned(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=4)
        assert b.sites[1].positions['A']['right_oh'] == b.overhangs[0]
        assert b.sites[1].positions['B']['left_oh'] == b.overhangs[0]
        assert b.sites[1].positions['B']['right_oh'] == b.overhangs[1]
        assert b.sites[1].positions['C']['left_oh'] == b.overhangs[1]

    def test_overhangs_shared_between_adjacent_positions(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=3)
        positions = b.sites[1].positions
        assert positions['A']['right_oh'] == positions['B']['left_oh']
        assert positions['B']['right_oh'] == positions['C']['left_oh']

    def test_resets_existing_barcodes(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        b.add_positions(n_per_site=3)
        assert b.sites[1].positions['A']['bc_only'] == []


class TestGenerateBarcodes:

    def test_correct_barcode_count_per_site(self):
        b = Barcodes(n_sites=2)
        b.generate_barcodes(30, 10)
        assert len(b.sites[1].positions['A']['bc_only']) == 10
        assert len(b.sites[2].positions['A']['bc_only']) == 10

    def test_correct_barcode_length(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        assert all(len(bc) == 30 for bc in b.sites[1].positions['A']['bc_only'])

    def test_bc_matches_bc_only_before_fixed_seq(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        pos = b.sites[1].positions['A']
        assert pos['bc'] == pos['bc_only']

    def test_long_barcode_over_60(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(90, 5)
        assert all(len(bc) == 90 for bc in b.sites[1].positions['A']['bc_only'])

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

    def test_pool_depleted_across_positions(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=4)
        pool_before = len(b._bc_pool)
        b.generate_barcodes(30, 5)  # 4 positions × 5 barcodes = 20 draws
        assert len(b._bc_pool) == pool_before - 20

    def test_barcodes_unique_across_sites(self):
        b = Barcodes(n_sites=2)
        b.generate_barcodes(30, 5)
        all_barcodes = (b.sites[1].positions['A']['bc_only'] +
                        b.sites[2].positions['A']['bc_only'])
        assert len(all_barcodes) == len(set(all_barcodes))

    def test_barcodes_generated_per_position(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=3)
        b.generate_barcodes(30, 5)
        for label in ['A', 'B', 'C']:
            assert len(b.sites[1].positions[label]['bc_only']) == 5


class TestAddFixedSequence:

    def test_fixed_left_prepended(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        raw = b.sites[1].positions['A']['bc_only'].copy()
        b.add_fixed_sequence('AAAA', site=1, side='left')
        for original, assembled in zip(raw, b.sites[1].positions['A']['bc']):
            assert assembled == 'AAAA' + original

    def test_fixed_right_appended(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        raw = b.sites[1].positions['A']['bc_only'].copy()
        b.add_fixed_sequence('TTTT', site=1, side='right')
        for original, assembled in zip(raw, b.sites[1].positions['A']['bc']):
            assert assembled == original + 'TTTT'

    def test_fixed_both_sides(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        raw = b.sites[1].positions['A']['bc_only'].copy()
        b.add_fixed_sequence('AAAA', site=1, side='left')
        b.add_fixed_sequence('TTTT', site=1, side='right')
        for original, assembled in zip(raw, b.sites[1].positions['A']['bc']):
            assert assembled == 'AAAA' + original + 'TTTT'

    def test_fixed_sequence_overwrite(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        b.add_fixed_sequence('AAAA', site=1, side='left')
        b.add_fixed_sequence('CCCC', site=1, side='left')
        bcs = b.sites[1].positions['A']['bc']
        assert all(bc.startswith('CCCC') for bc in bcs)
        assert not any(bc.startswith('AAAA') for bc in bcs)

    def test_bc_only_unchanged_after_fixed_seq(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        raw = b.sites[1].positions['A']['bc_only'].copy()
        b.add_fixed_sequence('AAAA', site=1, side='left')
        assert b.sites[1].positions['A']['bc_only'] == raw

    def test_fixed_seq_applies_to_all_positions(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=3)
        b.generate_barcodes(25, 5)
        b.add_fixed_sequence('AAAA', site=1, side='left')
        for label in ['A', 'B', 'C']:
            assert all(bc.startswith('AAAA') for bc in b.sites[1].positions[label]['bc'])

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
        assert list(df.columns) == ['site', 'position', 'barcode', 'forward_oligo', 'reverse_oligo']

    def test_row_count_single_position(self):
        b = Barcodes(n_sites=2)
        b.generate_barcodes(30, 5)
        df = b.view()
        assert len(df) == 10

    def test_row_count_multiple_positions(self):
        b = Barcodes(n_sites=2)
        b.add_positions(n_per_site=3)
        b.generate_barcodes(30, 5)
        df = b.view()
        assert len(df) == 30  # 2 sites × 3 positions × 5 barcodes

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

    def test_position_labels_reflect_add_positions(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=4)
        b.generate_barcodes(30, 3)
        df = b.view()
        assert set(df['position']) == {'A', 'B', 'C', 'D'}

    def test_barcode_matches_bc_only(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        df = b.view()
        assert list(df['barcode']) == b.sites[1].positions['A']['bc_only']

    def test_barcode_with_fixed_seq_reflects_fixed_seq(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(25, 5)
        b.add_fixed_sequence('GGGG', site=1, side='left')
        df = b.view()
        assert df['barcode_with_fixed_seq'].str.startswith('GGGG').all()
        assert not df['barcode'].str.startswith('GGGG').any()

    def test_returns_copy(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 5)
        df = b.view()
        df['barcode'] = 'MUTATED'
        assert b.sites[1].positions['A']['bc_only'][0] != 'MUTATED'

    def test_empty_before_generate(self):
        b = Barcodes(n_sites=2)
        df = b.view()
        assert list(df.columns) == ['site', 'position', 'barcode']
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

    def test_repr_multiple_positions(self):
        b = Barcodes(n_sites=1)
        b.add_positions(n_per_site=4)
        b.generate_barcodes(30, 3)
        r = repr(b.sites[1])
        assert '4 position(s)' in r
        assert '12 barcodes' in r


class TestValidation:

    def test_default_avoid_seqs_populated(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 3)
        b.validate()
        assert len(b.avoid_seqs) > 0

    def test_ignore_defaults_clears_avoid(self):
        b = Barcodes(n_sites=1)
        b.generate_barcodes(30, 3)
        b.validate(ignore_defaults=True)
        assert b.avoid_seqs == {}

    def test_contaminated_barcode_replaced(self, capsys):
        # Force contaminated barcode into pool for generation, then supply clean replacement
        contaminated = 'AAAGAATTCAAAAAAAAAAAAAAAAAAAAA'  # EcoRI site (GAATTC)
        clean = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        b = Barcodes(n_sites=1)
        b._bc_pool = [contaminated.ljust(60, 'A')]  # only option → always drawn
        b.generate_barcodes(30, 1)
        b._bc_pool = [clean.ljust(60, 'A')]         # clean replacement available
        b.validate(ignore_defaults=False)
        assert b.sites[1].positions['A']['bc_only'][0] == clean
        assert 'EcoRI' in capsys.readouterr().out

    def test_avoid_enzyme_name(self, capsys):
        contaminated = 'AAAGGATCCAAAAAAAAAAAAAAAAAAAAA'  # BamHI site (GGATCC)
        clean = 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        b = Barcodes(n_sites=1)
        b._bc_pool = [contaminated.ljust(60, 'A')]
        b.generate_barcodes(30, 1)
        b._bc_pool = [clean.ljust(60, 'A')]
        b.validate(ignore_defaults=True, motifs=['BamHI'])
        assert b.sites[1].positions['A']['bc_only'][0] == clean
        assert 'BamHI' in capsys.readouterr().out

    def test_contamination_from_fixed_seq(self, capsys):
        # fixed_left='AAAGAA' + bc_only starting with 'TTC...' → 'GAATTC' (EcoRI) at junction
        bc_only = 'TTCAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        clean   = 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
        b = Barcodes(n_sites=1)
        b._bc_pool = [bc_only.ljust(60, 'A')]
        b.generate_barcodes(30, 1)
        b.add_fixed_sequence('AAAGAA', site=1, side='left')
        b._bc_pool = [clean.ljust(60, 'A')]
        b.validate(ignore_defaults=False)
        assert b.sites[1].positions['A']['bc_only'][0] == clean
        assert 'EcoRI' in capsys.readouterr().out

    def test_pool_exhausted_raises(self):
        contaminated = 'AAAGAATTCAAAAAAAAAAAAAAAAAAAAA'
        b = Barcodes(n_sites=1)
        b._bc_pool = [contaminated.ljust(60, 'A')]
        b.generate_barcodes(30, 1)
        # pool is now empty — replacement impossible
        with pytest.raises(RuntimeError):
            b.validate(ignore_defaults=False)

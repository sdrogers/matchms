import os


from matchms.importing import load_from_mgf
from matchms.filtering import clean_inchis


def test_clean_inchis():
    """Draft for test.
    """
    module_root = os.path.join(os.path.dirname(__file__), '..')

    # Loading
    references_file = os.path.join(module_root, 'tests', 'testdata.mgf')

    reference_spectrums_raw = load_from_mgf(references_file)

    reference_spectrums = [clean_inchis(s) for s in reference_spectrums_raw]

    query_spectrum_raw = reference_spectrums_raw[0]

    query_spectrum = clean_inchis(query_spectrum_raw)

    assert query_spectrum_raw.get("inchi").startswith('InChI='), 'expected different InChI'
    assert query_spectrum.get("inchi").startswith('"InChI='), 'InChI style not as expected.'
    original_inchi = reference_spectrums_raw[2].get("inchi")
    assert original_inchi.startswith('"InChI=CCCCCCCCCCCCCCCC(=O)'), "expected misplaced smiles"
    modified_inchi = reference_spectrums[2].get("inchi")
    assert modified_inchi.startswith('"InChI=1S/C24H50NO7P/'), "inchi was not converted correctly"


if __name__ == '__main__':
    test_clean_inchis()
import os
import pytest
from matchms import calculate_scores, Scores
from matchms.filtering import add_parent_mass
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_mz
from matchms.filtering import select_by_relative_intensity
from matchms.importing import load_from_mgf
from matchms.similarity import CosineGreedy


def apply_my_filters(s):
    s = default_filters(s)
    s = add_parent_mass(s)
    s = normalize_intensities(s)
    s = select_by_relative_intensity(s, intensity_from=0.5, intensity_to=1.0)
    s = select_by_mz(s, mz_from=0, mz_to=1000)
    s = require_minimum_number_of_peaks(s, n_required=5)
    return s

module_root = os.path.join(os.path.dirname(__file__), "..")
spectrums_file = os.path.join("GNPS-LIBRARY.mgf")

# apply my filters to the data
spectrums = [apply_my_filters(s) for s in load_from_mgf(spectrums_file)]

# omit spectrums that didn't qualify for analysis
spectrums = [s for s in spectrums if s is not None]

# define similarity function
cosine_greedy = CosineGreedy()

# this will be a library grouping analysis, so queries = references = spectrums
queries = spectrums[:]
references = spectrums[:]

similarity_matrix = Scores(references=5*references,
                           queries=5*queries,
                           similarity_function=cosine_greedy).calculate().scores

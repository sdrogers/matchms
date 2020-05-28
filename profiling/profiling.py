import os
from timeit import default_timer as timer
from matchms import Scores
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


log = [(timer(), "start")]
module_root = os.path.join(os.path.dirname(__file__), "..")

spectrums_file = os.path.join(module_root, "tests", "pesticides.mgf")
log.append((timer(), "end load file"))

# apply my filters to the data
spectrums = [apply_my_filters(s) for s in load_from_mgf(spectrums_file)]
log.append((timer(), "end filtering"))

# omit spectrums that didn't qualify for analysis
spectrums = [s for s in spectrums if s is not None]
log.append((timer(), "end omit None"))

# define similarity function
cosine_greedy = CosineGreedy()

# this will be a library grouping analysis, so queries = references = spectrums
queries = references = spectrums

log.append((timer(), "start problem definition"))
problem_definition = Scores(references=references, queries=queries, similarity_function=cosine_greedy)
log.append((timer(), "end problem definition"))

log.append((timer(), "start calculating"))
problem_definition.calculate()
log.append((timer(), "end calculating"))

log.append((timer(), "end"))

for index, (t, msg) in enumerate(log):
    delta_t = t - log[index - 1][0] if index > 0 else 0
    print("{0:10.3f}: {1}".format(delta_t, msg))

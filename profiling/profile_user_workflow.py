import os
from timeit import default_timer as timer
from matchms import calculate_scores
from matchms.filtering import add_parent_mass
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_mz
from matchms.filtering import select_by_relative_intensity
from matchms.importing import load_from_mgf
from matchms.similarity import CosineGreedy


def test_user_workflow():

    def log_event(msg):
        log.append((timer(), msg))

    def apply_my_filters(s):
        s = default_filters(s)
        s = add_parent_mass(s)
        s = normalize_intensities(s)
        s = select_by_relative_intensity(s, intensity_from=0.0, intensity_to=1.0)
        s = select_by_mz(s, mz_from=0, mz_to=1000)
        s = require_minimum_number_of_peaks(s, n_required=5)
        return s

    log = []
    log_event("start")

    module_root = os.path.join(os.path.dirname(__file__), "..")
    spectrums_file = os.path.join(module_root, "tests", "pesticides.mgf")

    log_event("")

    # apply my filters to the data
    spectrums = [apply_my_filters(s) for s in load_from_mgf(spectrums_file)]

    # omit spectrums that didn't qualify for analysis
    spectrums = [s for s in spectrums if s is not None]
    log_event("loading and filtering")

    # this will be a library grouping analysis, so queries = references = spectrums
    queries = spectrums[:]
    references = spectrums[:]

    # define similarity function
    cosine_greedy = CosineGreedy()

    log_event("")

    # calculate_scores
    scores = list(calculate_scores(references,
                                   queries,
                                   cosine_greedy))

    log_event("calculating")

    # filter out self-comparisons, require at least 20 matching peaks:
    filtered = [(reference, query, score, n_matching) for (reference, query, score, n_matching) in scores
                if reference != query and n_matching >= 20]

    sorted_by_score = sorted(filtered, key=lambda elem: elem[2], reverse=True)

    log_event("sorting and filtering results")

    log_event("end")

    print("duration    description")
    print("----------- ------------------------------")
    for index, (t, msg) in enumerate(log):
        delta_t = t - log[index - 1][0] if index > 0 else 0
        print("{0:10.3f}: {1}".format(delta_t, msg))


if __name__ == "__main__":
    test_user_workflow()


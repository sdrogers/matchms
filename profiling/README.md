```
$ git clone https://github.com/matchms/matchms.git
$ git checkout bd588a2
$ cd matchms
$ conda env create --file conda/environment-dev.yml
$ conda activate matchms-dev
$ pip install pyprof2calltree
$ pip install --editable .
$ cd profiling
$ cat profile_user_workflow.py
import os
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

    def apply_my_filters(s):
        s = default_filters(s)
        s = add_parent_mass(s)
        s = normalize_intensities(s)
        s = select_by_relative_intensity(s, intensity_from=0.0, intensity_to=1.0)
        s = select_by_mz(s, mz_from=0, mz_to=1000)
        s = require_minimum_number_of_peaks(s, n_required=5)
        return s

    module_root = os.path.join(os.path.dirname(__file__), "..")
    spectrums_file = os.path.join(module_root, "tests", "pesticides.mgf")

    # apply my filters to the data
    spectrums = [apply_my_filters(s) for s in load_from_mgf(spectrums_file)]

    # omit spectrums that didn't qualify for analysis
    spectrums = [s for s in spectrums if s is not None]

    # this will be a library grouping analysis, so queries = references = spectrums
    queries = spectrums[:]
    references = spectrums[:]

    # define similarity function
    cosine_greedy = CosineGreedy()

    # calculate_scores
    scores = list(calculate_scores(references,
                                   queries,
                                   cosine_greedy))

    # filter out self-comparisons, require at least 20 matching peaks:
    filtered = [(reference, query, score, n_matching) for (reference, query, score, n_matching) in scores
                if reference != query and n_matching >= 20]

    sorted_by_score = sorted(filtered, key=lambda elem: elem[2], reverse=True)


if __name__ == "__main__":
    test_user_workflow()

$ python -m cProfile -o pprof.user_workflow.1 profile_user_workflow.py
$ pyprof2calltree -i pprof.user_workflow.1 -o cachegrind.out.user_workflow.1
$ sudo apt install kcachegrind

```

Use ``kcachegrind`` to open cahcegrind.out.user_workflow.1.

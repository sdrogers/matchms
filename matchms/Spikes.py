import numpy


class Spikes:
    """
    Stores arrays of intensities and M/z values, with some checks on their internal consistency.
    """
    def __init__(self, mz=None, intensities=None):

        assert isinstance(mz, numpy.ndarray), "Input argument 'mz' should be a numpy.array."
        assert isinstance(intensities, numpy.ndarray), "Input argument 'intensities' should be a numpy.array."
        assert mz.shape == intensities.shape, "Input arguments 'mz' and 'intensities' should be the same shape."
        assert mz.dtype == "float", "Input argument 'mz' should be an array of type float."
        assert intensities.dtype == "float", "Input argument 'intensities' should be an array of type float."

        self._mz = mz
        self._intensities = intensities
        self._make_immutable()

        assert self._is_sorted(), "mz values are out of order."

    def __eq__(self, other):
        return \
            self._mz.shape == other.mz.shape and \
            numpy.allclose(self._mz, other.mz) and \
            self.intensities.shape == other.intensities.shape and \
            numpy.allclose(self._intensities, other.intensities)

    def __len__(self):
        return self._mz.size

    def __getitem__(self, item):
        return [self._mz, self._intensities][item]

    def _is_sorted(self):
        return numpy.all(self._mz[:-1] <= self._mz[1:])

    def _make_immutable(self):
        self._mz.setflags(write=False)
        self._intensities.setflags(write=False)
        return self

    def clone(self):
        return Spikes(self._mz.copy(), self._intensities.copy())

    @property
    def mz(self):
        """getter method for mz private variable"""
        return self._mz

    @property
    def intensities(self):
        """getter method for intensities private variable"""
        return self._intensities

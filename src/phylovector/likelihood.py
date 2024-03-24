"""
Development method for computing likelihoods.
"""

# mantel from https://github.com/jwcarr/mantel

import itertools
import math

from itertools import permutations
import numpy as np
from scipy import spatial, stats


def mantel(X, Y, perms=10000, method="pearson", tail="two-tail", ignore_nans=False):
    """
    Takes two distance matrices (either redundant matrices or condensed vectors)
    and performs a Mantel test. The Mantel test is a significance test of the
    correlation between two distance matrices.
    Parameters
    ----------
    X : array_like
            First distance matrix (condensed or redundant).
    Y : array_like
            Second distance matrix (condensed or redundant), where the order of
            elements corresponds to the order of elements in the first matrix.
    perms : int, optional
            The number of permutations to perform (default: 10000). A larger
            number gives more reliable results but takes longer to run. If the
            number of possible permutations is smaller, all permutations will
            be tested. This can be forced by setting perms to 0.
    method : str, optional
            Type of correlation coefficient to use; either 'pearson' or 'spearman'
            (default: 'pearson').
    tail : str, optional
            Which tail to test in the calculation of the empirical p-value; either
            'upper', 'lower', or 'two-tail' (default: 'two-tail').
    ignore_nans : bool, optional
            Ignore NaN values in the Y matrix (default: False). This can be
            useful if you have missing values in one of the matrices.
    Returns
    -------
    r : float
            Veridical correlation
    p : float
            Empirical p-value
    z : float
            Standard score (z-score)
    """

    # Ensure that X and Y are represented as Numpy arrays.
    X = np.asarray(X)
    Y = np.asarray(Y)

    # Check that X and Y are valid distance matrices.
    if (
        spatial.distance.is_valid_dm(np.nan_to_num(X)) == False
        and spatial.distance.is_valid_y(X) == False
    ):
        raise ValueError("X is not a valid condensed or redundant distance matrix")
    if (
        spatial.distance.is_valid_dm(np.nan_to_num(Y)) == False
        and spatial.distance.is_valid_y(Y) == False
    ):
        raise ValueError("Y is not a valid condensed or redundant distance matrix")

    # If X or Y is a redundant distance matrix, reduce it to a condensed distance matrix.
    if len(X.shape) == 2:
        X = spatial.distance.squareform(X, force="tovector", checks=False)
    if len(Y.shape) == 2:
        Y = spatial.distance.squareform(Y, force="tovector", checks=False)

    # Check for size equality.
    if len(X) != len(Y):
        raise ValueError("X and Y are not of equal size")

    # Check for minimum size.
    if len(X) < 3:
        raise ValueError("X and Y should represent at least 3 objects")

    # Check finiteness of X and Y
    if not np.isfinite(X).all():
        raise ValueError(
            "X cannot contain NaNs (but Y may contain NaNs, so consider reordering X and Y)"
        )
    finite_Y = np.isfinite(Y)
    if not ignore_nans and not finite_Y.all():
        raise ValueError('Y may contain NaNs, but "ignore_nans" must be set to True')
    if ignore_nans and finite_Y.all():
        ignore_nans = False  # ignore_nans is True but Y contains no nans

    # If Spearman correlation is requested, convert X and Y to ranks.
    method = method.lower()
    if method == "spearman":
        X, Y = stats.rankdata(X), stats.rankdata(Y)
        Y[~finite_Y] = np.nan  # retain any nans, so that these can be ignored later

    # Check for valid method parameter.
    elif method != "pearson":
        raise ValueError('The method should be set to "pearson" or "spearman"')

    # Check for valid tail parameter.
    tail = tail.lower()
    if tail not in ["upper", "lower", "two-tail"]:
        raise ValueError('The tail should be set to "upper", "lower", or "two-tail"')

    # Now we're ready to start the Mantel test using a number of optimizations:
    #
    # 1. Rather than compute correlation coefficients, we'll just compute the
    #    covariances. This works because the denominator in the equation for the
    #    correlation coefficient will yield the same result however the objects
    #    are permuted, making it redundant. Removing the denominator leaves us
    #    with the covariance.
    #
    # 2. Rather than permute the Y distances and derive the residuals to calculate
    #    the covariance with the X distances, we'll represent the Y residuals in
    #    the matrix and shuffle those directly.
    #
    # 3. If the number of possible permutations is less than the number of
    #    permutations that were requested, we'll run a deterministic test where
    #    we try all possible permutations rather than sample the permutation
    #    space. This gives a faster, deterministic result.

    # Calculate the X and Y residuals, which will be used to compute the
    # covariance under each permutation.
    X_residuals = X - np.mean(X[finite_Y])
    Y_residuals = Y - np.mean(Y[finite_Y])

    # Expand the Y residuals to a redundant matrix.
    Y_residuals_as_matrix = spatial.distance.squareform(
        Y_residuals, force="tomatrix", checks=False
    )

    m = len(Y_residuals_as_matrix)  # number of objects
    n = np.math.factorial(m)  # number of possible matrix permutations

    # If the number of requested permutations is greater than the number of
    # possible permutations (m!) or the perms parameter is set to 0, then run a
    # deterministic Mantel test
    if perms >= n or perms == 0:
        if ignore_nans:
            correlations = deterministic_test_with_nans(m, n, X, Y_residuals_as_matrix)
        else:
            correlations = deterministic_test(m, n, X_residuals, Y_residuals_as_matrix)
        # correlations[0] is the veridical correlation

    else:
        if ignore_nans:
            correlations = stochastic_test_with_nans(m, perms, X, Y_residuals_as_matrix)
        else:
            correlations = stochastic_test(m, perms, X_residuals, Y_residuals_as_matrix)
        correlations[0] = sum(X_residuals[finite_Y] * Y_residuals[finite_Y]) / np.sqrt(
            sum(X_residuals[finite_Y] ** 2) * sum(Y_residuals[finite_Y] ** 2)
        )  # compute veridical correlation and place in positon 0

    r = correlations[0]

    if tail == "upper":
        p = sum(correlations >= r) / len(correlations)
    elif tail == "lower":
        p = sum(correlations <= r) / len(correlations)
    elif tail == "two-tail":
        p = sum(abs(correlations) >= abs(r)) / len(correlations)

    z = (r - np.mean(correlations)) / np.std(correlations)

    return r, p, z


def deterministic_test(m, n, X_residuals, Y_residuals_as_matrix):
    Y_residuals_permuted = np.zeros((m**2 - m) // 2)
    covariances = np.zeros(n)
    for i, order in enumerate(permutations(range(m))):
        Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]
        spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
            Y_residuals_as_matrix_permuted, Y_residuals_permuted
        )
        covariances[i] = (X_residuals * Y_residuals_permuted).sum()
    denominator = np.sqrt(sum(X_residuals**2) * sum(Y_residuals_permuted**2))
    return covariances / denominator


def deterministic_test_with_nans(m, n, X, Y_residuals_as_matrix):
    Y_residuals_permuted = np.zeros((m**2 - m) // 2)
    correlations = np.zeros(n)
    for i, order in enumerate(permutations(range(m))):
        Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]
        spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
            Y_residuals_as_matrix_permuted, Y_residuals_permuted
        )
        # Since on each permutation we will be ignoring different values in X,
        # the X_residuals need to be recomputed each time depending on which
        # values in permuted Y are finite.
        finite_Y_permuted = np.isfinite(Y_residuals_permuted)
        reduced_X = X[finite_Y_permuted]
        reduced_X_residuals = reduced_X - reduced_X.mean()
        reduced_Y_residuals = Y_residuals_permuted[finite_Y_permuted]
        covariance = (reduced_X_residuals * reduced_Y_residuals).sum()
        # The denominator will be different on each permutation
        denominator = np.sqrt(sum(reduced_X_residuals**2) * sum(reduced_Y_residuals**2))
        correlations[i] = covariance / denominator
    return correlations


def stochastic_test(m, n, X_residuals, Y_residuals_as_matrix):
    Y_residuals_permuted = np.zeros((m**2 - m) // 2)
    covariances = np.zeros(n)
    order = np.arange(m)
    for i in range(1, n):
        np.random.shuffle(order)
        Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]
        spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
            Y_residuals_as_matrix_permuted, Y_residuals_permuted
        )
        covariances[i] = (X_residuals * Y_residuals_permuted).sum()
    denominator = np.sqrt(sum(X_residuals**2) * sum(Y_residuals_permuted**2))
    return covariances / denominator


def stochastic_test_with_nans(m, n, X, Y_residuals_as_matrix):
    Y_residuals_permuted = np.zeros((m**2 - m) // 2)
    correlations = np.zeros(n)
    order = np.arange(m)
    for i in range(1, n):
        np.random.shuffle(order)
        Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]
        spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
            Y_residuals_as_matrix_permuted, Y_residuals_permuted
        )
        # Since on each permutation we will be ignoring different values in X,
        # the X_residuals need to be recomputed each time depending on which
        # values in permuted Y are finite.
        finite_Y_permuted = np.isfinite(Y_residuals_permuted)
        reduced_X = X[finite_Y_permuted]
        reduced_X_residuals = reduced_X - reduced_X.mean()
        reduced_Y_residuals = Y_residuals_permuted[finite_Y_permuted]
        covariance = (reduced_X_residuals * reduced_Y_residuals).sum()
        # The denominator will be different on each permutation
        denominator = np.sqrt(sum(reduced_X_residuals**2) * sum(reduced_Y_residuals**2))
        correlations[i] = covariance / denominator
    return correlations


def method1(vector, dm):
    # Collect information on number of taxa from vector length, collect the taxa,
    # and make sure there is a match
    n = int((math.sqrt(8 * len(vector) + 1) - 1) / 2)
    taxa = sorted(set(itertools.chain.from_iterable(dm.keys())))
    if n != len(taxa):
        raise ValueError(
            "Mismatch in number of taxa from vector length and distance matrix."
        )

    # Compute distances between pairs of taxa: note that this is different from the
    # information in the vector, which is the height of the pair's MRCA. We
    # also collect the height of the most recent common ancestor of each leaf,
    # so as to check in the step below
    dist = {}
    anc_height = {taxon: max(vector[n:]) for taxon in taxa}
    for idx, (leaf_i, leaf_j) in enumerate(itertools.combinations(range(n), 2)):
        com_anc_height = vector[n + idx]
        dist_i = vector[leaf_i]
        dist_j = vector[leaf_j]
        dist[taxa[leaf_i], taxa[leaf_j]] = (2 * com_anc_height) - dist_i - dist_j

        # Update height of most recent ancestor, if necessary
        if anc_height[taxa[leaf_i]] > com_anc_height:
            anc_height[taxa[leaf_i]] = com_anc_height
        if anc_height[taxa[leaf_j]] > com_anc_height:
            anc_height[taxa[leaf_j]] = com_anc_height

    # Reject trees where the (randomly selected) height of a leaf is before
    # its ancestor
    # TODO: make proportional
    penalty = 1.0
    for i, taxon in enumerate(taxa):
        if anc_height[taxon] <= vector[i]:
            penalty += 1.0

    # Compute differences between tree distances and matrix
    diff = sum(
        [
            abs(dm[taxa_i, taxa_j] - dist[taxa_i, taxa_j])
            for taxa_i, taxa_j in itertools.combinations(taxa, 2)
        ]
    )

    return diff * penalty


def primate(vector):
    dm = {
        ("bonobo", "chimpanzee"): 0.09696787751425995,
        ("bonobo", "gorilla"): 0.24977484238967274,
        ("bonobo", "homo"): 0.22095466826778742,
        ("bonobo", "orangutan"): 0.3161212848994296,
        ("chimpanzee", "gorilla"): 0.25788051636145304,
        ("chimpanzee", "homo"): 0.22695887120984692,
        ("chimpanzee", "orangutan"): 0.32182527769438607,
        ("gorilla", "homo"): 0.2668868207745422,
        ("gorilla", "orangutan"): 0.322425697988592,
        ("homo", "orangutan"): 0.31732212548784144,
    }

    # Collect information on number of taxa from vector length, collect the taxa,
    # and make sure there is a match
    n = int((math.sqrt(8 * len(vector) + 1) - 1) / 2)
    taxa = sorted(set(itertools.chain.from_iterable(dm.keys())))
    if n != len(taxa):
        raise ValueError(
            "Mismatch in number of taxa from vector length and distance matrix."
        )

    # Compute distances between pairs of taxa: note that this is different from the
    # information in the vector, which is the height of the pair's MRCA. We
    # also collect the height of the most recent common ancestor of each leaf,
    # so as to check in the step below
    dist = {}
    anc_height = {taxon: max(vector[n:]) for taxon in taxa}
    for idx, (leaf_i, leaf_j) in enumerate(itertools.combinations(range(n), 2)):
        com_anc_height = vector[n + idx]
        dist_i = vector[leaf_i]
        dist_j = vector[leaf_j]
        dist[taxa[leaf_i], taxa[leaf_j]] = (2 * com_anc_height) - dist_i - dist_j

        # Update height of most recent ancestor, if necessary
        if anc_height[taxa[leaf_i]] > com_anc_height:
            anc_height[taxa[leaf_i]] = com_anc_height
        if anc_height[taxa[leaf_j]] > com_anc_height:
            anc_height[taxa[leaf_j]] = com_anc_height

    # Reject trees where the (randomly selected) height of a leaf is before
    # its ancestor
    # TODO: make proportional? vector[i] - anc_height[taxon]  over vector?
    penalty = 1.0
    for i, taxon in enumerate(taxa):
        if anc_height[taxon] <= vector[i]:
            penalty += 1.0

    X, Y = [], []
    for taxa_i, taxa_j in itertools.combinations(taxa, 2):
        X.append(dm[taxa_i, taxa_j])
        Y.append(dist[taxa_i, taxa_j])

    # Compute differences between tree distances and matrix
    #    diff = sum(
    #        [
    #            abs(dm[taxa_i, taxa_j] - dist[taxa_i, taxa_j])
    #            for taxa_i, taxa_j in itertools.combinations(taxa, 2)
    #        ]
    #    )

    #    print(
    #        [
    #            (dm[taxa_i, taxa_j], dist[taxa_i, taxa_j])
    #            for taxa_i, taxa_j in itertools.combinations(taxa, 2)
    #        ]
    ##    )
    #   print(diff, penalty, diff * penalty)

    m = mantel(X, Y, method="pearson", tail="upper")
    print(m, penalty)

    return m[0] + penalty
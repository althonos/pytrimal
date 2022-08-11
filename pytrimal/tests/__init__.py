# noqa: D104

from . import (
    test_alignment,
    test_automatic_trimmer,
    test_doctest,
    test_manual_trimmer,
    test_overlap_trimmer,
    test_similarity_matrix,
)


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_alignment))
    suite.addTests(loader.loadTestsFromModule(test_automatic_trimmer))
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    suite.addTests(loader.loadTestsFromModule(test_manual_trimmer))
    suite.addTests(loader.loadTestsFromModule(test_overlap_trimmer))
    suite.addTests(loader.loadTestsFromModule(test_similarity_matrix))
    return suite

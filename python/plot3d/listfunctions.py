"""List utility functions for deduplication and pair enumeration."""


def unique_pairs(listOfItems=[]):
    """Yield unique unordered pairs from a list of 2-tuples.

    Filters out self-pairs ``(x, x)`` and duplicate reversed pairs
    so that ``(1, 2)`` and ``(2, 1)`` are treated as the same pair.

    Args:
        listOfItems: Iterable of ``(x, y)`` tuples.

    Yields:
        ``(x, y)`` pairs where ``x != y`` and the reversed pair has not
        already been yielded. Uses a set for O(1) duplicate lookups.
    """
    seen = set()
    for x, y in listOfItems:
        if x != y and (y, x) not in seen:
            seen.add((x, y))
            yield x, y

from __future__ import generators, nested_scopes

def itercat(*iterators):
    """Concatenate several iterators into one."""
    for i in iterators:
        for x in i:
            yield x

def iterwhile(func, iterator):
    """Iterate for as long as func(value) returns true."""
    iterator = iter(iterator)
    while 1:
        next = iterator.next()
        if not func(next):
            raise StopIteration
        yield next

def iterfirst(iterator, count=1):
    """Iterate through 'count' first values."""
    iterator = iter(iterator)
    for i in xrange(count):
        yield iterator.next()

def iterstep(iterator, n):
    """Iterate every nth value."""
    iterator = iter(iterator)
    while 1:
        yield iterator.next()
        # skip n-1 values
        for dummy in range(n-1):
            iterator.next()

def itergroup(iterator, count):
    """Iterate in groups of 'count' values. If there
    aren't enough values, the last result is padded with
    None."""
    iterator = iter(iterator)
    values_left = [1]
    def values():
        values_left[0] = 0
        for x in range(count):
            try:
                yield iterator.next()
                values_left[0] = 1
            except StopIteration:
                yield None
    while 1:
        value = tuple(values())
        if not values_left[0]:
            raise StopIteration
        yield value
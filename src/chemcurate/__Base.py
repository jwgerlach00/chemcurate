from typing import Iterable, Generator, List
from itertools import zip_longest


class __Base:
    _batch_size = 100 # NOTE: Modification not recommended due to URL length limit.
    
    def batch(iterable: Iterable[str], n: int = 1) -> Generator[List[str], None, None]:
        args = [iter(iterable)] * n
        yield from ([str(x) for x in tup if x is not None] for tup in zip_longest(*args))

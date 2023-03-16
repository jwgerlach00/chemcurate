from typing import Iterable, Generator


class __Base:
    _batch_size = 100 # NOTE: Modification not recommended due to URL length limit.
    
    @staticmethod
    def batch(iterable:Iterable, n:int=1) -> Generator:
        l = len(iterable)
        for ndx in range(0, l, n):
            yield iterable[ndx:min(ndx + n, l)]

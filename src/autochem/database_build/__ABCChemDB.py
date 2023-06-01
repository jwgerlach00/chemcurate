from abc import ABC, abstractmethod


class __ABCChemDB(ABC):
    @abstractmethod
    def build(self):
        pass


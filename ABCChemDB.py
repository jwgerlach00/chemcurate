from abc import ABC, abstractmethod


class ABCChemDB(ABC):
    @abstractmethod
    def build(self):
        pass


from abc import ABCMeta, abstractproperty

import rxncon.core.rxncon_system as sys


class RxnConInput(metaclass=ABCMeta):

    @abstractproperty
    def rxncon_system(self) -> sys.RxnConSystem:
        pass

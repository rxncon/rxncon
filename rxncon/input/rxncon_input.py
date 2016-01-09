from abc import ABCMeta, abstractproperty

import rxncon.core.rxncon_system as sys


class RxnConInput:
    __metaclass__ = ABCMeta

    @abstractproperty
    def rxncon_system(self) -> sys.RxnConSystem:
        pass

from typing import List, TypeVar
import inspect
from colorama import Fore


T = TypeVar('T')


def elems_eq(first_list: List[T], second_list: List[T]) -> bool:
    if all(isinstance(x, list) for x in first_list) and all(isinstance(x, list) for x in second_list):  # type: ignore
        uniq_first = [set(x) for x in first_list]    # type: ignore
        uniq_second = [set(x) for x in second_list]  # type: ignore

        return all(x in uniq_second for x in uniq_first) and all(x in uniq_first for x in uniq_second)

    else:
        return set(first_list) == set(second_list)


def current_function_name(colored: bool=True) -> str:
    name = inspect.stack()[1][3]
    if colored:
        return Fore.MAGENTA + name + Fore.RESET
    else:
        return name

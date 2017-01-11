import functools
from typing import List, TypeVar, Callable, Any
import inspect
from colorama import Fore


T = TypeVar('T')

def compose(*functions: Callable) -> Callable:
    return functools.reduce(lambda f, g: lambda x: f(g(x)), functions, lambda x: x)


def members(obj: Any) -> List[str]:
    return [x for x in dir(obj) if not x.startswith('__')]


def elems_eq(first_list: List[T], second_list: List[T]) -> bool:
    if all(isinstance(x, list) for x in first_list) and all(isinstance(x, list) for x in second_list):  # type: ignore
        uniq_first = [set(x) for x in first_list]    # type: ignore
        uniq_second = [set(x) for x in second_list]  # type: ignore

        return all(x in uniq_second for x in uniq_first) and all(x in uniq_first for x in uniq_second)

    else:
        return set(first_list) == set(second_list)

def all_type_eq(elems: List[T]) -> bool:
    assert elems
    the_type = type(elems[0])

    for elem in elems:
        if not isinstance(elem, the_type):
            return False

    return True

def all_eq(elems: List[T]) -> bool:
    assert elems
    the_val = elems[0]

    for elem in elems:
        if not elem == the_val:
            return False

    return True

def current_function_name(colored: bool=True) -> str:
    name = inspect.stack()[1][3]
    if colored:
        return Fore.MAGENTA + name + Fore.RESET
    else:
        return name

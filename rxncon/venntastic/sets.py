import functools
import itertools as itt
import typing as tg

METHOD_COMPLEMENTS_EXPANDED = '_complements_expanded'
METHOD_UNIONS_MOVED_TO_LEFT = '_unions_moved_to_left'
METHOD_SIMPLIFY_EMPTY_UNIVERSAL = '_simplify_empty_universal'
METHOD_SIMPLIFY_BINARY = '_simplify_binary'

METHOD_NAME_COMPLEMENTS = 'complements'

# @todo Make typechecking work in this module. Problem is the co/contravariance of function arguments and return types.
class Set:
    def simplified_form(self) -> 'Set':
        simplification_methods = [
            METHOD_COMPLEMENTS_EXPANDED,
            METHOD_UNIONS_MOVED_TO_LEFT,
            METHOD_SIMPLIFY_EMPTY_UNIVERSAL,
            METHOD_SIMPLIFY_BINARY
        ]

        return _call_method_list_until_stable(self, simplification_methods)

    def to_nested_list_form(self) -> tg.List[tg.List['Set']]:
        return _cleaned_nested_list_form(self.simplified_form()._to_nested_list())

    def to_boolean_function(self) -> 'BooleanFunction':
        return boolean_function_from_nested_list_form(self.to_nested_list_form())

    def to_union_list_form(self) -> tg.List['Set']:
        union_terms = []
        for term in self.to_nested_list_form():
            union_terms.append(nested_expression_from_list_and_binary_op(term, Intersection))

        return union_terms

    def to_full_simplified_form(self) -> 'Set':
        return nested_expression_from_list_and_binary_op(self.to_union_list_form(), Union)

    @property
    def cardinality(self):
        list_of_intersections = self.to_nested_list_form()

        cardinality = {}

        for i in range(len(list_of_intersections)):
            intersections = list(itt.combinations(list_of_intersections, i + 1))

            for intersection_parts in intersections:
                intersection = []

                for part in intersection_parts:
                    intersection += part

                partial_cardinality = _cardinality_of_intersection_term(_cleaned_intersection_term(intersection))

                if i % 2:
                    cardinality = _add_dicts(cardinality, _negate_dict(partial_cardinality))

                else:
                    cardinality = _add_dicts(cardinality, partial_cardinality)

        if () in cardinality.keys() and UniversalSet() in cardinality.keys():
            cardinality[(UniversalSet(),)] += cardinality.pop(())

        elif () in cardinality.keys():
            cardinality[(UniversalSet(),)] = cardinality.pop(())

        return {k: v for k, v in cardinality.items() if v != 0}

    def is_equivalent_to(self, other: 'Set') -> bool:
        return self.is_superset_of(other) and self.is_subset_of(other)

    def is_superset_of(self, other: 'Set') -> bool:
        return other.to_boolean_function().implies(self.to_boolean_function())

    def is_subset_of(self, other: 'Set') -> bool:
        return self.to_boolean_function().implies(other.to_boolean_function())

    def _complements_expanded(self) -> 'Set':
        return self

    def _unions_moved_to_left(self) -> 'Set':
        return self

    def _simplify_empty_universal(self):
        return self

    def _simplify_binary(self):
        return self

    def _to_nested_list(self) -> tg.List[tg.List['Set']]:
        pass


class UnarySet(Set):
    def _to_nested_list(self):
        return [[self]]


class PropertySet(UnarySet):
    def __init__(self, value):
        assert isinstance(hash(value), int)
        self.value = value

    def __eq__(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, PropertySet):
            if self.value is None and other.value is None:
                return True

            elif self.value is None or other.value is None:
                return False

            return self.value == other.value

        else:
            return False

    def __hash__(self) -> int:
        return hash('*property-set-{}*'.format(hash(self.value)))

    def __lt__(self, other: Set) -> bool:
        if isinstance(other, PropertySet):
            if self.value is None:
                return True

            return self.value < other.value

        elif isinstance(other, EmptySet):
            return False

        elif isinstance(other, Complement) or isinstance(other, BinarySet):
            return True

        else:
            raise AssertionError

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if self.value:
            return 'Property({})'.format(self.value)

        else:
            return 'UniversalSet'

    def _complements_expanded(self):
        return self


class EmptySet(UnarySet):
    def __eq__(self, other: Set) -> bool:
        return isinstance(other, EmptySet)

    def __hash__(self) -> int:
        return hash('*empty-set*')

    def __lt__(self, other: Set) -> bool:
        assert isinstance(other, Set)
        return False if isinstance(other, EmptySet) else True

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'EmptySet'

    def is_subset_of(self, other: 'Set') -> bool:
        return True


class Complement(UnarySet):
    def __init__(self, expr: Set):
        self.expr = expr

    def __eq__(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, Complement):
            return self.expr == other.expr

        else:
            return False

    def __hash__(self) -> int:
        return hash('*complement-{}*'.format(self.expr))

    def __lt__(self, other: Set) -> bool:
        if isinstance(other, Complement):
            return self.expr < other.expr

        elif isinstance(other, EmptySet) or isinstance(other, PropertySet):
            return False

        else:
            return True

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Complement({})'.format(self.expr)

    def _complements_expanded(self) -> Set:
        if isinstance(self.expr, Complement):
            return self.expr.expr._complements_expanded()

        elif isinstance(self.expr, EmptySet):
            return UniversalSet()

        elif self.expr == UniversalSet():
            return EmptySet()

        elif isinstance(self.expr, Union):
            return Intersection(Complement(self.expr.left_expr), Complement(self.expr.right_expr))._complements_expanded()

        elif isinstance(self.expr, Intersection):
            return Union(Complement(self.expr.left_expr), Complement(self.expr.right_expr))._complements_expanded()

        elif isinstance(self.expr, PropertySet) and hasattr(self.expr.value, METHOD_NAME_COMPLEMENTS):
            complementary_property_sets = [PropertySet(x) for x in getattr(self.expr.value, METHOD_NAME_COMPLEMENTS)()]
            return nested_expression_from_list_and_binary_op(complementary_property_sets, Union)

        else:
            return Complement(self.expr._complements_expanded())


class BinarySet(Set):
    def __init__(self, left_expr: Set, right_expr: Set):
        self.left_expr = left_expr
        self.right_expr = right_expr

    def __eq__(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, type(self)) and (self.left_expr == other.left_expr) and (self.right_expr == other.right_expr):
            return True
        else:
            return False

    def _complements_expanded(self) -> Set:
        return type(self)(self.left_expr._complements_expanded(), self.right_expr._complements_expanded())

    def _simplify_binary(self):
        if self.left_expr == self.right_expr:
            return self.left_expr
        else:
            return type(self)(self.left_expr._simplify_binary(), self.right_expr._simplify_binary())

class Intersection(BinarySet):
    def __lt__(self, other: Set) -> bool:
        if isinstance(other, EmptySet) or isinstance(other, PropertySet) or isinstance(other, Complement):
            return False

        elif isinstance(other, Intersection):
            if self.left_expr < other.left_expr:
                return True

            elif self.left_expr == other.left_expr:
                return self.right_expr < other.right_expr

            return False

        elif isinstance(other, Union):
            return True

        else:
            raise AssertionError

    def __hash__(self) -> int:
        return hash('*intersection-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Intersection({0}, {1})'.format(self.left_expr, self.right_expr)

    def _unions_moved_to_left(self) -> Set:
        if isinstance(self.left_expr, Union):
            # Distributivity of intersection with respect to union
            return Union(Intersection(self.left_expr.left_expr, self.right_expr), Intersection(self.left_expr.right_expr, self.right_expr))\
                ._unions_moved_to_left()

        elif isinstance(self.right_expr, Union):
            # Distributivity of intersection with respect to union
            return Union(Intersection(self.left_expr, self.right_expr.left_expr), Intersection(self.left_expr, self.right_expr.right_expr))\
                ._unions_moved_to_left()

        else:
            return Intersection(self.left_expr._unions_moved_to_left(), self.right_expr._unions_moved_to_left())

    def _simplify_empty_universal(self):
        if self.left_expr == UniversalSet():
            return self.right_expr._simplify_empty_universal()
        elif self.right_expr == UniversalSet():
            return self.left_expr._simplify_empty_universal()
        elif self.left_expr == EmptySet():
            return EmptySet()
        elif self.right_expr == EmptySet():
            return EmptySet()
        else:
            return Intersection(self.left_expr._simplify_empty_universal(), self.right_expr._simplify_empty_universal())

    def _to_nested_list(self) -> tg.List[tg.List[Set]]:
        return [self.left_expr._to_nested_list()[0] + self.right_expr._to_nested_list()[0]]


class Union(BinarySet):
    def __hash__(self) -> int:
        return hash('*union-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __lt__(self, other: Set) -> bool:
        if isinstance(other, EmptySet) or isinstance(other, PropertySet) or \
                isinstance(other, Complement) or isinstance(other, Intersection):
            return False

        elif isinstance(other, Union):
            if self.left_expr < other.left_expr:
                return True

            elif self.left_expr == other.left_expr:
                return self.right_expr < other.right_expr

            return False

        else:
            raise AssertionError

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Union({0}, {1})'.format(self.left_expr, self.right_expr)

    def _unions_moved_to_left(self) -> Set:
        if isinstance(self.right_expr, Union) and not isinstance(self.left_expr, Union):
            # Use commutativity to move all Union expressions to the left
            return Union(self.right_expr, self.left_expr)._unions_moved_to_left()

        elif isinstance(self.left_expr, Union) and isinstance(self.right_expr, Union):
            # Use associativity to move all Union expressions to the left
            return Union(Union(Union(self.left_expr.left_expr, self.left_expr.right_expr), self.right_expr.left_expr),
                         self.right_expr.right_expr)._unions_moved_to_left()

        else:
            return Union(self.left_expr._unions_moved_to_left(), self.right_expr._unions_moved_to_left())

    def _simplify_empty_universal(self):
        if self.left_expr == UniversalSet():
            return UniversalSet()
        elif self.right_expr == UniversalSet():
            return UniversalSet()
        elif self.left_expr == EmptySet():
            return self.right_expr._simplify_empty_universal()
        elif self.right_expr == EmptySet():
            return self.left_expr._simplify_empty_universal()
        else:
            return Union(self.left_expr._simplify_empty_universal(), self.right_expr._simplify_empty_universal())

    def _to_nested_list(self) -> tg.List[tg.List[Set]]:
        return self.left_expr._to_nested_list() + self.right_expr._to_nested_list()


class Difference(Set):
    def __new__(cls, *args, **kwargs) -> Set:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))


def UniversalSet() -> Set:
    return PropertySet(None)


### PUBLIC FUNCTIONS ###
def nested_expression_from_list_and_binary_op(xs: tg.List[Set], binary_op) -> Set:
    if binary_op == Intersection:
        unit = UniversalSet()
    elif binary_op == Union:
        unit = EmptySet()
    else:
        raise TypeError

    if len(xs) == 0:
        return unit
    elif len(xs) == 1:
        return xs[0]
    else:
        return functools.reduce(binary_op, xs[1:], xs[0])


def set_from_nested_list_form(xss: tg.List[tg.List[Set]]) -> Set:
    union_terms = [nested_expression_from_list_and_binary_op(xs, Intersection) for xs in xss]
    return nested_expression_from_list_and_binary_op(union_terms, Union)


def gram_schmidt_disjunctify(overlapping_sets: tg.List[Set]) -> tg.List[Set]:
    simplified_overlapping_sets = []

    for x in overlapping_sets:
        simplified_set = x.simplified_form()
        if simplified_set not in simplified_overlapping_sets:
            simplified_overlapping_sets.append(simplified_set)

    non_overlapping_sets = []

    complements = []

    for x in simplified_overlapping_sets:
        complements.append(Complement(x).simplified_form())

    for i, x in enumerate(overlapping_sets):
        assert len(x.to_union_list_form()) == 1

        for y in complements[:i]:
            x = Intersection(x, y)

        non_overlapping_sets.append(x)

    result = []

    for x in non_overlapping_sets:
        result += x.to_union_list_form()

    return result


### BOOLEAN FUNCTIONS FROM SETS ###
def boolean_function_from_nested_list_form(nested_list: tg.List[tg.List[Set]]) -> 'BooleanFunction':
    assert isinstance(nested_list, list)
    assert len(nested_list) > 0

    if all(isinstance(item, EmptySet) for term in nested_list for item in term):
        return BooleanFunctionAlwaysFalse()

    and_clauses = []

    for term in nested_list:
        assert isinstance(term, list)

        if all(item == UniversalSet() for item in term):
            return BooleanFunctionAlwaysTrue()

        if any(isinstance(item, EmptySet) for item in term):
            continue

        required_true = []
        required_false = []

        for item in term:
            if item == UniversalSet():
                continue

            elif isinstance(item, PropertySet):
                required_true.append(item)

            elif isinstance(item, Complement):
                required_false.append(item.expr)

            else:
                raise AssertionError

        if any(x in required_false for x in required_true) or any(x in required_true for x in required_false):
            continue

        and_clauses.append(BooleanAndClause(required_true, required_false))

    if not and_clauses:
        return BooleanFunctionAlwaysFalse()

    else:
        return BooleanFunction(and_clauses)


class BooleanFunction:
    def __init__(self, and_clauses):
        assert all([isinstance(x, BooleanAndClause) for x in and_clauses])
        self.and_clauses = and_clauses
        self.valid_statements = list(set([statement for and_clause in self.and_clauses for statement in and_clause.valid_statements]))

    def __call__(self, *args, **kwargs) -> bool:
        actual_true = kwargs.get('true', [])
        actual_false = kwargs.get('false', [])

        if any([statements not in self.valid_statements for statements in actual_true + actual_false]):
            raise NameError

        return any([and_clause(**{'true': actual_true, 'false': actual_false} ) for and_clause in self.and_clauses])

    def implies(self, other) -> bool:
        assert isinstance(other, BooleanFunction)

        for true_statements, false_statements in generate_boolean_value_lists(self.valid_statements):
            my_evaluation = self(true=true_statements, false=false_statements)
            other_evaluation = other(true=[statement for statement in true_statements if other.is_valid_statement(statement)],
                                     false=[statement for statement in false_statements if other.is_valid_statement(statement)])

            if my_evaluation and not other_evaluation:
                return False

        return True

    def is_valid_statement(self, statement) -> bool:
        return statement in self.valid_statements


class BooleanFunctionAlwaysFalse(BooleanFunction):
    def __init__(self):
        self.valid_statements = []
        self.and_clauses = []

    def __call__(self, *args, **kwargs) -> bool:
        return False

    def implies(self, other):
        # Since this 'function' has no 'valid_statements', we use the statement list of its counterpart.
        self.valid_statements = other.valid_statements
        implies = super().implies(other)
        self.valid_statements = []

        return implies

    def is_valid_statement(self, statement) -> bool:
        return True


class BooleanFunctionAlwaysTrue(BooleanFunction):
    def __init__(self):
        self.valid_statements = []
        self.and_clauses = []

    def __call__(self, *args, **kwargs) -> bool:
        return True

    def implies(self, other):
        # Since this 'function' has no 'valid_statements', we use the statement list of its counterpart.
        self.valid_statements = other.valid_statements
        implies = super().implies(other)
        self.valid_statements = []

        return implies

    def is_valid_statement(self, statement) -> bool:
        return True


class BooleanAndClause:
    def __init__(self, required_true=None, required_false=None):
        if required_true is None:
            required_true = []

        if required_false is None:
            required_false = []

        assert not any([x in required_false for x in required_true])
        assert not any([x in required_true for x in required_false])

        self.required_true = required_true
        self.required_false = required_false
        self.valid_statements = list(set(required_true + required_false))

    def __call__(self, *args, **kwargs) -> bool:
        if any([x in self.required_false for x in kwargs['true']]) or any([x in self.required_true for x in kwargs['false']]):
            return False

        else:
            return all([x in kwargs['true'] for x in self.required_true] + [x in kwargs['false'] for x in self.required_false])


def generate_boolean_value_lists(statements):
    for selector in itt.product([True, False], repeat=len(statements)):
        yield list(itt.compress(statements, selector)), list(itt.compress(statements, [not x for x in selector]))


### EXPRESSION SIMPLIFIERS ###
def _cleaned_nested_list_form(nested_list: tg.List[tg.List[Set]]) -> tg.List[tg.List[Set]]:
    clean_terms = []

    for term in nested_list:
        cleaned_intersection_term = _cleaned_intersection_term(term)

        if cleaned_intersection_term != [EmptySet()]:
            clean_terms.append(cleaned_intersection_term)

    if len(clean_terms) == 0:
        return [[EmptySet()]]

    clean_terms = _remove_subsets_from_union_terms(clean_terms)
    clean_terms = _remove_partial_complements_from_union_terms(clean_terms)
    clean_terms = _remove_subsets_from_union_terms(clean_terms)

    return clean_terms


def _cleaned_intersection_term(term: tg.List[Set]) -> tg.List[Set]:
    cleaned_term = []

    if all(item == UniversalSet() for item in term):
        return [UniversalSet()]

    for item in term:
        if item == EmptySet():
            cleaned_term = [EmptySet()]
            break

        elif item == UniversalSet():
            continue

        elif Complement(item) in cleaned_term or item in [Complement(x) for x in cleaned_term]:
            cleaned_term = [EmptySet()]
            break

        elif item not in cleaned_term:
            cleaned_term.append(item)

        else:
            continue

    return cleaned_term


def _remove_subsets_from_union_terms(union_terms: tg.List[tg.List[Set]]) -> tg.List[tg.List[Set]]:
    union_terms.sort(key=len)
    new_union_terms = [union_terms.pop(0)]

    for term in union_terms:
        is_subset = False
        for possible_superset in [x for x in new_union_terms]:
            if all(x in term for x in possible_superset):
                is_subset = True
                break

        if not is_subset:
            new_union_terms.append(term)

    return new_union_terms


def _remove_partial_complements_from_union_terms(union_terms: tg.List[tg.List[Set]]) -> tg.List[tg.List[Set]]:
    new_union_terms = []

    for term in union_terms:
        simplified_term = None
        for possible_counter_term in [x for x in union_terms if len(x) == len(term)]:
            possible_simplified_term = []
            for s in list(set(possible_counter_term).union(set(term))):
                if isinstance(s, Complement) and s.expr in possible_simplified_term:
                    possible_simplified_term.remove(s.expr)
                elif isinstance(s, PropertySet) and Complement(s) in possible_simplified_term:
                    possible_simplified_term.remove(Complement(s))
                else:
                    possible_simplified_term.append(s)

            if len(possible_simplified_term) < len(term):
                simplified_term = possible_simplified_term

        new_union_terms.append(simplified_term) if simplified_term else new_union_terms.append(term)

    return new_union_terms


### PROTECTED HELPERS ###
def _cardinality_of_intersection_term(term: tg.List[Set]) -> tg.Dict:
    if any(isinstance(x, EmptySet) for x in term):
        return {}

    if not any(isinstance(x, Complement) for x in term):
        return {tuple(sorted(term)): 1}

    term_copy = term.copy()

    head_of_list = []
    while term_copy:
        x = term_copy.pop()
        if isinstance(x, Complement):
            rest_of_list = head_of_list + term_copy
            return _add_dicts(_cardinality_of_intersection_term(rest_of_list),
                              _negate_dict(_cardinality_of_intersection_term(rest_of_list + [x.expr])))

        elif isinstance(x, PropertySet):
            head_of_list.append(x)

        else:
            raise AssertionError


def _add_dicts(x: tg.Dict, y: tg.Dict) -> tg.Dict:
    res = {}
    for k, v in x.items():
        res[k] = v

    for k, v in y.items():
        if k in res:
            res[k] += v

        else:
            res[k] = v

    return res


def _negate_dict(x: tg.Dict) -> tg.Dict:
    res = {}
    for k, v in x.items():
        res[k] = -1 * v

    return res


def _call_method_list_until_stable(expr: Set, methods: tg.List[str]):
    max_simplifications = 100
    i = 0

    previous_simplification = expr
    current_simplification = expr
    simplification_done = False

    while not simplification_done:
        i += 1
        if i > max_simplifications:
            raise RecursionError

        for method in methods:
            current_simplification = getattr(current_simplification, method)()

        if current_simplification == previous_simplification:
            simplification_done = True

        else:
            previous_simplification = current_simplification

    return current_simplification

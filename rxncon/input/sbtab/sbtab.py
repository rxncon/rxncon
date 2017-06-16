"""Module containing the class SBtabData, which can parse data in the tabular SBtab format. Includes classes
ValidatedSBtabData, which runs postprocessor functions to type-check the data in the table, and EntryBase
which serves as a parent class for the Entry class that gets dynamically generated based on the field names
that are presented. Also contains constructor function sbtab_data_from_file."""

from typing import List, Optional, Callable, Union, Dict, Any
import re


class SBtabData:
    def __init__(self, input: List[List[str]]) -> None:
        self.version = None  # type: Optional[str]
        self.entries = []  # type: List[EntryBase]
        self.document_name = None  # type: Optional[str]
        self.table_type = None  # type: Optional[str]
        self.table_name = None  # type: Optional[str]

        self._input = input
        self._column_names = []  # type: List[str]

        self._parse_header()
        self._parse_column_names()
        self._construct_entry_class()
        self._parse_entries()

    def _parse_header(self) -> None:
        REGEX_VERSION_A = 'SBtabVersion (\'|\").+?(\'|\")'
        REGEX_VERSION_B = 'SBtabVersion=(\'|\").+?(\'|\")'
        REGEX_DOCUMENT = 'Document=(\'|\").+?(\'|\")'
        REGEX_TABLE_TYPE = 'TableType=(\'|\").+?(\'|\")'
        REGEX_TABLE_NAME = 'TableName=(\'|\").+?(\'|\")'

        for col in self._input[0][1:]:
            assert not col.strip('\n')

        header = self._input[0][0]

        assert header.startswith('!!SBtab')

        match = re.search(REGEX_VERSION_A, header)
        if match:
            setattr(self, 'version', _unquote(match.group(0).split(' ')[1]))

        for regex, attr in zip([REGEX_VERSION_B, REGEX_DOCUMENT, REGEX_TABLE_TYPE, REGEX_TABLE_NAME],
                               ['version', 'document_name', 'table_type', 'table_name']):
            match = re.search(regex, header)
            if match:
                setattr(self, attr, _header_value(match.group(0)))

    def _parse_column_names(self) -> None:
        columns = self._input[1]
        assert len(columns) > 0

        self._column_names = [_cleaned_column_name(name) for name in columns]

    def _construct_entry_class(self) -> None:
        assert self.table_name
        assert self._column_names
        self._entry_class = type(_class_name_from_table_name(self.table_name), (EntryBase,),
                                 {'field_names': [_field_name_from_column_name(col) for col in self._column_names]})

    def _parse_entries(self) -> None:
        for row in self._input[2:]:
            entry = self._entry_class()

            for i, column_value in enumerate(row):
                setattr(entry, _field_name_from_column_name(self._column_names[i]), column_value.strip())

            self.entries.append(entry)


class ValidatedSBtabData(SBtabData):
    def __init__(self, input: List[List[str]], definition: SBtabData) -> None:
        super().__init__(input)

        self._definition = definition
        self._field_postprocessors = {}  # type: Dict[str, Callable[[str], Union[str, float, bool, int]]]
        self._construct_field_postprocessors()
        self._entry_class.field_postprocessors = {_field_name_from_column_name(col): func  # type: ignore
                                                  for col, func in self._field_postprocessors.items()}  # type: ignore
        self._postprocess_entries()

    def _construct_field_postprocessors(self) -> None:
        if hasattr(self._definition.entries[0], 'ComponentName'):
            type_definitions = {def_entry.ComponentName: def_entry.Format  # type: ignore
                                for def_entry in self._definition.entries
                                if def_entry.IsPartOf == self.table_type}  # type: ignore
        elif hasattr(self._definition.entries[0], 'Component'):
            type_definitions = {def_entry.Component: def_entry.Format      # type: ignore
                                for def_entry in self._definition.entries
                                if def_entry.IsPartOf == self.table_type}  # type: ignore
        else:
            raise AssertionError('Could not parse the definition file')

        for column in self._column_names:
            self._field_postprocessors[column] = _field_postprocessor_for_type_string(type_definitions[column])

    def _postprocess_entries(self) -> None:
        for entry in self.entries:
            entry.postprocess()


def sbtab_data_from_file(filename: str, separator: str = '\t', definitions: Optional[SBtabData] = None) -> SBtabData:
    sbtab_input = []

    with open(filename) as f:
        for row in f:
            if row.startswith('%'):
                continue
            else:
                sbtab_input.append(row.split(separator))

    if definitions:
        return ValidatedSBtabData(sbtab_input, definitions)
    else:
        return SBtabData(sbtab_input)


class EntryBase:
    def __getitem__(self, item: str) -> Any:
        return getattr(self, item)

    def postprocess(self) -> None:
        for field_name in self.field_names:  # type: ignore
            setattr(self, field_name, self.field_postprocessors[field_name](getattr(self, field_name)))  # type: ignore


def _unquote(x: str) -> str:
    return x.strip('\'\"')


def _header_value(header_definition: str) -> str:
    return _unquote(header_definition.split('=')[1])


def _cleaned_column_name(raw_name: str) -> str:
    assert raw_name.startswith('!') and not raw_name.startswith('!!')
    return raw_name[1:].strip()


def _class_name_from_table_name(table_name: str) -> str:
    for separator in [' ', '-', '_']:
        if separator in table_name:
            table_name = ''.join(x.capitalize() for x in table_name.split(separator))

    return table_name


def _field_name_from_column_name(column_name: str) -> str:
    return column_name.strip()


def _field_postprocessor_for_type_string(type_string: str) -> Callable[[str], Union[str, float, bool, int]]:
    if type_string is None or type_string == 'string':
        return lambda value: value
    elif type_string == 'float':
        return lambda value: float(value)
    elif type_string == 'Boolean':
        def str2bool(x: str) -> bool:
            if x in ['1', 'True']:
                return True
            elif x in ['0', 'False']:
                return False
            else:
                raise ValueError

        return str2bool
    elif type_string == 'integer':
        return lambda value: int(value)
    elif type_string == '{+,-,0}':
        def str2sign(x: str) -> str:
            if x not in ('+', '-', '0'):
                raise ValueError

            return x

        return str2sign
    else:
        raise TypeError

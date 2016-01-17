from typing import List


class SBtabData:
    def __init__(self, input: List[List[str]]):
        self.version = None
        self.entries = []
        self.document_name = None
        self.table_type = None
        self.table_name = None

        self._input = input
        self._column_names = []
        self._validation_functions = {}
        self._entry_class = None

        self._validate_input()
        self._parse_header()
        self._parse_column_names()
        self._construct_validation_functions()
        self._construct_entry_class()
        self._parse_entries()

        del self._input
        del self._column_names
        del self._validation_functions
        del self._entry_class

    def _validate_input(self):
        pass

    def _parse_header(self):
        assert len(self._input[0]) == 1
        header = self._input[0][0]

        assert header.startswith('!!SBtab')
        fields = header.split(' ')

        while fields:
            field = fields.pop()

            if field == 'SBtabVersion':
                self.version = _unquote(fields.pop())

            elif field.startswith('Document='):
                self.document_name = _header_value(field)

            elif field.startswith('TableType='):
                self.table_type = _header_value(field)

            elif field.startswith('TableName='):
                self.table_name = _header_value(field)

    def _parse_column_names(self):
        columns = self._input[1]
        assert len(columns) > 0

        self._column_names = [_cleaned_column_name(name) for name in columns]

    def _construct_validation_functions(self):
        for column in self._column_names:
            self._validation_functions[column] = lambda value: True

    def _construct_entry_class(self):
        self._entry_class = type(_class_name_from_table_name(self.table_name), (EntryBase,),
                                 {'validation_functions': {_field_name_from_column_name(col): func for col, func in self._validation_functions.items()},
                                  'field_names': [_field_name_from_column_name(col) for col in self._column_names]})

    def _parse_entries(self):
        for row in self._input[2:]:
            entry = self._entry_class()

            for i, column_value in enumerate(row):
                setattr(entry, _field_name_from_column_name(self._column_names[i]), column_value)

            entry.validate()
            self.entries.append(entry)


def sbtab_data_from_file(filename: str, separator='\t'):
    sbtab_input = []

    with open(filename) as f:
        for row in f:
            sbtab_input.append(row.split(separator))

    return SBtabData(sbtab_input)


class EntryBase:
    def validate(self):
        for field_name in self.field_names:
            assert self.validation_functions[field_name](getattr(self, field_name))


def _unquote(x: str):
    return x.strip('\'\"')


def _header_value(header_definition: str):
    return _unquote(header_definition.split('=')[1])


def _cleaned_column_name(raw_name: str):
    assert raw_name.startswith('!') and not raw_name.startswith('!!')
    return raw_name[1:]


def _class_name_from_table_name(table_name: str):
    for separator in [' ', '-', '_']:
        if separator in table_name:
            table_name = ''.join(x.capitalize() for x in table_name.split(separator))

    return table_name


def _field_name_from_column_name(column_name: str):
    return column_name.strip()

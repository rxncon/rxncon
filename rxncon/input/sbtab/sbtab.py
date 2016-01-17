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
        self._field_parsers = {}
        self._entry_class = None

        self._parse_header()
        self._parse_column_names()
        self._construct_field_parser()
        self._construct_entry_class()
        self._parse_entries()

        del self._input
        del self._column_names
        del self._field_parsers
        del self._entry_class

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

    def _construct_field_parser(self):
        for column in self._column_names:
            self._field_parsers[column] = lambda value: value

    def _construct_entry_class(self):
        self._entry_class = type(_class_name_from_table_name(self.table_name), (EntryBase,),
                                 {'field_parsers': {_field_name_from_column_name(col): func for col, func in self._field_parsers.items()},
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
            setattr(self, field_name, self.field_parsers[field_name](getattr(self, field_name)))


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


def _field_parser_for_type_string(type_string: str):
    if type_string is None or type_string == 'string':
        return lambda value: value

    elif type_string == 'float':
        return lambda value: float(value)

    elif type_string == 'Boolean':
        def str2bool(x: str):
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
        def str2sign(x: str):
            if x not in ['+', '-', '0']:
                raise ValueError

            return x
        return str2sign

    else:
        raise TypeError








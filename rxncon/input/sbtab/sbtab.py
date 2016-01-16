from typing import List


class SBtabFile:
    def __init__(self, input: List[List[str]]):
        self.version = None
        self.entries = []
        self.document_name = None
        self.table_type = None
        self.table_name = None

        self._input = input
        self._column_names = None

        self._validate_input()
        self._parse_header()
        self._parse_column_names()
        self._construct_entry_class()
        self._parse_entries()

        del self._input
        del self._header_info
        del self._column_names

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






def _unquote(x: str):
    return x.strip('\'\"')


def _header_value(header_definition: str):
    return _unquote(header_definition.split('=')[1])


def _cleaned_column_name(raw_name: str):
    assert raw_name.startswith('!') and not raw_name.startswith('!!')
    return raw_name[1:]

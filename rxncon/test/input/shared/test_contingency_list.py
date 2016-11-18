from rxncon.input.shared.contingency_list import contingency_list_entry_from_strs as cle_from_str, contingencies_from_contingency_list_entries

def test_nested_boolean():
    cles = [
        cle_from_str('<C1>', 'AND', 'A_[x]--B_[y]'),
        cle_from_str('<C1>', 'AND', 'A_[(r)]-{p}'),
        cle_from_str('<C2>', 'AND', 'B_[z]--D_[y]'),
        cle_from_str('<C2>', 'AND', 'B_[(r1)]-{p}'),
        cle_from_str('<C2>', 'AND', 'B_[(r2)]-{p}'),
        cle_from_str('<C1C2>', 'AND', '<C1>'),
        cle_from_str('<C1C2>', 'AND', '<C2>'),
        cle_from_str('A_[q]_ppi+_Q_[a]', '!', '<C1C2>')
    ]

    contingencies = contingencies_from_contingency_list_entries(cles)

    for x in contingencies:
        print(x)

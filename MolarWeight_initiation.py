import numpy as np
import json

molarmass = {}

def evaluate_term3_molemass_delta(term3):
    term3 = term3.strip("[").replace("]", "").split("(")
    mmass = term3[0]
    value = float(mmass)
    if len(term3) == 2:
        delta = term3[1].strip(")")
        nr_decimal = len(mmass.split(".")[-1])
        delta = int(delta) / (10 ** nr_decimal)
        return value, round(value-delta, nr_decimal), round(value+delta, nr_decimal)  # value, min, max
    else:
        return value, np.nan, np.nan

with open("MolarWeight.txt", "r") as f:
    for line in f.readlines():
        if line[0] != "$":
            terms = line.split()
            elem_symbol = terms[1]
            refNr = "-" if len(terms) < 5 else terms[-1]
            _val, _min, _max = evaluate_term3_molemass_delta(terms[3])
            molarmass[elem_symbol] = {"AtomNr": int(terms[0]), "Symbol": elem_symbol, "Name": terms[2],
                                      "MMass": _val, "MinMass": _min, "MaxMass": _max, "RefNr": refNr}

def evaluate_term3_molemass_minmax():
    _min, _max = float(terms[3]), float(terms[4])
    nr_decimal = len(terms[3].split(".")[-1])
    _val = round((_min + _max) / 2, nr_decimal)
    return _val, _min, _max

with open("MolarWeight_Table1.txt", "r") as f:
    for line in f.readlines():
        if line[0] != "$":
            terms = line.split()
            elem_symbol = terms[1]
            _val, _min, _max = evaluate_term3_molemass_minmax()
            molarmass[elem_symbol] = {"AtomNr": int(terms[0]), "Symbol": elem_symbol, "Name": terms[2],
                                      "MMass": _val, "MinMass": _min, "MaxMass": _max, "RefNr": "-"}

for value in molarmass.values():
    for val in value.values():
        print("{:>16}".format(val), end=" ")
    print()

with open("MolarWeight_Refs.txt", "r") as f:
    info_lines = f.readlines()[2:]   # get rid of comments lines
    molarmass["References"] = dict([term.split(":") for term in info_lines])

with open("molarmass_v2.json", "w") as f:
    json.dump(molarmass, f, indent=4)

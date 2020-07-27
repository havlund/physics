import json

class Elements:
    def __init__(self):
        pass

    def MolarMass(self):
        return MolarMass()

    def lookup_name_by_symbol(self, symbol):
        return self.MolarMass().get_molar_mass_elem_as_obj(symbol).get_name()

    def lookup_names_by_symbols(self, symbols):
        return [self.MolarMass().get_molar_mass_elem_as_obj(symbol).get_name() for symbol in symbols]


class MMelem:
    def __init__(self, AtomNr, Symbol, Name, MMass, RefNr, Ref, MinMass, MaxMass):
        self.AtomNr = AtomNr
        self.Symbol = Symbol
        self.Name = Name
        self.MMass = MMass
        self.Ref = Ref
        self.RefNo = RefNr
        self.MinMass = MinMass
        self.MaxMass = MaxMass

    def get_value(self):
        return self.MMass

    def get_reference(self):
        return self.Ref

    def get_boundries(self):
        return (self.MinMass, self.MaxMass)

    def get_name(self):
        return self.Name

class MolarMass:
    def __init__(self):
        self.json = "molarmass.json"
        with open(self.json, "r") as f:
            self.DB = json.load(f)

    def get_MolarMass_DB(self):
        return self.DB

    def _get_reference(self, refNr):
        refNr = str(refNr)
        refDct = self.DB["References"]
        try:
            return refNr if refNr == "-" else refDct[refNr]
        except ValueError as ex:
            raise ValueError("The RefNr '{}' cannot be identified!".format(refNr), ex)

    def get_molar_mass_elem_as_obj(self, elem):
        try:
            elemDct = self.DB[elem]
            refnr = elemDct["RefNr"]
            obj = MMelem(**elemDct, Ref=self._get_reference(refNr=refnr))
            return obj
        except KeyError as ex:
            raise ValueError(elem, "is not found in the DB.", ex)

    def get_molar_mass_elem(self, elem):
        try:
            return self.DB[elem]["MMass"]
        except KeyError as ex:
            print(elem, "is not found in the DB.\n", ex)
            return None

    def get_molar_mass_group(self, elems):
        return [self.get_molar_mass_elem(elem) for elem in elems]

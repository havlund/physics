import pandas as pd
import numpy as np

"""
version 1 is based on json
version 2 is based on pandas
"""

class Elements:
    def __init__(self):
        self.df = pd.read_excel("elements.xlsx", index_col="Symbol", sheet_name="data") # unique. as DB

    def get_molar_mass_elem(self, elem):
        return self.df.at[elem.title(), "MMass"]

    def get_molar_mass_group(self, elems):
        return self.df["MMass"][elems]

    def get_molar_mass_obj(self, elem):
        return MolarMass(self.df.loc[elem.title()])

    def lookup_name_by_symbol(self, symbol):
        return self.df.at[symbol.title(), "Name"]

    def lookup_names_by_symbols(self, symbols):
        symbols = [symb.title() for symb in symbols]
        return self.df["Name"][symbols]

    def lookup_symbol_by_name(self, name):
        arr = self.df["Name"].values
        idx = np.where(arr == name.title())
        return self.df.iloc[idx].index[0]

    def lookup_symbols_by_names(self, names):
        df_swt = self.df.set_index("Name")       # switch "Name" and "Symbol"; use "Name" as index
        df_swt["Symbol"] = self.df.index         # bring "Symbol" back as a column
        names = [name.title() for name in names]
        return df_swt["Symbol"][names]

class MolarMass:
    def __init__(self, dat_ser: pd.Series):
        self.dat_ser = dat_ser

    def get_reference(self):
        RefNr = self.dat_ser["RefNr"]
        ref_ser = pd.read_excel("elements.xlsx", index_col="quantity", sheet_name="ref").loc["MolarMass"]
        try:
            return RefNr if RefNr == "-" else ref_ser[RefNr]
        except ValueError as ex:
            raise ValueError("The RefNr '{}' cannot be identified!".format(RefNr), ex)

    def get_value(self):
        return self.dat_ser["MMass"]

    def get_boundries(self):
        MinMass = self.dat_ser["MinMass"]
        MaxMass = self.dat_ser["MaxMass"]
        return (MinMass, MaxMass)

if __name__ == "__main__":
    print(Elements().lookup_names_by_symbols(["Al", "Os"]))
    print(Elements().lookup_name_by_symbol("Cu"))
    print(Elements().get_molar_mass_group(["Al", "Os"]))
    print(Elements().get_molar_mass_elem("Cu"))
    obj = Elements().get_molar_mass_obj("Cu")
    print(obj.get_value())
    print(obj.get_boundries())
    print(obj.get_reference())
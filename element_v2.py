import pandas as pd
import numpy as np

"""
version 1 is based on json
version 2 is based on pandas
version 2.1 adds atom sizes
"""

class Elements:
    def __init__(self):
        self.dfs = pd.read_excel("elements.xlsx", sheet_name=None)
        self.dat_df = self.dfs["data"]
        self.dat_df.set_index("Symbol", inplace=True)
        self.ref_df = self.dfs["ref"]
        self.ref_df.set_index("quantity", inplace=True)

    def get_molar_mass_elem(self, elem):
        return self.dat_df.at[elem.title(), "MMass"]

    def get_molar_mass_group(self, elems):
        return self.dat_df["MMass"][elems]

    def get_atomic_radius(self, elem):
        return self.dat_df.at[elem.title(), "AtomicSize"]

    def get_atomic_radii_group(self, elems):
        return self.dat_df["AtomicSize"][elems]

    def get_covalent_radius(self, elem):
        return self.dat_df.at[elem.title(), "CovalentRadius"]

    def get_covalent_radii_group(self, elems):
        return self.dat_df["CovalentRadius"][elems]
    

    def get_molar_mass_obj(self, elem):
        return MolarMass(self.dat_df.loc[elem.title()])

    def lookup_name_by_symbol(self, symbol):
        return self.dat_df.at[symbol.title(), "Name"]

    def lookup_names_by_symbols(self, symbols):
        symbols = [symb.title() for symb in symbols]
        return self.dat_df["Name"][symbols]

    def lookup_symbol_by_name(self, name):
        arr = self.dat_df["Name"].values
        idx = np.where(arr == name.title())
        return self.dat_df.iloc[idx].index[0]

    def lookup_symbols_by_names(self, names):
        df_swt = self.dat_df.set_index("Name")       # switch "Name" and "Symbol"; use "Name" as index
        df_swt["Symbol"] = self.dat_df.index         # bring "Symbol" back as a column
        names = [name.title() for name in names]
        return df_swt["Symbol"][names]

    def update_DB_with_a_new_quantity(self, keyword: str, dat_ser: pd.Series):
        self.dat_df[keyword] = dat_ser
        self._write_dataframes()

    def _write_dataframes(self):
        writer = pd.ExcelWriter("elements.xlsx")
        for sht_name in self.dfs:
            self.dfs[sht_name].to_excel(writer, sht_name)
        writer.save()
        

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
    print(Elements().get_atomic_radius("Fe"))
    print(Elements().get_covalent_radii_group(["Al", "Cu", "Fr"]))

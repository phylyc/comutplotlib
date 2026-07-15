import pandas as pd


class MutationalSignatureSet(object):
    # The listed order determines the order in which the signatures are plotted
    signature_sets = {
        "clock-like": [
            "SBS1",
            "SBS5"
        ],
        "APOBEC": [
            "SBS2",
            "SBS13"
        ],
        "PolE/D/H": [
            "SBS9",
            "SBS10a",
            "SBS10b",
            "SBS10c",
            "SBS10d",
            "DBS3"
        ],
        "HRD": [
            "SBS3",
        ],
        "MMR": [
            "SBS6",
            "SBS8",
            "SBS14",
            "SBS15",
            "SBS20",
            "SBS21",
            "SBS26",
            "SBS30",
            "SBS36",
            "SBS44",
            "DBS7",
            "DBS10",
            "DBS13",
            "ID6",
            "ID7",
            "ID8",  # TOP2A
            "ID17",  # TOP2A
        ],
        "Other": [
            "SBS18",  # Oxog
            "SBS22a",  # Aristolochic acid
            "SBS22b",  # Aristolochic acid
            "SBS24",  # aflatoxin
            "SBS42",  # haloalkane
            "SBS85",  # ind eff of AID
            "SBS88",  # colibactin
            "SBS90",  # duocarmycin
            "SBS99",  # melphalan
            "DBS20",  # Aristolochic acid
            "ID23",  # Aristolochic acid
            "ID18",  # e.coli
        ],
        "Smoking": [
            "SBS4",
            "SBS29",
            "SBS92",
            "DBS2",
            "ID3"
        ],
        "UV": [
            "SBS7a",
            "SBS7b",
            "SBS7c",
            "SBS7d",
            "SBS38",
            "DBS1",
            "ID13"
        ],
        "Treatment": [
            "SBS17a", # gamma irradiation
            "SBS17b", # gamma irradiation
            "SBS40",  # gamma irradiation
            "SBS25",  # Chemotherapy
            "SBS31",  # Platinum chemotherapy
            "SBS35",  # Platinum chemotherapy
            "SBS86",  # Unknown chemotherapy
            "SBS87",  # Thiopurine chemotherapy
            "DBS5",  # Platinum chemotherapy
            "SBS11",  # Temozolomide
            "SBS32",  # Azathioprine
        ],
        "Error": [
            "SBS27",
            "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49",
            "SBS50", "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57", "SBS58", "SBS59",
            "SBS60",
            "SBS95",
            "DBS14"
        ],
        "Unknown": [
            "SBS41b",
            "SBS41c"
        ]
    }

    @classmethod
    def sort_signatures(cls, signatures: pd.Index) -> pd.Index:
        # Convert the index to a DataFrame for easier manipulation
        index_df = signatures.to_frame(name="signature")
        # Create a dictionary mapping each signature to its class and position within the class
        class_mapping = {}
        order_mapping = {}
        for c_idx, (class_name, sig_list) in enumerate(cls.signature_sets.items()):
            for s_idx, sig in enumerate(sig_list):
                class_mapping[sig] = c_idx
                order_mapping[sig] = s_idx
        # Add these mappings as columns to the DataFrame
        index_df["class"] = index_df["signature"].map(class_mapping)
        index_df["order"] = index_df["signature"].map(order_mapping)
        # Sort the DataFrame first by 'class' and then by 'order'
        sorted_df = index_df.sort_values(by=["class", "order"], ascending=[True, True])
        # Extract the sorted index
        sorted_index = sorted_df["signature"].to_list()
        return pd.Index(sorted_index, dtype=signatures.dtype, name=signatures.name)

    def __init__(self):
        pass

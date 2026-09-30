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
            "DBS13",
            "ID6",
        ],
        "MMR": [
            "SBS6",
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
            "SBS8",
            "SBS16",
            "SBS34",
            "SBS41",
            "SBS41b",
            "SBS41c"
        ]
    }

    @classmethod
    def get_signature_to_set_map(cls) -> dict[str, str]:
        """Map every known signature name onto the name of its etiology set.

        The etiology names themselves are included as identity entries so that
        grouping is idempotent: exposures that are already aggregated by
        etiology (or supplied that way) are passed through unchanged.
        """
        mapping = {set_name: set_name for set_name in cls.signature_sets}
        for set_name, sig_list in cls.signature_sets.items():
            for sig in sig_list:
                mapping[sig] = set_name
        return mapping

    @classmethod
    def sort_signature_sets(cls, signature_sets: pd.Index) -> pd.Index:
        """Sort etiology names by the order in which they are declared above.

        Unknown names are appended at the end, keeping their relative order
        (``sorted`` is stable), mirroring ``sort_signatures``.
        """
        order = {name: i for i, name in enumerate(cls.signature_sets)}
        sorted_sets = sorted(signature_sets, key=lambda s: order.get(s, len(order)))
        return pd.Index(sorted_sets, dtype=signature_sets.dtype, name=signature_sets.name)

    @classmethod
    def group_by_etiology(cls, signatures: pd.DataFrame | None) -> pd.DataFrame | None:
        """Aggregate signature exposures by their etiology (``signature_sets`` keys).

        Columns belonging to the same etiology are summed. Signatures that are
        not part of any set are kept as their own column (and sorted to the end)
        so no exposure is silently dropped.
        """
        if signatures is None:
            return None
        mapping = cls.get_signature_to_set_map()
        groups = pd.Index(
            [mapping.get(c, c) for c in signatures.columns],
            name=signatures.columns.name,
        )
        grouped = signatures.T.groupby(groups, sort=False).sum().T
        grouped.columns.name = signatures.columns.name
        return grouped.reindex(columns=cls.sort_signature_sets(grouped.columns))

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

"""Regression tests for hex-string colour support in :class:`Palette`.

Adding hex-string colours (e.g. ``"#cc66cc"``) as palette values used to break
``Palette.hash()`` because it iterated over the *characters* of the string
instead of RGB components, which then made ``from_hash`` raise
``ValueError: could not convert string to float: '#'``. Colours are now
normalized through ``matplotlib.colors.to_rgb`` so hex strings, named colours and
RGB tuples all serialize/round-trip uniformly.
"""

import matplotlib.colors as mc

from comutplotlib.palette import Palette


def test_to_rgb_normalizes_hex_and_named():
    assert Palette.to_rgb("#ff00ff") == mc.to_rgb("#ff00ff")
    assert Palette.to_rgb("black") == (0.0, 0.0, 0.0)
    assert Palette.to_rgb((0.1, 0.2, 0.3)) == (0.1, 0.2, 0.3)


def test_hash_roundtrip_with_hex_colors():
    palette = Palette({"a": "#cc66cc", "b": (0.0, 0.5, 1.0), "c": "black"})
    restored = Palette.from_hash(palette.hash())
    assert set(restored.keys()) == {"a", "b", "c"}
    for key in restored:
        assert list(restored[key]) == list(mc.to_rgb(palette[key]))


def test_condense_handles_hex_colors():
    # Two palettes with identical (hex) colours must condense to one entry.
    cmaps = {
        "X": Palette({"foo": "#cc66cc"}),
        "Y": Palette({"foo": "#cc66cc"}),
    }
    condensed = Palette.condense(cmaps)
    assert len(condensed) == 1
    (title,) = condensed.keys()
    assert set(title.split("\n")) == {"X", "Y"}


def test_mix_accepts_hex_strings():
    mixed = Palette.mix("#000000", "#ffffff", weight=0.5)
    assert all(abs(c - 0.5) < 1e-9 for c in mixed)


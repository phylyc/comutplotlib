import matplotlib.pyplot as plt

from comutplotlib.panel import Panel


class Layout(object):

    def __init__(self, xfigsize: float = 6.4, yfigsize: float = 4.8, panels: dict[str, Panel] = None, pad: int = 1, fig=None, gs=None, **kwargs):
        self.xfigsize = xfigsize
        self.yfigsize = yfigsize
        self.panels = panels if panels is not None else {}
        self.pad = pad
        self.fig = fig
        self.gs = gs
        self.kwargs = kwargs

    def add_panel(
        self, name, width, height, x=0, y=0,
        left_of: Panel = None, right_of: Panel = None, below: Panel = None, above: Panel = None,
        pad: int = None, align: str | float = "center",
        **kwargs
    ):
        """ Add a panel to the layout.
            The (x, y) coordinate is the top-left corner of the panel on the grid.
            x coordinate increases to the right.
            y coordinate increases downwards.
            The align parameter can be used to specify the alignment of the panel relative to the panel it is next to.
        """
        assert sum([left_of is not None, right_of is not None, below is not None, above is not None]) <= 1
        assert align in ["center", "left", "right", "top", "bottom"] or isinstance(align, float)

        pad = pad if pad is not None else self.pad

        if align == "center":
            shift = 0.5
        elif align in ["left", "top"]:
            shift = 0
        elif align in ["right", "bottom"]:
            shift = 1
        else:
            shift = align

        if left_of is not None:
            x, y = left_of.x - width - pad, left_of.y + int(shift * (left_of.height - height))
        elif right_of is not None:
            x, y = right_of.x + right_of.width + pad, right_of.y + int(shift * (right_of.height - height))
        elif below is not None:
            x, y = below.x + int(shift * (below.width - width)), below.y + below.height + pad
        elif above is not None:
            x, y = above.x + int(shift * (above.width - width)), above.y - height - pad

        panel = Panel(name, x, y, width, height, **kwargs)
        self.panels[name] = panel
        return panel

    def adjust_panel_coords(self):
        """ Adjust panel coordinates so that no panel has negative grid coordinates. """
        min_x = min([panel.x for panel in self.panels.values()])
        min_y = min([panel.y for panel in self.panels.values()])
        for panel in self.panels.values():
            panel.x -= min_x
            panel.y -= min_y

    def get_gridspec_shape(self):
        self.adjust_panel_coords()
        max_x = max([panel.x + panel.width for panel in self.panels.values()])
        max_y = max([panel.y + panel.height for panel in self.panels.values()])
        return max_y, max_x

    def place_panels_on_gridspec(self, autoscale_figsize: bool = False, scale: float = 1):
        max_y, max_x = self.get_gridspec_shape()
        if autoscale_figsize:
            self.xfigsize, self.yfigsize = max_x * scale, max_y * scale
        self.fig = plt.figure(figsize=(self.xfigsize, self.yfigsize))
        self.gs = self.fig.add_gridspec(nrows=max_y, ncols=max_x)
        for panel in self.panels.values():
            if panel.width and panel.height:
                try:
                    panel.ax = self.fig.add_subplot(self.gs[panel.y:panel.y + panel.height, panel.x:panel.x + panel.width])
                except Exception as e:
                    print(
                        f"Failed to place panel '{panel.name}' onto the grid "
                        f"of shape (ncols={self.gs.ncols}, nrows={self.gs.nrows}) "
                        f"at positions (cols={panel.y}:{panel.y + panel.height}, rows={panel.x}:{panel.x + panel.width})"
                    )
                    raise Exception(e)
        return self.fig, self.gs

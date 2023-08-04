import matplotlib.pyplot as plt
from matplotlib import rcParams, ticker
import numpy as np
import os
from typing import Iterable

from comutplotlib.palette import Palette


class Plotter(object):

    def __init__(self, output: str = "./plot.pdf", extra_palette: dict[str, tuple[float]] = None) -> None:
        super().__init__()
        self.out_dir, self.file_name = os.path.split(output)
        self.palette = Palette() | Palette(dict=extra_palette)

    def mk_out_path(self, recursive_folder_list: Iterable[str] = ()) -> str:
        """Make folder hierarchy from recursive_folder_list in self.out_dir
        and return full path.
        """
        path = os.path.join(self.out_dir, *recursive_folder_list)
        if len(path):
            os.makedirs(path, exist_ok=True)
        return path

    def save_figure(self, fig=None, ax=None, name=None, recursive_folder_list=(), **kwargs):
        _fig = fig if fig is not None else ax.get_figure()
        path = self.mk_out_path(recursive_folder_list=recursive_folder_list)
        name = name if name is not None else self.file_name
        if name.endswith(".png"):
            kwargs["dpi"] = kwargs.get("dpi", 600)
        _fig.savefig(os.path.join(path, name), **kwargs)

    @staticmethod
    def close_figure(fig=None, ax=None):
        _fig = ax.get_figure() if fig is None else fig
        plt.close(fig=_fig)

    # emulate grid using Axes.plot()
    @staticmethod
    def grid(
        ax=None,
        zorder=1.5,
        which=None,
        axis=None,
        alpha=None,
        color=None,
        linestyle=None,
        linewidth=None,
        **kwargs
    ):
        """From https://stackoverflow.com/a/66892353
        Caveats:
        - The x and y limits must be manually set before calling
          zorder_grid(). If specifying the ticks manually, that would also
          need to be done prior to plotting grid lines.
        - This will have issues if the axis limits and/or ticks change
          dynamically (e.g. animations). Possible solutions might be
          clearing the axis between frames or returning a Line2D list and
          toggling visibility.
        """
        # Honor rcParams values if keywords not specified
        if ax is None:
            ax = plt.gca()
        if which is None:
            which = rcParams["axes.grid.which"]
        if axis is None:
            axis = rcParams["axes.grid.axis"]
        if alpha is None:
            alpha = rcParams["grid.alpha"]
        if color is None:
            color = rcParams["grid.color"]
        if linestyle is None:
            linestyle = rcParams["grid.linestyle"]
        if linewidth is None:
            linewidth = rcParams["grid.linewidth"]

        # get coordinates for grid lines
        xlim = sorted(ax.get_xlim())
        ylim = sorted(ax.get_ylim())
        # ticks are sometimes assigned outside of limits, so remove them
        xticks = [t for t in ax.get_xticks() if xlim[0] < t < xlim[1]]
        yticks = [t for t in ax.get_yticks() if ylim[0] < t < ylim[1]]
        minor_xticks = [t for t in ax.get_xticks(minor=True) if xlim[0] < t < xlim[1]]
        minor_yticks = [t for t in ax.get_yticks(minor=True) if ylim[0] < t < ylim[1]]
        grid_xticks = []
        grid_yticks = []
        if which in ["major", "both"]:
            grid_xticks = np.concatenate((grid_xticks, xticks))
            grid_yticks = np.concatenate((grid_yticks, yticks))
        if which in ["minor", "both"]:
            grid_xticks = np.concatenate((grid_xticks, minor_xticks))
            grid_yticks = np.concatenate((grid_yticks, minor_yticks))

        # plot grid using Axes.plot()
        if axis in ["x", "both"]:
            for tick in grid_xticks:
                ax.axvline(
                    tick,
                    linestyle=linestyle,
                    color=color,
                    linewidth=linewidth,
                    alpha=alpha,
                    zorder=zorder,
                    **kwargs
                )
        if axis in ["y", "both"]:
            for tick in grid_yticks:
                ax.axhline(
                    tick,
                    linestyle=linestyle,
                    color=color,
                    linewidth=linewidth,
                    alpha=alpha,
                    zorder=zorder,
                    **kwargs
                )
        # ylim is being adjusted if something is being plotted until the edge.
        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())

    @staticmethod
    def set_integer_ticks(ax, xlim, xmin=None, n_major=3, n_minor=4):
        ax.set_xlim(xlim)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=n_major, steps=[1, 2, 3, 4, 5], integer=True))
        if xmin is not None:
            ax.set_xticks([t for t in ax.get_xticks() if t >= xmin])
            ax.set_xlim(xlim)
        if np.max(np.abs(xlim)) >= 2:
            # get minor tick interval as integer divisible of the major tick interval
            # that is closest to 4, first checking 5, then 3, then 6, then 2, etc.
            major_tick_interval = ax.get_xticks()[1]
            n = n_minor
            i = 0
            while True:
                if major_tick_interval % n == 0:
                    break
                i += 1
                n += i * (-1) ** (i + 1)
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=n))
            if xmin is not None:
                ax.xaxis.set_minor_locator(ticker.FixedLocator([t for t in ax.xaxis.get_minorticklocs() if t >= xmin]))
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    @staticmethod
    def set_spines(ax, **kwargs):
        for spine in ax.spines.values():
            spine.set(**kwargs)

    def no_spines(self, ax):
        self.set_spines(ax=ax, visible=False)

    def default_grid_layout(self, ax, axis=None, major=True, minor=True):
        if major:
            self.grid(ax=ax, axis=axis, which="major", zorder=0.5)
        if minor:
            self.grid(ax=ax, axis=axis, which="minor", zorder=0.1, color=self.palette.white)
        ax.set_facecolor(self.palette.backgroundgrey)
        self.no_spines(ax)

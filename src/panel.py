class Panel(object):

    def __init__(self, name, x, y, width, height, ax=None, plot_func=None, **kwargs):
        self.name = name
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.ax = ax
        self.plot_func = self.set_plot_func(plot_func, **kwargs) if plot_func is not None else None
        self.kwargs = kwargs

    def set_plot_func(self, plot_func, *args, **kwargs):
        self.plot_func = lambda ax: plot_func(ax, *args, **kwargs)

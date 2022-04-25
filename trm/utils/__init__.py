def style_mpl_axes(axes):
    """
    I don't like the default matplotlib axes as I prefer to
    have ticks marks going into the axes and for all 4 axes to
    have them (PGPLOT style). This routine applies fixes to accomplish
    this.

    Arguments::

      axes : :class:`~matplotlib.axes.Axes`
         the Axes object which will be updated.
    """
    axes.tick_params(axis="x", direction="in")
    axes.tick_params(axis="y", direction="in")
    axes.tick_params(bottom=True, top=True, left=True, right=True)

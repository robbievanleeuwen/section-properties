import matplotlib.pyplot as plt


def setup_plot(ax, pause):
    """Exectues code required to set up a matplotlib figure.

    :param ax: Axes object on which to plot
    :type ax: :class:`matplotlib.axes.Axes`
    :param bool pause: If set to true, the figure pauses the script until
        the window is closed. If set to false, the script continues
        immediately after the window is rendered.
    """

    if not pause:
        plt.ion()
        plt.show()

    ax.set_aspect("equal")


def finish_plot(pause):
    """Executes code required to finish a matplotlib figure.

    :param bool pause: If set to true, the figure pauses the script until
        the window is closed. If set to false, the script continues
        immediately after the window is rendered.
    """

    if pause:
        plt.show()
    else:
        plt.draw()
        plt.pause(0.001)

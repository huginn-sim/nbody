# -*- coding: utf-8 -*-
"""
.. module:: viz
   :synopsis: Provides a general interface for formatting plots of a multiple X,Y values.

.. moduleauthor:: Evin Özer
"""

#~ Modules
import  matplotlib.pyplot as plt, \
        numpy as np
#/~ Modules

#~ Functions
def configure(ax, title='', xlabel='', ylabel='', xbounds=None, ybounds=None, colors=('k','k'), suppress=False):
    """ Configures the plot with the specified parameters.

        :param ax: The current subplot to manipulate.
        :param title: The title of the current subplot.
        :param xlabel: The label of the x-axis.
        :param ylabel: The label of the y-axis.
        :param xbounds: The bounds of the x-axis.
        :param ybounds: The bounds of the y-axis.
    """
    if ax == None: return

    ax.set_axisbelow(True)
    ax.set_ylabel(ylabel, size=18)

    # Set the limits of the x,y axes.
    if xbounds == None:
        xbounds = ax.get_xlim()
    if ybounds == None:
        ybounds = ax.get_ylim()

    x_offset = abs(xbounds[1] - xbounds[0])/50.
    y_offset = abs(ybounds[1] - ybounds[0])/50.
    ax.set_xlim((xbounds[0] - x_offset, xbounds[1] + x_offset))
    ax.set_ylim((ybounds[0] - y_offset, ybounds[1] + y_offset))

    # Hide spines.
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


    # Set title, gridlines, and labels.
    if not suppress:
        ax.set_title(title, size=20)
        ax.grid(True, "major", alpha=.6)
        ax.grid(True, "minor", alpha=.4)

        ax.set_xlabel(xlabel, size=18)

        # X - Position and set major/minor ticks.
        ax.xaxis.set_ticks_position('bottom')
        offset = ax.get_xticks()[1] - ax.get_xticks()[0]
        ax.set_xticks(np.arange(ax.get_xticks()[1]+offset/2., ax.get_xticks()[-2], offset), minor=True)
        ax.xaxis.set_tick_params(width=4, length=5, color=colors[0])
        ax.xaxis.set_tick_params(which="minor", width=2, length=4, color=colors[0])

        # Draw x axis.
        x_edge = (xbounds[1] - xbounds[0]) / ((xbounds[1] - xbounds[0]) + x_offset)
        ax.axhline(ybounds[0] - y_offset, 1-x_edge, x_edge, lw=5, color=colors[0], alpha=1.)
    
    y_edge = (ybounds[1] - ybounds[0]) / ((ybounds[1] - ybounds[0]) + y_offset)
    # Y - Position and set major/minor ticks.
    if not suppress:
        ax.yaxis.set_ticks_position('left')
        ax.axvline(xbounds[0] - x_offset, 1-y_edge, y_edge, lw=4, color=colors[1], alpha=1.)
    else:
        ax.yaxis.set_ticks_position('right')
        ax.axvline(xbounds[1] + x_offset, 1-y_edge, y_edge, lw=4, color=colors[1], alpha=1.)

    offset = ax.get_yticks()[1] - ax.get_yticks()[0]
    ax.set_yticks(np.arange(ax.get_yticks()[1] + offset/2., ax.get_yticks()[-2], offset), minor=True)
    ax.yaxis.set_tick_params(width=4, length=5, color=colors[1])
    ax.yaxis.set_tick_params(which="minor", width=2, length=4, color=colors[1])
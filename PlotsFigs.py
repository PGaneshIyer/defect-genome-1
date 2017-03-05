import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# from matplotlib.patches import Ellipse, Polygon, Rectangle
import matplotlib.colors as colors
# from matplotlib import cm


class PlotsFigs:
    def plot_heatmap_of_matrix(self, arr1, figPath, xLabel=None, yLabel=None,
                               cmap='jet', figType='svg', dpi=300):
        maskedarr = np.ma.array(arr1,
                                mask=np.isnan(arr1))
        minVal = np.nanmin(maskedarr)
        maxVal = np.nanmax(maskedarr)
        fig, ax = plt.subplots()
        heatmap = ax.pcolor(maskedarr,
                            norm=colors.Normalize(vmin=minVal, vmax=maxVal),
                            cmap=cmap, edgecolor='black',
                            linestyle='-', lw=0.5,
                            vmin=minVal, vmax=maxVal)
        ax.patch.set(hatch='xx', edgecolor='black')
        p = patches.Rectangle((0, 0), 25, 15, hatch='xx',
                              fill=None, zorder=-10)
        ax.add_patch(p)
        fig.colorbar(heatmap)
        # put the major ticks at the middle of each cell
        ax.set_xticks(np.arange(maskedarr.shape[1])+0.5, minor=False)
        ax.set_yticks(np.arange(maskedarr.shape[0])+0.5, minor=False)

        # want a more natural, table-like display
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        if xLabel is not None:
            ax.set_xticklabels(xLabel, minor=False, rotation='vertical',
                               fontsize=8)
        if yLabel is not None:
            ax.set_yticklabels(yLabel, minor=False, fontsize=8)
        plt.tight_layout()
        plt.savefig(figPath, format=figType, dpi=dpi)
        plt.close()

    def plot_corr_matrix(self, corrMatrix, xLabels=None, yLabels=None,
                         filePath=None, saveFig=False, fileFormat='eps',
                         allow_masked=False, vmin=None, vmax=None,
                         cmap='jet', figsize=None,
                         rotateXLabel=0, rotateYLabel=0):
        fig, ax = plt.subplots(figsize=figsize)
        if allow_masked is True:
            if vmax is None:
                vmax = np.nanmax(corrMatrix)
            if vmin is None:
                vmin = np.nanmin(corrMatrix)
        else:
            if vmin is None:
                vmin = np.amin(corrMatrix)
            if vmax is None:
                vmax = np.amax(corrMatrix)
        heatmap = ax.pcolor(corrMatrix, cmap=cmap,
                            vmin=vmin, vmax=vmax)
        fig.colorbar(heatmap)
        # put the major ticks at the middle of each cell
        ax.set_xticks(np.arange(corrMatrix.shape[1])+0.5, minor=False)
        ax.set_yticks(np.arange(corrMatrix.shape[0])+0.5, minor=False)
        # want a more natural, table-like display
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        # row/column labels
        if (xLabels is not None):
            ax.set_xticklabels(xLabels, rotation=rotateXLabel, minor=False)
        if (yLabels is not None):
            ax.set_yticklabels(yLabels, rotation=rotateYLabel, minor=False)
        if (saveFig is True and filePath is not None):
            plt.savefig(filePath, format=fileFormat, dpi=1000)
            plt.close()
        elif (saveFig is True and filePath is None):
            print "WARNING: filePath not provided"
            print "WARNING: Returning plt"
            return plt
        else:
            return plt

    # OCT, 20, 2016 - FUNCTION NOT TESTED YET
    # FIXME: Make this general with *args and **kwargs
    # FIXME: The removeIndex must be a separate function
    #        (ex: get_value_in_range
    #        to handle multiple arrays and return both value and index)
    def plot2D(self, xVal, yVal, xminVal=None, xmaxVal=None,
               yminVal=None, ymaxVal=None, xLabel=None,
               yLabel=None, xLim=None, yLim=None,
               figOpt=[80, '^', 'g'], figType='eps',
               figFile=None, figShow=False):
        fig, ax = plt.subplots()
        plt.tick_params(axis='both', which='major',
                        length=8, width=2, labelsize=16)
        scatplot = ax.scatter(xVal, yVal, s=figOpt[0],
                              marker=figOpt[1], color=figOpt[2])
        ax.set_xlabel(xLabel, fontsize=20)
        ax.set_ylabel(yLabel, fontsize=20)
        ax.set_xlim(xLim)
        ax.set_ylim(yLim)
        if figFile:
            plt.savefig(figFile+'.'+figType,
                        bbox_inches='tight',
                        format=figType)
        if figShow:
            plt.show()
        else:
            plt.close()

    # OCT, 20, 2016 - FUNCTION NOT TESTED YET
    # FIXME: Make this general with *args and **kwargs
    # FIXME: The removeIndex must be a separate function
    #        (ex: get_value_in_range
    #        to handle multiple arrays and return both value and index)
    def plot2DWithColor(self, xVal, yVal, zVal,
                        xminVal=None, xmaxVal=None,
                        yminVal=None, ymaxVal=None,
                        zminVal=None, zmaxVal=None,
                        xLabel=None, yLabel=None, zLabel=None,
                        xLim=None, yLim=None, zLim=None,
                        figOpt=[200, 'o', 'jet'],
                        figType='eps', figFile=None, figShow=False):
        fig, ax = plt.subplots()
        plt.tick_params(axis='both', which='major',
                        length=8, width=3, labelsize=16)
        scatplot = ax.scatter(xVal, yVal, c=zVal, s=figOpt[0],
                              marker=figOpt[1], cmap=figOpt[2])
        ax.set_xlabel(xLabel, fontsize=20)
        ax.set_ylabel(yLabel, fontsize=20)
        ax.set_xlim(xLim)
        ax.set_ylim(yLim)
        # cbar = plt.colorbar(scatplot)
        # cbar.ax.set_ylabel(zLabel, rotation=90)
        # cbar.ax.tick_params(labelsize=14)
        if figFile:
            plt.savefig(figFile+'.'+figType,
                        bbox_inches='tight',
                        format=figType)
        if figShow:
            plt.show()
        else:
            plt.close()

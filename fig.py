import os
import json
import matplotlib.pyplot as plt

class PlotConfig():
    
    def __init__(self, config_name='slides_single', nrow=1, ncol=1, **kwargs):

        path = os.path.abspath(os.path.dirname(__file__))
        defaults = os.path.join(path, 'fig.json')

        with open(defaults, 'r') as f:
            configs = json.load(f)
            config = configs.get(config_name)

        for k, v in kwargs.items():
            config[k] = v

        self.plot_width = config.get('plot_width')
        self.margin_left = config.get('margin_left')
        self.margin_right = config.get('margin_right')
        self.margin_bottom = config.get('margin_bottom')
        self.margin_top = config.get('margin_top')
        self.space_width = config.get('space_width')
        self.space_height = config.get('space_height')
        self.subplot_ratio = config.get('subplot_ratio')
        self.ftsize = config.get('ftsize')

        self.nrow = nrow
        self.ncol = ncol

        self.subplot_width = ( self.plot_width 
                             - self.margin_left - self.margin_right 
                             - ( self.ncol - 1 ) * self.space_width
                             ) / self.ncol

        self.subplot_height = self.subplot_width * self.subplot_ratio

        self.plot_height = ( self.nrow * self.subplot_height
                           + self.margin_bottom + self.margin_top
                           + ( self.nrow -1 ) * self.space_height )

        font = {'family':'serif',
                'weight':'normal',
                'size':self.ftsize}

        # use TEX for interpreter
        plt.rc('text',usetex=True)

        plt.rc('text.latex', 
               preamble=[r'\usepackage{amsmath}',
                         r'\usepackage{bm}']
               )

        # use serif font
        plt.rc('font',**font)

    # cm inch transfer for matplotlib
    def __cm2inch(self, *tupl):
        inch = 2.54
        return tuple(i/inch for i in tupl)

    def get_fig(self, **kwargs):

        figsize = self.__cm2inch(self.plot_width, self.plot_height)

        fig = plt.figure(figsize=figsize, **kwargs)

        return fig

    def get_axes(self, fig, i=0, j=0, **kwargs):

        margin_height = ( self.margin_bottom
                         +(self.nrow-1-i)
                         *(self.space_height+self.subplot_height))
        margin_width = ( self.margin_left
                        +j*(self.space_width+self.subplot_width))

        rect = (margin_width/self.plot_width,
                margin_height/self.plot_height,
                self.subplot_width/self.plot_width,
                self.subplot_height/self.plot_height)

        ax = fig.add_axes(rect, **kwargs)

        return ax

    def get_simple(self):

        fig = self.get_fig()
        ax = self.get_axes(fig)

        return fig, ax

from .tools import *
import numpy as np


def Info():
    print("""
===== Thank you for using SiTePl =====
This library is being developped as a
       Simple Terminal Plotter.

Made by Me, for Me. 

Please Credit Leo BECHET if used.
""")



class Dataset:
    def __init__(self, x, y, marker, up_marker, color, label, type="curve", params={}):
        """
        if curve type, x,y coordinates.

        if histogram type, x=bins, y=bin edges
        
        
        """
        self.x = x
        self.y = y
        self.marker = marker
        self.up_marker = up_marker
        self.color = color
        self.label = label
        self.type = type
        self.params = params



class Color:
    # ANSI escape codes for colors
    RESET = "\033[0m"
    BLACK = "\033[30m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"







class Graph:
    def __init__(self, title=""):
        self.title = title

        self.datasets = []

        # Constructor

    def title(self, title):
        self.title = title
    

    def plot(self, x, y, marker="+", up_marker="|",color=Color.RESET, label="", type="curve"):
        self.datasets.append(Dataset(x, y, marker, up_marker, color, label, type))


    def show(self):
        """
        first 3 columns for numbers
        first rows for title
        last line for numbers        
        """
        # might support multiple lines one day, though it shouldn't be too hard
        rows, columns = get_terminal_size()
        # print(rows, columns)
        # graph_window = (rows-2-1, columns-3)
        graph_window = (rows-3, columns)
        # ^takes one line for the title, one is input, take another line for the numbers

        plot_screen = [[" " for i in range(graph_window[1])] for _ in range(graph_window[0])]


        # dynamic max/min size to respect scale inside graph
        graph_max = self.datasets[0].x[0]
        graph_min = self.datasets[0].x[0]
        for dataset in self.datasets:
            if graph_max< max(dataset.x):
                graph_max = max(dataset.x)
            if graph_min> min(dataset.x):
                graph_max = min(dataset.x)

        for dataset in self.datasets:
            # dataset = self.datasets[0]



            if dataset.type == "curve":
                # reduce curve, be careful of aliasing
                rx, ry = reduce_equidistant_2d(dataset.x, dataset.y, graph_window[1])



            elif dataset.type == "histogram":
                # y are bin edges. Not idiomatic ik
                # x_step = (dataset.y[-1]-dataset.y[0])/graph_window[0]
                rx = np.linspace(dataset.y[0], dataset.y[-1], graph_window[1])

                ry = []
                edge_index = 1
                for x in rx:
                    if x>dataset.y[edge_index]:
                        edge_index += 1

                    ry.append(dataset.x[edge_index-1])

                
            else:
                print(f'Invalid type, expected curve or histogram, got {dataset.type}')



            # in rx we have the step taken, average the thing and put | somewhere idk
            # the xlist used is a bunch of integers of the length of the window
            disp_x = [i for i in range(len(rx))]

            # need to scale ry to get the pixel on which to plot it
            foo = [i-graph_min for i in ry]
            foo = [i/graph_max*(graph_window[0]-1) for i in foo]
            disp_y = [   int(np.floor( i ))       for i in foo]



            for i in range(len(rx)):
                # print(i, disp_y[0], len(plot_screen))
                plot_screen[-disp_y[i]-1][disp_x[i]] = dataset.color+dataset.marker+dataset.color
                #                     ^ reverse to flip the screen, -1 or stuff in 0 is at the top
                
                # so we don't have cutout lines
                try:
                    a = disp_y[i]
                    b = disp_y[i+1]
                    start = min(a, b)
                    end = max(a, b)
                    
                    for j in range(start+1, end ):
                        plot_screen[-j-1][disp_x[i]] = dataset.color+dataset.up_marker+dataset.color
                        #            ^ reverse to flip the screen, -1 or stuff in 0 is at the top

                
                # just skip if we go out, handle case when we go out
                except IndexError:
                    # print("fail", i, len(disp_y))
                    pass


        # combines the screen
        lines = [''.join(row) for row in plot_screen]
        result = '\n'.join(lines)

        # adds the title
        padding = int(np.floor((columns - len(self.title)) /2)) * " "
        result = Color.RED + padding + self.title + padding + Color.RESET + "\n"+ result

        # adds start and end x values
        s_string = str(rx[0])
        e_string = str(rx[-1])
        padding = " " * (columns-(len(s_string) + len(e_string)))
        result = result + "\n" + s_string + padding + e_string


        print(result)



        
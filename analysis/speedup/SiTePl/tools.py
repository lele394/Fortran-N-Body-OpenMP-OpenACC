import os
import numpy as np

def get_terminal_size():
    rows, columns = os.popen('stty size', 'r').read().split()
    return int(rows), int(columns)





def reduce_equidistant_2d(x, y, m):
    """
    Reduces a list of 2D points described using separate x and y coordinates
    to m equidistant points while preserving the order.
    
    Parameters:
        x (list): List of x coordinates.
        y (list): List of y coordinates.
        m (int): Number of equidistant points to reduce to.
    
    Returns:
        reduced_x (list): List of reduced equidistant x coordinates.
        reduced_y (list): List of reduced equidistant y coordinates.
    """
    # Calculate total number of points
    num_points = len(x)
    
    # Calculate step size
    step = num_points / (m - 1)
    
    # Create index array for interpolation
    indices = np.arange(0, num_points, step=step)
    
    # Interpolate x and y coordinates separately
    reduced_x = np.interp(indices, np.arange(num_points), x)
    reduced_y = np.interp(indices, np.arange(num_points), y)
    
    return reduced_x.tolist(), reduced_y.tolist()
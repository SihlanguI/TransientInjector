import numpy as np


class TransientSources:
    
    def point_source(self, amplitude, x_size, y_size):
        """
        Create a delta function
        
        x_size : int
            The width of the 2D grid in pixels.
        y_size : int
            The height of the 2D grid in pixels.
        """
        # Define image size
        point_source = np.zeros((x_size, y_size))
        point_source[x_size//2, y_size//2] = amplitude  # Delta function at the center
        return point_source
    

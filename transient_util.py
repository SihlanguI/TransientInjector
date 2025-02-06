import numpy as np


class TransientSources:
    
    def point_source(self, amplitude, X_rot, Y_rot, sigma_major_pix, sigma_minor_pix):
        gaussian = amplitude * np.exp(-0.5 * ((X_rot / sigma_major_pix)**2 + (Y_rot / sigma_minor_pix)**2))
        return gaussian


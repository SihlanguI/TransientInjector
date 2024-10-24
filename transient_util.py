import numpy as np

def point_source(amplitude, X_rot, Y_rot, sigma_major_pix, sigma_minor_pix):
    gaussian = amplitude * np.exp(-0.5 * ((X_rot / sigma_major_pix)**2 + (Y_rot / sigma_minor_pix)**2))
    return gaussian

def create_fading_transient(amplitude, X_rot, Y_rot, sigma_major_pix, sigma_minor_pix, cube, scale_factor, t_indices):
    gaussian = elliptical_gaussian(amplitude, X_rot, Y_rot, sigma_major_pix, sigma_minor_pix)
    # Time decay: create a fading effect by reducing the amplitude at each time step
    for t in t_indices:
        decay_factor = np.exp(-scale_factor * (t / len(t_indices)))
        tran[t] = max_amplitude * decay_factor * gaussian

    return transient_cube



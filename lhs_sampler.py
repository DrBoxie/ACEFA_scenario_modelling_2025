import numpy as np

def lhs(params, N, rng = None):
    
# =============================================================================
# Latin Hypercube Sampling
# This is a classical LHS with the additional constraint that phi_V1 <= phi_S1
# and phi_V1 <= phi_V0. The paramater phi_V1 is therefore handled separate from
# the other parameters in the main loop.
# Input:
#       - params = Dictionary where the keys are the names of the parameters to 
#                  sample and the values are lists containing 2 elements
#                  representing the lower and upper bound of the range to 
#                  sample from.
#       - N = Number of particles
#       - rng = (optional) Random number generator
# Output:
#       - sample = Array of sample parameter values where each row corresponds
#                  to one such sample and each column corresponds to a 
#                  parameter, in the order in which they appeat in param_names.
#       - param_names = List of key names of sampled parameters (sorted 
#                       alphabetically).
# =============================================================================
    
    # If no random number generator has been given, create a new one
    if rng is None:
       rng = np.random.default_rng()
     
    # The sorted function ensures consistent ordering across versions of python, where .keys() is only order-preserving in newer versions of Python.
    param_names = sorted(params.keys())
    d = len(params)
    
    intervals = np.linspace(0, 1, N+1)
    samples = np.zeros((N, d))
    
    for j, name in enumerate(param_names):
        if name not in ["phi_V0", "phi_S1"]:                                                    # Skip phi_V0 and phi_S1 's for now, as they are treated in a separate loop further in the code
            low, high = params[name]
            # Sample one point per interval
            points = rng.uniform(intervals[:-1], intervals[1:], size=N)
            rng.shuffle(points)                                                     # Permute across this dimension. rng.shuffle automatically saves the output in the variable that was given to it, and there is therefore no need to assign it to a variable
            # Scale to parameter range
            samples[:, j] = low + points * (high - low)
    
    # Separate loop for phi_V0 and phi_S1    
    idx_S1, idx_V0, idx_V1 = param_names.index("phi_S1"), param_names.index("phi_V0"), param_names.index("phi_V1")
    for part in range(N):

        u_S1, u_V0 = rng.uniform(size = 2)
        samples[part, idx_S1] = samples[part, idx_V1] + (max(params["phi_S1"]) - samples[part, idx_V1]) * u_S1
        samples[part, idx_V0] = samples[part, idx_V1] + (max(params["phi_V0"]) - samples[part, idx_V1]) * u_V0
        
   
    return samples, param_names

if __name__ == "__main__":
    
    # These lines have no effect when calling the lhs function from outside 
    # this file, and are only meant here to allow for testing.
    N = 30
    beta_range = [0.01, 0.3]
    phi_S1_range = [0, 1]
    phi_V0_range = [0, 1]
    phi_V1_range = [0, 1]
    frac_S1_range = [0, 0.4]
    params = {
        "beta_range": beta_range,
        "phi_S1_range": phi_S1_range, 
        "phi_V0_range": phi_V0_range, 
        "phi_V1_range": phi_V1_range,
        "frac_S1_range":frac_S1_range
        }
    
    sample, param_names = lhs(params,  N)
    

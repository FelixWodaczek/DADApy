import numpy as np

def main():
    to_bohr = True
    mult = 1.
    if to_bohr:
        # Unit of eta is currently (1/A^2), this makes it (1/r_bohr^2)
        mult = (1./1.8897259886)**2.

    atom_types = ['H', 'O'] # Need to be in order of atomic numbers

    with open('n2p2_fitting/run_pot/input_mask.txt', 'r') as f: 
        input_header = f.read()
        f.close()

    g2_params = np.loadtxt('ice_in_water_data/g2_params_gridsearch_bohr_lambda.txt')
    g4_params = np.loadtxt('ice_in_water_data/g4_params_gridsearch_bohr_lambda.txt')

    select_gammas = True
    n_desc = 10
    lasso_gammas = np.load('trial_notebooks/lasso_gammas_hartbohr_lambda.npy')[n_desc]
    if not select_gammas:
        lasso_gammas = [1.0]*lasso_gammas.shape[0]
    if False: # For random generation
        lasso_gammas[:] = 0.
        lasso_gammas[np.random.choice(len(lasso_gammas), (n_desc), replace=False).astype(dtype=np.int16)] = 1.

    with open('n2p2_fitting/run_pot/input.nn', 'w') as f:
        f.write(input_header+'\n')

        counter = 0
        for to_at in atom_types:
            counter += 1 # add this for not included G1
            f.write(f'# radial to {to_at}\n')
            for eta, rshift in g2_params:
                if lasso_gammas[counter]:
                    for from_at in atom_types:
                        f.write(' '.join([
                            'symfunction_short',
                            from_at, '2', to_at, # from, type of function, to
                            f'{mult*eta:.3E}', f'{rshift:.3E}', '12.00', # eta, rshift, cutoff
                            '\n'
                        ]))
                counter+=1

        for ii_betw, betw_at in enumerate(atom_types):
            for ii_to, to_at in enumerate(atom_types):
                if ii_to >= ii_betw:
                    f.write(f'# angular to {betw_at} and {to_at}\n')
                    for eta, zeta, lambda_ in g4_params:
                        if lasso_gammas[counter]:
                            for from_at in atom_types:
                                f.write(' '.join([
                                    'symfunction_short',
                                    from_at, '3', betw_at, to_at, # from, type of function, middle atom, to
                                    f'{mult*eta:.3E}', f'{lambda_:.3E}', f'{zeta:.3E}', '12.00', # eta, lambda, zeta, cutoff
                                    '\n'
                                ]))
                        counter+=1
            
    print(f"Built input file with {counter} symmetry functions.")

if __name__ == '__main__':
    main()
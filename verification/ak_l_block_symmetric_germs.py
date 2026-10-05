"""Table 5: exact interval verification of all same-sign L-ray responses."""
def verify_ray_responses(d2_mass_resultants,d2_mass_rays):
    from smale_verify.d2_ray_response import replay
    return replay(d2_mass_resultants,d2_mass_rays)

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_l_block_symmetric_germs'))

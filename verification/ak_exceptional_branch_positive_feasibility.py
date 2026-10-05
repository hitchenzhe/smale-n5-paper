"""Table 5: global positive-mass L1/L8 ray and common-shape isolation."""
def verify_mass_rays(d2_mass_resultants):
    from smale_verify.d2_mass_rays import classify
    return classify(d2_mass_resultants)


if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_exceptional_branch_positive_feasibility'))

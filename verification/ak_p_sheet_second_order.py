"""Table 5: mass-parametric one-different fold exclusion (regular sheets pending)."""
def verify_general_fold():
    from smale_verify.d2_S5_blocks import verify_fold
    return verify_fold()

def verify_regular_X():
    from smale_verify.d2_regular_X import verify
    return verify()

def verify_regular_Euler(d2_mass_resultants):
    from smale_verify.d2_regular_euler import verify
    return verify(d2_mass_resultants)

def verify_regular_Lagrange():
    from smale_verify.d2_regular_lagrange import verify
    return verify()

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_p_sheet_second_order'))

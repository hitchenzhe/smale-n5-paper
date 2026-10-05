"""Table 5: exact ordinary and rank-19 fold jets at equal triangle masses."""
def verify_ordinary_seed():
    from smale_verify.d2_seed_jets import verify_ordinary
    return verify_ordinary()

def verify_fold_seed():
    from smale_verify.d2_seed_jets import verify_fold
    return verify_fold()

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_l5_obstruction'))

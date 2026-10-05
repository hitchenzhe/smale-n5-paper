"""Table 5 outer_census: full census and complete four-product coverage."""
def verify_elementary_four_products():
    from smale_verify.d5_elementary_four import verify
    return verify()

def verify_full_census():
    from smale_verify.d5_census import verify
    return verify()

def verify_dominant_products():
    from smale_verify.d5_dominant import verify
    return verify()

def verify_complete_four_products():
    from smale_verify.d5_four_complete import verify
    return verify()

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram5_outer_census'))

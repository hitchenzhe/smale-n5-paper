"""Table 5: exact response proof on the complete nonfold S5 mass sheet."""
def verify_S5_response():
    from smale_verify.d2_S5_blocks import verify_regular
    return verify_regular()

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_p_sheet_p5_exact'))

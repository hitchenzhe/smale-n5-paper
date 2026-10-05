"""Table 5 exceptional: D2 positive-mass polynomial exclusions."""
from smale_verify.d2_positivity import verify_positive_polynomial,save

def verify_K1_branches(d2_mass_resultants):
 from smale_verify.d2_K1_branches import verify
 return verify(d2_mass_resultants)

def verify_K4_arcs(d2_mass_resultants):
 from smale_verify.d2_K4_arcs import verify
 return verify(d2_mass_resultants)


def verify_K_and_L5_exclusions(d2_mass_resultants):
 result=d2_mass_resultants['resultants'];records=[]
 for j in (3,6,7,8):
  cover=verify_positive_polynomial(result[str(j)]['K_terms'])
  records.append(dict(polynomial='K'+str(j),sectors=cover))
 L5=result['5']['L_factors'][-1]['polynomial']
 if len(L5)!=91 or any(int(t['coefficient'])<=0 for t in L5):raise ArithmeticError('L5 prime coefficient positivity failed')
 equal_mass={j:sum(int(t['coefficient']) for t in r['L_resultant']) for j,r in result.items() if j!='5'}
 if any(not v for v in equal_mass.values()):raise ArithmeticError('Other shape meets Lambda=0 at equal masses')
 certificate=save('d2_positive_mass_exclusions.json',dict(K_covers=records,L5_prime_positive_coefficients=91,
   equal_mass_resultants=equal_mass,method='exact rational Bernstein coefficients with complete projective and corner-sector covers'))
 return [dict(check='D2.6.12.K3_K6_K7_K8_L5prime_and_equal_masses',status='PASS',items=12,
  detail=dict(certificate=certificate,K_polynomials=[3,6,7,8],
   coverage='three largest-mass charts; two order sectors each; full square covers with critical-corner blow-up',
   L5_prime_positive_coefficients=91,equal_mass_resultants=equal_mass))],{}


if __name__=='__main__':
 from smale_verify.runner import program_cli
 raise SystemExit(program_cli('ak_exceptional'))

"""Table 5: exact uniform-source resultant and Phi/Psi exclusions, Lemma 8.19."""
from fractions import Fraction as F
from smale_verify.exact import P,V,equal,det,substitute
from smale_verify.qpoly import gcd,sturm,variations


def verify_uniform_source_resultant():
    s3,s4,s5=map(V,('d5u_s3','d5u_s4','d5u_s5'))
    # Divide each cube difference by its nonzero difference of offsets.
    f=[P.coerce(3),-3*(s3+s4),s3**2+s3*s4+s4**2]
    g=[P.coerce(3),-3*(s3+s5),s3**2+s3*s5+s5**2]
    matrix=[f+[P()],[P()]+f,g+[P()],[P()]+g]
    expected=9*(s4-s5)**2*((s3-s4)**2+(s3-s5)**2+(s4-s5)**2)/2
    equal(det(matrix),expected,'D5 normalized forcing resultant')
    return [dict(check='D5.8.19.uniform_source_resultant',status='PASS',items=1,
        detail=dict(resultant='9(s4-s5)^2 * ((s3-s4)^2+(s3-s5)^2+(s4-s5)^2)/2',
            nonvanishing='three real pairwise distinct offsets',
            scope='part (a) single-source non-uniformity'))],{}


def verify_Phi_Psi_exclusions():
    b,r,omega=map(V,('d5u_b','d5u_r','d5u_omega'))
    mass=[b**2,b**2/r**2,r**-2]
    offsets=[P(),r,-r/b]
    records=[]
    for special in range(3):
        j,k=[j for j in range(3) if j!=special]
        for root in (omega,omega**2):
            # The two non-special forcings must agree. Their cube-root
            # ratio is root, with root^3=1 != root. Hence
            # theta=(s_j-root*s_k)/(1-root), and y_j/y_k=root.
            theta=(offsets[j]-root*offsets[k])*(2+root)/3
            theta=theta.reduce('d5u_omega',[1,1,1])
            equal((substitute(theta,{'d5u_omega':omega})*(1-root)
                  -offsets[j]+root*offsets[k]).reduce('d5u_omega',[1,1,1]),P(),
                  'D5 forcing coincidence theta')
            for kind in ('Phi','Psi'):
                if kind=='Phi':
                    obstruction=mass[special]**2*(theta-offsets[k])-(mass[j]*root+mass[k])**2*(theta-offsets[special])
                else:
                    obstruction=mass[special]**2*(theta-offsets[k])**3-(mass[j]+mass[k])**2*(theta-offsets[special])**3
                obstruction=obstruction.reduce('d5u_omega',[1,1,1]).clear('d5u_b').clear('d5u_r')
                degree=max(obstruction.exponents('d5u_r'),default=0)
                converted=P()
                for mon,coefficient in obstruction.terms.items():
                    exponents=dict(mon);power=exponents.pop('d5u_r',0)
                    converted+=P({tuple(sorted(exponents.items())):coefficient})*(b+1)**(degree-power)
                parts=[]
                for power in (0,1):
                    polynomial=converted.coefficient('d5u_omega',power)
                    if polynomial.names()-{'d5u_b'}:
                        raise AssertionError('unexpected gcd variable')
                    parts.append([polynomial.coefficient('d5u_b',n).terms.get((),F(0))
                                  for n in range(max(polynomial.exponents('d5u_b'),default=0)+1)])
                common=gcd(*parts)
                # Remove powers of b, since b>0, before using zero as a
                # Sturm endpoint. All denominator factors b,b+1 are units.
                while common and common[0]==0:common=common[1:]
                if not common:raise AssertionError('uniformity obstruction identically zero')
                chain=sturm(common)
                zero_variations=variations(chain,F(0))
                signs=[1 if row[-1]>0 else -1 for row in chain]
                infinity_variations=sum(x!=y for x,y in zip(signs,signs[1:]))
                if zero_variations!=infinity_variations:
                    raise AssertionError(('D5 uniformity has a positive-mass candidate',special,kind,common))
                records.append(dict(special_body=special+3,cube_root_power=1 if root is omega else 2,
                    kind=kind,component_coefficients_ascending=[[str(c) for c in part] for part in parts],
                    gcd_coefficients_ascending=[str(c) for c in common],
                    sturm_chain=[[str(c) for c in row] for row in chain],
                    positive_roots=0,variations_zero=zero_variations,variations_infinity=infinity_variations))
    if len(records)!=12:raise AssertionError('D5 uniformity case census')
    return [dict(check='D5.8.19.Phi_Psi_twelve_gcd_exclusions',status='PASS',items=12,
        detail=dict(records=records,method='necessary squared Phi/Psi identities; Q(omega) reduction; real/imaginary gcd; exact Sturm count',
            scope='part (c); the remaining parts of Lemma 8.19 follow from the explicit uniformity and order arguments in the manuscript'))],{}


if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram5_morse_bott_uniformity'))


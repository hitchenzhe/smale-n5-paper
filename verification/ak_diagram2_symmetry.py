"""Table 5: equal-mass D2 label counts, first obstruction and seed orbit."""
from itertools import product,permutations,combinations
from pathlib import Path
import json
from smale_verify.exact import V,substitute


def verify_equal_mass_labels():
    from smale_verify.d2_seed_jets import build
    from smale_verify.d2_fields import CF,s
    seed=build('fold');omega=CF.from_sympy(seed['omega']);one=CF.one
    def root_power(value):
        choices=[k for k in range(3)if not(value-omega**k)]
        if len(choices)!=1:raise ArithmeticError('Boundary coordinate is not a unique cube root')
        return choices[0]
    triangles=[]
    for exponents in permutations(range(3)):
        c=omega**((exponents[0]-exponents[2])%3)
        b=c+one-omega**((exponents[1]-exponents[2])%3)
        k=-b/CF.convert(2);w4=(k-one)/(k+CF.convert(2));w3=-w4-one
        phases=(root_power(w3),root_power(w4),0)
        if set(phases)!={0,1,2}:raise ArithmeticError('Fold shape is not the equal-mass S5 shape')
        triangles.append((tuple(exponents),phases))
    patterns=[]
    for signs in product((1,-1),repeat=6):
        S=sum((sgn*omega**phase for sgn,phase in zip(signs[:3],(1,2,0))),CF.zero)
        T=sum((sgn*omega**phase for sgn,phase in zip(signs[3:],(1,2,0))),CF.zero)
        zero=not S and not T
        if zero and not(signs[0]==signs[1]==signs[2]and signs[3]==signs[4]==signs[5]):raise ArithmeticError('Unexpected vanishing sign sums')
        if S and T:
            ratio=S/T
            choices=[(sign,power)for sign in(1,-1)for power in range(3)if not(ratio-sign*omega**power)]
            if len(choices)!=1:raise ArithmeticError('Mixed sign ratio not in +/- cube roots')
            sign,power=choices[0]
            # Re(omega**power) is 1 or -1/2, never zero. The two
            # positive real crossed amplitudes therefore cannot cancel.
            real=(omega**power+omega**((-power)%3))/CF.convert(2)
            if not real:raise ArithmeticError('First-obstruction phase argument failed')
        patterns.append(dict(signs=list(signs),Phi_zero_for_positive_pair=zero,
            proof='both sign sums zero'if zero else'one nonzero sum or nonzero real part of S/T'))
    if sum(p['Phi_zero_for_positive_pair']for p in patterns)!=4:raise ArithmeticError('Wrong first-obstruction sign census')
    compatible={(a,w,j,signs)for a,w in triangles for j in range(3)for signs in product((1,-1),repeat=2)}
    orbit=set();seed_edges={(3,4):0,(3,5):1,(4,5):2};seed_w={3:1,4:2,5:0}
    for perm,j,signs in product(permutations((3,4,5)),range(3),product((1,-1),repeat=2)):
        a=tuple((seed_edges[tuple(sorted((perm[i-3],perm[k-3])))]+j)%3 for i,k in combinations((3,4,5),2))
        w=tuple((seed_w[old]-seed_w[perm[2]])%3 for old in perm)
        orbit.add((a,w,j,signs))
    if orbit!=compatible or len(orbit)!=72:raise ArithmeticError('Rank-19 seed orbit does not exhaust compatible fold labels')
    # The paper's simultaneous cubic action has exact row covariance.
    w=V('orbit_omega');replacements={'t':w*V('t')}
    for j in (3,4,5):replacements['Z'+str(j)]=w*V('Z'+str(j))
    for j in(1,2):replacements['W'+str(j)]=w*V('W'+str(j))
    for p in((1,2),(3,4),(3,5),(4,5)):
        name='D'+''.join(map(str,p));replacements[name]=w*V(name)
    weights=[0,0,1,1,1,1,1,0,0,0]+[0]*12
    for row,power in zip(seed['raw'],weights):
        if (substitute(row,replacements)-w**power*row).reduce('orbit_omega',[1,1,1]).terms:
            raise ArithmeticError('Simultaneous cubic chart covariance failed')
    # The independent torus before gauge restoration: t'=r*t, Z_B'=r^-2 Z_B,
    # W_B'=r^2 W_B, crossed D'=r^-1 D. Coordinate rows have their matching
    # weights. A triangle permutation chooses r^2 W_perm(5)=1 locally.
    r=V('orbit_r');torus={'t':r*V('t')}
    for j in(3,4,5):
        torus['Z'+str(j)]=r**-2*V('Z'+str(j));torus['W'+str(j)]=r**2*V('W'+str(j))
    for i in(1,2):
        for j in(3,4,5):torus['D'+str(i)+str(j)]=r**-1*V('D'+str(i)+str(j))
    weights=[0,0,-2,-2,-2,0,0,2,2,2]+[0]*10
    for row,power in zip(seed['raw'][:20],weights):
        if (substitute(row,torus)-r**power*row).terms:raise ArithmeticError('Gauge-restoring chart torus covariance failed')
    payload=dict(check='D2.6.13.equal_mass_signs_and_72_fold_orbit',status='PASS',
        sign_patterns=patterns,fold_triangles=[dict(exponents=list(a),W_phases=list(w))for a,w in triangles],
        counts=dict(ordinary_Phi_zero=72,fold_total=1152,fold_Phi_zero=72,fold_Phi_nonzero=1080,total_in_lemma=1224),
        compatible_fold_orbit=[dict(exponents=list(a),W_phases=list(w),pair_phase=j,constant_row_signs=list(signs))for a,w,j,signs in sorted(orbit)],
        source_symmetries='Simultaneous cubic action and gauge-restoring torus checked row-by-row; permutations at equal triangle masses; reversal; complex conjugation followed by 3/4 swap. Gauges restore locally since 2(Z2-Z1)W5=2.',
        alternate_1080_exclusion='Phi is nonzero at all other 60 sign patterns for every positive pair mass. The exact torque tQ=0 excludes a non-boundary curve there. The separately verified finite boundary fibre excludes a boundary curve, hence the point is isolated. This replaces the manuscript rank-20/h^3 computation for those 1080 labels.',
        scope='Together with ordinary/general-fold seed checks and finite fibre reconstruction, covers every 1224 label of Lemma 6.13. No claim that the unused rank-20/h^3 identities were reproduced.')
    path=Path(__file__).resolve().parent/'certificates'/'d2_equal_mass_orbits.json';path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=1224,detail=dict(certificate=path.name,fold_orbit=72,first_order_exclusions=1080,ordinary_labels=72,scope=payload['scope']))],{}

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram2_symmetry'))

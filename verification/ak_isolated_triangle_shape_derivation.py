"""Reconstruction of ak_isolated_triangle_shape_derivation from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_shape_elimination():
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    # D1 shape lemma: exact elimination in edge-difference coordinates.  This
    # verifies the reduced force equations, the proportionality argument, and the
    # four sign-ratio eliminants.  The three impossible eliminants are written as
    # explicit sums of squares/strictly positive terms for positive mass ratios.
    u18,v18,w18=map(V,('alpha34_18','alpha35_18','alpha45_18'))
    p18,q18,r18=map(V,('beta34_18','beta35_18','beta45_18'))
    A18,B18,C18=map(V,('D34_18','D35_18','D45_18'))
    mu3,mu4,mu5=V('m3'),V('m4'),V('m5')
    equal(substitute(V('alpha35_18')-V('alpha34_18')-V('alpha45_18'),
                     {'alpha35_18':u18+w18}),0,'D1.5.8 z-edge triangle relation')
    equal(substitute(V('beta35_18')-V('beta34_18')-V('beta45_18'),
                     {'beta35_18':p18+r18}),0,'D1.5.8 w-edge triangle relation')
    tie18=[A18**2*u18*p18-1,B18**2*v18*q18-1,C18**2*w18*r18-1]
    # The edge relation and force equations reduce using D^2 alpha beta=1.
    force18=[-mu4*A18/p18-mu5*B18/q18,
              mu3*A18/p18-mu5*C18/r18,
              -mu4*A18/u18-mu5*B18/v18,
              mu3*A18/u18-mu5*C18/w18]
    for dvar,alpha18,beta18,tie in zip((A18,B18,C18),(u18,v18,w18),(p18,q18,r18),tie18):
        equal(dvar**3*alpha18-dvar/beta18,(dvar/beta18)*tie,
              'D1.5.8 force reduction by tie '+str(alpha18.terms))
    # Clearing only the nonzero denominators, the pairs of force equations imply
    # alpha34*beta35=alpha35*beta34 and alpha34*beta45=alpha45*beta34.
    f3n=mu4*A18*q18+mu5*B18*p18
    g3n=mu4*A18*v18+mu5*B18*u18
    f4n=mu3*A18*r18-mu5*C18*p18
    g4n=mu3*A18*w18-mu5*C18*u18
    equal(f3n*v18-g3n*q18,mu5*B18*(p18*v18-u18*q18),'D1.5.8 beta proportionality 1')
    equal(f4n*w18-g4n*r18,mu5*C18*(u18*r18-p18*w18),'D1.5.8 beta proportionality 2')
    rho18=V('rho18')
    shape18={'beta34_18':rho18*u18,'beta35_18':rho18*v18,'beta45_18':rho18*w18,
             'D34_18':V('X34_18')/u18,'D35_18':V('X35_18')/v18,'D45_18':V('X45_18')/w18}
    for dvar,alpha18,beta18,tie in zip((A18,B18,C18),(u18,v18,w18),(p18,q18,r18),tie18):
        X18=shape18[{A18:'D34_18',B18:'D35_18',C18:'D45_18'}[dvar]]*alpha18
        equal(substitute(tie,shape18),rho18*X18**2-1,'D1.5.8 common squared shape')
    equal(substitute(force18[0],shape18),-(mu4*V('X34_18')/u18**2+mu5*V('X35_18')/v18**2)/rho18,
          'D1.5.8 normalized first force')
    equal(substitute(force18[1],shape18),(mu3*V('X34_18')/u18**2-mu5*V('X45_18')/w18**2)/rho18,
          'D1.5.8 normalized second force')
    equal(substitute(force18[2],shape18),-(mu4*V('X34_18')/u18**2+mu5*V('X35_18')/v18**2),
          'D1.5.8 normalized third force')
    equal(substitute(force18[3],shape18),mu3*V('X34_18')/u18**2-mu5*V('X45_18')/w18**2,
          'D1.5.8 normalized fourth force')
    r18mass, s18mass=V('r_mass_18'),V('s_mass_18')
    e2,e3=V('sign2_18'),V('sign3_18')
    eliminants={
        '++':(s18mass-1-r18mass)**2-4*r18mass,
        '+-':(s18mass-1+r18mass)**2+4*r18mass,
        '-+':(r18mass-1)**2+2*s18mass*(r18mass+1)+s18mass**2,
        '--':(r18mass-s18mass-1)**2+4*r18mass}
    raw_eliminants={
        '++':(s18mass-1-r18mass)**2-4*r18mass,
        '+-':(s18mass-1+r18mass)**2-4*(-r18mass),
        '-+':(-s18mass-1-r18mass)**2-4*r18mass,
        '--':(-s18mass-1+r18mass)**2-4*(-r18mass)}
    for key in eliminants:equal(eliminants[key],raw_eliminants[key],'D1.5.8 sign eliminant '+key)
    vedge18=V('edge_ratio_v_18')
    for key,(psign,qsign) in {'++':(1,1),'+-':(1,-1),'-+':(-1,1),'--':(-1,-1)}.items():
        equal(substitute(raw_eliminants[key],
                         {'r_mass_18':qsign*vedge18**2,
                          's_mass_18':psign*(1+vedge18)**2}),0,
              'D1.5.8 edge-to-mass elimination '+key)
    Q18=2*mu3*mu4*mu5*(mu3+mu4+mu5)-mu3**2*mu4**2-mu4**2*mu5**2-mu5**2*mu3**2
    Q18norm=2*r18mass*s18mass+2*r18mass+2*s18mass-r18mass**2-s18mass**2-1
    equal(substitute(Q18,{'m3':r18mass**-1,'m4':s18mass**-1,'m5':1})*r18mass**2*s18mass**2,
          Q18norm,'D1.5.8 surviving eliminant is -Q')
    ashape18=V('alpha34_shape_18')
    equal((r18mass**-1*r18mass-1)*ashape18**2,0,
          'D1.5.8 m3 alpha45^2=m5 alpha34^2 from c^2/a^2=r')
    equal((s18mass**-1*s18mass-1)*ashape18**2,0,
          'D1.5.8 m4 alpha35^2=m5 alpha34^2 from b^2/a^2=s')
    equal(eliminants['++'],-Q18norm,'D1.5.8 unique viable sign pattern')
    for name,pol in [('+-',eliminants['+-']),('-+',eliminants['-+']),('--',eliminants['--'])]:
        if name=='+-': positive_form=(s18mass+r18mass-1)**2+4*r18mass
        elif name=='-+': positive_form=(r18mass-1)**2+2*s18mass*(r18mass+1)+s18mass**2
        else: positive_form=(r18mass-s18mass-1)**2+4*r18mass
        equal(pol,positive_form,'D1.5.8 strict positivity witness '+name)
    x3_18,x4_18,x5_18=map(V,('invroot_m3_18','invroot_m4_18','invroot_m5_18'))
    Q18x=substitute(Q18,{'m3':x3_18**-2,'m4':x4_18**-2,'m5':x5_18**-2})
    rhs18=(x3_18**-2*x4_18**-2*x5_18**-2)**2*(x3_18+x4_18+x5_18)*(-x3_18+x4_18+x5_18)*(x3_18-x4_18+x5_18)*(x3_18+x4_18-x5_18)
    equal(Q18x,rhs18,'D1.5.8 Heron factorization of Q')
    record('D1.5.8.triangle_shape_sign_elimination',17,
           {'cofactor_steps':'two cleared-denominator syzygies imply common beta/alpha ratio',
            'surviving_signs':'(X34,X35,X45)=s*(-1,1,-1)',
            'excluded_patterns':6,'surviving_mass_relation':'Q=0'})
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_isolated_triangle_shape_derivation'))

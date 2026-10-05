"""Finite state action in Diagram 6; labels' cube-root scales are separate."""
from itertools import product

SHAPE_SWAP=(0,1,5,4,3,2)
SHAPE_CONJUGATE=(0,1,4,5,2,3)


def generators(state):
    z,w,a,b,c,d=state
    return (
        (SHAPE_SWAP[z],w,-c,d,a,-b),
        (z,w,-a,-b,-c,-d),
        (SHAPE_CONJUGATE[z],SHAPE_CONJUGATE[w],-a,b,c,-d),
        (w,z,a,c,b,d),
    )


def census():
    states=set(product(range(6),range(6),(1,-1),(1,-1),(1,-1),(1,-1)))
    remaining=set(states)
    orbits=[]
    while remaining:
        seed=min(remaining);orbit={seed};pending=[seed]
        while pending:
            for image in generators(pending.pop()):
                if image not in states:raise AssertionError('D6 state action leaves state space')
                if image not in orbit:orbit.add(image);pending.append(image)
        remaining-=orbit
        kinds={b==-c and a==-d for z,w,a,b,c,d in orbit}
        if len(kinds)!=1:raise AssertionError('D6 corank sign condition not invariant')
        orbits.append(dict(seed=list(seed),states=[list(s) for s in sorted(orbit)],
            corank_one=kinds.pop(),size=len(orbit),generic_labels=9*len(orbit)))
    if (len(states),len(orbits),sum(o['corank_one'] for o in orbits))!=(576,29,8):
        raise AssertionError('D6 orbit census mismatch')
    return orbits

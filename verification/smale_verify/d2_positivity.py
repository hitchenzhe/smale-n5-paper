"""Exact positive-mass exclusions for D2 Lemma 6.12."""
from pathlib import Path
import hashlib,json
from .bernstein import coefficients,corner_sectors,blowup_critical_corner,positive_cover


def verify_positive_polynomial(terms):
 records=[]
 for sector in corner_sectors(terms):
  power=sector.pop('power_terms');matrix=coefficients(power)
  if min(v for row in matrix for v in row)>=0:
   cover=positive_cover(matrix,1)
   if cover['unresolved'] or cover['sign']!=1:raise ArithmeticError('Direct positive-mass sector not certified')
   records.append(dict(sector=sector,mode='direct',cover=cover));continue
  for mode in (0,1):
   transformed,low=blowup_critical_corner(power,mode);matrix=coefficients(transformed)
   cover=positive_cover(matrix,10)
   if cover['unresolved'] or cover['sign']!=1:raise ArithmeticError('Critical-corner positive-mass sector not certified')
   # mode 1, k=0 is a genuine positive-mass boundary (r=1,s=h).
   if mode==1 and min(row[0] for row in matrix)<=0:
    raise ArithmeticError('Uncertified k=0 positive-mass face')
   records.append(dict(sector=sector,mode=mode,removed_h_power=low,cover=cover,
    critical_substitution='r=1-h,s=h*k' if mode==0 else 'r=1-h*k,s=h',
    first_column_minimum=str(min(row[0] for row in matrix))))
 return records


def save(name,payload):
 path=Path(__file__).resolve().parents[1]/'certificates'/name
 path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
 return path.name

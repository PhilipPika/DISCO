# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/simpson.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
# This code is taken from: http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method#Python  on 15 november 2011
def simpsons_rule(f,a,b):
    c = (a+b) / 2.0
    h3 = abs(b-a) / 6.0
    return h3*(f(a) + 4.0*f(c) + f(b))
 
def recursive_asr(f,a,b,eps,whole):
    "Recursive implementation of adaptive Simpson's rule."
    c = (a+b) / 2.0
    left = simpsons_rule(f,a,c)
    right = simpsons_rule(f,c,b)
    if abs(left + right - whole) <= 15*eps:
        return left + right + (left + right - whole)/15.0
    return recursive_asr(f,a,c,eps/2.0,left) + recursive_asr(f,c,b,eps/2.0,right)
 
def adaptive_simpsons_rule(f,a,b,eps):
    "Calculate integral of f from a to b with max error of eps."
    return recursive_asr(f,a,b,eps,simpsons_rule(f,a,b))

if __name__ == "__main__":
    from math import sin,cos,log,exp
    print("Estim solution: ",adaptive_simpsons_rule(sin,0,1,.000000001))
    print("Exact solution: ",-cos(1) + cos(0)) 
    print("Estim solution: ",adaptive_simpsons_rule(sin,0,1,.001)) 
    print("Exact solution: ",-cos(1) + cos(0)) 


from sympy import *
a12, a13, a14, a23, a24, a34  = symbols("a12 a13 a14 a23 a24 a34")
A = Matrix([[0,a12,a13,a14], [-a12,0,a23,a24],[-a13,-a23,0,a34],[-a14,-a24,-a34,0]])
#B = Matrix([[0,1,0,0], [-1,0,0,0],[0,0,0,0],[0,0,0,0]])
B = A.diff(a14)
#print(B)
x = A.inv() 
y = x*B
z = y.trace()
print(simplify(z))

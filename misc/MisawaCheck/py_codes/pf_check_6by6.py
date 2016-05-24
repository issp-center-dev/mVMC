from sympy import *
a12, a13, a14, a15, a16 = symbols("a12 a13 a14 a15 a16")
a23, a24, a25, a26  = symbols("a23 a24 a25 a26")
a34, a35, a36  = symbols("a34 a35 a36")
a45, a46  = symbols("a45 a46")
a56 = symbols("a56")
A = Matrix([[0,a12,a13,a14,a15,a16], [-a12,0,a23,a24,a25,a26],[-a13,-a23,0,a34,a35,a36],[-a14,-a24,-a34,0,a45,a46],[-a15,-a25,-a35,-a45,0,a56],[-a16,-a26,-a36,-a46,-a56,0]])
B = A.diff(a14)
print(B)
x = A.inv()
y = x*B
z = y.trace()
print(simplify(z))

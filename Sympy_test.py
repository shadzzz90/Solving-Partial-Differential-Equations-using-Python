from sympy import *
import numpy as np
import matplotlib.pyplot as plt

x= symbols('x')
f = exp(-5*x**2)*cos(pi*x)

diff_f = diff(f,x)
diff_f_x = diff_f.evalf(subs={x:0.2})

# Backward Difference for O(h)
h = np.logspace(-5,1, 50)
x_prev = 0.2-h
x_prev_prev = 0.2 - (2*h)
x_next_next = 0.2 + (2*h)
x_next = 0.2+h
x_curr = np.full(len(h),0.2)
f_prev = np.zeros(len(h))
f_prev_prev = np.zeros(len(h))
f_curr = np.zeros(len(h))
f_next = np.zeros(len(h))
f_next_next = np.zeros(len(h))

for i in range(0,len(h)):
     f_prev[i] = f.evalf(subs ={x:x_prev[i]})
     f_curr[i] = f.evalf(subs = {x:x_curr[i]})
     f_next[i] = f.evalf(subs= {x:x_next[i]})
     f_prev_prev[i] = f.evalf(subs ={x:x_prev_prev[i]})
     f_next_next[i] = f.evalf(subs={x: x_next_next[i]})

diff_f_approx_first_accu = (f_curr-f_prev)/h
diff_f_approx_second_accu = (f_next-f_prev)/(2*h)
diff_f_approx_fourth_accu = (-f_next_next + 8*f_next-8*f_prev+f_prev_prev)/(12*h)


e_first_order_accu = np.abs(diff_f_x - diff_f_approx_first_accu)
e_second_order_accu = np.abs(diff_f_x - diff_f_approx_second_accu)
e_fourth_order_accu = np.abs(diff_f_x - diff_f_approx_fourth_accu)


# np.reshape(h,(10,1))
# np.reshape(e,(10,1))
plt.loglog(h, e_first_order_accu,'r', h, e_second_order_accu, 'b', h, e_fourth_order_accu,'g')
plt.title('Plot of $e(x)$ vs $\Delta x$ ($f^\prime (x)$ at $x=0.2$)')
plt.grid()
plt.xlim(h[0], h[len(h)-1])
plt.xlabel('$\Delta x$' )
plt.ylabel('$e(x)$' )
plt.gca().legend(('$ O (\Delta x)$', '$ O (\Delta x^2) $', '$ O (\Delta x^4) $'))

plt.show()
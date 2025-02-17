from temp import temp 

y0,t0,tf,npoints=1,0,1,10
coeff=1.0/10.0
equations = "double_system_of_equations"
params = [coeff]
y0 = [1,2]
t,yt=temp.rk4(equations, [t0,tf], y0, params, npoints)

print(yt)


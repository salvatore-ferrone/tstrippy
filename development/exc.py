from temp import temp 

y0,t0,tf,npoints=10,0,1,100
coeff=1.0/100.0
equations = "exponential_decay"
params = [coeff]
t,yt=temp.rk4(equations, y0,t0,tf,params, npoints)


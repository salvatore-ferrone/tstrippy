import matplotlib.pyplot as plt
from temp import temp 

npoints=1000
r,w,dwdr=temp.solve_king_potential_profile(10.,npoints)
r_break = r[npoints//2]
fig,axis=plt.subplots()
axis.plot(r,w,'.')
axis.stem([r_break],[w[npoints//2]],'r')
# axis.plot(r,dwdr,'.')
axis.set_xlabel('r')
axis.set_title('King Potential Profile')
plt.show()
fig.savefig('king_potential_profile.png')
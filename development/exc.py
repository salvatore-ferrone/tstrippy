import matplotlib.pyplot as plt
from temp import temp 

npoints=1001
r,w,dwdr=temp.solve_king_potential_profile(10.,npoints)

fig,axis=plt.subplots()
axis.plot(r,w,'.')
axis.plot(r,dwdr,'.')
axis.set_xlabel('r')
axis.legend(['w','dwdr'])
axis.set_title('King Potential Profile')
plt.show()
fig.savefig('king_potential_profile.png')
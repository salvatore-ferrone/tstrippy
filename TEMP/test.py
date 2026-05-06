import gravitymini


g = gravitymini.gravity


g.addgravitycomponent('exponentialoblatehalo', [1,1,2])
g.finalizegravity()
print("ok")
ax,ay,az=g.evaluategravityforces([1],[3],[1])
print(ax,ay,ax)
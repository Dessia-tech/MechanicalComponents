import numpy as npy
ligne1='56.097561,477.97196 81.707319,-1.21951 81.70732,1.21951 82.92682,0 79.2683,-1.21951 82.92683,0 81.70731,-1.21951 81.70732,0 81.70732,1.21951'
ligne2='57.317073,476.75245 -1.219512,-82.92683 0,-82.92683 1.219512,-84.14634 0,-82.92683 1.219512,-82.92683'
ligne3='56.097561,58.459766 76.829269,59.756094 71.95122,43.90244 89.02439,45.12195 78.04878,32.92683 75.60976,19.5122 74.39024,7.31707 90.2439,1.21951 97.56098,1.21952'

axe_x=[0,5,10,15,20,25,30,35,40]
axe_y=[0.5,0.6,0.7,0.8,0.9,1]

def analyze(ligne):
    data=[]
    sol=ligne.split(' ')
    for item in sol:
        temp=item.split(',')
        data.append([float(temp[0]),float(temp[1])])
        
    export=[data[0]]
    for item in data[1::]:
        export.append([item[0]+export[-1][0],item[1]+export[-1][1]])
    return export

export1=npy.array(analyze(ligne1))
export2=npy.array(analyze(ligne2))
export3=npy.array(analyze(ligne3))

vect_x=npy.linspace(export1[0,0],export1[-1,0],len(axe_x))
vect_y=npy.linspace(export2[0,1],export2[-1,1],len(axe_y))

ax=(axe_x[0]-axe_x[-1])/(vect_x[0]-vect_x[-1])
bx=axe_x[0]-ax*vect_x[0]
ay=(axe_y[0]-axe_y[-1])/(vect_y[0]-vect_y[-1])
by=axe_y[0]-ay*vect_y[0]    

export=[]
for item in export3:
    export.append([item[0]*ax+bx,item[1]*ay+by])

fichier=open('coeff_yb_iso.txt','w')
fichier.write('helix_angle,coeff_yb_iso\n')
for item in export:
    fichier.write('{},{}\n'.format(item[0]/180*npy.pi,item[1]))
fichier.close()


import volmdlr as vm

fichier=open("coeff_yb_iso.txt",'r')
lines=fichier.readlines()
fichier.close()
tab_coeff_yb_iso={}
list_col=[]
for i in lines[0].split('\n')[0].split(','):
    tab_coeff_yb_iso[i]=[]
    list_col.append(i)
list_coeff_yb_iso=[]
for i in lines[1::]:
    temp=i.split('\n')[0].split(',')
    tab_coeff_yb_iso[list_col[0]].append(float(temp[0]))
    tab_coeff_yb_iso[list_col[1]].append(float(temp[1]))
    list_coeff_yb_iso.append(vm.Point2D((float(temp[0]),float(temp[1]))))
            
    

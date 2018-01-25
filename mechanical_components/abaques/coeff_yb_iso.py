#import volmdlr as vm
#
#fichier=open("coeff_yb_iso.txt",'r')
#lines=fichier.readlines()
#fichier.close()
#tab_coeff_yb_iso={}
#list_col=[]
#for i in lines[0].split('\n')[0].split(','):
#    tab_coeff_yb_iso[i]=[]
#    list_col.append(i)
#list_coeff_yb_iso=[]
#for i in lines[1::]:
#    temp=i.split('\n')[0].split(',')
#    tab_coeff_yb_iso[list_col[0]].append(float(temp[0]))
#    tab_coeff_yb_iso[list_col[1]].append(float(temp[1]))
#    list_coeff_yb_iso.append(vm.Point2D((float(temp[0]),float(temp[1]))))
            
    
tab_coeff_yb_iso={'coeff_yb_iso': [1.0029325508401201,
  0.9310850480431024,
  0.8782991233021732,
  0.8240469255458759,
  0.784457481990179,
  0.7609970656504502,
  0.7521994155347822,
  0.7507331425194141,
  0.7492668574805859],
 'helix_angle': [0.0,
  0.08205652498638352,
  0.15890311169191082,
  0.25398448077837904,
  0.33734348899993005,
  0.4180975328872207,
  0.49754908392700825,
  0.5939329351806066,
  0.6981317007977318]}
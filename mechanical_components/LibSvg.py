import numpy as npy

class SVGTrace:
    def __init__(self,scale=1,swift=0,angle=0,decal_x=0,decal_y=0):
        
        self.scale=scale
        self.swift=swift
        self.angle=angle
        self.decal_x=decal_x
        self.decal_y=decal_y
        self.data=[]
        self.box=[]
        self.zpd=[]
        self.begin=[]
        self.end=[]
        self.list_name={}
        self.indice=1
        self.InitBorne()
        
    def InitBorne(self):
        
        self.boxX_min=[npy.inf]
        self.boxX_max=[-npy.inf]
        self.boxY_min=[npy.inf]
        self.boxY_max=[-npy.inf]
        
    def UpdateBox(self,x1,y1,r1=None):
        
        X1=x1
        Y1=y1
        if not r1==None:
            X1=[x1[0]-r1,x1[0],x1[0]+r1,x1[0]]
            Y1=[y1[0],y1[0]+r1,y1[0],y1[0]-r1]
        self.boxX_min=min(self.boxX_min,[min(X1)])
        self.boxX_max=max(self.boxX_max,[max(X1)])
        self.boxY_min=min(self.boxY_min,[min(Y1)])
        self.boxY_max=max(self.boxY_max,[max(Y1)])
        
    def CreateBegin(self):
        
        scale_html=100-100/(1500/(1000/(self.boxX_max[0]-self.boxX_min[0])*(self.boxY_max[0]-self.boxY_min[0])))
        self.begin.append('''<html>
<head>
<meta charset="UTF-8"> 
<script src="https://cdn.dessia.tech/snap.svg/v0.5.1/snap.svg-min.js"></script>
    <script src="https://cdn.dessia.tech/snap.svg.zpd/v0.0.11/snap.svg.zpd.js"></script>

</script>
</head>

<body>
<svg id="Gear" width="100%" height="100%">

</svg>

<script>''')
        
    def CreateEnd(self,animate):
        
        if not animate==None:
            master=list(animate.keys())[0]
            i=0
            Temp=[]
            for (j,k) in animate.items():
                self.end.append('var myMatrix'+str(int(i))+' = new Snap.Matrix();')
                self.end.append('myMatrix'+str(int(i))+'.add(1,0,0,1,'+str(k['R'][1])+', '+str(k['R'][2])+');')
                self.end.append(j+'.transform(myMatrix'+str(int(i))+');')
                Temp.append(j+'.attr({ transform: myMatrix'+str(int(i))+'});')
                i+=1
                self.end.append('var myMatrix'+str(int(i))+' = new Snap.Matrix();')
                self.end.append('myMatrix'+str(int(i))+'.add('+str(npy.cos(k['R'][0]))+','+str(npy.sin(k['R'][0]))+','+str(-npy.sin(k['R'][0]))+','+str(npy.cos(k['R'][0]))+','+str(k['R'][1])+', '+str(k['R'][2])+');')
                i+=1
                
            Temp2=[]
            i=0
            for (j,k) in animate.items():
                if not j==master:
                    Temp2.append(str(j)+'''.animate(
			{ transform: myMatrix'''+str(int(i+1))+'''}, 
			1000)''')
                if j==master:
                    indice=i
                i+=2
                
            self.end.append('''function anim(){
		'''+str(master)+'''.animate(
			{ transform: myMatrix'''+str(int(indice+1))+'''}, 
			1000,
			function(){ 
				'''+Temp[0]+'''
                  '''+Temp[1]+'''
				anim();
			});
             '''+Temp2[0]+'''
	}''')
    
            self.end.append('s.click( anim);') 
        self.end.append('</script></body></html>')
        
    def CreateBox(self):
        
        self.view_x=(self.boxX_max[0]-self.boxX_min[0])
        self.view_y=(self.boxY_max[0]-self.boxY_min[0])
        self.box.append('var s = Snap("#Gear");')
#        self.box.append('var s = Snap('+str(self.view_x)+','+str(self.view_y)+');')
#        self.box.append('s.attr({ viewBox: "'+str(self.boxX_min[0]+self.decal_x)+' '+str(self.boxY_min[0]+self.decal_y)+' '+str(self.boxX_max[0]-self.boxX_min[0]+selfe.decal_x)+' '+str(self.boxY_max[0]-self.boxY_min[0]+self.decal_y)+'" });')

    def Transform(self,data,init=None):
        
        alpha=self.angle
        MatRot=npy.array([[npy.cos(alpha),-npy.sin(alpha),0],[npy.sin(alpha),npy.cos(alpha),0],[0,0,1]])
        if self.swift==1:
            MatPermut=npy.array([[0,1,0],[1,0,0],[0,0,1]])
            Mat=npy.matmul(MatRot,MatPermut)
        else:
            Mat=MatRot
        dataS=npy.matmul(Mat,data)
        
        if not init==None:
            self.InitBorne()
        self.UpdateBox(list(dataS[0,:]),list(dataS[1,:]))
        
        scaleX=self.scale/(self.boxX_max[0]-self.boxX_min[0])
        scaleY=self.scale/(self.boxY_max[0]-self.boxY_min[0])
        scaleZPD=min(scaleX,scaleY)
        MatScale=npy.array([[scaleZPD,0,-scaleZPD*self.boxX_min[0]],[0,scaleZPD,-scaleZPD*self.boxY_min[0]],[0,0,1]])
        
#        Mat=npy.matmul(Mat1,npy.array([[1,0,-self.boxX_min[0]],[0,1,-5],[0,0,1]]))
        Mat=npy.matmul(MatScale,Mat)
        return Mat
    
    def CreateZPD(self):
        
        data=npy.array([[self.boxX_min[0],self.boxX_max[0],self.boxX_max[0],self.boxX_min[0]],[self.boxY_max[0],self.boxY_max[0],self.boxY_min[0],self.boxY_min[0]],[1,1,1,1]])
        Mat=self.Transform(data,'ok')
        
#        self.zpd.append('s.zpd({ load: {a:'+str(scaleZPD)+',b:0,c:0,d:'+str(scaleZPD)+',e:'+str(-scaleZPD*self.boxX_min[0])+',f:'+str(-scaleZPD*self.boxY_min[0])+'}});')
        
        for n,m in self.list_name.items():
            Temp='var '+n+' = s.group('
            for i in m:
                Temp+=i+','
            self.zpd.append(Temp[0:-1]+');')
            
        self.zpd.append('s.zpd({ load: {a:'+str(Mat[0,0])+',b:'+str(Mat[1,0])+',c:'+str(Mat[0,1])+',d:'+str(Mat[1,1])+',e:'+str(Mat[0,2])+',f:'+str(Mat[1,2])+'}});')
        self.zpd.append('s.zpd({ drag: false}, function (err, paper) {console.log(paper)});')
                
    def ConvertPrimitive2D(self,L,group_name,stroke,strokeWidth,impact_box,strokeDasharray="none"):
        
        for i,j in enumerate(L):
            if 'primitives2D' in str(j.__class__):
                for n,m in enumerate(j.primitives):
                    if 'Line2D' in str(m.__class__):
                        temp='var primitive2D'+str(self.indice)+' = s.polyline(['
                        self.list_name[group_name].append('primitive2D'+str(self.indice))
                        self.indice+=1
                        for k in m.points[0:-1]:
                            if impact_box==0:
                                self.UpdateBox([k.vector[0]],[k.vector[1]])
                            temp+=str(k.vector[0])+','+str(k.vector[1])+','
                        temp+=str(m.points[-1].vector[0])+','+str(m.points[-1].vector[1])
                        temp+=']).attr({stroke: \''+stroke+'\',strokeWidth: '+str(strokeWidth)+',fill:"none" , strokeDasharray: "'+strokeDasharray+'"}).addClass("the-class");'
                        self.data.append(temp)
                    if 'Arc2D' in str(m.__class__):
                        temp='var path'+str(self.indice)+' = s.path('
                        self.list_name[group_name].append('path'+str(self.indice))
                        self.indice+=1
                        temp+='\'M '+str(m.start.vector[0])+','+str(m.start.vector[1])
                        temp+=' A '+str(m.radius)+','+str(m.radius)+' '+str((m.angle1+m.angle2)/2/npy.pi*180+90)+' 0 '+' 1 '+str(m.end.vector[0])+' '+str(m.end.vector[1])
                        temp+='\').attr({stroke: \''+stroke+'\',strokeWidth: '+str(strokeWidth)+',fill:"none" , strokeDasharray: "'+strokeDasharray+'"}).addClass("the-class");'
                        self.data.append(temp)
                        
            
    def ConvertCircle2D(self,L,group_name,stroke,strokeWidth,impact_box,strokeDasharray="none"):
        
        for i,j in enumerate(L):
            if 'Circle2D' in str(j.__class__):
                if impact_box==0:
                    self.UpdateBox([j.center.vector[0]],[j.center.vector[1]],j.radius)
                temp='var circle2D'+str(self.indice)+' = s.circle('
                self.list_name[group_name].append('circle2D'+str(self.indice))
                self.indice+=1
                temp+=str(j.center.vector[0])+','+str(j.center.vector[1])+','
                temp+=str(j.radius)
                temp+=').attr({stroke: \''+stroke+'\',strokeWidth: '+str(strokeWidth)+',fill:"none" , strokeDasharray: "'+strokeDasharray+'"}).addClass("the-class");'
                self.data.append(temp)
                
    def ConvertLine2D(self,L,group_name,stroke,strokeWidth,impact_box,strokeDasharray="none"):
        
        for i,j in enumerate(L):
            if 'Line2D' in str(j.__class__):
                temp='var primitive2D'+str(self.indice)+' = s.polyline(['
                self.list_name[group_name].append('primitive2D'+str(self.indice))
                self.indice+=1
                if impact_box==0:
                    self.UpdateBox([j.center.vector[0]],[j.center.vector[1]])
                for k in j.points[0:-1]:
                    if impact_box==0:
                        self.UpdateBox([k.vector[0]],[k.vector[1]])
                    temp+=str(k.vector[0])+','+str(k.vector[1])+','
                temp+=str(j.points[-1].vector[0])+','+str(j.points[-1].vector[1])
                temp+=']).attr({stroke: \''+stroke+'\',strokeWidth: '+str(strokeWidth)+',fill:"none" , strokeDasharray: "'+strokeDasharray+'"}).addClass("the-class");'
                self.data.append(temp)
    
    def Convert(self,L,group_name,stroke,strokeWidth,impact_box=0,strokeDasharray="none"):
        
        self.list_name[group_name]=[]
        self.ConvertPrimitive2D(L,group_name,stroke,strokeWidth,impact_box,strokeDasharray)
        self.ConvertCircle2D(L,group_name,stroke,strokeWidth,impact_box,strokeDasharray)
        self.ConvertLine2D(L,group_name,stroke,strokeWidth,impact_box,strokeDasharray)
    
    def Export(self,name='export.html',animate=None):
        fichier=open(name,'w')
        self.CreateBegin()
        for i in self.begin:
            fichier.write(i+'\n')
        self.CreateZPD()
        self.CreateBox()
        for i in self.box:
            fichier.write(i+'\n')
        for i in self.data:
            fichier.write(i+'\n')
        for i in self.zpd:
            fichier.write(i+'\n')
        self.CreateEnd(animate)
        for i in self.end:
            fichier.write(i+'\n')
        fichier.close()
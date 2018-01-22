import numpy as npy
from jinja2 import Environment, PackageLoader, select_autoescape

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
        self.primitives=[]
        self.name=[]
        
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
        
    def CreateBox(self):
        
        self.view_x=(self.boxX_max[0]-self.boxX_min[0])
        self.view_y=(self.boxY_max[0]-self.boxY_min[0])
        self.box.append('var s = Snap("#Gear");')
                                      
        self.pt1=self.boxX_min[0]+self.decal_x
        self.pt2=self.boxY_min[0]+self.decal_y
        self.pt3=self.boxX_max[0]-self.boxX_min[0]+self.decal_x
        self.pt4=self.boxY_max[0]-self.boxY_min[0]+self.decal_y 
                
    def ConvertPrimitive2D(self,L,group_name,stroke,strokeWidth,
                           impact_box,strokeDasharray="none"):
        
        for i,j in enumerate(L):
            if 'primitives2D' in str(j.__class__):
                for n,m in enumerate(j.primitives):
                    if 'Line2D' in str(m.__class__):
                        temp='var {} = svg.selectAll("#{}").append("path").attr("d","M{},{}L{},{}").attr("stroke-width", {}).attr("stroke", "{}").attr("stroke-dasharray","{}").attr("fill", "none");'.format(group_name,group_name,m.points[0].vector[0],m.points[0].vector[1],m.points[-1].vector[0],m.points[-1].vector[1],strokeWidth,stroke,strokeDasharray)
                        self.primitives.append(temp)
                        if impact_box==0:
                            self.UpdateBox([m.points[0].vector[0]],[m.points[0].vector[1]])
                            
                    if 'Arc2D' in str(m.__class__):
                        temp='var {} = svg.selectAll("#{}").append("path").attr("d", "M{},{}A{},{} {} 0,1 {},{}").attr("stroke-width", {}).attr("stroke", "{}").attr("stroke-dasharray","{}").attr("fill", "none");'.format(group_name,group_name,m.start.vector[0],m.start.vector[1],m.radius,m.radius,(m.angle1+m.angle2)/2/npy.pi*180+90,m.end.vector[0],m.end.vector[1],strokeWidth,stroke,strokeDasharray)
                        self.primitives.append(temp)                        
            
    def ConvertCircle2D(self,L,group_name,stroke,strokeWidth,impact_box,strokeDasharray="none"):
        
        for i,j in enumerate(L):
            if 'Circle2D' in str(j.__class__):
                if impact_box==0:
                    self.UpdateBox([j.center.vector[0]],[j.center.vector[1]],j.radius)
                temp='var {} = svg.selectAll("#{}").append("circle").attr("stroke-width", {}).attr("stroke", "{}").attr("stroke-dasharray","{}").attr("fill", "none").attr("cx",{}).attr("cy",{}).attr("r",{});'.format(group_name,group_name,strokeWidth,stroke,strokeDasharray,j.center.vector[0],j.center.vector[1],j.radius)
                self.primitives.append(temp)
                
    def ConvertLine2D(self,L,group_name,stroke,strokeWidth,impact_box,strokeDasharray="none"):
        
        for i,j in enumerate(L):
            if 'Line2D' in str(j.__class__):
                temp='var {} = svg.selectAll("#{}").append("line").attr("stroke-width", {}).attr("stroke", "{}").attr("stroke-dasharray","{}").attr("fill", "none").attr("x1",{}).attr("y1",{}).attr("x2",{}).attr("y2",{});'.format(group_name,group_name,strokeWidth,stroke,strokeDasharray,j.points[0].vector[0],j.points[0].vector[1],j.points[-1].vector[0],j.points[-1].vector[1])
                self.primitives.append(temp)
                if impact_box==0:
                    self.UpdateBox([j.center.vector[0]],[j.center.vector[1]])
                for k in j.points[0:-1]:
                    if impact_box==0:
                        self.UpdateBox([k.vector[0]],[k.vector[1]])
    
    def Convert(self,L,group_name,stroke,strokeWidth,impact_box=0,strokeDasharray="none"):
        
        self.name.append('{} "group":"{}"{},'.format(chr(123),group_name,chr(125)))
        
        self.ConvertPrimitive2D(L,group_name,stroke,strokeWidth,impact_box,strokeDasharray)
        self.ConvertCircle2D(L,group_name,stroke,strokeWidth,impact_box,strokeDasharray)
        self.ConvertLine2D(L,group_name,stroke,strokeWidth,impact_box,strokeDasharray)
    
    def ExportAnimate(self,name='export.html',animate=None):
        
        env = Environment(loader=PackageLoader('mechanical_components', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))
        
        template = env.get_template('template_animate.html')
        
        #box animate
        self.view_x=(self.boxX_max[0]-self.boxX_min[0])
        self.view_y=(self.boxY_max[0]-self.boxY_min[0])
        width=700
        scale=width/self.view_x
        height=scale*self.view_y
                                      
        self.vb1=self.boxX_min[0]
        self.vb2=self.boxY_min[0]
        self.vb3=self.boxX_max[0]-self.boxX_min[0]
        self.vb4=self.boxY_max[0]-self.boxY_min[0]
        
        return template.render(list_polyline=self.primitives,rot1=(animate['gear1']['R'][0]/npy.pi*180),pos1_x=animate['gear1']['R'][1],pos1_y=animate['gear1']['R'][2],rot2=(animate['gear2']['R'][0]/npy.pi*180),pos2_x=animate['gear2']['R'][1],pos2_y=animate['gear2']['R'][2],width=width,height=height,vb1=self.vb1,vb2=self.vb2,vb3=self.vb3,vb4=self.vb4)
    
    def Export(self,name='export.html'):
        
        env = Environment(loader=PackageLoader('mechanical_components', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))
        
        template = env.get_template('template.html')
        
        #box animate
        self.view_x=(self.boxX_max[0]-self.boxX_min[0])
        self.view_y=(self.boxY_max[0]-self.boxY_min[0])
        width=700
        scale=width/self.view_x
        height=scale*self.view_y
                                      
        self.vb1=self.boxX_min[0]
        self.vb2=self.boxY_min[0]
        self.vb3=self.boxX_max[0]-self.boxX_min[0]
        self.vb4=self.boxY_max[0]-self.boxY_min[0]
        
        return template.render(list_name=self.name,list_polyline=self.primitives,width=width,height=height,vb1=self.vb1,vb2=self.vb2,vb3=self.vb3,vb4=self.vb4)
    
    def Show(self,name='export.html',animate=None):
#        print(self.BabylonScript())
        if not animate==None:
            with open(name,'w') as file:
                file.write(self.ExportAnimate(name,animate))
        else:
            with open(name,'w') as file:
                file.write(self.Export(name))


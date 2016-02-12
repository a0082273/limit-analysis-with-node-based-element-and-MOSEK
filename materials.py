material_type=None
param=[]
cohesion=None
inner_friction_angle=None
alpha=None
k=None

class Material():
    # def __init__(self):

    def get_material_type(self,infile):
        global material_type,param

        infile.f.seek(0)
        for iline in range(infile.nline):
            line=infile.f.readline()
            if line.find('material')!=-1: break
        material_type=infile.f.readline().strip()
        print('Material type :',material_type)

#get_param() wo tukuru
        param=[]
        if material_type=='druckerprager' or material_type=='mohrcoulomb':
            para=infile.f.readline().strip()
            param.append(float(para.split()[0]))
            param.append(float(para.split()[1]))
            print('Material property :',param[0],param[1])
        else:
            print('material type must be druckerprager or mohrcoulomb')


    def get_cohesion(self):
        global cohesion
        cohesion=param[0]#/3**0.5   #for drucker->von


    def get_inner_friction_angle(self):
        global inner_friction_angle
        from math import pi
        inner_friction_angle=(param[1]/180.)*pi


    def get_alpha(self):   #yf=J_2**0.5+alpha*tr(stress)-k
        from math import sin
        global alpha
        phi=inner_friction_angle
        alpha=sin(phi)/(9.+3*sin(phi)**2)**0.5

    def get_k(self):   #yf=J_2**0.5+alpha*tr(stress)-k
        from math import sin,cos
        global k
        phi=inner_friction_angle
        c=cohesion
        k=3**0.5*c*cos(phi)/(3.+sin(phi)**2)**0.5


    def get_yield_function(self,stress):
        if material_type=='druckerprager':   #yf=J_2**0.5+alpha*tr(stress)-k
            from mesh import Point
            poi=Point()
            m=poi.get_mean_stress(stress)
            J2=poi.get_second_invariant_of_deviatoric_stress(stress)
            yf=J2**0.5+3*alpha*m-k
        return yf


    def get_derivative_of_yield_function(self,stress):
        import numpy as np
        if material_type=='druckerprager':
            from mesh import Point
            poi=Point()
            J2=poi.get_second_invariant_of_deviatoric_stress(stress)
            J2+=1.e-10
            dJ2=np.empty(3,dtype=float)
            for i in range(3):
                j=(i+1)%3
                k=(i+2)%3
                dJ2[i]=(2*stress[i]-stress[j]-stress[k])/3
            dyf=np.empty(6,dtype=float)
            dyf[0]=0.5*dJ2[0]*J2**(-0.5)+alpha
            dyf[1]=0.5*dJ2[1]*J2**(-0.5)+alpha
            dyf[2]=0.5*dJ2[2]*J2**(-0.5)+alpha
            dyf[3]=stress[3]*J2**(-0.5)
            dyf[4]=stress[4]*J2**(-0.5)
            dyf[5]=stress[5]*J2**(-0.5)
        return dyf


    def get_material_parameters(self):
        self.get_cohesion()
        self.get_inner_friction_angle()
        self.get_alpha()
        self.get_k()

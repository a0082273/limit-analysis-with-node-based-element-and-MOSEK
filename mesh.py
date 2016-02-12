import numpy as np

elem_type=None
nelem=None
elem=[]
npoint=None
point=[]

class Element():
    def __init__(self):
        self.node=np.empty(4,dtype=int)
        self.volume=None



    def get_element_type(self,infile):
        infile.f.seek(0)
        for iline in range(infile.nline):
            line=infile.f.readline()
            if line.find('element_type')!=-1: break
        elem_type=infile.f.readline().strip()
        print('Element type :',elem_type)



    def get_number_of_element(self,infile):
        global nelem
        infile.f.seek(0)
        for iline in range(infile.nline):
            line=infile.f.readline()
            if line.find('elements')!=-1: break
        nelem=int(line.strip().split()[1])
        print('Number of elements :',nelem)



    def get_node_of_element(self,infile):
        infile.f.seek(0)
        for iline in range(infile.nline):
            line=infile.f.readline()
            if line.find('elements')!=-1: break
        for iline in range(nelem):
            line=infile.f.readline()
            for i in range(4):
                elem.append(Element())
                elem[iline].node[i]=int(line.strip().split()[i+1])-1
            # for i in range(3):
            #     for j in range(i+1,4):
            #         if(elem[iline].node[i]>elem[iline].node[j]):
            #             elem[iline].node[i],elem[iline].node[j]=elem[iline].node[j],elem[iline].node[i]




    def get_derivative_of_shape_function_T4(self,ielem):   #volume is multiplied
        import sys
        dNdx=np.empty(4,dtype=float)
        dNdy=np.empty(4,dtype=float)
        dNdz=np.empty(4,dtype=float)
        xe=np.empty(4,dtype=float)
        ye=np.empty(4,dtype=float)
        ze=np.empty(4,dtype=float)

        def get_local_node_coordinate(ielem,elem,point):
            for ilocal in range(4):
                iglobal=elem[ielem].node[ilocal]
                xe[ilocal]=point[iglobal].x[0]
                ye[ilocal]=point[iglobal].x[1]
                ze[ilocal]=point[iglobal].x[2]
            return xe,ye,ze

        xe,ye,ze=get_local_node_coordinate(ielem,elem,point)
        dNdx[0]=((ye[2]-ye[3])*(ze[1]-ze[3])-(ye[1]-ye[3])*(ze[2]-ze[3]))/6.
        dNdx[1]=((ye[0]-ye[2])*(ze[0]-ze[3])-(ye[0]-ye[3])*(ze[0]-ze[2]))/6.
        dNdx[2]=((ye[0]-ye[3])*(ze[0]-ze[1])-(ye[0]-ye[1])*(ze[0]-ze[3]))/6.
        dNdx[3]=((ye[0]-ye[1])*(ze[0]-ze[2])-(ye[0]-ye[2])*(ze[0]-ze[1]))/6.
        dNdy[0]=((ze[2]-ze[3])*(xe[1]-xe[3])-(ze[1]-ze[3])*(xe[2]-xe[3]))/6.
        dNdy[1]=((ze[0]-ze[2])*(xe[0]-xe[3])-(ze[0]-ze[3])*(xe[0]-xe[2]))/6.
        dNdy[2]=((ze[0]-ze[3])*(xe[0]-xe[1])-(ze[0]-ze[1])*(xe[0]-xe[3]))/6.
        dNdy[3]=((ze[0]-ze[1])*(xe[0]-xe[2])-(ze[0]-ze[2])*(xe[0]-xe[1]))/6.
        dNdz[0]=((xe[2]-xe[3])*(ye[1]-ye[3])-(xe[1]-xe[3])*(ye[2]-ye[3]))/6.
        dNdz[1]=((xe[0]-xe[2])*(ye[0]-ye[3])-(xe[0]-xe[3])*(ye[0]-ye[2]))/6.
        dNdz[2]=((xe[0]-xe[3])*(ye[0]-ye[1])-(xe[0]-xe[1])*(ye[0]-ye[3]))/6.
        dNdz[3]=((xe[0]-xe[1])*(ye[0]-ye[2])-(xe[0]-xe[2])*(ye[0]-ye[1]))/6.
        elem[ielem].volume=xe[0]*dNdx[0]+xe[1]*dNdx[1]+xe[2]*dNdx[2]+xe[3]*dNdx[3]
        return dNdx,dNdy,dNdz






class Point():
    def __init__(self):
        self.x=np.empty(3,dtype=float)
        self.b0=np.zeros(3,dtype=float)
        self.b1=np.zeros(3,dtype=float)
        self.bc=['NEUMANN','NEUMANN','NEUMANN']
        self.nnode=1
        self.node=[]
        self.volume=0.
        # self.dNdx=[] #koreni touitsh sitai
        self.dNdx=[]
        self.dNdy=[]
        self.dNdz=[]
        self.stress=np.empty(6,dtype=float)
        self.u=np.empty(3,dtype=float)
        self.s=np.empty(6,dtype=float)



    def get_number_of_point(self,infile):
        global npoint
        infile.f.seek(0)
        for iline in range(infile.nline):
            line=infile.f.readline()
            if line.find('node_coordinates')!=-1: break
        npoint=int(line.strip().split()[1])
        print('Number of nodes :',npoint)



    def get_coordinate(self,infile):
        infile.f.seek(0)
        for iline in range(infile.nline):
            line=infile.f.readline()
            if line.find('node_coordinates')!=-1: break
        for iline in range(npoint):
            line=infile.f.readline()
            for i in range(3):
                point.append(Point())
                point[iline].x[i]=float(line.strip().split()[i+1])



    def get_prescribed_velocity(self,infile):
        infile.f.seek(0)
        for iline in range(infile.nline):
            line=infile.f.readline()
            if line.find('prescribed_velocities')!=-1: break
        nprevel=int(line.strip().split()[1])
        print('Number of prescribed velocities :',nprevel)
        for iline in range(nprevel):
            line=infile.f.readline()
            ipoint=int(line.strip().split()[0])-1
            idir=int(line.strip().split()[1])-1
            prevel=float(line.strip().split()[2])
            if prevel==0.:
                # point[ipoint].bc[idir]='HOMOGENEOUSDIRICHLET'
                point[ipoint].bc[idir]='HOMOGEN'
                point[ipoint].bc[idir]+='EOUSDIRICHLET'
            else:
                # point[ipoint].bc[idir]='NONHOMOGENEOUSDIRICHLET'
                point[ipoint].bc[idir]='NONHOMO'
                point[ipoint].bc[idir]+='GENEOUSDIRICHLET'



    def get_prescribed_force(self,infile):
        infile.f.seek(0)
        for iline in range(infile.nline):
            line=infile.f.readline()
            if line.find('prescribed_forces')!=-1: break
        npreforce=int(line.strip().split()[1])
        print('Number of prescribed forces :',npreforce)
        for iline in range(npreforce):
            line=infile.f.readline()
            ipoint=int(line.strip().split()[0])-1
            idir=int(line.strip().split()[1])-1
            preforce=float(line.strip().split()[2])
            point[ipoint].b0[idir]=preforce



    def get_applied_force(self,infile):
        infile.f.seek(0)
        for iline in range(infile.nline):
            line=infile.f.readline()
            if line.find('applied_forces')!=-1: break
        nforce=int(line.strip().split()[1])
        print('Number of applied forces :',nforce)
        for iline in range(nforce):
            line=infile.f.readline()
            ipoint=int(line.strip().split()[0])-1
            idir=int(line.strip().split()[1])-1
            force=float(line.strip().split()[2])
            point[ipoint].b1[idir]=force



    def get_mean_stress(self,stress):
        m=(stress[0]+stress[1]+stress[2])/3.
        return m



    def get_deviatoric_stress(self,stress):
        dev=np.empty(6,dtype=float)
        m=self.get_mean_stress(stress)
        for i in range(3):
            dev[i]=stress[i]-m
        for i in range(3,6):
            dev[i]=stress[i]
        return dev



    def get_second_invariant_of_deviatoric_stress(self,stress):
        dev=self.get_deviatoric_stress(stress)
        J2=(dev[0]**2+dev[1]**2+dev[2]**2)*0.5+dev[3]**2+dev[4]**2+dev[5]**2
        return J2






    def get_node_of_point(self):
        for ipoint in range(npoint):
            point[ipoint].node.append(ipoint)
        for ielem in range(nelem):
            for inode in range(4):
                jnode=(inode+1)%4
                knode=(inode+2)%4
                lnode=(inode+3)%4
                ipoint=elem[ielem].node[inode]
                jpoint=elem[ielem].node[jnode]
                kpoint=elem[ielem].node[knode]
                lpoint=elem[ielem].node[lnode]
                point[ipoint].nnode,point[ipoint].node=self.add_node(point[ipoint].nnode,point[ipoint].node,jpoint) #check
                point[ipoint].nnode,point[ipoint].node=self.add_node(point[ipoint].nnode,point[ipoint].node,kpoint)
                point[ipoint].nnode,point[ipoint].node=self.add_node(point[ipoint].nnode,point[ipoint].node,lpoint)




    def add_node(self,nnode,node,ipoint):
        for inode in range(nnode):
            if node[inode]==ipoint:
                return nnode,node
        nnode+=1
        node.append(ipoint)
        return nnode,node



    def get_derivative_of_shape_function_NS(self):   #volume is multiplied
        import sys

        self.get_node_of_point()
        from mesh import Element
        el=Element()
        for ipoint in range(npoint):
            for inode in range(point[ipoint].nnode):
                point[ipoint].dNdx.append(0.)
                point[ipoint].dNdy.append(0.)
                point[ipoint].dNdz.append(0.)

        for ielem in range(nelem):
            dNdxT4,dNdyT4,dNdzT4=el.get_derivative_of_shape_function_T4(ielem)
            for inode in range(4):
                jnode=(inode+1)%4
                knode=(inode+2)%4
                lnode=(inode+3)%4
                ipoint=elem[ielem].node[inode]
                jpoint=elem[ielem].node[jnode]
                kpoint=elem[ielem].node[knode]
                lpoint=elem[ielem].node[lnode]
                nnode=point[ipoint].nnode
                point[ipoint].dNdx[0]+=dNdxT4[inode]/4.
                point[ipoint].dNdy[0]+=dNdyT4[inode]/4.
                point[ipoint].dNdz[0]+=dNdzT4[inode]/4.
                j=self.searchnode(nnode,point[ipoint].node,jpoint)
                point[ipoint].dNdx[j]+=dNdxT4[jnode]/4.
                point[ipoint].dNdy[j]+=dNdyT4[jnode]/4.
                point[ipoint].dNdz[j]+=dNdzT4[jnode]/4.
                k=self.searchnode(nnode,point[ipoint].node,kpoint)
                point[ipoint].dNdx[k]+=dNdxT4[knode]/4.
                point[ipoint].dNdy[k]+=dNdyT4[knode]/4.
                point[ipoint].dNdz[k]+=dNdzT4[knode]/4.
                l=self.searchnode(nnode,point[ipoint].node,lpoint)
                point[ipoint].dNdx[l]+=dNdxT4[lnode]/4.
                point[ipoint].dNdy[l]+=dNdyT4[lnode]/4.
                point[ipoint].dNdz[l]+=dNdzT4[lnode]/4.
                point[ipoint].volume+=elem[ielem].volume/4.




    def searchnode(self,nnode,node,ipoint):
        import sys
        for inode in range(nnode):
            if node[inode]==ipoint:
                i=inode
                return i
        print ipoint,'searchnode miss'
        sys.exit()



    def get_angle(self,ipoint):
        from math import atan2
        if point[ipoint].x[0]!=0.:
            angle=-atan2(point[ipoint].x[1],point[ipoint].x[0])
        else:
            angle=0.
        return angle

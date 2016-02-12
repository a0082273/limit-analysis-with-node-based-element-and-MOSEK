import numpy as np

class Primal():
    def __init__(self,mesh,mate):#mate iranai?
        self.nelem=mesh.nelem
        self.elem=mesh.elem
        self.poi=mesh.Point()
        self.npoint=mesh.npoint
        self.point=mesh.point
        self.mate=mate
        self.mat=mate.Material()
        self.mattype=mate.material_type
        self.reduction=removeHOMOpoint()
        self.asub=[]
        self.aval=[]
        self.xm=[]






    def assemble(self):
        from numpy import linalg as la
        for i in range(6*self.npoint+1):
            self.asub.append([])
            self.aval.append([])
        D=Dmatrix()
        invD=la.inv(D)
        self.poi.get_derivative_of_shape_function_NS()
        for ipoint in range(self.npoint):
            for ilocal in range(self.point[ipoint].nnode):
                B=Bmatrix(self.point,ipoint,ilocal)*self.point[ipoint].volume
                iglobal=self.point[ipoint].node[ilocal]
                R=Rmatrix(iglobal)
                RTBTinvD=R.T.dot(B.T.dot(invD))
                for idir in range(3):
                    if self.point[iglobal].bc[idir]!='HOMOGENEOUSDIRICHLET':
                        for i in range(6):
                            self.asub[6*ipoint+i].\
                                append(self.reduction[3*iglobal+idir])
                            self.aval[6*ipoint+i].append(RTBTinvD[idir,i])

        for ipoint in range(self.npoint):
            for idir in range(3):
                if self.point[ipoint].bc[idir]!='HOMOGENEOUSDIRICHLET':
                    if self.point[ipoint].b1[idir]!=0.:
                        self.asub[6*self.npoint].\
                            append(self.reduction[3*ipoint+idir])
                        self.aval[6*self.npoint].\
                            append(-self.point[ipoint].b1[idir])



    def solve(self):
#
#  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File:    cqo1.py
#
#  Purpose: Demonstrates how to solve small linear
#           optimization problem using the MOSEK Python API.
##
        import sys
        import mosek
        from numpy import array,zeros,ones


# Since the actual value of Infinity is ignores, we define it solely
# for symbolic purposes:
        inf=0.0

# Define a stream printer to grab output from MOSEK
        def streamprinter(text):
            sys.stdout.write(text)
            sys.stdout.flush()
# We might write everything directly as a script, but it looks nicer
# to create a function.
        def mosek_main():
# Make a MOSEK environment
            env=mosek.Env()
# Attach a printer to the environment
            env.set_Stream(mosek.streamtype.log,streamprinter)
  
# Create a task
            task=env.Task(0,0)
# Attach a printer to the task
            task.set_Stream(mosek.streamtype.log,streamprinter)

            numcon=self.reduction[3*self.npoint]

            bkc=[]
            for num in range(numcon):
                bkc.append(mosek.boundkey.fx)

            blc=[]
            for ipoint in range(self.npoint):
                for idir in range(3):
                    if self.point[ipoint].bc[idir]!='HOMOGENEOUSDIRICHLET':
                        blc.append(self.point[ipoint].b0[idir])
            d=dvector()
            D=Dmatrix()
            from numpy import linalg as la
            invD=la.inv(D)
            for ipoint in range(self.npoint):
                for ilocal in range(self.point[ipoint].nnode):
                    B=Bmatrix(self.point,ipoint,ilocal)*self.point[ipoint].volume
                    iglobal=self.point[ipoint].node[ilocal]
                    R=Rmatrix(iglobal)
                    RTBTinvDd=R.T.dot(B.T.dot(invD.dot(d)))
                    for idir in range(3):
                        if self.point[iglobal].bc[idir]!='HOMOGENEOUSDIRICHLET':
                            blc[self.reduction[3*iglobal+idir]]+=RTBTinvDd[idir]

            buc=blc

            c=[]
            for num in range(6*self.npoint):
                c.append(0.)
            c.append(-1.)

            numvar=6*self.npoint+1

            bkx=[]
            for num in range(6*self.npoint):
                bkx.append(mosek.boundkey.fr)
            bkx.append(mosek.boundkey.lo)

            blx=[]
            for num in range(6*self.npoint):
                blx.append(-inf)
            blx.append(0.)

            bux=[]
            for num in range(6*self.npoint):
                bux.append(inf)
            bux.append(inf)
            # bkx=[]
            # for ipoint in range(self.npoint):
            #     for i in range(6):
            #         bkx.append(mosek.boundkey.fr)
            #     bkx.append(mosek.boundkey.up)
            # bkx.append(mosek.boundkey.lo)

            # blx=[]
            # for ipoint in range(self.npoint):
            #     for i in range(6):
            #         blx.append(-inf)
            #     blx.append(-inf)
            # blx.append(1.e-8)

            # bux=[]
            # for ipoint in range(self.npoint):
            #     for i in range(6):
            #         bux.append(inf)
            #     bux.append(0.)
            # bux.append(inf)

            asub=self.asub

            aval=self.aval

            numanz=0
            for row in range(len(self.asub)):
                numanz+=len(self.asub[row])
            numanz+=1

# Append 'numcon' empty constraints.
# The constraints will initially have no bounds. 
            task.appendcons(numcon)
     
#Append 'numvar' variables.
# The variables will initially be fixed at zero (x=0). 
            task.appendvars(numvar)

            for j in range(numvar):
# Set the linear term c_j in the objective.
                task.putcj(j,c[j])
# Set the bounds on variable j
# blx[j] <= x_j <= bux[j] 
                task.putbound(mosek.accmode.var,j,bkx[j],blx[j],bux[j])

            for j in range(len(aval)):
# Input column j of A 
                task.putacol(j,            # Variable (column) index.
                             asub[j],      # Row index of non-zeros in column j.
                             aval[j])      # Non-zero Values of column j. 
            for i in range(numcon):
                task.putbound(mosek.accmode.con,i,bkc[i],blc[i],buc[i])

# Input the cones
            numcone=self.npoint
            if self.mattype=='druckerprager':
                for num in range(numcone):
                    task.appendconeseq(mosek.conetype.quad,0.0,6,6*num)
            # elif self.mattype=='mohrcoulomb':
            #     for num in range(numcone):
            #         task.appendconeseq(mosek.conetype.quad,0.0,6,6*num)

# Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)

# Optimize the task
            task.optimize()
# Print a summary containing information
# about the solution for debugging purposes
            task.solutionsummary(mosek.streamtype.msg)
            prosta=task.getprosta(mosek.soltype.itr)
            solsta=task.getsolsta(mosek.soltype.itr)

# Output a solution
            self.xm=zeros(numvar,float)
            task.getxx(mosek.soltype.itr,
                                self.xm)

            if(solsta==mosek.solsta.optimal or \
               solsta==mosek.solsta.near_optimal):
                print("Optimal solution: %s" % self.xm)
            elif(solsta==mosek.solsta.dual_infeas_cer): 
                print("Primal or dual infeasibility.\n")
            elif(solsta==mosek.solsta.prim_infeas_cer):
                print("Primal or dual infeasibility.\n")
            elif(solsta==mosek.solsta.near_dual_infeas_cer):
                print("Primal or dual infeasibility.\n")
            elif(solsta==mosek.solsta.near_prim_infeas_cer):
                print("Primal or dual infeasibility.\n")
            elif(mosek.solsta.unknown):
                print("Unknown solution status")
            else:
                print("Other solution status")




        def extract_stress():
            import mesh
            from numpy import linalg as la
            D=Dmatrix()
            invD=la.inv(D)
            d=dvector()
            Dstressplusd=np.empty(6,dtype=float)
            for ipoint in range(self.npoint):
                for i in range(6):
                    Dstressplusd[i]=self.xm[6*ipoint+i]
                mesh.point[ipoint].stress=invD.dot(Dstressplusd-d)




        mosek_main()
        extract_stress()









class Dual(Primal):
    def __init__(self,mesh,mate):
        Primal.__init__(self,mesh,mate)
        self.asub=[]
        self.aval=[]
        self.xm=[]




    def assemble(self):
        from math import sin,cos

        for i in range(self.reduction[3*self.npoint]):
            self.asub.append([])
            self.aval.append([])
        for i in range(6*self.npoint+1):
            self.asub.append([])
            self.aval.append([])

        D=Dmatrix()
        from numpy import linalg as la
        invD=la.inv(D)
        for ipoint in range(self.npoint):
            for ilocal in range(self.point[ipoint].nnode):
                B=Bmatrix(self.point,ipoint,ilocal)
                iglobal=self.point[ipoint].node[ilocal]
                R=Rmatrix(iglobal)
                invDTBR=invD.T.dot(B.dot(R))
                for idir in range(3):
                    if self.point[iglobal].bc[idir]!='HOMOGENEOUSDIRICHLET':
                        for i in range(6):
                            self.asub[self.reduction[3*iglobal+idir]].\
                                append(6*ipoint+i)
                            self.aval[self.reduction[3*iglobal+idir]].\
                                append(invDTBR[i,idir])

            for idir in range(3):
# kore nakatta 
                if self.point[ipoint].bc[idir]!='HOMOGENEOUSDIRICHLET':
                    if self.point[ipoint].b1[idir]!=0.:
                        self.asub[self.reduction[3*ipoint+idir]].\
                            append(6*self.npoint)
                        self.aval[self.reduction[3*ipoint+idir]].\
                            append(-self.point[ipoint].b1[idir])

        for num in range(6*self.npoint+1):
            self.asub[self.reduction[3*self.npoint]+num].append(num)
            self.aval[self.reduction[3*self.npoint]+num].append(1.)





    def solve(self):
#
#  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File:    cqo1.py
#
#  Purpose: Demonstrates how to solve small linear
#           optimization problem using the MOSEK Python API.
##
        import sys
        import mosek
        from numpy import array,zeros,ones


# Since the actual value of Infinity is ignores, we define it solely
# for symbolic purposes:
        inf=0.0

# Define a stream printer to grab output from MOSEK
        def streamprinter(text):
            sys.stdout.write(text)
            sys.stdout.flush()
# We might write everything directly as a script, but it looks nicer
# to create a function.
        def mosek_main():
# Make a MOSEK environment
            env=mosek.Env()
# Attach a printer to the environment
            env.set_Stream(mosek.streamtype.log,streamprinter)
  
# Create a task
            task=env.Task(0,0)
# Attach a printer to the task
            task.set_Stream(mosek.streamtype.log,streamprinter)

            numcon=6*self.npoint+1

            bkc=[]
            for num in range(6*self.npoint+1):
                bkc.append(mosek.boundkey.fx);

# Bound values for constraints
            blc=[]
            for num in range(6*self.npoint):
                blc.append(0.)
            blc.append(-1.)

            buc=blc

            c=[]
            for ipoint in range(self.npoint):
                for idir in range(3):
                    if self.point[ipoint].bc[idir]!='HOMOGENEOUSDIRICHLET':
                        c.append(self.point[ipoint].b0[idir])
            d=dvector()
            D=Dmatrix()
            from numpy import linalg as la
            invD=la.inv(D)
            for ipoint in range(self.npoint):
                for ilocal in range(self.point[ipoint].nnode):
                    B=Bmatrix(self.point,ipoint,ilocal)*self.point[ipoint].volume
                    iglobal=self.point[ipoint].node[ilocal]
                    R=Rmatrix(iglobal)
                    RTBTinvDd=R.T.dot(B.T.dot(invD.dot(d)))
                    for idir in range(3):
                        if self.point[iglobal].bc[idir]!='HOMOGENEOUSDIRICHLET':
                            c[self.reduction[3*iglobal+idir]]+=RTBTinvDd[idir]
            for num in range(6*self.npoint+1):
                c.append(0.)

# Bound keys for variables
            numvar=len(self.asub)

            bkx=[]
            for ipoint in range(self.npoint):
                for idir in range(3):
                    if self.point[ipoint].bc[idir]!="HOMOGENEOUSDIRICHLET":
                        bkx.append(mosek.boundkey.fr);
            for num in range(6*self.npoint):
                bkx.append(mosek.boundkey.fr);
            bkx.append(mosek.boundkey.lo);

# Bound values for variables
            blx=[]
            for ipoint in range(self.npoint):
                for idir in range(3):
                    if self.point[ipoint].bc[idir]!="HOMOGENEOUSDIRICHLET":
                        blx.append(-inf)
            for num in range(6*self.npoint):
                blx.append(-inf)
            blx.append(0.)
                        
            bux=[]
            for ipoint in range(self.npoint):
                for idir in range(3):
                    if self.point[ipoint].bc[idir]!="HOMOGENEOUSDIRICHLET":
                        bux.append(inf)
            for num in range(6*self.npoint):
                bux.append(inf)
            bux.append(inf)

# Below is the sparse representation of the A
# matrix stored by column. 
            asub=self.asub

            aval=self.aval

            numanz=0
            for row in range(len(self.asub)):
                numanz+=len(self.asub[row])
            numanz+=1

# Append 'numcon' empty constraints.
# The constraints will initially have no bounds. 
            task.appendcons(numcon)
     
#Append 'numvar' variables.
# The variables will initially be fixed at zero (x=0). 
            task.appendvars(numvar)

            for j in range(numvar):
# Set the linear term c_j in the objective.
                task.putcj(j,c[j])
# Set the bounds on variable j
# blx[j] <= x_j <= bux[j] 
                task.putbound(mosek.accmode.var,j,bkx[j],blx[j],bux[j])

            for j in range(len(aval)):
# Input column j of A 
                task.putacol(j,            # Variable (column) index.
                             asub[j],      # Row index of non-zeros in column j.
                             aval[j])      # Non-zero Values of column j. 
            for i in range(numcon):
                task.putbound(mosek.accmode.con,i,bkc[i],blc[i],buc[i])

# Input the cones
            numcone=self.npoint
            if self.mattype=='druckerprager':
                for num in range(numcone):
                    task.appendconeseq(mosek.conetype.quad,0.0,6,self.reduction[3*self.npoint]+6*num)
            # elif self.mattype=='mohrcoulomb':
            #     for num in range(numcone):
            #         task.appendconeseq(mosek.conetype.quad,0.0,6,6*num)

# Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.maximize)

# Optimize the task
            task.optimize()
# Print a summary containing information
# about the solution for debugging purposes
            task.solutionsummary(mosek.streamtype.msg)
            prosta=task.getprosta(mosek.soltype.itr)
            solsta=task.getsolsta(mosek.soltype.itr)

# Output a solution
            self.xm=zeros(numvar,float)
            task.getxx(mosek.soltype.itr,
                                self.xm)

            if(solsta==mosek.solsta.optimal or \
               solsta==mosek.solsta.near_optimal):
                print("Optimal solution: %s" % self.xm)
            elif(solsta==mosek.solsta.dual_infeas_cer): 
                print("Primal or dual infeasibility.\n")
            elif(solsta==mosek.solsta.prim_infeas_cer):
                print("Primal or dual infeasibility.\n")
            elif(solsta==mosek.solsta.near_dual_infeas_cer):
                print("Primal or dual infeasibility.\n")
            elif(solsta==mosek.solsta.near_prim_infeas_cer):
                print("Primal or dual infeasibility.\n")
            elif(mosek.solsta.unknown):
                print("Unknown solution status")
            else:
                print("Other solution status")






        def extract_velocity_and_plastic_multiplier():
            invRu=np.empty(3,dtype=float)
            for ipoint in range(self.npoint):
                for idir in range(3):
                    if self.point[ipoint].bc[idir]=='HOMOGENEOUSDIRICHLET':
                        invRu[idir]=0.
                    else:
                        invRu[idir]=self.xm[self.reduction[3*ipoint+idir]]
                R=Rmatrix(ipoint)
                self.point[ipoint].u=R.dot(invRu)

            for ipoint in range(self.npoint):
                for i in range(6):
                    self.point[ipoint].s[i]=self.xm[self.reduction\
                                                    [3*self.npoint]+6*ipoint+i]
            sl=self.xm[self.reduction[3*self.npoint]+6*self.npoint]
            print 'sl=',sl



        mosek_main()
        extract_velocity_and_plastic_multiplier()







def removeHOMOpoint():
    import mesh
    reduction=np.empty(3*mesh.npoint+1,dtype=int)
    i=0
    for ipoint in range(mesh.npoint):
        for idir in range(3):
            if mesh.point[ipoint].bc[idir]!='HOMOGENEOUSDIRICHLET':
                reduction[3*ipoint+idir]=i
                i+=1
            else:
                reduction[3*ipoint+idir]=-111
    reduction[3*mesh.npoint]=i
    return reduction



def Dmatrix():
    import materials
    D=np.zeros([6,6],dtype=float)
    alpha=materials.alpha
    D[0,0]=-2.*alpha
    D[0,1]=-2.*alpha
    D[0,2]=-2.*alpha
    D[1,0]=1./3.**0.5
    D[1,1]=1./3.**0.5
    D[1,2]=-2./3.**0.5
    D[2,0]=1.
    D[2,1]=-1.
    D[3,3]=2.
    D[4,4]=2.
    D[5,5]=2.
    return D




def dvector():
    import materials
    d=np.zeros(6,dtype=float)
    k=materials.k
    d[0]=2.*k
    return d



def Rmatrix(ipoint):
    from mesh import Point
    from math import sin,cos
    poi=Point()
    R=np.zeros([3,3],dtype=float)
    angle=poi.get_angle(ipoint)
    # angle=0.
    R[0,0]=cos(angle)
    R[0,1]=sin(angle)
    R[1,0]=-sin(angle)
    R[1,1]=cos(angle)
    R[2,2]=1.
    return R



def Bmatrix(point,ipoint,ilocal):
    B=np.zeros([6,3],dtype=float)
    B[0,0]=point[ipoint].dNdx[ilocal]
    B[1,1]=point[ipoint].dNdy[ilocal]
    B[2,2]=point[ipoint].dNdz[ilocal]
    B[3,0]=point[ipoint].dNdy[ilocal]
    B[3,1]=point[ipoint].dNdx[ilocal]
    B[4,1]=point[ipoint].dNdz[ilocal]
    B[4,2]=point[ipoint].dNdy[ilocal]
    B[5,0]=point[ipoint].dNdz[ilocal]
    B[5,2]=point[ipoint].dNdx[ilocal]
    B/=point[ipoint].volume
    return B




def output(infile,mesh,mat):
    i=len(infile.f.name)
    f=open(infile.f.name[0:i-3]+'vtk','w')
    f.write('# vtk DataFile Version 2.0\n')
    f.write('title\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')

    f.write('POINTS %s double \n' %str(mesh.npoint))
    for ipoint in range(mesh.npoint):
        f.write(str(mesh.point[ipoint].x[0])+' '+\
                str(mesh.point[ipoint].x[1])+' '+\
                str(mesh.point[ipoint].x[2])+'\n')

    f.write('CELLS %s %s \n' %(str(mesh.nelem),str(5*mesh.nelem)))
    for ielem in range(mesh.nelem):
        f.write(str(4)+' '+\
                str(mesh.elem[ielem].node[0])+' '+\
                str(mesh.elem[ielem].node[1])+' '+\
                str(mesh.elem[ielem].node[2])+' '+\
                str(mesh.elem[ielem].node[3])+'\n')

    f.write('CELL_TYPES %s \n' %str(mesh.nelem))
    for ielem in range(mesh.nelem):
        f.write(str(10)+'\n')

    f.write('POINT_DATA %s \n' %str(mesh.npoint))
    f.write('VECTORS displacement_rate double\n')
    for ipoint in range(mesh.npoint):
        for idir in range(3):
            f.write(str(mesh.point[ipoint].u[idir])+' ')
        f.write('\n')
    f.write('VECTORS force double\n')
    for ipoint in range(mesh.npoint):
        for idir in range(3):
            f.write(str(mesh.point[ipoint].b1[idir])+' ')
        f.write('\n')
    f.write("TENSORS stress double\n")
    for ipoint in range(mesh.npoint):
        f.write(str(mesh.point[ipoint].stress[0])+' '+\
                str(mesh.point[ipoint].stress[3])+' '+\
                str(mesh.point[ipoint].stress[4])+'\n')
        f.write(str(mesh.point[ipoint].stress[3])+' '+\
                str(mesh.point[ipoint].stress[1])+' '+\
                str(mesh.point[ipoint].stress[5])+'\n')
        f.write(str(mesh.point[ipoint].stress[4])+' '+\
                str(mesh.point[ipoint].stress[5])+' '+\
                str(mesh.point[ipoint].stress[2])+'\n')
    f.write("TENSORS s double\n")
    for ipoint in range(mesh.npoint):
        f.write(str(mesh.point[ipoint].s[0])+' '+\
                str(mesh.point[ipoint].s[3])+' '+\
                str(mesh.point[ipoint].s[4])+'\n')
        f.write(str(mesh.point[ipoint].s[3])+' '+\
                str(mesh.point[ipoint].s[1])+' '+\
                str(mesh.point[ipoint].s[5])+'\n')
        f.write(str(mesh.point[ipoint].s[4])+' '+\
                str(mesh.point[ipoint].s[5])+' '+\
                str(mesh.point[ipoint].s[2])+'\n')
    f.write("TENSORS epsilon double\n")
    for ipoint in range(mesh.npoint):
        strain=np.zeros(6,dtype=float)
        for ilocal in range(mesh.point[ipoint].nnode):
            B=Bmatrix(mesh.point,ipoint,ilocal)
            iglobal=mesh.point[ipoint].node[ilocal]
            strain+=np.dot(B,mesh.point[iglobal].u)
        strain[3:6]/=2.
        f.write(str(strain[0])+' '+\
                str(strain[3])+' '+\
                str(strain[4])+'\n')
        f.write(str(strain[3])+' '+\
                str(strain[1])+' '+\
                str(strain[5])+'\n')
        f.write(str(strain[4])+' '+\
                str(strain[5])+' '+\
                str(strain[2])+'\n')
    f.write("TENSORS sbyvolume double\n")
    for ipoint in range(mesh.npoint):
        f.write(str(mesh.point[ipoint].s[0]/mesh.point[ipoint].volume)+' '+\
                str(mesh.point[ipoint].s[3]/mesh.point[ipoint].volume)+' '+\
                str(mesh.point[ipoint].s[4]/mesh.point[ipoint].volume)+'\n')
        f.write(str(mesh.point[ipoint].s[3]/mesh.point[ipoint].volume)+' '+\
                str(mesh.point[ipoint].s[1]/mesh.point[ipoint].volume)+' '+\
                str(mesh.point[ipoint].s[5]/mesh.point[ipoint].volume)+'\n')
        f.write(str(mesh.point[ipoint].s[4]/mesh.point[ipoint].volume)+' '+\
                str(mesh.point[ipoint].s[5]/mesh.point[ipoint].volume)+' '+\
                str(mesh.point[ipoint].s[2]/mesh.point[ipoint].volume)+'\n')
    f.write('VECTORS boundary_condition double\n')
    for ipoint in range(mesh.npoint):
        for idir in range(3):
            if mesh.point[ipoint].bc[idir]=='HOMOGENEOUSDIRICHLET':
                f.write(str(1.)+' ')
            else:
                f.write(str(0.)+' ')
        f.write('\n')


    f.close()
    return

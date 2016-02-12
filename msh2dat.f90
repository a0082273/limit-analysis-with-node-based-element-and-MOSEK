program msh2dat
  implicit none
  integer i,j,k,l,in,ie,nelem_tmp,inum_tag,ielem_type,idir,ij
  integer ielem(4),elm_type(4),tag(2),bc(3)
  integer,allocatable::elem1d_tmp(:,:),elem_tmp(:,:)
  real(8) aa,a(3),b(3),c(3),xl(3,4),bj(3,4)
  character(len=20) inputfile,outputfile

  integer,parameter::ndim=3
  integer nprevel,npreforce,nforce,nelem,npoint
  integer,parameter::HOMOGENEOUSDIRICHLET=0,NEUMANN=2,forced=4
  type elem_type
     real(8) v
     integer nd(4)
  end type elem_type
  type(elem_type),allocatable::elem(:)
  type point_type
     real(8) x(ndim),b1(ndim)
     integer bc(ndim),exist
  end type point_type
  type(point_type),allocatable::point(:)

  elm_type(1)=2 ! line element
  elm_type(2)=3 ! triangular element
  elm_type(4)=4 ! quadric element
  read(*,*) inputfile
  ! inputfile='arai.msh'
  ! inputfile='psp.msh'
  ! inputfile='koma.msh'
  ! inputfile='cube.msh'
  i=len_trim(inputfile)
  outputfile=inputfile(1:i-3)//'dat'
  open(10,file=inputfile)
  read(10,*) ! $MeshFormat
  read(10,*) aa,i,j
!!!  if (a/=2.2) stop 'gmsh version is not 2nd'
  if (i/=0) stop 'gmsh version is not 2nd'
  if (j/=8) stop 'imput data is not double precision'
  read(10,*) !$EndMeshFormat
  read(10,*) !$Nodes
  read(10,*) npoint
  allocate(point(npoint))
  do in=1,npoint
     point(in)%bc(:)=NEUMANN
  end do
  do in=1,npoint
     point(in)%b1(:)=0.d0
  end do
  do in=1,npoint
     read(10,*) k,point(in)%x(1),point(in)%x(2),point(in)%x(3)
  end do
  read(10,*) !$EndNodes
  read(10,*) !$Elements
  read(10,*) nelem_tmp
  allocate(elem_tmp(4,nelem_tmp))
  elem_tmp(:,:)=0
  ielem(:)=0
  nelem=0
  do ie=1,nelem_tmp
     read(10,*)  i,ielem_type,inum_tag,(tag(j),j=1,inum_tag),&
          & (ielem(k),k=1,elm_type(ielem_type))
     if(ielem_type==2) then ! triangle
!!!boundary condition!!!
!bc=2->free,bc=0->fixed,bc=4->forced
        do l=1,3
           if(point(ielem(l))%bc(1)==NEUMANN)then
              point(ielem(l))%bc(1)=mod(tag(1)/100,10)
           end if
           if(point(ielem(l))%bc(2)==NEUMANN)then
              point(ielem(l))%bc(2)=mod(tag(1)/10,10)
           end if
           if(point(ielem(l))%bc(3)==NEUMANN)then
              point(ielem(l))%bc(3)=mod(tag(1),10)
           end if
        end do
     else if(ielem_type==4) then ! tetrahedron
        nelem=nelem+1
        elem_tmp(:,nelem)=ielem(1:4)
     end if
  end do
  allocate(elem(nelem))
  do ie=1,nelem
     elem(ie)%nd(:)=elem_tmp(:,ie)
  enddo
  deallocate(elem_tmp)

!!!node existance check!!!
  point(:)%exist=0
  do ie=1,nelem
     point(elem(ie)%nd(:))%exist=1
  end do
  do in=1,npoint
     if(point(in)%exist==0) then
        write(*,*) in,'existance_error'
        stop
     end if
  end do
!!!


  nprevel=0
  do in=1,npoint
     do idir=1,ndim
        if(point(in)%bc(idir)==HOMOGENEOUSDIRICHLET)then
           nprevel=nprevel+1
        end if
     end do
  end do
  npreforce=0




  do ie=1,nelem
     do i=1,3
        do k=1,4
           xl(i,k)=point(elem(ie)%nd(k))%x(i)
        end do
     end do
     bj(1,1)=((xl(2,3)-xl(2,4))*(xl(3,2)-xl(3,4))-(xl(2,2)-xl(2,4))*(xl(3,3)-xl(3,4)))/6.d0
     bj(1,2)=((xl(2,1)-xl(2,3))*(xl(3,1)-xl(3,4))-(xl(2,1)-xl(2,4))*(xl(3,1)-xl(3,3)))/6.d0
     bj(1,3)=((xl(2,1)-xl(2,4))*(xl(3,1)-xl(3,2))-(xl(2,1)-xl(2,2))*(xl(3,1)-xl(3,4)))/6.d0
     bj(1,4)=((xl(2,1)-xl(2,2))*(xl(3,1)-xl(3,3))-(xl(2,1)-xl(2,3))*(xl(3,1)-xl(3,2)))/6.d0
     bj(2,1)=((xl(3,3)-xl(3,4))*(xl(1,2)-xl(1,4))-(xl(3,2)-xl(3,4))*(xl(1,3)-xl(1,4)))/6.d0
     bj(2,2)=((xl(3,1)-xl(3,3))*(xl(1,1)-xl(1,4))-(xl(3,1)-xl(3,4))*(xl(1,1)-xl(1,3)))/6.d0
     bj(2,3)=((xl(3,1)-xl(3,4))*(xl(1,1)-xl(1,2))-(xl(3,1)-xl(3,2))*(xl(1,1)-xl(1,4)))/6.d0
     bj(2,4)=((xl(3,1)-xl(3,2))*(xl(1,1)-xl(1,3))-(xl(3,1)-xl(3,3))*(xl(1,1)-xl(1,2)))/6.d0
     bj(3,1)=((xl(1,3)-xl(1,4))*(xl(2,2)-xl(2,4))-(xl(1,2)-xl(1,4))*(xl(2,3)-xl(2,4)))/6.d0
     bj(3,2)=((xl(1,1)-xl(1,3))*(xl(2,1)-xl(2,4))-(xl(1,1)-xl(1,4))*(xl(2,1)-xl(2,3)))/6.d0
     bj(3,3)=((xl(1,1)-xl(1,4))*(xl(2,1)-xl(2,2))-(xl(1,1)-xl(1,2))*(xl(2,1)-xl(2,4)))/6.d0
     bj(3,4)=((xl(1,1)-xl(1,2))*(xl(2,1)-xl(2,3))-(xl(1,1)-xl(1,3))*(xl(2,1)-xl(2,2)))/6.d0

     elem(ie)%v=xl(1,1)*bj(1,1)+xl(1,2)*bj(1,2)+xl(1,3)*bj(1,3)+xl(1,4)*bj(1,4)
     if(elem(ie)%v<=0) then
        write(*,*) ie,elem(ie)%v,'elem%v<=0'
        stop
     end if

     do j=1,4
        ij=elem(ie)%nd(j)
        point(ij)%b1(3)=point(ij)%b1(3)-elem(ie)%v/4.d0
     end do
  end do

  nforce=0
  do in=1,npoint
     do idir=1,ndim
        if(point(in)%b1(idir)/=0.d0)then
           nforce=nforce+1
        end if
     end do
  end do





  open(50,file=outputfile)
  write(50,'(a)') 'title'
  write(50,'(a)') 'created_by_gmsh'
  write(50,'(a)')
  write(50,'(a)') 'problem'
  write(50,'(a)') 'loadfactor'
  write(50,'(a)')
  write(50,'(a)') 'element_type'
  write(50,'(a)') 'T4'
  write(50,'(a)')
  write(50,'(a)') 'material'
  write(50,'(a)') 'druckerprager'
  write(50,'(a)') '1.0 20.0'
  write(50,'(a)')
  write(50,'(a, i7)') 'elements',nelem
  do ie=1,nelem
     write(50,'(7(i7))') ie,(elem(ie)%nd(i),i=1,4)
  end do
  write(50,'(a)')
  write(50,'(a,i7)') 'node_coordinates',npoint
  do in=1,npoint
     write(50,'(i7,3(e24.17))') in,(point(in)%x(i),i=1,3)
  end do
  write(50,'(a)')
  write(50,'(a,i7)') 'prescribed_velocities',nprevel
  do in=1,npoint
     do idir=1,ndim
        if(point(in)%bc(idir)==HOMOGENEOUSDIRICHLET)then
           write(50,'(2(i7),e24.17)') in,idir,0.d0
        end if
     end do
  end do
  write(50,'(a)')
  write(50,'(a,i7)') 'prescribed_forces',npreforce
  write(50,'(a)')
  write(50,'(a,i7)') 'applied_forces',nforce
  do in=1,npoint
     do idir=1,ndim
        if(point(in)%b1(idir)/=0.d0)then
           write(50,'(2(i7),x1e24.17)') in,idir,point(in)%b1(idir)
        end if
     end do
  end do
  write(50,'(a)')
  write(50,'(a)') 'output 1'
  write(50,'(a)') 'vtk'

end program msh2dat

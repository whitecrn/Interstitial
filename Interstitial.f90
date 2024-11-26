program Interstitial
implicit none
type :: atom
        integer :: atom_type
        real :: x(3)
end type atom
integer :: i,j,k,l,m,n
integer :: alive
integer :: step
integer,dimension(5) :: atom_number
real,parameter :: cutoff_r=2.9
real,dimension(6) :: x_b
real :: dx,dy,dz,dr,time,a,b,c
real,dimension(6) :: vector
type(atom),allocatable :: R(:,:)
integer,allocatable :: co(:,:)
step=0
atom_number(:)=0
!-----------------------------------------------------------------------|
!                                                                       |
!                               read file                               |
!                                                                       |
!-----------------------------------------------------------------------|
open(unit=100,file='dump.atom',status='old',action='read')
open(unit=200,file='data/d-900k.dat',status='replace',action='write')
open(unit=300,file='data/m-900.dat',status='replace',action='write')
open(unit=400,file='o.atom',status='old',action='read')
!-----------------------------------------------------------------------|
!                                                                       |
!                 confirm atom number and crystal constant              |
!                                                                       |
!-----------------------------------------------------------------------| 
do i=1,3
        read(100,*)
end do
read(100,*) atom_number(1)
read(100,*)
read(100,*) x_b(1),x_b(2)
read(100,*) x_b(3),x_b(4)
read(100,*) x_b(5),x_b(6)
a=x_b(2)-x_b(1)
b=x_b(4)-x_b(3)
c=x_b(6)-x_b(5)
rewind(100)
!-----------------------------------------------------------------------|
!                                                                       |
!                               confirm step                            |
!                                                                       |
!-----------------------------------------------------------------------| 
m=1
do while (.true.)
        do i=1,9
                read(100,*,iostat=alive)
                if (alive/=0) then
                        goto 55
                end if
        end do
        
        do i=1,atom_number(1)
                read(100,*) 
        end do
        step=step+1
 end do
!-----------------------------------------------------------------------|
!                                                                       |
!                               read Li atom                            |
!                                                                       |
!-----------------------------------------------------------------------| 
55 continue
rewind(100)
allocate(R(atom_number(1),step))
do i=1,step
        do j=1,9
                read(100,*)
        end do

        do j=1,atom_number(1)
                read(100,*) n,R(n,i)%atom_type,R(n,i)%x(1),R(n,i)%x(2),R(n,i)%x(3)
                R(n,i)%x(1)=R(n,i)%x(1)-x_b(1)
                R(n,i)%x(2)=R(n,i)%x(2)-x_b(3)
                R(n,i)%x(3)=R(n,i)%x(3)-x_b(5)
                if (i==1) then
                        atom_number(1+R(n,i)%atom_type)=atom_number(1+R(n,i)%atom_type)+1
                end if
        end do
end do
!-----------------------------------------------------------------------|
!                                                                       |
!                                O atom                                 |
!                                                                       |
!-----------------------------------------------------------------------|
do i=atom_number(2)+1,atom_number(1)
        write(*,*) 'enter loop:',i
        read(400,*) n,R(n,1)%atom_type,R(n,1)%x(1),R(n,1)%x(2),R(n,1)%x(3)
        R(n,:)%atom_type=R(n,1)%atom_type
        R(n,:)%x(1)=R(n,1)%x(1)
        R(n,:)%x(2)=R(n,1)%x(2)
        R(n,:)%x(3)=R(n,1)%x(3)
end do
!-----------------------------------------------------------------------|
!                                                                       |
!                               statistic                               |
!                                                                       |
!-----------------------------------------------------------------------| 
allocate(co(atom_number(1),step))
co(:,:)=0
do i=1,step  !change here
        m=0
        do j=1,atom_number(2)
                do k=atom_number(2)+1,atom_number(1)
                        dx=(R(j,i)%x(1)-R(k,i)%x(1))
                        if (dx>=0.5*a) then
                                dx=a-dx
                        else if (dx<=-0.5*a) then
                                dx=a+dx
                        end if
                        dy=(R(j,i)%x(2)-R(k,i)%x(2))
                        if (dy>=0.5*b) then
                                dy=b-dy
                        else if (dy<=-0.5*b) then
                                dy=b+dy
                        end if
                        dz=(R(j,i)%x(3)-R(k,i)%x(3))
                        if (dz>=0.5*c) then
                                dz=c-dz
                        else if (dz<=-0.5*c) then
                                dz=c+dz
                        end if
                        dr=sqrt(dx**2+dy**2+dz**2)

                        if (dr<=cutoff_r) then
                                co(j,i)=co(j,i)+1
                        end if
                end do
                !write(200,*) j,i,co(j,i)
                if (co(j,i)>=5) then
                        write(200,*) j,i,co(j,i)
                        m=m+1
                end if
        end do
        if (mod(i,10)==0) then
                call cpu_time(time)
                write(*,*) i,'is done,cpu time: ',time
        end if
        write(300,*) m
end do
!-----------------------------------------------------------------------|
!                                                                       |
!                               Vector calculation                      |
!                                                                       |
!-----------------------------------------------------------------------|

stop
end program
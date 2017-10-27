Module TPSAMod
  use MathFuncMod
  use PolyMapMod
  
  implicit none
  integer, parameter ::MAX_TPS_ORDERS = 6
  integer :: Maximum_TPS_Degree,TPS_Dim
  type (polymap) :: pmap

  type dctps
    integer, private          :: degree
    integer(kind=4), private  :: terms
    double precision, allocatable, dimension(:) :: map
    integer(kind=4),  allocatable, dimension(:) :: index
    
  contains
    PROCEDURE, PASS :: redegree => TpsaRedegree
    PROCEDURE, PASS :: output => TpsaPrint
    
    PROCEDURE :: pTpsaAdd   =>  TpsaAdd
    GENERIC :: OPERATOR(+)  => pTpsaAdd
    PROCEDURE :: pTpsaDec   =>  TpsaDec
    GENERIC :: OPERATOR(-)  => pTpsaDec
    PROCEDURE :: pTpsaMul   =>  TpsaMul
    GENERIC :: OPERATOR(*)  => pTpsaMul
    PROCEDURE :: pTpsaDiv   =>  TpsaDiv
    GENERIC :: OPERATOR(/)  => pTpsaDiv
    
    !PROCEDURE :: pTpsaEq      =>  TpsaEq
    PROCEDURE :: pTpsaEqN     =>  TpsaEqN
    GENERIC :: assignment(=)  => pTpsaEqN
    
    PROCEDURE :: pTpsaIsEq    =>  TpsaIsEq
    GENERIC :: OPERATOR(.eq.) => pTpsaIsEq
    GENERIC :: OPERATOR(==)   => pTpsaIsEq
    
    PROCEDURE :: pTpsaNotEq   =>  TpsaNotEq
    GENERIC :: OPERATOR(.ne.) => pTpsaNotEq
  end type dctps
  
  interface dctps
    module procedure constructor
  end interface
  
  interface assign
    module procedure assignV1,assignV2
  end interface
  
  interface operator (+)
    module procedure TpsaAddN1,TpsaAddN2,TpsaAddN3
  end interface
  
  interface operator (-)
    module procedure TpsaDecN1,TpsaDecN2,TpsaDecN3
  end interface
  
  interface operator (*)
    module procedure TpsaMulN1,TpsaMulN2
  end interface
  
  interface operator (/)
    module procedure TpsaDivN1,TpsaDivN2
  end interface
  
  interface exp
    module procedure TpsaExp
  end interface
  
  interface log
    module procedure TpsaLog
  end interface
  
  interface sqrt
    module procedure TpsaSqrt
  end interface
  
  interface pow
    module procedure TpsaPow
  end interface
  
  interface sin
    module procedure TpsaSin
  end interface
  
  interface cos
    module procedure TpsaCos
  end interface
  
  interface tan
    module procedure TpsaTan
  end interface
  
  interface sinh
    module procedure TpsaSinh
  end interface
  
  interface cosh
    module procedure TpsaCosh
  end interface
  
contains

    subroutine dctps_Initialize(dim,order)
    
        implicit none
        integer, intent(in) :: dim
        integer, optional, intent(in) :: order
        
        !write(*,*) 'dctps_Initialize works'
        
        TPS_Dim   = dim
        if(.not. present(order)) then
          Maximum_TPS_Degree = MAX_TPS_ORDERS
        else
          Maximum_TPS_Degree = order
        endif
        !write(*,*) TPS_Dim,Maximum_TPS_Degree
        pmap = polymap(TPS_Dim,Maximum_TPS_Degree)
        
    end subroutine
    
    subroutine TpsaRedegree(this,degree)
        implicit none
        class(dctps) :: this
        integer  :: degree,i
        double precision, allocatable, dimension(:) :: temp
        
        this%degree = degree
        if (this%degree > Maximum_TPS_Degree) this%degree=Maximum_TPS_Degree;
        this%terms=binomial(TPS_Dim+degree, degree);
        
        if (ALLOCATED(this % map)) then
        
          allocate( temp(size(this%map)) )
          temp = this%map
          
          deallocate(this % map)
          allocate(this % map(this%terms))
        
          do i=1,min(size(this%map),size(temp))
            this % map(i) = temp(i)
          enddo
          do i=min(size(this%map),size(temp))+1,size(this%map)
            this % map(i) = 0
          enddo
        
          deallocate(temp)
        else
          allocate(this % map(this%terms))
          this%map=0
        endif
    end subroutine
    
    subroutine assignV1(this,a)
        implicit none
        class(dctps) :: this
        double precision :: a
        
        this%degree=0
        this%terms=1
        if (ALLOCATED(this % map)) deallocate(this % map)
        allocate(this % map(1))
        this%map(1)     = a

    end subroutine assignV1
    
    subroutine assignV2(this,a,n_var)
        implicit none
        class(dctps) :: this
        integer          :: n_var,i
        double precision :: a
        
        if (n_var <= TPS_Dim .and. n_var>0) then 
          this%degree=1
          this%terms=TPS_Dim+1;
          if (ALLOCATED(this % map)) deallocate(this % map)
          allocate(this % map(binomial(TPS_Dim+Maximum_TPS_Degree,TPS_Dim)))
          do i=1,this%terms
            this%map(i)   = 0.0
          enddo
          this%map(n_var+1) = 1.0
          this%map(1)     = a
        else
          print*, "Error: Num of var out of range in CTPS"
        endif
    end subroutine assignV2
    
    double precision Function TpsaElement(this,inivalue,siz) result(res)
      class(dctps) ,Intent(In) :: this
      integer      ,Intent(In) :: siz
      integer      ,dimension(siz), Intent(In) :: inivalue
      integer(kind=4) :: ind
      
      ind = TpsaFindindex(inivalue,siz)
      res = this%map(ind+1)
    End Function TpsaElement
    
    double precision Function TpsaEvaluate(this,inivalue,siz) result(res)
      class(dctps) ,Intent(In) :: this
      integer      ,Intent(In) :: siz
      double precision ,dimension(siz), Intent(In) :: inivalue
      integer :: i,j
      double precision :: product
      if (siz .ne. TPS_Dim) then
        write(*,*) "Inconsistent dimension to evaluate CTPS"
        stop
      endif
      
      res = this%map(1)
      do i=2,this%terms
        product = 1.0d0
        do j=1,this%terms
          product = product * inivalue(j)**pmap%map(j+1,i)
        enddo
        res = res + product * this%map(i)
      enddo
    End Function TpsaEvaluate
    
    Type(dctps) Function TpsaDerivative(this,ndim,order) result(res)
      implicit none
      class(dctps) ,Intent(In) :: this
      integer      ,Intent(In) :: ndim, order
      integer(kind=4) :: new_i
      integer :: i,new_max_order,thisdim
      integer, allocatable, dimension(:) :: indexlist
      
      if(order <= 0) then 
        res = this
        return
      elseif(ndim <= TPS_Dim .and. ndim>0) then
        res          = this * 0.0d0
        !print*, size(this%map),size(res%map)
        new_max_order = 0
        allocate(indexlist(size(pmap%map,1)))
        do i=1,this%terms-1
          indexlist = pmap%map(:,i+1)
          if(indexlist(ndim+1) >= order) then
            thisdim             = indexlist(ndim+1)
            indexlist(ndim+1)   = indexlist(ndim+1) - order
            indexlist(     1)   = indexlist(     1) - order
            if (new_max_order<indexlist(1)) new_max_order=indexlist(1)
            new_i               = TpsaFindindex(indexlist,size(indexlist))
            print*, new_i+1,i+1,thisdim,order,binomial(thisdim,order)
            res%map(new_i+1) = this%map(i+1) * DBLE(binomial(thisdim,order))
          endif
        enddo
        deallocate(indexlist)
      else
        write(*,*) "Inconsistent dimension to take derivative"
        call assign(res,1.0d0)
      endif
    End Function TpsaDerivative

    Type(dctps) Function TpsaAdd(M,N)
      implicit none
      class(dctps) ,Intent(In) :: M,N
      integer :: i
      
      if(M%degree > N%degree) then
        TpsaAdd = M
        do i= 1,N%terms
          TpsaAdd%map(i) = TpsaAdd%map(i) + N%map(i)
        enddo
      else
        TpsaAdd = N
        do i= 1,M%terms
          TpsaAdd%map(i) = TpsaAdd%map(i) + M%map(i)
        enddo
      endif
    End Function TpsaAdd
    
    Type(dctps) Function TpsaAddN1(M,N) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      double precision, Intent(in) :: N
      type(dctps)  :: temp
      
      call assign(temp,N)
      this = TpsaAdd(M,temp)
    End Function TpsaAddN1
    
    Type(dctps) Function TpsaAddN2(N,M) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      double precision, Intent(in) :: N
      type(dctps)  :: temp
      
      call assign(temp,N)
      this = TpsaAdd(M,temp)
    End Function TpsaAddN2
    
    Type(dctps) Function TpsaAddN3(M) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      
      this = M
    End Function TpsaAddN3
    
    Type(dctps) Function TpsaDec(M,N)
      implicit none
      class(dctps) ,Intent(In) :: M,N
      integer :: i
      
      if(M%degree > N%degree) then
        TpsaDec = M
        do i= 1,N%terms
          TpsaDec%map(i) = TpsaDec%map(i) - N%map(i)
        enddo
      else
        TpsaDec = -N+M
      endif
    End Function TpsaDec
    
    Type(dctps) Function TpsaDecN1(M,N) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      double precision, Intent(in) :: N
      type(dctps)  :: temp
      
      call assign(temp,N)
      this = TpsaDec(M,temp)
    End Function TpsaDecN1
    
    Type(dctps) Function TpsaDecN2(N,M) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      double precision, Intent(in) :: N
      type(dctps)  :: temp
      
      call assign(temp,N)
      this = TpsaDec(temp,M)
    End Function TpsaDecN2
    
    Type(dctps) Function TpsaDecN3(M) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      type(dctps)  :: temp
      
      call assign(temp,-1.0d0)
      this = M * temp
    End Function TpsaDecN3
    
    
    Type(dctps) Function TpsaMul(M,N) result(this)
      implicit none
      class(dctps) ,Intent(In) :: M,N
      
      if (M%degree > N%degree) then
        this = TpsaMulM(M,N)
      else
        this = TpsaMulM(N,M)
      endif
      
    End Function TpsaMul
    Type(dctps) Function TpsaMulM(M,N) result(this)
      implicit none
      class(dctps) ,Intent(In) :: M,N
      type(dctps) :: temp
      integer :: i, j, k, newdegree
      integer(kind=4) :: j_max, target_ind
      integer, allocatable, dimension(:) :: indexmap
      
      if(M%degree >= N%degree) then
        this = M
        if(N%degree == 0) then
          do i=1,this%terms
            this%map(i) = this%map(i) * N%map(1)
          enddo
          return
        else
          newdegree=this%degree+N%degree
          temp = this
          deallocate(this%map)
          call this%redegree(min(Maximum_TPS_Degree, newdegree))
          !print*, size(this%map),size(temp%map),size(N%map)
          
          allocate(indexmap(size(pmap%map,1)))
          do i=1,size(temp%map)
            if(abs(temp%map(i))<1e-14) cycle
            
            j_max=min(size(N%map),binomial(TPS_Dim+Maximum_TPS_Degree-pmap%map(1,i),TPS_Dim))
            do j=1,j_max
              if(abs(N%map(j))<1e-14) cycle
              
              do k=1,size(pmap%map,1)
                indexmap(k) =  pmap%map(k,i) + pmap%map(k,j)
              enddo
              target_ind = TpsaFindindex(indexmap,size(indexmap))
              !print*,target_ind,this%map(2)
              this%map(target_ind+1) = this%map(target_ind+1) + temp%map(i) * N%map(j)
            enddo
          enddo
          deallocate(indexmap)
          return
        endif
      else
      endif
    End Function TpsaMulM
    
    Type(dctps) Function TpsaMulN1(M,N) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      double precision, Intent(in) :: N
      type(dctps)  :: temp
      
      call assign(temp,N)
      this = TpsaMul(M,temp)
    End Function TpsaMulN1
    
    Type(dctps) Function TpsaMulN2(N,M) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      double precision, Intent(in) :: N
      type(dctps)  :: temp
      
      call assign(temp,N)
      this = TpsaMul(temp,M)
    End Function TpsaMulN2
    
    Type(dctps) Function TpsaDiv(M,N) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M,N
      integer ::i
      
      if(abs(N%map(1)) <1.0e-14) then
        write(*,*) "Error: Divide by zero, in CTPS"
        return
      endif
      
      if(N%degree==0) then
        this = M
        do i=1,this%terms
          this % map(i) = this%map(i) / N%map(1)
        enddo
      else
        this = M * TpsaInv(N)
      endif
      
    End Function TpsaDiv
    
    Type(dctps) Function TpsaDivN1(M,N) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      double precision, Intent(in) :: N
      type(dctps)  :: temp
      
      call assign(temp,N)
      this = TpsaDiv(M,temp)
    End Function TpsaDivN1
    
    Type(dctps) Function TpsaDivN2(N,M) result(this)
      implicit none
      class(dctps)    , Intent(in) :: M
      double precision, Intent(in) :: N
      type(dctps)  :: temp
      
      call assign(temp,N)
      this = TpsaDiv(temp,M)
    End Function TpsaDivN2
    
    Subroutine TpsaEq(M,N)
      implicit none
      class(dctps) ,Intent(out) :: M
      class(dctps) ,Intent(In)  :: N
      
      M%degree = N%degree
      M%terms  = N%terms
      M%map    = N%map
    End Subroutine TpsaEq
    
    Subroutine TpsaEqN(M,N)
      implicit none
      class(dctps) ,Intent(out) :: M
      double precision ,Intent(In)  :: N
      
      call assign(M,N)
    End Subroutine TpsaEqN
    
        
    Logical Function TpsaIsEq(M,N) result(res)
      implicit none
      class(dctps) ,Intent(in) :: M,N
      integer :: i,si
      
      si = size(M%map)
      if(si .ne. size(N%map)) then 
        res = .false.
        return
      endif
      
      do i=1,si
        if(abs(M%map(i)) > 1.0d-30 .or. abs(N%map(i)) > 1.0d-30) then
          if(abs(M%map(i) / N%map(i)- 1.0d0) > 1.0d-15) then
            res = .false.
            print*, "BBB",M%map(i),N%map(i)
            return
          endif
        endif
      enddo
      
      res = .true.
    End Function TpsaIsEq
    
    Logical Function TpsaNotEq(M,N) result(res)
      class(dctps) ,Intent(in) :: M,N
      
      res = .not. M.eq.N
    End Function TpsaNotEq
    
    Type(dctps) Function TpsaInv(M) result(this)
      implicit none
      class(dctps), Intent(in) :: M
      type(dctps)  :: temp, temp2, term_by_order
      integer      :: i
      
      if (abs(M%map(1)) < 1e-14) then 
        write(*,*) "Error: Divide by zero, in CTPS"
      endif
      
      temp = -(M - M%map(1))  / M%map(1)
      call assign(temp2,1.0d0)
      term_by_order = temp2 / M%map(1)
      
      call assign(this,0.0d0)
      this = this + term_by_order
      
      do i=1,Maximum_TPS_Degree
        term_by_order = term_by_order * temp
        this = this + term_by_order
      enddo
    End Function TpsaInv
    
    
    Type(dctps) Function TpsaExp(M) result(this)
      implicit none
      class(dctps), Intent(in) :: M
      type(dctps)  :: temp, term_by_order
      integer      :: i
      double precision :: index

      temp = M - M%map(1)
      call assign(term_by_order,1.0d0)
      call assign(this,0.0d0)
      this = this + term_by_order
      
      do i=1,Maximum_TPS_Degree
        index = 1.0d0 / factorial(i)
        term_by_order = term_by_order * temp
        this = this + term_by_order * index
      enddo
      this = this * exp(M%map(1))
    End Function TpsaExp
    
    
    Type(dctps) Function TpsaLog(M) result(this)
      implicit none
      class(dctps), Intent(in) :: M
      type(dctps)  :: temp, term_by_order
      integer      :: i

      temp = -(M - M%map(1))/M%map(1)
      term_by_order = -temp
      call assign(this,0.0d0)
      this = this + term_by_order
      
      do i=2,Maximum_TPS_Degree
        term_by_order = term_by_order * temp
        this = this + term_by_order / DBLE(i)
      enddo
      this = this + log(M%map(1))
    End Function TpsaLog
    
    Type(dctps) Function TpsaSqrt(M) result(this)
      implicit none
      class(dctps), Intent(in) :: M
      type(dctps)  :: temp, term_by_order
      integer      :: i
      double precision :: a0, index

      a0   = sqrt(M%map(1))
      temp = -(M - M%map(1))/ a0
      term_by_order = -temp 
      call assign(this,0.0d0)
      this = this + term_by_order / 2.0d0
      
      do i=2,Maximum_TPS_Degree
        index = 1.0d0 *  doublefactorial(2*i - 3) / doublefactorial(2 * i)
        term_by_order = term_by_order * temp
        this = this + term_by_order * index
      enddo
      this = this + a0
    End Function TpsaSqrt
    
    Type(dctps) Function TpsaPow(M,b) result(this)
      implicit none
      class(dctps),     Intent(in) :: M
      double precision, Intent(in) :: b
      type(dctps)  :: temp, term_by_order
      integer      :: i
      double precision :: index, factor
      
      if(b-1.0d0<1e-14) then
        this = M
        return
      else if (abs(b)<1e-14) then
        call assign(this,1.0d0)
        return
      else
        index = b
        temp = M - M%map(1)
        call assign(term_by_order,1.0d0)
        factor = M%map(1)**b
        call assign(this,factor)
        do i=1,Maximum_TPS_Degree
          factor = factor / M%map(1) * DBLE(index / i)
          index = index -1
          term_by_order = term_by_order * temp
          this = this + term_by_order * factor
          if(abs(index) < 1.0e-14) exit 
        enddo
      endif
    End Function TpsaPow
    
    Type(dctps) Function TpsaSin(M) result(this)
      implicit none
      class(dctps), Intent(in) :: M
      type(dctps)  :: temp, term_by_order
      integer      :: i
      double precision :: a0, sin_a0, cos_a0, index

      a0     = M%map(1)
      sin_a0 = sin(a0)
      cos_a0 = cos(a0)
      temp = M - a0
      call assign(term_by_order,1.0d0)
      call assign(this,0.0d0)
      
      do i=1,Maximum_TPS_Degree
        if( mod(i,2) ==1) then
          index = cos_a0 * DBLE((-1.0)**((i - 1) / 2) / factorial(i))
        else
          index = sin_a0 * DBLE((-1.0)**( i      / 2) / factorial(i))
        endif
        term_by_order = term_by_order * temp
        this = this + term_by_order * index
      enddo
      this = this + sin_a0
    End Function TpsaSin
    
    Type(dctps) Function TpsaCos(M) result(this)
      implicit none
      class(dctps), Intent(in) :: M
      type(dctps)  :: temp, term_by_order
      integer      :: i
      double precision :: a0, sin_a0, cos_a0, index

      a0     = M%map(1)
      sin_a0 = sin(a0)
      cos_a0 = cos(a0)
      temp = M - a0
      call assign(term_by_order,1.0d0)
      call assign(this,0.0d0)
      
      do i=1,Maximum_TPS_Degree
        if( mod(i,2) ==1) then
          index = sin_a0 * DBLE((-1.0)**((i + 1) / 2) / factorial(i))
        else
          index = cos_a0 * DBLE((-1.0)**( i      / 2) / factorial(i))
        endif
        term_by_order = term_by_order * temp
        this = this + term_by_order * index
      enddo
      this = this + cos_a0
    End Function TpsaCos
    
    Type(dctps) Function TpsaTan(M) result(this)
      implicit none
      class(dctps), Intent(in) :: M
      this = sin(M)/cos(M)
    End Function TpsaTan
    
    Type(dctps) Function TpsaSinh(M) result(this)
      implicit none
      class(dctps), Intent(in) :: M
      type(dctps)  :: temp, term_by_order
      integer      :: i
      double precision :: a0, sinh_a0, cosh_a0, index

      a0     = M%map(1)
      sinh_a0 = sinh(a0)
      cosh_a0 = cosh(a0)
      temp = M - a0
      call assign(term_by_order,1.0d0)
      call assign(this,0.0d0)
      
      do i=1,Maximum_TPS_Degree
        if( mod(i,2) ==1) then
          index = cosh_a0 / DBLE(factorial(i))
        else
          index = sinh_a0 / DBLE(factorial(i))
        endif
        term_by_order = term_by_order * temp
        this = this + term_by_order * index
      enddo
      this = this + sinh_a0
    End Function TpsaSinh
    
    Type(dctps) Function TpsaCosh(M) result(this)
      implicit none
      class(dctps), Intent(in) :: M
      type(dctps)  :: temp, term_by_order
      integer      :: i
      double precision :: a0, sinh_a0, cosh_a0, index

      a0     = M%map(1)
      sinh_a0 = sinh(a0)
      cosh_a0 = cosh(a0)
      temp = M - a0
      call assign(term_by_order,1.0d0)
      call assign(this,0.0d0)
      
      do i=1,Maximum_TPS_Degree
        if( mod(i,2) ==1) then
          index = sinh_a0 / DBLE(1.0*factorial(i))
        else
          index = cosh_a0 / DBLE(1.0*factorial(i))
        endif
        !print*,index
        term_by_order = term_by_order * temp
        this = this + term_by_order * index
      enddo
      this = this + cosh_a0
    End Function TpsaCosh
    
    Function TpsaFindindex(indexmap,si) result(res)
      integer :: si,i
      integer, dimension(si) :: indexmap,summ
      integer(kind=4) :: res
      
      if(si .ne. TPS_Dim+1) then 
        write(*,*) "Index map does not have correction length"
        stop
      endif
      
      summ(1) = indexmap(1)

      do i = 2,TPS_Dim+1
         if (indexmap(i)<0) then 
           stop "The index map has invalid component"
         endif
         summ(i)=summ(i-1)-indexmap(i);
      enddo
      
      res = 0
      
      do i = TPS_Dim+1,2,-1
         if (summ(TPS_Dim-i+2)==0) exit 
         res = res + binomial(summ(TPS_Dim-i+2)-2+i,i-1)
      enddo
    End Function
    
    subroutine TpsaPrint(this)
      class(dctps) ,Intent(In) :: this
      integer :: current_order, acc_return, just_endl
      integer :: i
      
      current_order = 0
      acc_return    = 0
      just_endl     = 0
      write(*,*) "Order", 0,":"

      do i = 1,this%terms
        if(abs(this%map(i)) <1.0E-14 .and. i>1) cycle
        if(pmap%map(1,i)>current_order) then
          write(*,*) "Order", pmap%map(1,i), ":"
          current_order = pmap%map(1,i)
          acc_return = 0
        endif
        write(*,*) pmap%map(2:,i),this%map(i)
        acc_return = acc_return + 1
      enddo
    end subroutine

    subroutine pushPtc_TPSA(ptc,tpsa)
      !assume it is (x,px,y,py,z,pz)
      integer                            :: nPtc
      double precision, dimension(:,:),    Intent(Inout) :: ptc
      type(dctps), dimension(:),      Intent(In)    :: tpsa
      double precision,      dimension(6)                     :: ptcTemp
      double precision :: sumTemp
      integer :: i,j,k,m
      
      nPtc = size(ptc,2)
      do k = 1,nPtc
        ptcTemp=0.d0
        do j = 1,6
          do i = 1,tpsa(j)%terms
            if(abs(tpsa(j)%map(i)) <1.0E-14 .and. i>1) cycle
            sumTemp = tpsa(j)%map(i)
            do m = 1,6
              if(pmap%map(m+1,i) ==0 ) cycle
              sumTemp = sumTemp * ptc(m,k) ** pmap%map(m+1,i)
            enddo
            ptcTemp(j) = ptcTemp(j) + sumTemp
            enddo
        enddo
        ptc(:,k) = ptcTemp
      enddo
    end subroutine pushPtc_TPSA

End Module TPSAMod

subroutine NpPotential(x_t, y_t, z_t, pier)  

  use M_grid
  use M_bas
  !aggiunta per gestione file e scritte da Leo_e_Pier
  use KF
  !fine
  use Vartypes

  implicit none

  type(grid)     :: G     ! integration grid 
  type(type_bas) :: BAS   ! basis functions

  integer(KINT) :: iblock, i, mu, nu, iu15     ! loop indices
  real(KREAL), allocatable :: bas0(:,:)   ! basis functions
  real(KREAL), allocatable :: matrix(:)   ! matrix in basisrepresentation, triangular
  real(KREAL), allocatable :: op(:)       ! operator
  real(KREAL), intent(in) :: x_t, y_t, z_t
  integer(KINT), intent(in) :: pier
  character(5000) :: tessera_attuale

  call create (G)
  call create (BAS, G)
  allocate ( bas0 (G%npoints,BAS%num) )
  allocate ( matrix (BAS%num_tri) )
  allocate (op(G%npoints))

  matrix  = 0.0_KREAL

  call timers ('pp')

  iblock_: do iblock = 1, G%nblocks

     if (skip_block(G, iblock)) cycle iblock_

     call get_block(G)
     call calc_bas_block (BAS, G, bas0)
     do i=1,G%npoints
        op(i) = 1/sqrt((x_t - G%coord(i,1))**2 + (y_t - G%coord(i,2))**2 + (z_t - G%coord(i,3))**2)
     enddo
     op = op * G%w

     call calc_bas0_op_bas0_block (BAS, G%npoints, bas0, op, matrix)
     
  enddo iblock_

  call ppcbnr (matrix, size(matrix), 'matrix')
  call timere ('pp')
  if (pier.lt.10) then 
     write(tessera_attuale, '(i1.1)') pier
  elseif (pier.ge.10.and.pier.lt.100) then
     write(tessera_attuale, '(i2.2)') pier
  elseif (pier.ge.100.and.pier.lt.1000) then
     write(tessera_attuale, '(i3.3)') pier
  elseif (pier.ge.1000.and.pier.lt.10000) then
     write(tessera_attuale, '(i4.4)') pier
  elseif (pier.ge.10000) then
     write(iuout,*) 'ATTENTION! Too much tesserae in the input'
  endif 
 


  call kfopfl(iu15, 'TAPE15')
  call kfwrnr (iu15, 'Matrices%potential_tessera'//trim(tessera_attuale), matrix, BAS%num_tri, 1)
  call kfclfl (iu15)

  deallocate (bas0, matrix)
  call delete (BAS)
  call delete (G)
endsubroutine NpPotential

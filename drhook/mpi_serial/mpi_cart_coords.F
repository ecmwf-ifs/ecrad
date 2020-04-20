      subroutine mpi_cart_coords(comm,rank,maxdims,coords,ierr)
      integer comm, rank, maxdims
      integer coords(maxdims)

      coords(1:maxdims) = 0
      ierr = 0
      return
      end

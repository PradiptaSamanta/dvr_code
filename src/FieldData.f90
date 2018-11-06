module FieldData

  use constants

  implicit none

  public

  complex(idp), target, allocatable :: PrimFieldInts(:,:,:,:)
  complex(idp), target, allocatable :: FieldInts(:,:,:)

end module FieldData

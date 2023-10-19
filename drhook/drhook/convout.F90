! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

INTEGER FUNCTION CONVOUT(KCOUNT,KTYPE)
REAL ZCONV(0:5)
INTEGER KCOUNT,KTYPE
ZCONV=(/0.25,1.0,2.0,2.0,1.0,1.0/)
CONVOUT=CEILING(KCOUNT/ZCONV(KTYPE))
RETURN
END FUNCTION CONVOUT

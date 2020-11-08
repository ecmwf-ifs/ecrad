/**
 * (C) Copyright 2014- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

INTERFACE
SUBROUTINE EC_MEMINFO(KU,CDSTRING,KCOMM,KBARR,KIOTASK,KCALL)
USE PARKIND1, ONLY : JPIM
INTEGER(KIND=JPIM), INTENT(IN) :: KU,KCOMM,KBARR,KIOTASK,KCALL
CHARACTER(LEN=*), INTENT(IN) :: CDSTRING
END SUBROUTINE EC_MEMINFO
END INTERFACE

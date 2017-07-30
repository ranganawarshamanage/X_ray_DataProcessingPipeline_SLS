MODULE resolution_mod

REAL :: a=0.,b,c,alf,bet,gam
REAL(kind(1.d0)), PRIVATE :: as,bs,cs,cosgs,cosbs,cosas
LOGICAL, PRIVATE :: firstcall=.TRUE.

contains

SUBROUTINE getres(ih,ik,il,s2,resol)
IMPLICIT NONE
INTEGER, INTENT(IN):: ih,ik,il
REAL, INTENT(OUT)  :: s2,resol
REAL(kind(1.d0)) :: vol,sina,sinb,sing,cosa,cosb,cosg

IF (firstcall) THEN
  firstcall = .FALSE.
  if (a==0.) stop 'a,b,c,alf,bet,gam unknown'
  cosa=COS(alf/57.295779513082320877d0)
  cosb=COS(bet/57.295779513082320877d0)
  cosg=COS(gam/57.295779513082320877d0)
  sina=SIN(alf/57.295779513082320877d0)
  sinb=SIN(bet/57.295779513082320877d0)
  sing=SIN(gam/57.295779513082320877d0)
  cosas=(cosb*cosg-cosa)/(sinb*sing)
  cosbs=(cosa*cosg-cosb)/(sina*sing)
  cosgs=(cosa*cosb-cosg)/(sina*sinb)
  vol=a*b*c*SQRT(1.d0-cosa**2-cosb**2-cosg**2+2.d0*cosa*cosb*cosg)
  as=b*c*sina/vol
  bs=a*c*sinb/vol
  cs=a*b*sing/vol
END IF
s2=((ih*as)**2+(ik*bs)**2+(il*cs)**2+2.d0*ih*ik*as*bs*cosgs  &
      +2.d0*ih*il*as*cs*cosbs+2.d0*ik*il*bs*cs*cosas)/4.d0
         resol = SQRT(S2)
         RESOL=1/(2.*resol)
END SUBROUTINE getres
END MODULE resolution_mod

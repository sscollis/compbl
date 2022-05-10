c------------------------------------------------------------------------------
c.... Routines from Numerical Recipes
c------------------------------------------------------------------------------
      SUBROUTINE NR_RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
C     implicit double precision (a-h,o-z)
      PARAMETER (NMAX=5)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      external derivs
      HH=H*0.5d0
      H6=H/6.0d0
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.0d0*DYM(I))
14    CONTINUE
      RETURN
      END
c
      SUBROUTINE NR_RUNGE(VSTART,NVAR,X1,X2,NSTEP,DERIVS)
C     implicit double precision (a-h,o-z)
      PARAMETER (NMAX=5,NSTPMX=5000)
      COMMON /NR_RUNGE_PATH/ XX(NSTPMX),Y(NMAX,NSTPMX)
      DIMENSION VSTART(NVAR),V(NMAX),DV(NMAX)
      external derivs
      DO 11 I=1,NVAR
        V(I)=VSTART(I)
        Y(I,1)=V(I)
11    CONTINUE
      XX(1)=X1
      X=X1
      H=(X2-X1)/NSTEP
      DO 13 K=1,NSTEP
        CALL DERIVS(X,V,DV)
        CALL NR_RK4(V,DV,NVAR,X,H,V,DERIVS)
        IF(X+H.EQ.X)PAUSE 'Stepsize not significant in RKDUMB.'
        X=X+H
        XX(K+1)=X
        DO 12 I=1,NVAR
          Y(I,K+1)=V(I)
12      CONTINUE
13    CONTINUE
      RETURN
      END

      SUBROUTINE NR_ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,
     &                     RKQC)
C     implicit double precision (a-h,o-z)
      PARAMETER (MAXSTP=5000,NMAX=5,TWO=2.0,ZERO=0.0,TINY=1.E-30)
      COMMON /NR_ODEINT_PATH/ KMAX,KOUNT,DXSAV,XP(MAXSTP),YP(NMAX,MAXSTP)
      DIMENSION YSTART(NVAR),YSCAL(NMAX),Y(NMAX),DYDX(NMAX)
      external derivs, rkqc
      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=1
      XP(KOUNT)=X
      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
        YP(I,KOUNT)=Y(I)
11    CONTINUE
      IF (KMAX.GT.0) XSAV=X-DXSAV*TWO
      DO 16 NSTP=1,MAXSTP
        CALL DERIVS(X,Y,DYDX)
        DO 12 I=1,NVAR
          YSCAL(I)=ABS(Y(I))+ABS(H*DYDX(I))+TINY
12      CONTINUE
        IF(KMAX.GT.0)THEN
          IF(ABS(X-XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX-1)THEN
              KOUNT=KOUNT+1
              XP(KOUNT)=X
              DO 13 I=1,NVAR
                YP(I,KOUNT)=Y(I)
13            CONTINUE
              XSAV=X
            ENDIF
          ENDIF
        ENDIF
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        ENDIF
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
c          DO 14 I=1,NVAR
c            YSTART(I)=Y(I)
c14        CONTINUE
          IF(KMAX.NE.0)THEN
            KOUNT=KOUNT+1
            XP(KOUNT)=X
            DO 15 I=1,NVAR
              YP(I,KOUNT)=Y(I)
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) PAUSE 'Stepsize smaller than minimum.'
        H=HNEXT
16    CONTINUE
      PAUSE 'Too many steps.'
      RETURN
      END

      SUBROUTINE NR_RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
C     implicit double precision (a-h,o-z)
      PARAMETER (NMAX=5,FCOR=.0666666667,ONE=1.,SAFETY=0.9,ERRCON=6.E-4)
      EXTERNAL DERIVS
      DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      PGROW=-0.20
      PSHRNK=-0.25
      XSAV=X
      DO 11 I=1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
11    CONTINUE
      H=HTRY
1     HH=0.5*H
      CALL NR_RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL NR_RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV)PAUSE 'Stepsize not significant in RKQC.'
      CALL NR_RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.
      DO 12 I=1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX=MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
12    CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.*H
        ENDIF
      ENDIF
      DO 13 I=1,N
        Y(I)=Y(I)+YTEMP(I)*FCOR
13    CONTINUE
      RETURN
      END

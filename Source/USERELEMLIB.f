      
*deck,UserElem     USERDISTRIB                                     jxw
c Copyright ANSYS.  All Rights Reserved.      
   
      subroutine UserElem (elId, matId, keyMtx, lumpm, nDim, nNodes,
     &                     Nodes, nIntPnts, nUsrDof, kEStress, 
     &                     keyAnsMat, keySym, nKeyOpt, KeyOpt,
     &                     temper, temperB, tRef, kTherm, 
     &                     nPress, Press, kPress, nReal, RealConst, 
     &                     nSaveVars, saveVars, xRef, xCur, 
     &                     TotValDofs, IncValDofs, ItrValDofs,
     &                     VelValDofs, AccValDofs,
     &                     kfstps, nlgeom, nrkey, outkey, elPrint, iott,
     &                     keyHisUpd, ldstep, isubst, ieqitr, timval, 
     &                     keyEleErr, keyEleCnv,
     &                     eStiff, eMass, eDamp, eSStiff,
     &                     fExt, fInt, elVol, elMass, elCG, 
     &                     nRsltBsc, RsltBsc, nRsltVar, RsltVar, 
     &                     nElEng, elEnergy)
     &      
*deck,usermat      USERDISTRIB  parallel                                gal
  
      
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERELEM"::UserElem

 
c*************************************************************************
c
c *** Primary function: General User Element Subroutine
c *** Note:
c       This routine is completed with an example, see more details later.
c
c
c     PROGRAMMER SHOULD NOT CHANGE ANY PURE INPUT ARGUMENTS (marked by ....,in)!
c
c     elId      (int,sc,in)        element number
c     matId     (int,sc,in)        material number of this element
c     keyMtx    (int,ar(10),in)    matrix and load vector form requests
c                                     0 = not requested, 1 = requested
c                                     see below for more details
c     lumpm     (int,sc,in)        mass matrix format
c                                    = 0 no lumped mass matrix
c                                    = 1 lumped mass matrix
c     nDim      (int,sc,in)        number of dimensions of the problem
c                                       (defined on USRELEM command as NDIM)
c                                    = 2 2D
c                                    = 3 3D
c     nNodes    (int,sc,in)        number of nodes of the element
c                                       (defined on USRELEM command as NNODES)
c     Nodes     (int,ar(nNodes),in)node list of this element 
c     nIntPnts  (int,sc,in)        maximum number of integration points
c                                       (defined on USRELEM command as NINTPNTS)
c     nUsrDof   (int,sc,in)        number of DOFs of this element (matrix and 
c                                     load vector size)
c     kEStress  (int,sc,in)        kEStress 
c                                       (defined on USRELEM command as KESTRESS)
c     keyAnsMat (int,sc,in)        key to indicate if ANSYS material
c                                     routine is going to be called
c                                     (defined on USRELEM command as KEYANSMAT)
c                                     = 0, No
c                                     = 1, Yes
c     keySym    (int,sc,in)        key to indicate if element matrices
c                                     is symmetric
c                                       (defined on USRELEM command as KEYSYM)
c                                     = 0, symmetric
c                                     = 1, unsymmetric
c     nKeyOpt   (int,sc,in)        number of element key options able to be
c                                     used in this routine
c     KeyOpt    (int,ar(nKeyOpt),in) values of element key option defined
c                                     by et or keyopt command for the
c                                     user elements, only the first
c                                     nKeyOpt values are passed in and can
c                                     be used to branch the routine for
c                                     different formulations
c     temper    (dp,ar(nNodes),in) nodal temperatures at current time
c     temperB   (dp,ar(nNodes),in) nodal temperatures at the beginning of this
c                                     incremental step (substep)
c     tRef      (dp,sc,in)         reference temperature
c     kTherm    (int,sc,inout)     input:  flag for thermal loading 
c                                      = 1, Temperatures at nodes are different 
c                                      from the reference temperature, 
c                                      thermal loading might be needed.
c                                      = 0, Temperatures at nodes are the same
c                                      as the reference temperature, 
c                                      thermal loading is not needed.
c                                  output:  flag for thermal strains
c     nPress    (int,sc,in)        number of pressure values for this element
c     Press     (dp,ar(nPress),in) applied elemental face load (pressure)
c     kPress    (int,sc,in)        flag for pressure loading 
c                                      = 1, pressure load is applied and 
c                                      equivalent nodal forces should be 
c                                      calculated
c                                      = 0, no pressure loading
c     nReal     (int,sc,in)        number of real constants
c                                       (defined on USRELEM command as NREAL)
c     RealConst (dp,ar(nReal),in)  user defined real constants 
c     nSaveVars (int,sc,in)        number of saved variables
c                                      (defined on USRELEM command as NSAVEVARS)
c     saveVars  (dp,ar(nSaveVars),inout) user saved variables
c     xRef      (dp,ar(nDim,nNodes),in)
c                                  nodal coordinates in initial configuration
c     xCur      (dp,ar(nDim,nNodes),in)
c                                  nodal coordinates in current configuration
c     TotValDofs (dp,ar(nUsrDof),in) total values of DOFs (displacements) 
c                                     from time = 0
c     IncValDofs (dp,ar(nUsrDof),in) incremental values of DOFs (displacements) 
c                                     for the current step
c     ItrValDofs (dp,ar(nUsrDof),in) iterative values of DOFs (displacements)
c                                     for the current iteration
c                                     (normally needed for debug only)
c     VelValDofs (dp,ar(nUsrDof),in) first time derivatives of DOFs 
c                                             (velocities) (normally not needed)
c     AccValDofs (dp,ar(nUsrDof),in) second time derivatives of DOFs 
c                                          (accelerations) (normally not needed)
c     kfstps    (int,sc,in)        key for the first iteration of first 
c                                     substep of the first load step
c                                     = 1 yes
c                                     = 0 no
c     nlgeom    (int,sc,in)        large deformation key [from nlgeom command]
c                                     = 0 NLGEOM,OFF
c                                     = 1 NLGEOM, ON
c     nrkey     (int,sc,in)        key to indicate a newton-raphson
c                                     (incremental) procedure
c                                     = 0 No
c                                     = 1 Yes
c     outkey    (int,sc,in)        key to indicate if any element output is
c                                     to be placed on the print file or the 
c                                     result file
c                                     = 0 No
c                                     = 1 Yes
c     elPrint   (int,sc,in)        key to indicate if any element output is 
c                                     to be placed on the print file
c                                     = 0 No
c                                     = 1 Yes
c     iott      (int,sc,in)        print output file unit number
c     keyHisUpd (int,sc,in)        key to indicate if history-dependent
c                                    variables need to be updated, like
c                                    equivalent plastic strain, back stress
c                                    etc. since the iteration is already
c                                    converged
c                                     = 0 not converged, don't need to update
c                                         history dependent variables
c                                     = 1 yes, converged, need to update
c                                         history dependent variables
c
c     --- The following 7 variable group can usually be ignored.
c     --- The variables are used for debug, timing, and convergence control.
c     ldstep    (int,sc,in)        current load step number
c     isubst    (int,sc,in)        current substep number
c     ieqitr    (int,sc,in)        current equilibium iteration  number
c     timval    (int,sc,in)        current time value
c     keyEleErr (int,sc,inout)     key to indicate if there is any element 
c                                     formulation error, like negative Jacobian.
c                                     The error could be caused by too
c                                     large incremental step, illegal model.
c                                     = 0 no error (preset value before calling)
c                                     = 1 some error happens. ANSYS will
c                                     decide to stop the analysis or cutback
c                                     the substep (bi-section) based on other
c                                     user input and information at higher
c                                     level.
c     keyEleCnv (int,sc,inout)     key to flag if this element satisfies
c                                     the user defined element convergence
c                                     criterion. 
c                                     = 1, yes, the criterion is satisfied
c                                       or don't have any criterion at all
c                                       it is preset value before calling
c                                     = 0, no, the element doesn't satisfy
c                                       element convergence criterion. If
c                                       this is the case, the iteration will
c                                       not converge even when both force
c                                       and displacement converge 
c       ---- end of 7 variable group -----
c
c                                                                  requested if
c     eStiff(dp,ar(nUsrDof,nUsrDof),inout) stiffness matrix         keyMtx(1)=1
c     eMass (dp,ar(nUsrDof,nUsrDof),inout) mass matrix              keyMtx(2)=1
c     eDamp (dp,ar(nUsrDof,nUsrDof),inout) damping matrix           keyMtx(3)=1
c     eSStiff(dp,ar(nUsrDof,nUsrDof),inout)stress stiffness matrix  keyMtx(4)=1
c     fExt      (dp,ar(nUsrDof),out)       applied f vector         keyMtx(5)=1
c     fInt      (dp,ar(nUsrDof),out)       internal force vector    keyMtx(6)=1

c     elVol     (dp,sc,out)        element volume
c     elMass    (dp,sc,out)        element mass
c     elCG      (dp,ar(3),out)     element centroid coordinates in current
c                                     configuration
c     nRsltBsc  (dp,sc,in)         number of basic elemental results saved in
c                                   result files 
c     RsltBsc   (dp,ar(nRsltBsc),out) basic elemental results 
c                                       (see EXPLANATION below)
c     nRsltVar  (int,sc,in)        number of elemental results saved in 
c                                     result file as non-summable miscellaneous
c                                     variables 
c                                       (defined on USRELEM command as NRSLTVAR)
c     RsltVar   (dp,ar(nRsltVar),out) variables to saved in result files as
c                                      non-summable miscellaneous variables 
c                                      requested when outkey = 1
c
c     nElEng    (int,sc,in)        number of energies (fixed at 3)
c
c     elEnergy  (dp,ar(nElEng),out) elemental energy
c                                     elEnergy(1), element strain energy
c                                     elEnergy(2), element plastic strain energy
c                                     elEnergy(3), element creep strain energy
c
c     EXPLANATION OF RsltBsc
c     
c       Basic element results are saved and total number of result 
c     quantities is nRsltBsc, where:
c            nRsltBsc = (7+7)* number of corner nodes at one element.
c       To process the quantites by post processing properly, the results 
c     must be in the following order:
c       1.) Stresses: Sx Sy Sz Sxy Syz Sxz Seqv at all corner points,
c     followed by:
c       2.) Strains : Ex Ey Ez Exy Eyz Exz Eeqv at all corner points
c     where Seqv and Eeqv = equivalent stress and strain respectively
c
c
************************************************************************
c

#include "impcom.inc"
c
      EXTERNAL         ElemGetMat  ! this routine may be user-programmed

      INTEGER          elId, matId, keyMtx(10), lumpm,nDim, nNodes,
     &                 Nodes(nNodes), nIntPnts, nUsrDof, kEStress, 
     &                 keyAnsMat, keySym, nKeyOpt, KeyOpt(nKeyOpt),
     *                 kTherm, nPress, kPress, nReal, nSaveVars, 
     &                 kfstps, nlgeom, nrkey, outkey, jdim, 
     &                 elPrint, iott, keyHisUpd, l, inod,
     &                 ldstep, isubst, ieqitr, keyEleErr, keyEleCnv,
     &                 nRsltBsc, nRsltVar, nElEng, gggggg, intPnttt


      DOUBLE PRECISION temper(nNodes), temperB(nNodes), tRef, 
     &                 Press(nPress), RealConst(nReal),
     &                 saveVars(nSaveVars), klk(30),
     &                 xRef(nDim,nNodes), xCur(nDim,nNodes),
     &                 TotValDofs(nUsrDof), IncValDofs(nUsrDof), 
     &                 ItrValDofs(nUsrDof), VelValDofs(nUsrDof),
     &                 AccValDofs(nUsrDof),      timval,
     &                 eStiff(nUsrDof,nUsrDof), eMass(nUsrDof,nUsrDof), 
     &                 eDamp(nUsrDof,nUsrDof), eSStiff(nUsrDof,nUsrDof), 
     &                 fExt(nUsrDof), fInt(nUsrDof), 
     &                 elVol, elMass, elCG(3),
     &                 RsltBsc(nRsltBsc), RsltVar(nRsltVar), 
     &                 elEnergy(nElEng)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c *** CODE EXAMPLE ***
c
c --- The element code is only to show how to use the routine to create user 
c     elements.  Two element types are coded.   Only the stiffness matrix, mass 
c     matrix and internal load vector are shown.
c
c       When KeyOpt(1) = 0, it is a structural 2D plane strain element 
c                                          with 4 nodes and 4 integration points
c       When KeyOpt(1) = 1, it is a structural 3D solid elements 
c                                         with 20 nodes and 8 integration points
c       No advanced element technology is employed,
c                         and they are only coded for geometric linear analysis.
c
c --- This demonstration code only shows how to create eStiff, eMass, and fInt.
c --- Other matrices and/or vectors should be created similarly. 
c --- This is coded for good readability.
c
c --- Included decks, functions, variables defined for the example

#include "locknm.inc"
      EXTERNAL         vzero, vmove, vmult, vdot, vidot,
     &                 maxv, matxb, matba, maat, matsym, getMatProp,
     &                 erhandler, equivStrain, ElemJac, ElemMass,
     &                 ElemRsltNode, ElemShpFn, pplock, ppunlock

      DOUBLE PRECISION vdot, vidot

      INTEGER          nUsrDof2, intPnt, iNode, nTens, flgSingular,
     &                 k1, k2, k3, nComp, iDim, iDim1, iComp,numU,
     &                 nNodesCorner, nDirect, kThermIP, i, j, kFlagF,
     &                 kFlagS, typemodel, linear, drivingforce,
     &                 func_name, LSTR, keyHmin   
      DOUBLE PRECISION BMat(nDim*2,nUsrDof), Ex, nu, density, G, workDb,
     &                 con1, con2, cMat(nDim*2,nDim*2), shIsoC(nNodes),
     &                 shIso(nNodes), shDerIso(nDim,nNodes), wtIP(1),
     &                 workArr(360), elJac(nDim*nDim), detJac, dperr(2),
     &                 shDerEl(nDim,nNodes), dVol, Strain(nDim*2), 
     &                 Stress(nDim*2), wStrain(48), wStress(48),
     &                 nStrain(48), nStress(48), sigm, tem, prop(3),
     &                 IncStrain(nDim*2),  defG(3,3), 
     &                 defG0(3,3), xCurIP(nDim), TemperIP, 
     &                 TemperIPB, StressTh(nDim*2), MatProp(5),
     &                 StrainPl(nDim*2), StrainCr(nDim*2), 
     &                 StrainTh(nDim*2), StrainSw, StressBk(nDim*2),
     &                 MatRotGlb(3,3), wStrainTh(48), wStrainPl(48),
     &                 wStrainCr(48), eMassb(nNodes,nNodes), EnergyD(3),
     &                 phik(nNodes), phi, dNdx(nDim,nNodes) , phin, H,               !
     &                 du(nDim*nNodes), Gc, xlc, xk, dstran(nDim*2,1),             !
     &                 oldStrain(nDim*2), oldStress(nDim*2), Hn, Psi,                    !     
     &                 Stress2(nDim*2,1), dN(nNodes,1), w, dw, ddw, cw,                             !
     &                 b(nDim*2,nUsrDof-nNodes), rhs(nUsrdof), trE,                      !
     &                 amatrx(nUsrDof,nUsrDof),eStiff1(nUsrDof,nUsrDof),             !   
     &                 u(nUsrDof-nNodes), pl, pln, plast, Psi1, sedd,
     &                 alph, alphn, alphBn, Fdeg, Ac, alphT, seddn,                ! 
     &                 alphB, coordx, pi, bulk, Hmin, e, d, gNum, dgNum, 
     &                 ddgNum, gDen, dgDen, ddgDen, a, StrainEl(6), 
     &                 gf, dgf, sigc, ddgf, mm, EdevS, Edev(3), phinn,
     &                 psip, psin, eg, trEp2, trEn2, trEp, trEn, ggf,
     &                 STR(3), Stress3(6), psipPl, psipn
      
      
      CHARACTER*4      label(3)

c --- temporary debug key
      INTEGER debug, ix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c --- B E G I N   E X E C U T A B L E  C O D I N G
c
c --- initialization
      xlc= RealConst(1)
      Gc=  RealConst(2)
      kFlagF = RealConst(3)
      typemodel = RealConst(4)
      drivingforce = RealConst(5)
      linear = RealConst(6)
      kFlagS = RealConst(7)
      sigc = RealConst(8)
      keyHmin = RealConst(9)
      xk = 1e-07
      Ac = 0.d0
      nTens = nDim*2
      nComp = nDim*nDim
      nDirect = 3     
      amatrx = 0.d0
      rhs = 0.d0
      numU = nUsrDof/nNodes
      pi = 3.14159
      
      
      IF (numU.EQ.3.d0) THEN
          KeyOpt(1) = 0.0d0
      ELSE
          KeyOpt(1) = 1.0d0
      END IF    
      
      if (kFlagF.eq.1) then
        alphT=Gc/(2.d0*6.d0*xlc)
      else
        alphT=1.d10
      endif   
      
      keyMtx = 0.d0
      keyMtx(1) = 1.d0
      keyMtx(6) = 1.d0
      nUsrDof2 = nUsrDof*nUsrDof
      CALL vzero (BMat(1,1),nUsrDof*nTens)
      IF (keyMtx(1).EQ.1) CALL vzero (eStiff(1,1),nUsrDof2)
      IF (keyMtx(2).EQ.1) CALL vzero (eMass(1,1) ,nUsrDof2)
      IF (keyMtx(5).EQ.1) CALL vzero (fExt(1)    ,nUsrDof)
      IF (keyMtx(6).EQ.1) CALL vzero (fInt(1)    ,nUsrDof)

*!*!*!*!*!*!*!*!*!*!*!*****start integration loop*!*!*!*!*!*!*!*!*!*!*!*!

      DO intPnt = 1, nIntPnts

c --- obtain shape functions and derivatives of shape functions

         IF (KeyOpt(1).EQ.0) THEN
            CALL ElemShpFn (1, intPnt, 1, shIso(1), nNodes)
            CALL ElemShpFn (1, intPnt, 2, shDerIso(1,1), nUsrDof)
            CALL ElemShpFn (1, intPnt, 3, wtIP(1), 1)
         ELSE
            CALL ElemShpFn (2, intPnt, 1, shIso(1), nNodes)
            CALL ElemShpFn (2, intPnt, 2, shDerIso(1,1), nUsrDof)
            CALL ElemShpFn (2, intPnt, 3, wtIP(1), 1)
         END IF
         
         dN = 0.d0
         do i=1,nNodes
           dN(i,1) = shIso(i)
         end do  

c --- coordinates at integration points

         DO iDim = 1, nDim
            xCurIP(iDim) =  vidot(shIso(1), 1, xCur(iDim,1), nDim,
     &                 nNodes)
         END DO

c --- temperatures at integration points

         TemperIP = vdot(shIso(1), temper(1), nNodes)
         TemperIPB = vdot(shIso(1), temperB(1), nNodes)
         

c******derivatives of shape functions****************************************
         CALL vzero (workArr(1), nComp)
         iComp = 1
         DO iDim = 1, nDim
            DO iDim1 = 1, nDim
               DO iNode = 1, nNodes
                  workArr(iComp) = workArr(iComp) 
     &                + shDerIso(iDim1,iNode)*xCur(iDim,iNode)
               END DO
               iComp = iComp + 1
            END DO
         END DO
         CALL ElemJac (workArr(1), elJac(1), nDim, detJac, 
     &                   flgSingular)
         IF (flgSingular.LE.0) THEN
            dperr(1) = detJac
            dperr(2) = elId
            CALL erhandler ('UserElem', 1100, 3, 'Negative element
     &                       Jacobian value %I at element %I. This
     &                       is due to wrong element order or bad  
     &                       mesh.',dperr(1),' ') 
            GOTO 990
         END IF
         DO iNode = 1, nNodes
            IF (KeyOpt(1).EQ.0) THEN
               shDerEl(1,iNode) =  elJac(1)*shDerIso(1,iNode)
     &                           + elJac(3)*shDerIso(2,iNode)
               shDerEl(2,iNode) =  elJac(2)*shDerIso(1,iNode)
     &                           + elJac(4)*shDerIso(2,iNode)
            ELSE
               shDerEl(1,iNode) =  elJac(1)*shDerIso(1,iNode)
     &                           + elJac(4)*shDerIso(2,iNode)
     &                           + elJac(7)*shDerIso(3,iNode)
               shDerEl(2,iNode) =  elJac(2)*shDerIso(1,iNode)
     &                           + elJac(5)*shDerIso(2,iNode)
     &                           + elJac(8)*shDerIso(3,iNode)
               shDerEl(3,iNode) =  elJac(3)*shDerIso(1,iNode)
     &                           + elJac(6)*shDerIso(2,iNode)
     &                           + elJac(9)*shDerIso(3,iNode)
            END IF
         END DO
         
         dNdx = 0.d0
         do i=1,nDim
             do j=1,nNodes
                 dNdx(i,j) = shDerEl(i,j) 
             end do
         end do    
         dVol = detJac*wtIP(1)
c************************************************************************
         
c******creatr B-matrix*************************************************** 
          b = 0.d0
          k1 = 1
          DO iNode = 1,nNodes
            k2  = k1 + 1
            b(1,k1) =  shDerEl(1,iNode)
            b(2,k2) =  shDerEl(2,iNode)
            b(4,k1) =  shDerEl(2,iNode)
            b(4,k2) =  shDerEl(1,iNode)
            IF (KeyOpt(1).EQ.1) THEN
               k3 = k2 + 1
               b(3,k3) =  shDerEl(3,iNode)
               b(5,k2) =  shDerEl(3,iNode)
               b(5,k3) =  shDerEl(2,iNode)
               b(6,k3) =  shDerEl(1,iNode)
               b(6,k1) =  shDerEl(3,iNode)
            END IF         
            k1 = k1 + nDim
          END DO
c************************************************************************
          
c******calculate damage parameter**************************************** 
          
          DO iNode = 1,nNodes
              phik(iNode) = TotValDofs(numU*iNode)
          end do    
          phi = 0.d0
          do inod=1,nNodes
              phi=phi+shIso(inod)*phik(inod)
          end do
          if (phi.gt.1.d0) then
              phi=1.d0
          end if
          if (TotValDofs(2).ne.0.d0) then
              cMat(4,4) = G
          end if
c*************************************************************************           

c******reading data from the previous step********************************          
          do i=1,ntens
              oldStress(i) = saveVars(20*(intPnt-1) + i)
              oldStrain(i) = saveVars(20*(intPnt-1) + i + 6)
          end do    
          phin = saveVars(20*(intPnt-1) + 13)
          Hn = saveVars(20*(intPnt-1) + 14)
          alphn = saveVars(20*(intPnt-1) + 15)
          alphBn = saveVars(20*(intPnt-1) + 16) 
          psip = saveVars(20*(intPnt-1) + 18)
          psipn = psip 
          phinn = phin
c*************************************************************************  
       
c******calculate crack length*********************************************           
          coordx=0.d0
          if (phi.ge.0.95d0) then
            do i=1,nNodes
                coordx=coordx+dN(i,1)*xCur(1,i)
            enddo
            if (coordx.gt.Ac) then
                Ac=coordx
            endif
          endif 
c*************************************************************************            

c******calculate strains and stress***************************************
         u = 0.d0
         do i=1,nUsrDof
             IF (MOD(i,numU).EQ.0.d0) THEN
                 cMat(4,4) = G
             ELSE 
                 du(i-INT(i/numU)) = ItrValDofs(i)
                 u(i-INT(i/numU)) = TotValDofs(i)
             END IF    
         end do 


         CALL maxv (b(1,1), u(1), IncStrain(1), nTens, 
     &              nUsrDof-nNodes)
          
         nlgeom = 0
         IF (nlgeom.EQ.0) THEN
             DO iDim = 1, 3
                DO iDim1 = 1, 3
                    defG0(iDim, iDim1) = 0.0D0
                END DO
                defG0(iDim, iDim) = 1.0D0
             END DO
             CALL vmove (defG0(1,1),defG(1,1),9)
         ELSE
             CALL straininc(ntens,ndim,nNodes,nUsrDof,b,u,defG0)
             CALL vmove (defG0(1,1),defG(1,1),9)
         END IF
         
         IF (keyAnsMat.EQ.1) THEN
c           ---- Use standard ANSYS material (METHOD 1)
c          ---- USERMAT is called from here
            CALL ElemGetMat (elId, matId, nDim, nTens, nDirect,
     &                         intPnt, xCurIP(1), TemperIP,
     &                         TemperIPB, kThermIP, IncStrain(1),
     &                         defG0(1,1), defG(1,1),
     &                         cMat(1,1), MatProp(1), Stress(1), 
     &                         Strain(1), StressTh(1), StrainTh(1),
     &                         StrainPl(1), StrainCr(1), 
     &                         StressBk(1), StrainSw, EnergyD(1),
     &                         MatRotGlb(1,1))
            if (kThermIP .eq. 1) kTherm = 1

            nu = MatProp(5)            
            Ex = MatProp(1)*2.d0*(1.d0+nu)
            bulk = Ex/3.d0/(1.d0-2.d0*nu)
         ELSE
c          ---- Make up your own material (METHOD 2)
            IF (nlgeom.EQ.0) CALL vmove (IncStrain(1), Strain(1),
     &                                   nTens)
            Stress(1) = con1*Strain(1) + con2*(Strain(2)+Strain(3))
            Stress(2) = con1*Strain(2) + con2*(Strain(3)+Strain(1))
            Stress(3) = con1*Strain(3) + con2*(Strain(1)+Strain(2))
            Stress(4) = G*Strain(4)
            IF (KeyOpt(1).EQ.1) THEN
               Stress(5) = G*Strain(5)
               Stress(6) = G*Strain(6)
            END IF
         END IF 

c*************************************************************************
         
******Choice type model***************************************************
         if (typemodel .EQ. 1) then
             gf = (1.d0-phi)**2 + xk
             ggf = (1.d0-phinn)**2 + xk
             dgf = 2*(phi -1)
             ddgf = 2.d0         
             w = phi**2
             dw = 2*phi
             ddw = 2.d0
             cw = 0.5
             if (keyHmin.eq.1) then
                 Hmin = 3.d0*Gc/(16.d0*xlc)
             else
                 Hmin = 0.d0
             endif   
         else if (typemodel .EQ. 2) then
             gf = (1.d0-phi)**2 + xk
             ggf = (1.d0-phinn)**2 + xk
             dgf = 2*(phi -1)
             ddgf = 2.d0
             w = phi
             if (phi.EQ.0) then
                 dw = 0.d0
             else    
                 dw = 1.d0
             end if    
             ddw = 0.d0
             cw = 2.d0/3.d0
             Hmin = 0.d0 
         else
            mm =4.d0*Ex*Gc/(xlc*sigc**2) 
            w =2.d0*phi-phi**2
            dw = 2.d0 - 2.d0*phi
            ddw = -2.d0
            cw = pi/4.d0
            if (linear .EQ. 1) then
                e = -0.5d0
                d = 2.d0
            else
                e = 2.d0**(5.d0/3.d0)-3.d0
                d = 2.5d0
            end if
            gNum=(1.d0-phi)**d
            dgNum=-d*(1.d0-phi)**(d-1.d0);
            ddgNum=d*(d-1.d0)*(1.d0-phi)**(d-2.d0)
            gDen=gNum+mm*phi+mm*e*phi**2.d0
            dgDen=dgNum+mm+2.d0*mm*e*phi
            ddgDen=ddgNum+2.d0*mm*e

            gf=gNum/gDen + xk
            ggf=(1.d0-phinn)**d/(gNum+mm*phinn+mm*e*phinn**2.d0) + xk
            dgf=(dgNum*gDen-gNum*dgDen)/(gDen**2.d0)
            ddgf=((ddgNum*gDen-gNum*ddgDen)*gDen-2.d0*
     1 (dgNum*gDen-gNum*dgDen)*dgDen)/(gDen**3.d0)
            Hmin = 0.5d0*sigc**2/Ex           
         end if
c*************************************************************************
         
c********calculate driving force fracture*********************************
         psin = 0.d0
         StrainEl(1:4) = Strain-StrainPl
         IF (numU.EQ.3.d0) THEN   
             StrainEl(5) = 0.d0
             StrainEl(6) = 0.d0
         ELSE
             StrainEl(5) = Strain(5)-StrainPl(5)
             StrainEl(6) = Strain(6)-StrainPl(6)             
         END IF
         CALL SPRIND (StrainEl, STR, 2)
         trE=STR(1)+STR(2)+STR(3)
         trEp=0.5d0*(trE+abs(trE))
         trEn=0.5d0*(trE-abs(trE))
         Edev(1:3)=STR(1:3)-trE/3.d0
         EdevS=Edev(1)**2+Edev(2)**2+Edev(3)**2 
         eg = Ex/(1.d0+nu)/2.d0
         if (drivingforce.eq.2) then ! Amor et al.
             psip = 0.d0
             psip = 0.5d0*bulk*trEp**2+eg*EdevS
             psin = 0.5d0*bulk*trEn**2
             psipPl = EnergyD(2)
         elseif (drivingforce.eq.3) then ! Miehe et al.
             trEp2 = 0.d0
             trEn2 = 0.d0
             do i=1,3
                 trEp2 = trEp2+(STR(i)+abs(STR(i)))**2.d0/4.d0
                 trEn2 = trEn2+(STR(i)-abs(STR(i)))**2.d0/4.d0
             end do
             psip = 0.d0
             psip = nu*eg/(1d0-2d0*nu)*trEp**2d0+eg*trEp2
             psin = nu*eg/(1d0-2d0*nu)*trEn**2d0+eg*trEn2
             psipPl = EnergyD(2)
         elseif  (drivingforce.eq.1) then! no split
              psip = EnergyD(1)
              psipPl = EnergyD(2)
         elseif  (drivingforce.eq.4) then! maxstress
            Stress3(1:4) = Stress
            IF (numU.EQ.3.d0) THEN   
                Stress3(5) = 0.d0
                Stress3(6) = 0.d0
            ELSE
                Stress3(5) = Stress(5)
                Stress3(6) = Stress(6)             
            END IF
            CALL SPRIND (Stress3, STR, 1) 
            psip = 0.5d0*max(STR(1),STR(2),STR(3))**2/Ex
            psipPl = EnergyD(2)
         else ! no split
            do i =1, ntens 
                psip = psip + (oldStress(i) + 
     1          Stress(i))*(Strain(i) - oldStrain(i))*0.5d0
            end do
            psipPl = 0.d0
         endif            

         if(kflagS.eq.0) then
             H=max(Hmin, Hn, psip, psipn)
         else
             H=max(Hmin, psipn, Hn)
         end if      
         
c*************************************************************************
          
c********calculate fatigue degradation function***************************
         alph=(psip+psipPl)*gf

         if (alph.ge.alphn) then
             alphB = alphBn+abs(alph-alphn)
         else
             alphB=alphBn
         endif            
          
          if (alphB.lt.alphT) then
              Fdeg= 1.d0
          else
              Fdeg=(2.d0*alphT/(alphB+alphT))**2.d0
          endif      
c*************************************************************************
          
c*********data recording**************************************************
          DO i=1,ntens
              saveVars(20*(intPnt-1) + i) = Stress(i)
              saveVars(20*(intPnt-1) + i + 6) = Strain(i)
          END DO    
          saveVars(20*(intPnt-1) + 13) = phi
          saveVars(20*(intPnt-1) + 14) = H    
          saveVars(20*(intPnt-1) + 15) = alph 
          saveVars(20*(intPnt-1) + 16) = alphB
          saveVars(20*(intPnt-1) + 17) = Ac
          saveVars(20*(intPnt-1) + 18) = psip
          saveVars(20*(intPnt-1) + 19) = H + EnergyD(2) 
c*************************************************************************
          
c*********form and assemble stiffness matrix and internal force vector****         
           IF (keyMtx(1).EQ.1) THEN
              amatrx(1:ndim*nNodes,1:ndim*nNodes)=
     1      amatrx(1:ndim*nNodes,1:ndim*nNodes)+
     1      dvol*(ggf*
     1      matmul(matmul(transpose(b),cMat),b)) 
            
          
              amatrx((ndim*nNodes+1):nUsrDof,(ndim*nNodes+1):nUsrDof)=
     1  amatrx((ndim*nNodes+1):nUsrDof,(ndim*nNodes+1):nUsrDof)
     1 +dvol*(matmul(transpose(dNdx),dNdx)*Gc*xlc*Fdeg/2.d0/cw
     2 +matmul(dN,transpose(dN))*(Gc/xlc*Fdeg/4.d0/cw*ddw
     2  +ddgf*(H+psipPl)))
              
          END IF  
          
          IF (keyMtx(6).EQ.1) THEN 
              rhs(1:ndim*nNodes)=rhs(1:ndim*nNodes)+
     1 dvol*(ggf*matmul(transpose(b),Stress)) 
              
              rhs((ndim*nNodes+1):nUsrDof)=
     1 +rhs((ndim*nNodes+1):nUsrDof)
     1 +dvol*(matmul(transpose(dNdx),matmul(dNdx,phik(1:nNodes)))*Fdeg
     2 *Gc*xlc/2.d0/cw+dN(1:nNodes,1)*(Gc/xlc*Fdeg/4.d0/cw*dw
     2  +dgf*(H+psipPl)))
            
          END IF  
      END DO
*!*!*!*!*!*!*!*!*!*!*!*!*!*end integration loop*!*!*!*!*!*!*!*!*!*!*!*!
         
          do i=1,nUsrDof
               IF (MOD(i,numU).EQ.0.d0) THEN
                  fint(i) = rhs(ndim*nNodes + INT(i/numU))
              ELSE 
                  fint(i) = rhs(i-INT(i/numU))
              END IF    
          end do
          
          do i=1,nUsrDof
               IF (MOD(i,numU).EQ.0.d0) THEN
                   do j=1,nUsrDof
                       IF (MOD(j,numU).EQ.0.d0) THEN
                          eStiff(i,j) = 
     &    amatrx(ndim*nNodes + INT(i/numU),ndim*nNodes + INT(j/numU))
                       ELSE
                          eStiff(i,j) = 0.d0
                       END IF    
                   end do    
              ELSE 
                   do j=1,nUsrDof
                       IF (MOD(j,numU).EQ.0.d0) THEN
                          eStiff(i,j) = 0.d0 
                       ELSE
                          eStiff(i,j) = 
     &                        amatrx(i-INT(i/numU),j-INT(j/numU))     
                       END IF    
                   end do   
              END IF    
          end do
c*************************************************************************
          
 990  CONTINUE
      RsltVar = saveVars 
      RETURN
      END
            
      SUBROUTINE SPRIND(S, PS, LSTR)

      INTEGER LSTR, NDI, NSHR
      DOUBLE PRECISION S(6), PS(3), A(3), B(3), C(3), V(3), AN(3,3)
      DOUBLE PRECISION I1, I2, I3, J1, J2, J3, R, T, Q, ALP, SCALC,
     1     PRINC1, PRINC2, PRINC3, L(3), M(3), N(3)
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0, TWOSEVEN = 27.D0,
     1     THIRD = ONE/THREE, TWENTYSEVENTH= ONE/TWOSEVEN,
     2     TWOTWENTYSEVENTH = TWO/TWOSEVEN, ONETWENTYDEG=2.094395102D0,
     3     TWOFORTYDEG=4.188790205D0)
C

C======================================================================+
C-----------
C  INPUT :
C-----------
C  S     	: STRESS OR STRAIN TENSOR
C  LSTR 	: FLAG DETERMINING STRESS OR STRAIN CALCULATION
C  NDI  	: NUMBER OF DIRECT STRESS/STRAIN COMPONENTS
C  NSHR         : NUMBER OF SHEAR COMPONENTS
C-----------
C  OUTPUT :
C-----------
C  PS(I), I=1,2,3       : THE THREE PRINCIPAL VALUES
C  AN(K1,I), I=1,2,3	: THE DIRECTION COSINES OF THE PRINCIPAL DIRECTIONS CORRESPONDING TO PS(K1)
C----------------------------------------------------------------------+
C=======================================================================
C
C     Calculate stress or strain invariants based on LSTR value
      IF (LSTR.EQ.1) THEN
         I1 = S(1)+S(2)+S(3)
         I2 = (S(1)*S(2))+(S(2)*S(3))+(S(1)*S(3))
     1        -(S(4)**2)-(S(5)**2)-(S(6)**2)
         I3 = (S(1)*S(2)*S(3))+(2*S(4)*S(5)*S(6))-(S(1)*S(5)**2)
     1        -(S(2)*S(6)**2)-(S(3)*S(4)**2)
         R = (THIRD*I1**2)-I2
         T = SQRT(TWENTYSEVENTH*R**3)
         Q = (THIRD*I1*I2)-I3-(TWOTWENTYSEVENTH*I1**3)
         IF (Q.EQ.0.d0 .and. T.EQ.0.d0) THEN
             ARG = 0.d0
         ELSE    
             ARG = -Q/(TWO*T)
         ENDIF          
         ARG = -Q/(TWO*T)
         IF (ARG.GT.1) THEN
            ARG = ARG - 1.E-10
         END IF
         ALP = ACOS(ARG)
         SCALC = SQRT(THIRD*R)
         PRINC1 = (2*SCALC*COS(ALP/THREE))+(THIRD*I1)
         PRINC2 = (2*SCALC*COS((ALP/THREE)+TWOFORTYDEG))+(THIRD*I1)
         PRINC3 = (2*SCALC*COS((ALP/THREE)+ONETWENTYDEG))+(THIRD*I1)
      ELSE
         J1 = S(1)+S(2)+S(3)
         J2 = S(1)*S(2)+S(2)*S(3)+S(1)*S(3)
     1        -(S(4)**2)-(S(5)**2)-(S(6)**2)
         J3 = S(1)*S(2)*S(3)+(2*S(4)*S(5)*S(6))-(S(1)*S(5)**2)
     1        -(S(2)*S(6)**2)-(S(3)*S(4)**2)
         R = (THIRD*J1**2)-J2
         T = SQRT(TWENTYSEVENTH*R**3)
         Q = (THIRD*J1*J2)-J3-(TWOTWENTYSEVENTH*J1**3)
         IF (Q.EQ.0.d0 .and. T.EQ.0.d0) THEN
             ARG = 0.d0
         ELSE    
             ARG = -Q/(TWO*T)
         ENDIF    
         IF (ARG.GT.1) THEN
            ARG = ARG - 1.E-10
         END IF
         ALP = ACOS(ARG)
         SCALC = SQRT(THIRD*R)
         PRINC1 = (2*SCALC*COS(ALP/THREE))+(THIRD*J1)
         PRINC2 = (2*SCALC*COS((ALP/THREE)+TWOFORTYDEG))+(THIRD*J1)
         PRINC3 = (2*SCALC*COS((ALP/THREE)+ONETWENTYDEG))+(THIRD*J1)
      END IF
C     Assign Principal Stress/Strains values to array
C
      PS(1) = PRINC1
      PS(2) = PRINC2
      PS(3) = PRINC3
C
C     Calculate cofactors and factor
C
      DO K=1, 3
         A(K)=((S(2)-S(K))*(S(3)-S(K)))-(S(5)**2)
         B(K)=-(S(4)*(S(3)-S(K)))-(S(5)*S(6))
         C(K)=(S(4)*S(5))-((S(2)-S(K))*S(6))
      END DO
C      PRINT *, A
C      PRINT *, B
C      PRINT *, C
      DO K=1, 3
         V(K)=1/SQRT(A(K)**2+B(K)**2+C(K)**2)
      END DO
C      PRINT *, V
C
C     Calculate Direction Cosines
C
      DO K=1, 3
         L(K)=A(K)*V(K)
         M(K)=B(K)*V(K)
         N(K)=C(K)*V(K)
      END DO
C      PRINT *, L
C      PRINT *, M
C      PRINT *, N
C
C     Assign Direction Cosines to array locations
C
      AN(1,1)=L(1)
      AN(1,2)=M(1)
      AN(1,3)=N(1)
      AN(2,1)=L(3)
      AN(2,2)=M(3)
      AN(2,3)=N(3)
      AN(3,1)=L(2)
      AN(3,2)=M(2)
      AN(3,3)=N(2)
C      PRINT *, AN
      RETURN
      END
      
       subroutine straininc(ntens,ndim,nNodes,nUsrDof,bmat,utmp,xx1)
c
c     Notation:  
c       dstran(i)  incremental strain component 
c       note:      i = 1   xx direction 
c                    = 2   yy direction 
c                    = 3   zz direction
c                    = 4   xy direction
c    u() - displacement
c   
#include "impcom.inc"

      DOUBLE PRECISION bmat(ndim,nnodes),
     & xx1(3,3),utmp(nUsrDof-nNodes), utt(ndim)
      
      INTEGER k1,i, ntens,ndim,nNodes,nUsrDof, nodi, incr_row
      
      DOUBLE PRECISION dNidx, dNidy

      ! set xx1 to Identity matrix
      xx1=0.d0
      do k1=1,3
       xx1(k1,k1)=1.d0       
      end do

c************************************
c    Compute incremental strains
c************************************
c
      do nodi=1,nNodes
       dNidx=bmat(1,nodi)
       dNidy=bmat(2,nodi)
       incr_row=(nodi-1)*ndim
       do i=1,ndim
        utt(i)=utmp(i+incr_row)
       end do       
c     deformation gradient
       xx1(1,1)=xx1(1,1)+dNidx*utt(1)
       xx1(1,2)=xx1(1,2)+dNidy*utt(1)
       xx1(2,1)=xx1(2,1)+dNidx*utt(2)
       xx1(2,2)=xx1(2,2)+dNidy*utt(2)
c
      end do          



!! function to calculate the row index of the nonzero elements (oidx,midx) given a column index (gidx)
!! for a given graph matrix S (fixed in this case and known)
SUBROUTINE rowidx(gidx,oidx,midx,nvert)
    IMPLICIT NONE

    INTEGER:: gidx, nvert, blknbr, idxcol, idxrow, oidx, midx
    
    ! calculate row indices of non zero elements in S matrix (non-symmetric case)
    blknbr = CEILING(REAL(gidx) / REAL((nvert - 1)))
    idxcol = gidx - (blknbr-1) * (nvert-1)
    idxrow  = idxcol + (blknbr-1) * nvert
    IF(idxcol >= blknbr) idxrow = idxrow + 1
    oidx = (blknbr-1) * nvert + blknbr              ! index for 1
    midx = idxrow                                   ! index for -1

    RETURN 
END SUBROUTINE rowidx



!! function to calculate the column index (cidx1,cidx2) for the non-symmetric S from the column index (cidx) of symmetric S 
SUBROUTINE colidxSymToNonsym(cidx,cidx1,cidx2,nvert)
    IMPLICIT NONE

    INTEGER:: cidx, cidx1, cidx2, nvert, tmp1, tmp2
    
    tmp1 = CEILING(REAL((2*nvert-1)-SQRT(REAL((2*nvert-1)**2-8*cidx)))/REAL(2))
    tmp2 = cidx - ((tmp1 - 1) * (2 * nvert - tmp1)) / 2
    
    cidx1 = (tmp1-1) + tmp2 + (tmp1 - 1) * (nvert - 1)
    cidx2 = (tmp1 + tmp2 - 1) * (nvert - 1) + tmp1

    RETURN 
END SUBROUTINE colidxSymToNonsym


!!  --------------------- non-symmetric

!! lasso for graphical input matrix 
SUBROUTINE admmGraph(y, coeffs, res, postv, normp, nobs, nvars, nvert, &
                        pen, rho, maxinit, epsin, errin, maxoutit, epsout, errout, errcode)
    IMPLICIT NONE
    
    
    ! data specific params
    INTEGER:: nobs, nvars

    ! penalty specific params
    INTEGER::gidx, postv
    DOUBLE PRECISION:: pen, rho, normp

    ! lasso loop specific params 
    INTEGER:: jidx, iter, nnzgr, maxinit, errcode
    INTEGER:: oidx, midx, nvert
    DOUBLE PRECISION:: epsin
    DOUBLE PRECISION:: inter, uabs, tem, errin, diffloc, bloc, ustp
    DOUBLE PRECISION:: y(nobs), bvec(1:nvars), bvecloc(1:nvars), res(1:nobs)
    DOUBLE PRECISION:: nzind(1:nvars), nzidx(1:nvars), nzge(1:nvars), coeffs(nvars)

    ! admm loop
    INTEGER:: liter, maxoutit
    INTEGER:: subidx(nvert)
    DOUBLE PRECISION:: epsout, errout
    DOUBLE PRECISION:: uvec(nobs), etavec(nobs), etavecold(nobs), uvecold(nobs)

    ! set some variables
    errcode = -1
    liter = 0
    nnzgr = 0
    nzind = 0.0D0
    nzidx = 0.0D0
    nzge  = 1  ! so that we start with all 
    
    ! definine admm starting parameters
    uvec      = 0.0D0
    uvecold   = 0.0D0
    etavec    = 0.0D0
    etavecold = 0.0D0
    bvec      = coeffs

    ! weight with rho
    res = res + y * (1/rho-1)
    y   = y / rho
    pen = pen / rho
    
    
    ! get subidx for projection 
    subidx = 0
    DO gidx = 1, nvert
            subidx(gidx) = (gidx-1) * nvert + gidx
    ENDDO
    
    
    
    ! check which ones can stay at zero 
    DO gidx = 1, nvars
            CALL rowidx(gidx,oidx,midx,nvert)
            ustp = (res(oidx) - res(midx)) / nobs
            tem = ABS(ustp)
            IF(tem >= pen) nzge(gidx) = 1
            IF(tem < pen)  nzge(gidx) = 0
    ENDDO

    
    !! --- start admm
    DO 
            liter = liter + 1

            !! --- first step
            !! --- lasso start --------------------------------------------------------------------------
            iter    = 0 
            diffloc = 0
            bloc    = 0

            ! check which ones can stay at zero 
            !DO gidx = 1, nvars
            !    CALL rowidx(gidx,oidx,midx,nvert)
            !    ustp = (res(oidx) - res(midx)) / nobs
            !    tem = ABS(ustp)
            !    IF(tem >= pen) nzge(gidx) = 1
            !    IF(tem < pen)  nzge(gidx) = 0
            !ENDDO


            DO
            ! --- iterate over all coefficients - if one coefficient was found to be zero it jumps to next iterate --------
                    iter  = iter + 1
                    errin = 0
                    DO gidx = 1,nvars
                            IF(nzge(gidx) == 0) CYCLE   ! check if coefficient reenters         
                            bloc = bvec(gidx)
                  
                            !! --- calculate gradient specific of this problem
                            CALL rowidx(gidx,oidx,midx,nvert)
                            ustp =  (res(oidx) - res(midx)) / nobs
                            !! --- gradient end
                            
                            ustp = bvec(gidx)+ustp
                            IF(postv == 1) ustp = max(ustp,REAL(0.0D0))               ! if only positiv wanted project
                            uabs = ABS(ustp)
                            tem=uabs-pen
                            IF(tem >= 0.0D0) THEN
                               bvec(gidx)=ustp*tem/uabs
                            ELSE
                               bvec(gidx)=0.0D0
                            ENDIF
                            diffloc = bvec(gidx) - bloc
                            IF(diffloc/=0.0D0) THEN
                                    errin=max(errin,(diffloc)**2)
                                    res(oidx) = res(oidx) - diffloc 
                                    res(midx) = res(midx) + diffloc  
                                    IF(nzind(gidx) == 0) THEN
                                           nnzgr = nnzgr+1
                                           nzind(gidx) = 1
                                           nzidx(nnzgr) = gidx 
                                    ENDIF
                            ENDIF
                    ENDDO

                    ! --- check errors 
                    IF((errin < epsin) .AND. (iter > 1)) THEN
                        errcode = 0
                        EXIT
                    ENDIF    
                    IF(iter > maxinit) THEN
                        errcode = 1
                        EXIT
                    ENDIF  
            ! --- iterate over reduced set of groups ----------------------------------------------------------
                    iter = iter + 1
                    DO jidx = 1,nnzgr
                            gidx  = nzidx(jidx)
                            bloc = bvec(gidx)
                            
                            !! ---- calculate gradient specific of this problem
                            CALL rowidx(gidx,oidx,midx,nvert)
                            ustp = (res(oidx) - res(midx)) / nobs
                            !! ---- gradient end
                            
                            ustp = bvec(gidx)+ustp
                            IF(postv == 1) ustp = max(ustp,REAL(0.0D0))               ! if only positiv wanted project
                            uabs = ABS(ustp)
                            tem=uabs-pen
                            IF(tem>=0.0D0) THEN
                               bvec(gidx)=ustp*tem/uabs
                            ELSE
                               bvec(gidx)=0.0D0
                            ENDIF
                            diffloc = bvec(gidx) - bloc
                            IF(diffloc/=0.0D0) THEN
                                    errin=max(errin,diffloc**2)
                                    res(oidx) = res(oidx) - diffloc 
                                    res(midx) = res(midx) + diffloc  
                                    
                            ENDIF
                    ENDDO
                    
                    ! --- check errors 
                    IF((errin < epsin) .AND. (iter > 1)) THEN
                        errcode = 0
                        EXIT
                    ENDIF    
                    IF(iter > maxinit) THEN
                        errcode = 1
                        EXIT
                    ENDIF  
                    
            ! --- check which coeffs can be jumped for next iteration -----------------------------------------
                  
                    
                    ! calculate  bvecloc = matmul(res, Z)/nobs
                    DO gidx = 1, nvars
                            CALL rowidx(gidx,oidx,midx,nvert)
                            bvecloc(gidx) = (res(oidx) - res(midx)) / nobs
                    ENDDO
                    
                    
                    DO gidx = 1, nvars
                        IF(nzge(gidx) == 1) CYCLE
                        ustp = bvecloc(gidx)
                        tem = ABS(ustp)
                        IF(tem >= pen)THEN
                            nzge(gidx) = 1
                        ENDIF
                    ENDDO
                    
            ENDDO
            !! --- lasso end ----------------------------------------------------------------------------
            
            !! --- projection step
            etavecold = etavec
            etavec    = etavec + y - res
            IF(normp <= 0) THEN
                    tem = MAXVAL(ABS(etavec(subidx)))
            ELSE 
                    tem = SUM(ABS(etavec(subidx))**(REAL(normp)))**(1.0 / REAL(normp)) !NORM2(etavec)
            END IF
            IF(tem >= 1) etavec = etavec / tem
            res = res + (etavec - etavecold)
            
            
            !! --- third step 
            uvecold = uvec
            uvec    = res - y
            
            !! -- update error and res 
            errout = NORM2(uvec - uvecold)
            res    = res + (uvec - uvecold)
            
            ! --- check errors
            IF(errout < epsout) THEN
                errcode = 0
                EXIT
            ENDIF
            IF(liter > maxoutit) THEN
                        errcode = 1
                        EXIT
            ENDIF  
            
   ENDDO
   ! final preparation of data
   res    = res - uvec - etavec
   coeffs = bvec
   pen = pen * rho
   res = res + y * (rho - 1)
   y   = rho * y
    
   RETURN
END SUBROUTINE admmGraph



! subroutine to calculate vec(L*covariance matrix ) 
! out: result is saved in y
! in: L is given as wcurr (columnwise), Sigma is usually the covariance matrix 
SUBROUTINE leftvec(ynew, Sigma, wcurr, nvert)

        IMPLICIT NONE
      
        ! parameters 
        INTEGER:: nvert, iidx, jidx, llidx1, llidx2
        DOUBLE PRECISION:: wsums(nvert), wcurr(nvert*(nvert-1)), ynew(nvert*nvert), Sigma(nvert,nvert)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: wSprod
        
        
        DO iidx  = 1,nvert
                llidx1       = (iidx-1)*(nvert-1)+1
                llidx2       = iidx*(nvert-1)
                wsums(iidx)  = SUM(wcurr(llidx1:llidx2))
        END DO 
        
        DO jidx = 1,nvert
                ALLOCATE(wSprod(1:nvert))
                DO iidx = 1,nvert
                    llidx1       = (iidx-1)*(nvert-1)+1
                    llidx2       = iidx*(nvert-1)
                    IF(iidx == 1) THEN
                        wSprod(iidx) = DOT_PRODUCT(wcurr(llidx1:llidx2), Sigma(2:nvert,jidx))
                    ELSE IF(iidx == nvert) THEN
                        wSprod(iidx) = DOT_PRODUCT(wcurr(llidx1:llidx2), Sigma(1:(nvert-1),jidx))
                    ELSE 
                        wSprod(iidx) = DOT_PRODUCT(wcurr(llidx1:(llidx1+iidx-2)), Sigma(1:(iidx-1),jidx))
                        wSprod(iidx) = wSprod(iidx) +  DOT_PRODUCT(wcurr((llidx1 + iidx - 1):llidx2), Sigma((iidx+1):nvert,jidx))
                    ENDIF
                END DO
                ynew(((jidx-1)*nvert+1):(jidx*nvert)) =  Sigma(:,jidx) * wsums  - wSprod        
                DEALLOCATE(wSprod)
        END DO

END SUBROUTINE leftvec




! subroutine to calculate graph structure on the basis of total variance
SUBROUTINE varGraph(Sigma, Ws, pens, npens, nvert, postv, normp, &
                    maxoutit, epsout, errout, &
                    lassoMaxinit, LassoEpsin, lassoErrin, &
                    admmMaxoutit, admmEpsout, admmErrout, errcode)
    IMPLICIT NONE
    
    ! define
    INTEGER:: npens, nvert, postv, nobs, nvars
    DOUBLE PRECISION:: pen, pens(npens), Sigma(nvert,nvert)
    DOUBLE PRECISION:: normp
    
    !admm loop
    INTEGER:: lassoMaxinit, admmMaxoutit, errcode
    DOUBLE PRECISION:: rho, LassoEpsin, lassoErrin, admmEpsout, admmErrout
    
    ! outer loop parameters 
    INTEGER:: lamidx, lidx, iidx, jidx, cnt, maxoutit
    DOUBLE PRECISION:: epsout, errout
    DOUBLE PRECISION:: ynew(nvert*nvert), yold(nvert*nvert), coeffs(nvert*(nvert-1)), coeffsOld(nvert*(nvert-1))
    DOUBLE PRECISION:: start(nvert*nvert), res(nvert*nvert), Ws(nvert,nvert,npens)
    
    

    ! init loop parameters
    nobs    = nvert * nvert
    nvars   = nvert * (nvert-1)
    errcode = 1
    
    
    ! init parameters
    start = 1
    ynew = 0.0D0 
    yold = 0.0D0
    coeffs    = 1
    coeffsOld = 0.0D0
    
    ! get starting res (Z*coeffs)
    res = -1
    DO iidx = 1,nvert
            res((iidx-1)*(nvert) + iidx) = nvert - 1
    END DO
    res = - res
    
    DO lamidx = 1, npens
            
            
            ! start with res: y = vec(LX'X) with L coda and Z*coeffs, coeffs from before
            yold = ynew 
            CALL leftvec(ynew, Sigma, start, nvert)    ! we start with coeffs = start (usually all ones)
            res = res + (ynew - yold)

            ! get new lambda 
            cnt = 1
            rho = pens(lamidx) * (nvert*nvert)
            pen = pens(lamidx) 
            

            ! run optimization for fixed lambda 
            DO lidx = 1,maxoutit
                    
                    ! save coeffs
                    coeffsOld = coeffs
                    
                    ! run admm to find next best solution coeffs
                    CALL admmGraph(ynew, coeffs, res, postv, normp, nobs, nvars, nvert, &
                                         pen, rho, lassoMaxinit, LassoEpsin, lassoErrin, &
                                         admmMaxoutit, admmEpsout, admmErrout, errcode)
                    !save old y
                    yold = ynew 
                    
                    !calculate vec(L*X'X)
                    CALL leftvec(ynew, Sigma, coeffs, nvert)
                     
                    ! update residual and save coeffs
                    res = res + (ynew - yold)
                    
            
                    ! calculate error
                    errout = NORM2(coeffs - coeffsOld) / (1 + NORM2(coeffsOld))
                    
                    ! check conditions 
                    IF(errout < epsout .AND. lidx > 1) THEN
                        errcode = 0
                        EXIT
                    END IF
                    IF(lidx == maxoutit) errcode = 1
            END DO
            
            ! return the solution
            DO iidx = 1,nvert
                DO jidx = 1,nvert
                    IF(iidx /= jidx) THEN
                        Ws(iidx,jidx,lamidx) = coeffs(cnt)
                        cnt = cnt + 1
                    END IF    
                END DO
            END DO
    END DO 

    
   RETURN
END SUBROUTINE varGraph





!!  --------------------- symmetric


!! lasso for graphical input matrix 
SUBROUTINE admmGraphSym(y, coeffs, res, postv, normp, nobs, nvars, nvert, &
                        pen, rho, maxinit, epsin, errin, maxoutit, epsout, errout, errcode)
    IMPLICIT NONE
    
    ! data specific params
    INTEGER:: nobs, nvars
    DOUBLE PRECISION:: normp

    ! penalty specific params
    INTEGER::gidx, postv
    DOUBLE PRECISION:: pen

    ! lasso loop specific params 
    INTEGER:: jidx, iter, nnzgr, maxinit, errcode
    INTEGER:: oidx1, oidx2, midx1, midx2, gidx1, gidx2, nvert
    DOUBLE PRECISION:: epsin
    DOUBLE PRECISION:: inter, uabs, tem, errin, diffloc, bloc, ustp
    DOUBLE PRECISION:: y(nobs) , bvec(1:nvars), bvecloc(1:nvars), res(1:nobs)
    DOUBLE PRECISION:: nzind(1:nvars), nzidx(1:nvars), nzge(1:nvars), coeffs(nvars)

    ! admm loop
    INTEGER:: liter, maxoutit
    DOUBLE PRECISION:: epsout, errout
    DOUBLE PRECISION:: rho, yloc(nobs), uvec(nobs), etavec(nobs), etavecold(nobs), uvecold(nobs)

  
    ! initialize 
    errcode = -1
    liter   = 0

    nnzgr   = 0
    nzind   = 0.0D0
    nzidx   = 0.0D0
    nzge    = 1  ! so that we start with all 
    
    ! definine admm starting parameters
    uvec    = 0.0D0
    uvecold = 0.0D0
    etavec  = 0.0D0
    etavecold = 0.0D0
    bvec      = coeffs
    
    ! weight with rho
    res = res + y * (1/rho-1)
    y   = y / rho
    pen = pen / rho
    

    ! check which ones can stay at zero 
    DO gidx = 1, nvars
            CALL colidxSymToNonsym(gidx,gidx1,gidx2,nvert)    ! column indices (gidx1,gidx2) in full matrix 
            CALL rowidx(gidx1,oidx1,midx1,nvert)              ! row index (idx1,midx1) from col index gidx1
            CALL rowidx(gidx2,oidx2,midx2,nvert)              ! row index (idx2,midx2) from col index gidx2
            ustp =  (res(oidx1) + res(oidx2) - res(midx1) - res(midx2)) / nobs
            tem = ABS(ustp)
            IF(tem >= pen) nzge(gidx) = 1
            IF(tem < pen) nzge(gidx) = 0
    ENDDO


    !! --- start admm
    DO 
            liter = liter + 1
            
            ! check which ones can stay at zero 
            !DO gidx = 1, nvars
            !    CALL colidxSymToNonsym(gidx,gidx1,gidx2,nvert)    ! column indices (gidx1,gidx2) in full matrix 
            !    CALL rowidx(gidx1,oidx1,midx1,nvert)              ! row index (idx1,midx1) from col index gidx1
            !    CALL rowidx(gidx2,oidx2,midx2,nvert)              ! row index (idx2,midx2) from col index gidx2
            !    ustp =  (res(oidx1) + res(oidx2) - res(midx1) - res(midx2)) / nobs
            !    tem = ABS(ustp)
            !    IF(tem >= pen) nzge(gidx) = 1
            !    IF(tem < pen) nzge(gidx) = 0
            !ENDDO

            !! --- first step
            !! ---  lasso start --------------------------------------------------------------------------
            DO
            ! --- iterate over all coefficients - if one coefficient was found to be zero it jumps to next iterate --------
                    iter  = iter + 1
                    errin = 0
                    DO gidx = 1,nvars
                            IF(nzge(gidx) == 0) CYCLE   ! check if coefficient reenters         
                            bloc = bvec(gidx)
                  
                            !! ---- calculate gradient specific of this problem
                            CALL colidxSymToNonsym(gidx,gidx1,gidx2,nvert)    ! column indices (gidx1,gidx2) in full matrix 
                            CALL rowidx(gidx1,oidx1,midx1,nvert)              ! row index (idx1,midx1) from col index gidx1
                            CALL rowidx(gidx2,oidx2,midx2,nvert)              ! row index (idx2,midx2) from col index gidx2
                            ustp =  (res(oidx1) + res(oidx2) - res(midx1) - res(midx2)) / nobs
                            !! ---- gradient end
                            
                            ustp = bvec(gidx)+ustp
                            IF(postv == 1) ustp = max(ustp,REAL(0.0D0))               ! if only positiv wanted project
                            uabs = ABS(ustp)
                            tem=uabs-pen
                            IF(tem>=0.0D0) THEN
                               bvec(gidx)=ustp*tem/uabs
                            ELSE
                               bvec(gidx)=0.0D0
                            ENDIF
                            diffloc = bvec(gidx) - bloc
                            IF(diffloc/=0.0D0) THEN
                                    errin=max(errin,(diffloc)**2)
                                    res(oidx1) = res(oidx1) - diffloc 
                                    res(midx1) = res(midx1) + diffloc 
                                    res(oidx2) = res(oidx2) - diffloc 
                                    res(midx2) = res(midx2) + diffloc 
                                    IF(nzind(gidx) == 0) THEN
                                           nnzgr = nnzgr+1
                                           nzind(gidx) = 1
                                           nzidx(nnzgr) = gidx 
                                    ENDIF
                            ENDIF
                    ENDDO
                    ! --- check errors 
                    IF((errin < epsin) .AND. (iter > 1)) THEN
                        errcode = 0
                        EXIT
                    ENDIF    
                    IF(iter > maxinit) THEN
                        errcode = 1
                        EXIT
                    ENDIF  
            ! --- iterate over reduced set of groups ----------------------------------------------------------
                    iter = iter + 1
                    DO jidx = 1,nnzgr
                            gidx  = nzidx(jidx)
                            bloc = bvec(gidx)
                            
                            !! ---- calculate gradient specific of this problem
                            CALL colidxSymToNonsym(gidx,gidx1,gidx2,nvert)    ! column indices (gidx1,gidx2) in full matrix 
                            CALL rowidx(gidx1,oidx1,midx1,nvert)              ! row index (idx1,midx1) from col index gidx1
                            CALL rowidx(gidx2,oidx2,midx2,nvert)              ! row index (idx2,midx2) from col index gidx2
                            ustp =  (res(oidx1) + res(oidx2) - res(midx1) - res(midx2)) / nobs
                            !! ---- gradient end
                            
                            ustp = bvec(gidx)+ustp
                            IF(postv == 1) ustp = max(ustp,REAL(0.0D0))               ! if only positiv wanted project
                            uabs = ABS(ustp)
                            tem=uabs-pen
                            IF(tem>=0.0D0) THEN
                               bvec(gidx)=ustp*tem/uabs
                            ELSE
                               bvec(gidx)=0.0D0
                            ENDIF
                            diffloc = bvec(gidx) - bloc
                            IF(diffloc/=0.0D0) THEN
                                    errin=max(errin,diffloc**2)
                                    res(oidx1) = res(oidx1) - diffloc 
                                    res(midx1) = res(midx1) + diffloc 
                                    res(oidx2) = res(oidx2) - diffloc 
                                    res(midx2) = res(midx2) + diffloc 
                            ENDIF
                    ENDDO
                    ! --- check errors 
                    IF((errin < epsin) .AND. (iter > 1)) THEN
                        errcode = 0
                        EXIT
                    ENDIF    
                    IF(iter > maxinit) THEN
                        errcode = 1
                        EXIT
                    ENDIF  
                    
            ! --- check which coeffs can be jumped for next iteration -----------------------------------------
                  
                    ! calculate  bvecloc = matmul(res, Z)/nobs
                    DO gidx = 1, nvars
                            CALL colidxSymToNonsym(gidx,gidx1,gidx2,nvert)    ! column indices (gidx1,gidx2) in full matrix 
                            CALL rowidx(gidx1,oidx1,midx1,nvert)              ! row index (idx1,midx1) from col index gidx1
                            CALL rowidx(gidx2,oidx2,midx2,nvert)              ! row index (idx2,midx2) from col index gidx2
                            bvecloc(gidx) = (res(oidx1) + res(oidx2) - res(midx1) - res(midx2)) / nobs
                    ENDDO
                    
                    
                    DO gidx = 1, nvars
                        IF(nzge(gidx) == 1) CYCLE
                        ustp = bvecloc(gidx)
                        tem = ABS(ustp)
                        IF(tem >= pen)THEN
                            nzge(gidx) = 1
                        ENDIF
                    ENDDO
                    
            ENDDO
            !! --- lasso end ----------------------------------------------------------------------------

            !! --- projection step
            etavecold = etavec
            etavec    = etavec + y - res
            IF(normp <= 0) THEN
                    tem = MAXVAL(ABS(etavec))
            ELSE 
                    tem = SUM(ABS(etavec)**(REAL(normp)))**(1.0 / REAL(normp)) !NORM2(etavec)
            END IF
            IF(tem >= 1) etavec = etavec / tem
            res = res + (etavec - etavecold)
            
            
            !! --- third step 
            uvecold = uvec
            uvec    = res - y 
            
            !! -- update error and res 
            errout = NORM2(uvec - uvecold)
            res    = res + (uvec - uvecold)
            
            ! --- check errors
            IF(errout < epsout) THEN
                errcode = 0
                EXIT
            ENDIF
            IF(liter > maxoutit) THEN
                        errcode = 1
                        EXIT
            ENDIF  
            
   ENDDO
   ! final preparation of data
   coeffs = bvec
   pen = pen * rho
   res = res - uvec - etavec + y * (rho - 1)
   y   = rho * y



   RETURN
END SUBROUTINE admmGraphSym



! subroutine to calculate vec(L*covariance matrix ) 
! out: result is saved in y
! in: L is given as wcurr (columnwise), Sigma is usually the covariance matrix 
SUBROUTINE leftvecSym(ynew, Sigma, wcurr, nvert)

        IMPLICIT NONE
      
        ! parameters 
        INTEGER:: nvert, iidx, jidx, cidx, cidx1, cidx2, llidx1, llidx2
        DOUBLE PRECISION:: wsums(nvert), wcurr((nvert*(nvert-1))/2), ynew(nvert*nvert), Sigma(nvert,nvert)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: wcurrFull, wSprod
  
        !map from smaller wcurr to bigger wcurrFull
        ALLOCATE(wcurrFull(1:(nvert*(nvert-1))))
        
        DO cidx = 1,((nvert*(nvert-1))/2)
               ! get the indices cidx1, cidx2
               CALL colidxSymToNonsym(cidx,cidx1,cidx2,nvert)
               
               ! set so symmetric 
               wCurrFull(cidx1) = wcurr(cidx)
               wCurrFull(cidx2) = wcurr(cidx)
        END DO
        
        
        ! now do vec(L*covariance matrix ) 
        DO iidx  = 1,nvert
                llidx1       = (iidx-1)*(nvert-1)+1
                llidx2       = iidx*(nvert-1)
                wsums(iidx)  = SUM(wcurrFull(llidx1:llidx2))
        END DO 
        
        DO jidx = 1,nvert
                ALLOCATE(wSprod(1:nvert))
                DO iidx = 1,nvert
                    llidx1       = (iidx-1)*(nvert-1)+1
                    llidx2       = iidx*(nvert-1)
                    IF(iidx == 1) THEN
                        wSprod(iidx) = DOT_PRODUCT(wcurrFull(llidx1:llidx2), Sigma(2:nvert,jidx))
                    ELSE IF(iidx == nvert) THEN
                        wSprod(iidx) = DOT_PRODUCT(wcurrFull(llidx1:llidx2), Sigma(1:(nvert-1),jidx))
                    ELSE 
                        wSprod(iidx)=DOT_PRODUCT(wcurrFull(llidx1:(llidx1+iidx-2)), Sigma(1:(iidx-1),jidx))
                        wSprod(iidx)=wSprod(iidx)+DOT_PRODUCT(wcurrFull((llidx1 + iidx - 1):llidx2), Sigma((iidx+1):nvert,jidx))
                    ENDIF
                END DO
                ynew(((jidx-1)*nvert+1):(jidx*nvert)) =  Sigma(:,jidx) * wsums  - wSprod        
                DEALLOCATE(wSprod)
        END DO
        DEALLOCATE(wcurrFull)
        
        
END SUBROUTINE leftvecSym


! subroutine to calculate graph structure on the basis of total variance if symmetric is assumed
SUBROUTINE varGraphSym(Sigma, Ws, pens, npens, nvert, postv, normp, &
                    maxoutit, epsout, errout, &
                    lassoMaxinit, LassoEpsin, lassoErrin, &
                    admmMaxoutit, admmEpsout, admmErrout, errcode)
    IMPLICIT NONE
    
    ! define
    INTEGER:: npens, nvert, postv, nobs, nvars
    DOUBLE PRECISION:: pens(npens), pen, Sigma(nvert,nvert)
    DOUBLE PRECISION:: normp
    
    !admm loop
    INTEGER:: lassoMaxinit, admmMaxoutit, errcode
    DOUBLE PRECISION:: rho, LassoEpsin, lassoErrin, admmEpsout, admmErrout
    
    ! outer loop parameters 
    INTEGER:: lamidx, lidx, iidx, jidx, cnt, maxoutit
    DOUBLE PRECISION:: epsout, errout
    DOUBLE PRECISION::  ynew(nvert*nvert), yold(nvert*nvert), res(nvert*nvert), Ws(nvert,nvert,npens)
    DOUBLE PRECISION::  start(nvert*nvert), coeffs((nvert * (nvert-1)) / 2), coeffsOld((nvert * (nvert-1)) / 2)


    ! init loop parameters
    nobs    = nvert * nvert
    nvars   = (nvert * (nvert-1)) / 2
    errcode = 1
    
    
    ! init parameters
    start = 1
    ynew = 0.0D0 
    yold = 0.0D0
    coeffs  = 1
    coeffsOld = 0.0D0
    
    ! get starting res (Z*coeffs)
    res = -1
    DO iidx = 1,nvert
            res((iidx-1)*(nvert) + iidx) = nvert - 1
    END DO
    res = -res
    
    DO lamidx = 1, npens
    
            ! start with res: y = vec(LX'X) with L coda and Z*coeffs, coeffs from before
            yold = ynew 
            CALL leftvecSym(ynew, Sigma, start, nvert)
            res = res + (ynew - yold)
    
            ! get new lambda 
            cnt = 1
            rho = pens(lamidx)  * (nvert*nvert)
            pen = pens(lamidx) 
    
            ! run optimization for fixed lambda 
            DO lidx = 1,maxoutit
                    
                    ! save coeffs
                    coeffsOld = coeffs
                    
                    ! run admm to find next best solution coeffs
                    CALL admmGraphSym(ynew, coeffs, res, postv, normp, nobs, nvars, nvert, &
                                         pen, rho, lassoMaxinit, LassoEpsin, lassoErrin, &
                                        admmMaxoutit, admmEpsout, admmErrout, errcode)
            
            
                    !save old y
                    yold = ynew 
                    
                    !calculate vec(L*X'X)
                    CALL leftvecSym(ynew, Sigma, coeffs, nvert)
                     
                    ! update residual and save coeffs
                    res = res + (ynew - yold)
            
                    ! calculate error
                    errout = NORM2(coeffs - coeffsOld) / (1 + NORM2(coeffsOld))
                    
                    ! check conditions 
                    IF(errout < epsout .AND. lidx > 1) THEN
                        errcode = 0
                        EXIT
                    END IF
                    IF(lidx == maxoutit) errcode = 1
            END DO
            
            ! return the solution
            DO iidx = 1,nvert
               DO jidx = 1,nvert
                   IF(iidx < jidx) THEN
                        Ws(iidx,jidx,lamidx) = coeffs(cnt)
                        Ws(jidx,iidx,lamidx) = coeffs(cnt)
                        cnt = cnt + 1
                    END IF    
                END DO
            END DO
    END DO 
    


   RETURN
END SUBROUTINE varGraphSym



!! ----- simple single L which is symmetric

SUBROUTINE singleVar(Sigma, Ws, pens, npens, nvert, postv, normp, &
                    lassoMaxinit, LassoEpsin, lassoErrin, &
                    admmMaxoutit, admmEpsout, admmErrout, errcode)
        IMPLICIT NONE
        
        ! define
        INTEGER:: npens, nvert, postv, nobs, nvars
        DOUBLE PRECISION:: pens(npens), pen, Sigma(nvert,nvert)
        DOUBLE PRECISION:: normp
        
        !admm loop
        INTEGER:: lassoMaxinit, admmMaxoutit, errcode
        DOUBLE PRECISION:: rho, LassoEpsin, lassoErrin, admmEpsout, admmErrout
        
        ! outer loop parameters 
        INTEGER:: lamidx, iidx, jidx, cnt
        DOUBLE PRECISION:: ynew(nvert*nvert), res(nvert*nvert), Ws(nvert,nvert,npens)
        DOUBLE PRECISION:: coeffs((nvert * (nvert-1)) / 2)
    
        ! init loop parameters
        nobs    = nvert * nvert
        nvars   = (nvert * (nvert-1)) / 2
        errcode = 1
            
        ! initialize parameters
        ynew   = 0.0D0 
        coeffs = 0.0D0

        ! calculate left vector vec(Sigma) to start 
        DO iidx = 1,nvert
                ynew(((iidx-1)*nvert+1):(iidx*nvert)) = Sigma(:,iidx)
        END DO 
        res = ynew
        
        DO lamidx = 1,npens
        
        
              ! get new lambda 
              cnt = 1
              rho = pens(lamidx) 
              pen = pens(lamidx) 
        
              ! run admm to find next best solution coeffs - res is automatically updated
              CALL admmGraphSym(ynew, coeffs, res, postv, normp, nobs, nvars, nvert, &
                                   pen, rho, lassoMaxinit, LassoEpsin, lassoErrin, &
                                   admmMaxoutit, admmEpsout, admmErrout, errcode)
                       

             ! return the solution
              DO iidx = 1,nvert
                 DO jidx = 1,nvert
                     IF(iidx < jidx) THEN
                          Ws(iidx,jidx,lamidx) = coeffs(cnt)
                          Ws(jidx,iidx,lamidx) = coeffs(cnt)
                          cnt = cnt + 1
                     END IF    
                  END DO
              END DO
              
        END DO      
        

        RETURN
END SUBROUTINE singleVar



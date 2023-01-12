        chbra => this%chbra
        chket => this%chket
        jbra = chbra%GetJ()
        tbra = chbra%GetT()
        jket = chket%GetJ()
        tket = chket%GetT()
        !$omp parallel
        !$omp do private(bra, mbra, pbra, lbra, sbra, ket, mket, pket, lket, sket, me )
        do bra = 1, chbra%GetNumberStates()
          mbra => chbra%getp(bra)
          pbra = mbra%GetP()*hc
          lbra = mbra%GetL()
          sbra = mbra%GetS()
          do ket = 1, chket%GetNumberStates()
            mket => chket%getp(ket)
            pket = mket%GetP()*hc
            lket = mket%GetL()
            sket = mket%GetS()
            me = op_class%calc_matrix_element(pbra, lbra, sbra, jbra, tbra, pket, lket, sket, jket, tket, .false.)
            this%m(bra,ket) = me * hc**3 ! in unit of MeV fm3
          end do
        end do
        !$omp end do
        !$omp end parallel


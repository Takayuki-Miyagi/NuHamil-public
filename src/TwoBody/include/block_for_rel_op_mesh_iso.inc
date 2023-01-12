        do ichbra = 1, this%ms%GetNumberChannels()
          do ichket = 1, this%ms%GetNumberChannels()
            if( this%MatCh(ichbra,ichket)%is_zero ) cycle
            op_ch => this%MatCh(ichbra,ichket)
            chbra => op_ch%chbra
            chket => op_ch%chket
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
                op_ch%m(bra,ket) = me * hc**3 ! in unit of MeV fm3
              end do
            end do
            !$omp end do
            !$omp end parallel
          end do
        end do


/**
* \brief Applies the controlled-gate \a A to the part \a subsys
* of the multi-partite state vector or density matrix \a state
* \see qpp::Gates::CTRL()
*
* \note The dimension of the gate \a A must match
* the dimension of \a subsys.
* Also, all control subsystems in \a ctrl must have the same dimension.
*
* \param state Eigen expression
* \param A Eigen expression
* \param ctrl Control subsystem indexes
* \param subsys Subsystem indexes where the gate \a A is applied
* \param dims Dimensions of the multi-partite system
* \return CTRL-A gate applied to the part \a subsys of \a state
*/
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar>
applyCTRL(const Eigen::MatrixBase<Derived1>& state,
          const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
          const std::vector<idx>& subsys, const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived1>::EvalReturnType& rstate =
        state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::applyCTRL()");

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::applyCTRL()");

    // check zero sizes
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::applyCTRL()");

    // check square matrix for the gate
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::applyCTRL()");

    // check that all control subsystems have the same dimension
    idx d = ctrl.size() > 0 ? dims[ctrl[0]] : 1;
    for (idx i = 1; i < ctrl.size(); ++i)
        if (dims[ctrl[i]] != d)
            throw exception::DimsNotEqual("qpp::applyCTRL()");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::applyCTRL()");

    // check subsys is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(subsys, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL()");

    // check that gate matches the dimensions of the subsys
    std::vector<idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];
    if (!internal::check_dims_match_mat(subsys_dims, rA))
        throw exception::MatrixMismatchSubsys("qpp::applyCTRL()");

    std::vector<idx> ctrlgate = ctrl; // ctrl + gate subsystem vector
    ctrlgate.insert(std::end(ctrlgate), std::begin(subsys), std::end(subsys));
    std::sort(std::begin(ctrlgate), std::end(ctrlgate));

    // check that ctrl + gate subsystem is valid
    // with respect to local dimensions
    if (!internal::check_subsys_match_dims(ctrlgate, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL()");
    // END EXCEPTION CHECKS

    // construct the table of A^i and (A^dagger)^i
    std::vector<dyn_mat<typename Derived1::Scalar>> Ai;
    std::vector<dyn_mat<typename Derived1::Scalar>> Aidagger;
    for (idx i = 0; i < std::max(d, static_cast<idx>(2)); ++i) {
        Ai.push_back(powm(rA, i));
        Aidagger.push_back(powm(adjoint(rA), i));
    }

    idx D = static_cast<idx>(rstate.rows()); // total dimension
    idx N = dims.size();                     // total number of subsystems
    idx ctrlsize = ctrl.size();              // number of ctrl subsystem
    idx ctrlgatesize = ctrlgate.size();      // number of ctrl+gate subsystems
    idx subsyssize = subsys.size(); // number of subsystems of the target
    // dimension of ctrl subsystem
    idx Dctrl = static_cast<idx>(std::llround(std::pow(d, ctrlsize)));
    idx DA = static_cast<idx>(rA.rows()); // dimension of gate subsystem

    idx Cdims[maxn];          // local dimensions
    idx CdimsA[maxn];         // local dimensions
    idx CdimsCTRL[maxn];      // local dimensions
    idx CdimsCTRLA_bar[maxn]; // local dimensions

    // compute the complementary subsystem of ctrlgate w.r.t. dims
    std::vector<idx> ctrlgate_bar = complement(ctrlgate, N);
    // number of subsystems that are complementary to the ctrl+gate
    idx ctrlgate_barsize = ctrlgate_bar.size();

    idx DCTRLA_bar = 1; // dimension of the rest
    for (idx i = 0; i < ctrlgate_barsize; ++i)
        DCTRLA_bar *= dims[ctrlgate_bar[i]];

    for (idx k = 0; k < N; ++k)
        Cdims[k] = dims[k];
    for (idx k = 0; k < subsyssize; ++k)
        CdimsA[k] = dims[subsys[k]];
    for (idx k = 0; k < ctrlsize; ++k)
        CdimsCTRL[k] = d;
    for (idx k = 0; k < ctrlgate_barsize; ++k)
        CdimsCTRLA_bar[k] = dims[ctrlgate_bar[k]];

    // worker, computes the coefficient and the index for the ket case
    // used in #pragma omp parallel for collapse
    auto coeff_idx_ket = [&](idx i_, idx m_, idx r_) noexcept
                             ->std::pair<typename Derived1::Scalar, idx> {
        idx indx = 0;
        typename Derived1::Scalar coeff = 0;

        idx Cmidx[maxn];          // the total multi-index
        idx CmidxA[maxn];         // the gate part multi-index
        idx CmidxCTRLA_bar[maxn]; // the rest multi-index

        // compute the index

        // set the CTRL part
        for (idx k = 0; k < ctrlsize; ++k) {
            Cmidx[ctrl[k]] = i_;
        }

        // set the rest
        internal::n2multiidx(r_, N - ctrlgatesize, CdimsCTRLA_bar,
                             CmidxCTRLA_bar);
        for (idx k = 0; k < N - ctrlgatesize; ++k) {
            Cmidx[ctrlgate_bar[k]] = CmidxCTRLA_bar[k];
        }

        // set the A part
        internal::n2multiidx(m_, subsyssize, CdimsA, CmidxA);
        for (idx k = 0; k < subsyssize; ++k) {
            Cmidx[subsys[k]] = CmidxA[k];
        }

        // we now got the total index
        indx = internal::multiidx2n(Cmidx, N, Cdims);

        // compute the coefficient
        for (idx n_ = 0; n_ < DA; ++n_) {
            internal::n2multiidx(n_, subsyssize, CdimsA, CmidxA);
            for (idx k = 0; k < subsyssize; ++k) {
                Cmidx[subsys[k]] = CmidxA[k];
            }
            coeff +=
                Ai[i_](m_, n_) * rstate(internal::multiidx2n(Cmidx, N, Cdims));
        }

        return std::make_pair(coeff, indx);
    }; /* end coeff_idx_ket */

    // worker, computes the coefficient and the index
    // for the density matrix case
    // used in #pragma omp parallel for collapse
    auto coeff_idx_rho = [&](idx i1_, idx m1_, idx r1_, idx i2_, idx m2_,
                             idx r2_) noexcept
                             ->std::tuple<typename Derived1::Scalar, idx, idx> {
        idx idxrow = 0;
        idx idxcol = 0;
        typename Derived1::Scalar coeff = 0, lhs = 1, rhs = 1;

        idx Cmidxrow[maxn];          // the total row multi-index
        idx Cmidxcol[maxn];          // the total col multi-index
        idx CmidxArow[maxn];         // the gate part row multi-index
        idx CmidxAcol[maxn];         // the gate part col multi-index
        idx CmidxCTRLrow[maxn];      // the control row multi-index
        idx CmidxCTRLcol[maxn];      // the control col multi-index
        idx CmidxCTRLA_barrow[maxn]; // the rest row multi-index
        idx CmidxCTRLA_barcol[maxn]; // the rest col multi-index

        // compute the ket/bra indexes

        // set the CTRL part
        internal::n2multiidx(i1_, ctrlsize, CdimsCTRL, CmidxCTRLrow);
        internal::n2multiidx(i2_, ctrlsize, CdimsCTRL, CmidxCTRLcol);

        for (idx k = 0; k < ctrlsize; ++k) {
            Cmidxrow[ctrl[k]] = CmidxCTRLrow[k];
            Cmidxcol[ctrl[k]] = CmidxCTRLcol[k];
        }

        // set the rest
        internal::n2multiidx(r1_, N - ctrlgatesize, CdimsCTRLA_bar,
                             CmidxCTRLA_barrow);
        internal::n2multiidx(r2_, N - ctrlgatesize, CdimsCTRLA_bar,
                             CmidxCTRLA_barcol);
        for (idx k = 0; k < N - ctrlgatesize; ++k) {
            Cmidxrow[ctrlgate_bar[k]] = CmidxCTRLA_barrow[k];
            Cmidxcol[ctrlgate_bar[k]] = CmidxCTRLA_barcol[k];
        }

        // set the A part
        internal::n2multiidx(m1_, subsyssize, CdimsA, CmidxArow);
        internal::n2multiidx(m2_, subsyssize, CdimsA, CmidxAcol);
        for (idx k = 0; k < subsys.size(); ++k) {
            Cmidxrow[subsys[k]] = CmidxArow[k];
            Cmidxcol[subsys[k]] = CmidxAcol[k];
        }

        // we now got the total row/col indexes
        idxrow = internal::multiidx2n(Cmidxrow, N, Cdims);
        idxcol = internal::multiidx2n(Cmidxcol, N, Cdims);

        // check whether all CTRL row and col multi indexes are equal
        bool all_ctrl_rows_equal = true;
        bool all_ctrl_cols_equal = true;

        idx first_ctrl_row, first_ctrl_col;
        if (ctrlsize > 0) {
            first_ctrl_row = CmidxCTRLrow[0];
            first_ctrl_col = CmidxCTRLcol[0];
        } else {
            first_ctrl_row = first_ctrl_col = 1;
        }

        for (idx k = 1; k < ctrlsize; ++k) {
            if (CmidxCTRLrow[k] != first_ctrl_row) {
                all_ctrl_rows_equal = false;
                break;
            }
        }
        for (idx k = 1; k < ctrlsize; ++k) {
            if (CmidxCTRLcol[k] != first_ctrl_col) {
                all_ctrl_cols_equal = false;
                break;
            }
        }

        // at least one control activated, compute the coefficient
        for (idx n1_ = 0; n1_ < DA; ++n1_) {
            internal::n2multiidx(n1_, subsyssize, CdimsA, CmidxArow);
            for (idx k = 0; k < subsyssize; ++k) {
                Cmidxrow[subsys[k]] = CmidxArow[k];
            }
            idx idxrowtmp = internal::multiidx2n(Cmidxrow, N, Cdims);

            if (all_ctrl_rows_equal) {
                lhs = Ai[first_ctrl_row](m1_, n1_);
            } else {
                lhs = (m1_ == n1_) ? 1 : 0; // identity matrix
            }

            for (idx n2_ = 0; n2_ < DA; ++n2_) {
                internal::n2multiidx(n2_, subsyssize, CdimsA, CmidxAcol);
                for (idx k = 0; k < subsyssize; ++k) {
                    Cmidxcol[subsys[k]] = CmidxAcol[k];
                }

                if (all_ctrl_cols_equal) {
                    rhs = Aidagger[first_ctrl_col](n2_, m2_);
                } else {
                    rhs = (n2_ == m2_) ? 1 : 0; // identity matrix
                }

                idx idxcoltmp = internal::multiidx2n(Cmidxcol, N, Cdims);

                coeff += lhs * rstate(idxrowtmp, idxcoltmp) * rhs;
            }
        }

        return std::make_tuple(coeff, idxrow, idxcol);
    }; /* end coeff_idx_rho */

    //************ ket ************//
    if (internal::check_cvector(rstate)) // we have a ket
    {
        // check that dims match state vector
        if (!internal::check_dims_match_cvect(dims, rstate))
            throw exception::DimsMismatchCvector("qpp::applyCTRL()");
        if (D == 1)
            return rstate;

        dyn_mat<typename Derived1::Scalar> result = rstate;

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif // WITH_OPENMP_
        for (idx m = 0; m < DA; ++m)
            for (idx r = 0; r < DCTRLA_bar; ++r) {
                if (ctrlsize == 0) // no control
                {
                    result(coeff_idx_ket(1, m, r).second) =
                        coeff_idx_ket(1, m, r).first;
                } else
                    for (idx i = 0; i < d; ++i) {
                        result(coeff_idx_ket(i, m, r).second) =
                            coeff_idx_ket(i, m, r).first;
                    }
            }

        return result;
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rstate)) // we have a density operator
    {
        // check that dims match state matrix
        if (!internal::check_dims_match_mat(dims, rstate))
            throw exception::DimsMismatchMatrix("qpp::applyCTRL()");

        if (D == 1)
            return rstate;

        dyn_mat<typename Derived1::Scalar> result = rstate;

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(4)
#endif // WITH_OPENMP_
        for (idx m1 = 0; m1 < DA; ++m1)
            for (idx r1 = 0; r1 < DCTRLA_bar; ++r1)
                for (idx m2 = 0; m2 < DA; ++m2)
                    for (idx r2 = 0; r2 < DCTRLA_bar; ++r2)
                        if (ctrlsize == 0) // no control
                        {
                            auto coeff_idxes =
                                coeff_idx_rho(1, m1, r1, 1, m2, r2);
                            result(std::get<1>(coeff_idxes),
                                   std::get<2>(coeff_idxes)) =
                                std::get<0>(coeff_idxes);
                        } else {
                            for (idx i1 = 0; i1 < Dctrl; ++i1)
                                for (idx i2 = 0; i2 < Dctrl; ++i2) {
                                    auto coeff_idxes =
                                        coeff_idx_rho(i1, m1, r1, i2, m2, r2);
                                    result(std::get<1>(coeff_idxes),
                                           std::get<2>(coeff_idxes)) =
                                        std::get<0>(coeff_idxes);
                                }
                        }

        return result;
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::applyCTRL()");
}

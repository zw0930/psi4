  if (params_.print_mp2_amps) amp_write();

    tau_build();
    taut_build();
    outfile->Printf("                Solving CC Amplitude Equations\n");
    outfile->Printf("                ------------------------------\n");
    outfile->Printf("  Iter             Energy              RMS        T1Diag      D1Diag    New D1Diag    D2Diag\n");
    outfile->Printf("  ----     ---------------------    ---------   ----------  ----------  ----------   --------\n");
    moinfo_.ecc = energy();
    pair_energies(&emp2_aa, &emp2_ab);
    double last_energy = 0;

    moinfo_.t1diag = diagnostic();
    moinfo_.d1diag = d1diag();
    moinfo_.new_d1diag = new_d1diag();

    moinfo_.d2diag = d2diag();
    update();
    checkpoint();
    for (moinfo_.iter = 1; moinfo_.iter <= params_.maxiter; moinfo_.iter++) {
        sort_amps_sp();
        
        timer_on("F build");
        Fme_build();
        Fae_build();
        Fmi_build();
        if (params_.print & 2) status("F intermediates", "outfile");
        timer_off("F build");

        t1_build();
        if (params_.print & 2) status("T1 amplitudes", "outfile");

        if (params_.wfn == "CC2" || params_.wfn == "EOM_CC2") {
            cc2_Wmnij_build();
            if (params_.print & 2) status("Wmnij", "outfile");

            timer_on("Wmbij build");
            cc2_Wmbij_build();
            if (params_.print & 2) status("Wmbij", "outfile");
            timer_off("Wmbij build");

            timer_on("Wabei build");
            cc2_Wabei_build();
            if (params_.print & 2) status("Wabei", "outfile");
            timer_off("Wabei build");

            timer_on("T2 Build");
            cc2_t2_build();
            if (params_.print & 2) status("T2 amplitudes", "outfile");
            timer_off("T2 Build");
        } else {
            timer_on("Wmbej build");
            Wmbej_build();
            if (params_.print & 2) status("Wmbej", "outfile");
            timer_off("Wmbej build");

            Z_build();
            if (params_.print & 2) status("Z", "outfile");
            Wmnij_build();
            if (params_.print & 2) status("Wmnij", "outfile");

            timer_on("T2 Build");
            t2_build();
            if (params_.print & 2) status("T2 amplitudes", "outfile");
            timer_off("T2 Build");

            if (params_.wfn == "CC3" || params_.wfn == "EOM_CC3") {
                /* step1: build cc3 intermediates, Wabei, Wmnie, Wmbij, Wamef */
                cc3_Wmnij();
                cc3_Wmbij();
                cc3_Wmnie();
                cc3_Wamef();
                cc3_Wabei();

                /* step2: loop over T3's and add contributions to T1 and T2 as you go */
                cc3();
            }
        }
        
        if (!params_.just_residuals) denom(); /* apply denominators to T1 and T2 */

        if (converged(last_energy - moinfo_.ecc)) {
            done = 1;

            tsave();
            tau_build();
            taut_build();
            last_energy = moinfo_.ecc;
            moinfo_.ecc = energy();
            moinfo_.t1diag = diagnostic();
            moinfo_.d1diag = d1diag();
            moinfo_.new_d1diag = new_d1diag();
            moinfo_.d2diag = d2diag();
            sort_amps();
            update();
            outfile->Printf("\n    Iterations converged.\n");

            outfile->Printf("\n");
            amp_write();
            if (params_.analyze != 0) analyze();
            break;
        }
        if (params_.diis) diis(moinfo_.iter);
        tsave();
        tau_build();
        taut_build();
        last_energy = moinfo_.ecc;
        moinfo_.ecc = energy();
        moinfo_.t1diag = diagnostic();
        moinfo_.d1diag = d1diag();
        moinfo_.new_d1diag = new_d1diag();
        moinfo_.d2diag = d2diag();
        update();
        checkpoint();
        // Cast the t1 and t2 back to single-precision for the next iteration if not converged
        cast_dtof();
    }  // end loop over iterations



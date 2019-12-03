psi4/adc/compute_energy.cc|59 col 27| if (options_.get_str("REFERENCE") == "RHF") {
psi4/cc/amps.cc|42 col 34| auto ref = options_.get_str("REFERENCE");
psi4/cc/ccdensity/get_params.cc|59 col 32| // junk = options.get_str("REFERENCE");
psi4/cc/ccdensity/get_params.cc|64 col 48| //  printf("Invalid value of input keyword REFERENCE: %s\n", junk.c_str());
psi4/cc/ccenergy/get_params.cc|72 col 29| junk = options.get_str("REFERENCE");
psi4/cc/ccenergy/get_params.cc|86 col 60| throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);
psi4/cc/cceom/get_params.cc|67 col 45| std::string read_ref = options.get_str("REFERENCE");
psi4/cc/cceom/get_params.cc|78 col 53| std::string read_eom_ref = options.get_str("EOM_REFERENCE");
psi4/cc/cceom/get_params.cc|96 col 61| params.eom_ref = 2; /* run in UHF mode - ignore EOM_REFERENCE */
psi4/cc/cchbar/get_moinfo.cc|79 col 53| std::string read_eom_ref = options.get_str("EOM_REFERENCE");
psi4/cc/cchbar/get_moinfo.cc|80 col 33| //  errcod = ip_string("EOM_REFERENCE", &(read_eom_ref),0);
psi4/cc/ccresponse/get_params.cc|77 col 29| junk = options.get_str("REFERENCE");
psi4/cc/ccresponse/get_params.cc|86 col 60| throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);
psi4/cc/ccresponse/get_params.cc|91 col 35| outfile->Printf("Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n", ref, params.ref);
psi4/cc/cctriples/get_moinfo.cc|98 col 29| junk = options.get_str("REFERENCE");
psi4/cc/cctriples/get_moinfo.cc|110 col 60| throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);
psi4/cctransort/cctransort.cc|88 col 26| if (options.get_str("REFERENCE") == "RHF")
psi4/cctransort/cctransort.cc|90 col 31| else if (options.get_str("REFERENCE") == "ROHF" &&
psi4/cctransort/cctransort.cc|98 col 33| } else if (options.get_str("REFERENCE") == "ROHF")
psi4/cctransort/cctransort.cc|100 col 31| else if (options.get_str("REFERENCE") == "UHF")
psi4/cctransort/cctransort.cc|103 col 57| outfile->Printf("Invalid value of input keyword REFERENCE: %s\n", options.get_str("REFERENCE").c_str());
psi4/cctransort/cctransort.cc|103 col 92| outfile->Printf("Invalid value of input keyword REFERENCE: %s\n", options.get_str("REFERENCE").c_str());
psi4/cctransort/cctransort.cc|396 col 26| if (options.get_str("REFERENCE") == "RHF")
psi4/cctransort/cctransort.cc|399 col 31| else if (options.get_str("REFERENCE") == "ROHF") {
psi4/cctransort/cctransort.cc|408 col 33| } else if (options.get_str("REFERENCE") == "UHF")
psi4/dct/dct_compute.cc|93 col 31| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|329 col 27| if (options_.get_str("REFERENCE") == "RHF") {
psi4/dct/dct_df_tensor.cc|480 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|580 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|678 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|777 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|843 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|908 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|959 col 27| if (options_.get_str("REFERENCE") == "RHF") {
psi4/dct/dct_df_tensor.cc|1135 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|1275 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|1364 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|2869 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_df_tensor.cc|2969 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/dct/dct_gradient.cc|51 col 27| if (options_.get_str("REFERENCE") == "RHF") {
psi4/dct/dct_oo_UHF.cc|212 col 27| if (options_.get_str("REFERENCE") == "RHF")
psi4/detci/compute_mpn.cc|478 col 92| set_scalar_variable("CURRENT CORRELATION ENERGY", Empn2 - scalar_variable("CURRENT REFERENCE ENERGY"));
psi4/detci/compute_mpn.cc|487 col 93| set_scalar_variable("CURRENT CORRELATION ENERGY", Empn2a - scalar_variable("CURRENT REFERENCE ENERGY"));
psi4/detci/compute_mpn.cc|495 col 91| set_scalar_variable("CURRENT CORRELATION ENERGY", Empn - scalar_variable("CURRENT REFERENCE ENERGY"));
psi4/detci/detci.1|25 col 5| .SH REFERENCES
psi4/detci/detci.1|132 col 9| .IP "\fBREFERENCE =\fP \fIstring\fP"
psi4/detci/diag_h.cc|478 col 34| set_scalar_variable("CURRENT REFERENCE ENERGY", CalcInfo_->escf);
psi4/detci/opdm.cc|101 col 38| set_scalar_variable("CURRENT REFERENCE ENERGY", CalcInfo_->escf);
psi4/detci/params.cc|98 col 45| Parameters_->ref_sym = options.get_int("REFERENCE_SYM");
psi4/detci/params.cc|191 col 41| Parameters_->ref = options.get_str("REFERENCE");
psi4/detci/params.cc|205 col 78| throw InputException("Invalid DETCI reference " + Parameters_->ref, "REFERENCE", __FILE__, __LINE__);
psi4/dfmp2/mp2.cc|170 col 30| // if (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "CUHF")
psi4/dfmp2/mp2.cc|170 col 73| // if (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "CUHF")
psi4/dfmp2/wrapper.cc|56 col 26| if (options.get_str("REFERENCE") == "RHF" || options.get_str("REFERENCE") == "RKS") {
psi4/dfmp2/wrapper.cc|56 col 67| if (options.get_str("REFERENCE") == "RHF" || options.get_str("REFERENCE") == "RKS") {
psi4/dfmp2/wrapper.cc|58 col 33| } else if (options.get_str("REFERENCE") == "UHF" || options.get_str("REFERENCE") == "UKS") {
psi4/dfmp2/wrapper.cc|58 col 74| } else if (options.get_str("REFERENCE") == "UHF" || options.get_str("REFERENCE") == "UKS") {
psi4/dfmp2/wrapper.cc|60 col 33| } else if (options.get_str("REFERENCE") == "ROHF") {
psi4/dfocc/dfocc.cc|92 col 35| reference = options_.get_str("REFERENCE");
psi4/dfocc/manager.cc|188 col 25| variables_["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|297 col 29| variables_["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|514 col 25| variables_["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|817 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|1201 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|1591 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|1861 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|2189 col 47| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|2409 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|2722 col 47| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|2923 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|3208 col 47| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager.cc|3407 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|104 col 25| variables_["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|200 col 29| variables_["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|313 col 25| variables_["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|581 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|951 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|1315 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|1573 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|1894 col 47| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|2101 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|2389 col 47| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|2577 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|2829 col 47| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/manager_cd.cc|3018 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/dfocc/qchf.cc|512 col 37| //========================= UHF REFERENCE ==================================================
psi4/dfocc/semi_canonic.cc|132 col 37| //========================= UHF REFERENCE ==================================================
psi4/fnocc/frozen_natural_orbitals.cc|82 col 27| if (options_.get_str("REFERENCE") != "RHF") {
psi4/libscf_solver/frac.cc|75 col 33| if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS"))
psi4/libscf_solver/frac.cc|75 col 75| if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS"))
psi4/libscf_solver/frac.cc|226 col 29| if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS" ||
psi4/libscf_solver/frac.cc|226 col 71| if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS" ||
psi4/libscf_solver/frac.cc|227 col 29| options_.get_str("REFERENCE") == "CUHF"))
psi4/libscf_solver/hf.cc|251 col 81| VBase::build_V(basisset_, functional_, options_, (options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));
psi4/libscf_solver/hf.cc|314 col 47| std::string reference = options_.get_str("REFERENCE");
psi4/libscf_solver/hf.cc|490 col 87| outfile->Printf("                             %4s Reference\n", options_.get_str("REFERENCE").c_str());
psi4/libscf_solver/hf.cc|870 col 47| std::string reference = options_.get_str("REFERENCE");
psi4/libscf_solver/hf.cc|1001 col 51| std::string reference = options_.get_str("REFERENCE");
psi4/libscf_solver/hf.cc|1224 col 47| std::string reference = options_.get_str("REFERENCE");
psi4/libscf_solver/mom.cc|142 col 31| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/libscf_solver/mom.cc|142 col 73| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/libscf_solver/mom.cc|240 col 38| } else if (options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS") {
psi4/libscf_solver/mom.cc|240 col 80| } else if (options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS") {
psi4/libscf_solver/mom.cc|584 col 38| } else if (options_.get_str("REFERENCE") == "ROHF") {
psi4/libscf_solver/rhf.cc|620 col 89| outfile->Printf("   ==> Coupled-Perturbed %s Solver <==\n\n", options_.get_str("REFERENCE").c_str());
psi4/libscf_solver/uhf.cc|603 col 89| outfile->Printf("   ==> Coupled-Perturbed %s Solver <==\n\n", options_.get_str("REFERENCE").c_str());
psi4/mcscf/mcscf.cc|81 col 26| if (options.get_str("REFERENCE") == "RHF" || options.get_str("REFERENCE") == "ROHF" ||
psi4/mcscf/mcscf.cc|81 col 67| if (options.get_str("REFERENCE") == "RHF" || options.get_str("REFERENCE") == "ROHF" ||
psi4/mcscf/mcscf.cc|82 col 26| options.get_str("REFERENCE") == "UHF" || options.get_str("REFERENCE") == "TWOCON") {
psi4/mcscf/mcscf.cc|82 col 67| options.get_str("REFERENCE") == "UHF" || options.get_str("REFERENCE") == "TWOCON") {
psi4/mcscf/mcscf.cc|94 col 43| wfn->set_scalar_variable("CURRENT REFERENCE ENERGY", wfn->energy());
psi4/mcscf/mcscf.cc|97 col 33| } else if (options.get_str("REFERENCE") == "MCSCF") {
psi4/mcscf/mcscf.cc|98 col 29| throw PSIEXCEPTION("REFERENCE = MCSCF not implemented yet");
psi4/mcscf/scf.cc|82 col 27| if (options_.get_str("REFERENCE") == "RHF") {
psi4/mcscf/scf.cc|86 col 34| } else if (options_.get_str("REFERENCE") == "ROHF") {
psi4/mcscf/scf.cc|90 col 34| } else if (options_.get_str("REFERENCE") == "UHF") {
psi4/mcscf/scf.cc|94 col 34| } else if (options_.get_str("REFERENCE") == "TWOCON") {
psi4/mcscf/scf_save_info.cc|104 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = total_energy;
psi4/mrcc/mrcc.cc|599 col 26| if (options.get_str("REFERENCE") != "RHF")
psi4/mrcc/mrcc.cc|605 col 26| if (options.get_str("REFERENCE") == "UHF") restricted = false;
psi4/mrcc/mrcc.cc|607 col 36| if (pertcc && options.get_str("REFERENCE") == "ROHF") {
psi4/mrcc/mrcc.cc|726 col 26| if (options.get_str("REFERENCE") == "ROHF") canonical = false;
psi4/mrcc/mrcc.cc|728 col 36| if (pertcc && options.get_str("REFERENCE") == "ROHF") {
psi4/mrcc/mrcc.cc|746 col 26| if (options.get_str("REFERENCE") == "ROHF") {
psi4/occ/manager.cc|231 col 29| variables_["CURRENT REFERENCE ENERGY"] = Escf;
psi4/occ/manager.cc|355 col 25| variables_["CURRENT REFERENCE ENERGY"] = Escf;
psi4/occ/manager.cc|707 col 47| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/occ/manager.cc|855 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/occ/manager.cc|1076 col 47| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/occ/manager.cc|1187 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
psi4/occ/manager.cc|1432 col 47| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
psi4/occ/manager.cc|1536 col 43| Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
psi4/occ/occwave.cc|95 col 35| reference = options_.get_str("REFERENCE");
psi4/occ/semi_canonic.cc|176 col 12| // UHF REFERENCE
psi4/sapt/wrapper.cc|65 col 30| if (options.get_str("REFERENCE") == "RHF") {
psi4/sapt/wrapper.cc|76 col 34| if (options.get_str("REFERENCE") == "ROHF") {
psi4/scfgrad/scf_grad.cc|127 col 27| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|127 col 69| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|144 col 27| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|144 col 69| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|154 col 31| if (options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|331 col 27| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|331 col 69| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|348 col 27| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|348 col 69| if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|362 col 34| // if (options_.get_str("REFERENCE") == "RKS") {
psi4/scfgrad/scf_grad.cc|1044 col 27| if (options_.get_str("REFERENCE") == "RHF") {
read_options.cc|279 col 26| options.add_str("REFERENCE", "RHF", "RHF ROHF");
read_options.cc|551 col 26| options.add_int("REFERENCE_SYM", -1);
read_options.cc|1063 col 26| options.add_str("REFERENCE", "RHF", "UHF RHF ROHF");
read_options.cc|1219 col 26| options.add_str("REFERENCE", "RHF", "RHF ROHF UHF CUHF RKS UKS");
read_options.cc|1629 col 26| options.add_str("REFERENCE", "RHF");
read_options.cc|1647 col 26| options.add_str("REFERENCE", "RHF");
read_options.cc|1659 col 26| options.add_str("REFERENCE", "RHF");
read_options.cc|1768 col 26| options.add_str("REFERENCE", "RHF", "RHF");
read_options.cc|1794 col 30| options.add_str("EOM_REFERENCE", "RHF");
read_options.cc|1813 col 26| options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
read_options.cc|1815 col 30| options.add_str("EOM_REFERENCE", "RHF", "RHF ROHF UHF");
read_options.cc|1936 col 26| options.add_str("REFERENCE", "RHF");
read_options.cc|2002 col 28| //    options.add_str("REFERENCE", "RHF");
read_options.cc|2024 col 26| options.add_str("REFERENCE", "RHF", "RHF ROHF UHF TWOCON MCSCF GENERAL");
read_options.cc|2084 col 26| options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
read_options.cc|4104 col 32| options.add_str("CFOUR_REFERENCE", "RHF", "RHF UHF ROHF TCSCF CASSCF");
psi4/libmints/bessel.h|35 col 9| REFERENCES:
psi4/libmints/ecpint.h|36 col 9| REFERENCES:
psi4/libmints/gaussquad.h|34 col 9| REFERENCES:

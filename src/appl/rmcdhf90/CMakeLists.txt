add_executable(rmcdhf
    rscfvu.f90

    cofpot.f90 cofpot_I.f90
    consis.f90 consis_I.f90
    csfwgt.f90 csfwgt_I.f90
    dacon.f90 dacon_I.f90
    dampck.f90 dampck_I.f90
    dampor.f90 dampor_I.f90
    defcor.f90 defcor_I.f90
    dsubrs.f90 dsubrs_I.f90
    eigen.f90 eigen_I.f90
    endsum.f90 endsum_I.f90
    engoutgg.f90 engoutgg_I.f90
    estim.f90 estim_I.f90
    fco.f90 fco_I.f90
    gco.f90 gco_I.f90
    getald.f90 getald_I.f90
    getaldwt.f90 getaldwt_I.f90
    getold.f90 getold_I.f90
    getoldwt.f90 getoldwt_I.f90
    getscd.f90 getscd_I.f90
    hmout.f90 hmout_I.f90
    improv.f90 improv_I.f90
    in.f90 in_I.f90
    ispar.f90 ispar_I.f90
    itjpo.f90 itjpo_I.f90
    lagcon.f90 lagcon_I.f90
    lodcsh2GG.f90 lodcsh2GG_I.f90
    maneig.f90 maneig_I.f90
    matrix.f90 matrix_I.f90
    maxarr.f90 maxarr_I.f90
    newco.f90 newco_I.f90
    newe.f90 newe_I.f90
    orbout.f90 orbout_I.f90
    orthor.f90 orthor_I.f90
    orthy.f90 orthy_I.f90
    out.f90 out_I.f90
    outbnd.f90 outbnd_I.f90
    prtrsl.f90 prtrsl_I.f90
    prwf.f90 prwf_I.f90
    scf.f90 scf_I.f90
    setcof.f90 setcof_I.f90
    setcsl.f90 setcsl_I.f90
    setdbg.f90 setdbg_I.f90
    setham.f90 setham_I.f90
    setlag.f90 setlag_I.f90
    setmcp.f90 setmcp_I.f90
    setmix.f90 setmix_I.f90
    setsum.f90 setsum_I.f90
    setxuv.f90 setxuv_I.f90
    setxv.f90 setxv_I.f90
    setxz.f90 setxz_I.f90
    solve.f90 solve_I.f90
    strsum.f90 strsum_I.f90
    xpot.f90 xpot_I.f90
    ypot.f90 ypot_I.f90

    mpi_s.f90
)
target_link_libraries_Fortran(rmcdhf PUBLIC mod 9290 dvd90)
install(TARGETS rmcdhf DESTINATION bin/)

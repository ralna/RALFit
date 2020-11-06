add_library(ral_nlls_lapack SHARED
         dbdsqr.f dgebd2.f dgebrd.f dgebrd.f
         dgebrd.f dgelq2.f dgelqf.f dgels.f
         dgeqr2.f dgeqrf.f dgesvd.f dgesv.f
         dgetrf2.f dgetrf.f dgetrs.f dggbak.f
         dggbal.f dggev.f dggev.f dgghrd.f
         dhgeqz.f disnan.f dlabad.f dlabrd.f
         dlacpy.f dladiv.f dlae2.f dlaebz.f
         dlaev2.f dlag2.f dlagtf.f dlagts.f
         dlaisnan.f dlaln2.f dlange.f dlanhs.f
         dlanst.f dlansy.f dlapy2.f dlapy3.f
         dlarfb.f dlarf.f dlarfg.f dlarft.f
         dlarnv.f dlartg.f dlaruv.f dlas2.f
         dlascl.f dlaset.f dlasq1.f dlasq2.f
         dlasq3.f dlasq4.f dlasq5.f dlasq6.f
         dlasr.f dlasrt.f dlassq.f dlasv2.f
         dlaswp.f dlatrd.f dlauu2.f dlauum.f
	 dorg2l.f dorg2r.f
         dorgbr.f dorgl2.f dorglq.f dorgql.f
         dorgqr.f dorgtr.f dorm2l.f dorm2r.f
         dormbr.f dorml2.f dormlq.f dormql.f 
         dormqr.f dormtr.f dposv.f dpotrf2.f
         dpotrf.f dpotri.f dpotrs.f dstebz.f
	 dstein.f
         dsteqr.f dsterf.f dsyev.f dsyevx.f
         dsytd2.f dsytrd.f dtgevc.f dtrtri.f
	 dtrti2.f dtrtrs.f
         ieeeck.f iladlc.f iladlr.f ilaenv.f
         iparmq.f dlamch.f)

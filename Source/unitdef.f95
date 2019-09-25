!   unit definitions.



module unitdef
implicit none


!	unitalpha is used to save alpha and cholnhess.
!	unitbeta is used to print beta parameters.
!	unitbetacov is the unit number for the estimated asymptotic covariance matrix of the beta estimates.
!	unitbetacov_complex is the unit number for the complex-sampling estimated covariance matrix of the beta estimates.
!	unitbetacov_louis is the unit number for the Louis estimated covariance matrix of the beta estimates.
!	unitbetacov_sandwich is the unit number for the sandwich estimated covariance matrix of the beta estimates.
!	uniteap is used for output of eap data on latent vectors.
!	uniteapscale is used for output of eap data on scaled scores.
!	uniteapskill is used for output of eap data for transformed latent vectors.
!	uniteaptrans is used for output of eap data on transformations of the latent vectors.
!	uniteapwt is used for output of eap data on weighted sums.
!	unitfreq is used for output of distractor frequencies.
!	unitgrad is used for output of gradients.
!   unitguess is the unit for output of guessing tests.
!	unitinfo is the unit number for the information summary.
!	unitinfodist is the unit number for the information summary for multiple choices.
!	unitirf is the unit for item response functions.
!	unititeration is the unit number for iterations.
!	unititerationstart is the unit number for starting iterations.
!	unitmargin is the unit number for marginal distributions.
!	unitmargins2 is the unit number for two-way cross-products of item scores.
!	unitmarginwtsum is the unit number of marginal distributions of weighted sums.
!	unitmargin2 is the unit number for two-way marginal distributions.
!	unitmeans is the unit number for means and covariances of underlying latent vectors.
!   unitml is the unit number for maximum likelihood.
!	unitmp is the unit number for maximum posterior likelihood.
!   unitobservedscale is the unit number for individual observed scale scores.
!   unitobsscale is the unit number for the observed scale score summaries.
!   unitobsscalerel is the unit number for the observed scale score reliability.
!	unitoutdata is the unit number for the data summary.
!	unitparam is the unit number for the gamma parameters.
!	unitparamcov is the unit number for the estimated asymptotic covariance matrix of the gamma estimates.
!	unitparamcov_complex is the unit number for the complex-sampling estimated covariance matrix of the gamma estimates.
!	unitparamcov_louis is the unit number for the Louis estimated covariance matrix of the gamma estimates.
!	unitparamcov_sandwich is the unit number for the sandwich estimated covariance matrix of the gamma estimates.
!	unitpost is the unit number for posterior distributions.
!	unitpreditem is the unit number for the sums and averages of products of category indicators and weighted sums.
!	unitprob is the unit number for marginal probabilities of observed response vectors.
!	unitprobdist is the unit number for marginal probabilities of observed distractor vectors.
!   unitpwt is for output of individual significance levels based on weighted sums.
!   unitpwtdist is for output of individiual significance levels based on weighted distractor sums.
!	unitrel is for output of reliabilities of eap's for latent vectors.
!	unitrelskill is for output of reliabilities of eap's for transformed latent vectors.
!	unitrelwt is for output of reliabilities of eap's for weighted sums.
!   unitreswt is for output of individual residuals for weighted sums.
!   unitreswtdist is for output of individual residuals for weighted distractor sums.
!	unitscalerel is for reliability of scaled scores.
!	unitthetaitem is the unit number for the sums and averages of products of category indicators and transformed latent vectors.
!	unittitle is the unit number for the title.
!	unitwtitem is the unit number for the sums and averages of products of category indicators and predictors.
!	unitwtitemdist is the unit number for the sums and averages of products of category indicators for distractors and predictors.
!	unitcounts counts print specifications tied to physical units.

integer::unitalpha=10
integer::unitbeta=10,unitbetacov=10,unitbetacov_complex=10,unitbetacov_louis=10,unitbetacov_sandwich=10
integer::uniteap=10,uniteapscale=10,uniteapskill=10,uniteaptrans=10,uniteapwt=10
integer::unitfreq=10
integer::unitgrad=10,unitguess=10
integer::unitinfo=10,unitinfodist=10,unitirf=10,unititeration=10,unititerationstart=10
integer::unitmargin=10,unitmargins2=10,unitmarginwtsum=10,unitmargin2=10,unitmeans=10,unitml=10
integer::unitmp=10
integer::unitobservedscale=10,unitobsscale=10,unitobsscalerel=10
integer::unitoutdata=10
integer::unitparam=10,unitparamcov=10,unitparamcov_complex=10,unitparamcov_louis=10,unitparamcov_sandwich=10
integer::unitpost=10,unitpreditem=10,unitprob=10,unitprobdist=10,unitpwt=10,unitpwtdist=10
integer::unitrel=10,unitrelskill=10,unitreltrans=10,unitrelwt=10,unitreswt=10,unitreswtdist=10
integer::unitscalerel=10
integer::unitthetaitem=10,unittitle=10
integer::unitwtitem=10,unitwtitemdist=10
!   Counter for number of output runs on each unit number.
integer::unitcount(10:99)=0

!	printalpha is used to save alpha and cholnhess.
!	printbeta indicates if the beta statistics are to be provided.
!	printbetacov indicates if the estimated asymptotic covariance matrix is provided for beta estimates.
!	printbetacov_complex indicates if the complex estimated asymptotic covariance matrix is provided for beta estimates.
!	printbetacov_louis indicates if the Louis estimated asymptotic covariance matrix is provided for beta estimates.
!	printbetacov_sandwich indicates if the sandwich estimated asymptotic covariance matrix is provided for beta estimates.
!	printeap indicates if eap data are desired for latent vectors.
!	printeapscale is .true. if eap data are desired for scaled scores.
!	printeapskill is .true. if eap data are sought for transformed latent vectors.
!	printeaptrans is .true. if eap data are sought for transformations of latent vectors.
!	printeapwt is .true. if eap data for weighted sums is desired.
!	printent indicates if  the entropy statistics are to be provided.
!	printentdist indicates if  the entropy statistics are to be provided for distractors.
!	printfreq indicates if distractor frequencies are to be provided.
!	printgrad indicates that gradients are to be printed.
!   printguess indicates that guessing tests are to be supplied.
!	printirf indicates that item response functions are to be printed.
!	printirfres indicates that residuals for item response functions are to be printed.
!	printmargin is .true. if marginal distributions are printed.
!	prinmarginres is .true. if residuals are printed for marginal distributions.
!	printmargins2 is .true. if cross products are printed for item scors.
!	prinmargins2res is .true. if residuals are printed for cross products of item scores.
!	printmarginwtsum is .true. if marginal distributions of weighted sums are to be printed.
!	printmarginwtsumres is .true. if residuals for marginal distributions of weighted sums are to be printed.
!	printmargin2 is .true. if two-way marginal distributions are printed.
!	prinmargin2res is .true. if residuals are printed for two-way marginal distributions.
!	printmeans is .true. if means and covariance matrices of underlying latent vectors are requested.
!   printml is .true. if maximum likelihood proficiency estimates are provided.
!	printmp is .true. if maximum posterior likelihood estimates are provided.
!   printobservedscale is .true. if individual observed scale scores are desired.
!   printobsscale is .true. if observed scale score summeries are desired.
!   printobsscalerel is .true. if reliability of observed scale score is provided.
!   printobsscaleres is .true. if residuals for observed scale scores are provided.
!	printparam indicates if the gamma statistics are to be provided.
!	printparamcov indicates if the estimated asymptotic covariance matrix is provided for gamma estimates.
!	printparamcov_complex indicates if the complex estimated asymptotic covariance matrix is provided for gamma estimates.
!	printparamcov_louis indicates if the Louis estimated asymptotic covariance matrix is provided for gamma estimates.
!	printparamcov_sandwich indicates if the sandwich estimated asymptotic covariance matrix is provided for gamma estimates.
!	printpost is for printing posterior distributions.
!	printpreditem is .true. if sums and averages are sought for products of category indicators and predictors.
!	printpreditemres is .true. if residuals are sought for sums and averages of products of category indicators
!		and predictors.
!	printprob leads to printing of probabilities for each individual.
!	printprobdist leads to printing of distractor probabilities for each individual.
!	printprog indicates that regular iteration progress is to be printed to unititeration.
!	printprogstart indicates that starting iteration progress is to be printed to unititerationstart.
!	printprogstartstd indicates that starting iteration progress is to go to standard output.
!	printprogstd indicates that regular iteration progress is to go to standard output.
!   printpwt indicates output of individual significance levels based on weighted sums.
!   printpwtdist indicates output of individiual significance levels based on weighted distractor sums.
!	printrel indicates if reliabilities of eap's are to be  obtained for latent vectors.
!	printrelskill indicates if reliabilities of eap's are to be obtained for transformed latent vectors.
!	printreltrans indicates if reliabilities of eap's are to be obtained for transformations of latent vectors.
!	printrelwt indicates if reliabilities of eap's for weighted sums are to be  obtained.
!   printreswt indicates if individual residuals for weighted sums are to be obtained.
!   printreswtdist indicates if individual residuals for weighted distractor sums are to be obtained.
!	printscalerel indicates if reliabilities of scale scores are to be obtained.
!	printthetaitem is .true. if totals and averages are sought for category indicators and transformed latent vectors.
!	printthetaitemres is .true. if residuals are sought totals and averages  category indicators and transformed latent predictors.
!	printwtitem is .true. if totals and averages are sought for category indicators and weighted sums.
!	printwtitemres is .true. if residuals are sought for totals and averages of category indicators and weighted sums.
!	printwtitemdist is .true. if totals and averages are saught for distractor indicators and weighted sums.
!	printwtitemdistres is .true. if residuals are sought for totals and averages of distractor indicators and weighted sums.
logical::printalpha=.FALSE.
logical::printbeta=.FALSE.,printbetacov=.FALSE.,printbetacov_complex=.FALSE.,printbetacov_louis=.FALSE.,printbetacov_sandwich=.FALSE.
logical::printeap=.FALSE.,printeapscale=.FALSE.,printeapskill=.FALSE.,printeaptrans=.FALSE.,printeapwt=.FALSE.,printent=.TRUE.,printentdist=.FALSE.
logical::printfreq=.FALSE.
logical::printgrad=.FALSE.,printguess=.FALSE.
logical::printirf=.FALSE.,printirfres=.FALSE.
logical::printmargin=.FALSE.,printmarginres=.FALSE.,printmargins2=.FALSE.,printmargins2res=.FALSE.
logical::printmarginwtsum=.FALSE.,printmarginwtsumres=.FALSE.
logical::printmargin2=.FALSE.,printmargin2res=.FALSE.,printmeans=.FALSE.,printml=.FALSE.
logical::printmp=.FALSE.
logical::printobservedscale=.FALSE.,printobsscale=.FALSE.,printobsscalerel=.FALSE.,printobsscaleres=.FALSE.
logical::printparam=.TRUE.,printparamcov=.FALSE.,printparamcov_louis=.FALSE.,printparamcov_sandwich=.FALSE.,printparamcov_complex=.FALSE.
logical::printpost=.FALSE.,printpreditem=.FALSE.,printpreditemres=.FALSE.,printprob=.FALSE.,printprobdist=.FALSE.
logical::printprog=.TRUE.,printprogstart=.TRUE.,printprogstartstd=.TRUE.,printprogstd=.TRUE.,printpwt=.FALSE.,printpwtdist=.FALSE.
logical::printrel=.FALSE.,printrelskill=.FALSE.,printreltrans=.FALSE.,printrelwt=.FALSE.,printreswt=.FALSE.,printreswtdist=.FALSE.
logical::printscalerel=.FALSE.
logical::printthetaitem=.FALSE.,printthetaitemres=.FALSE.
logical::printwtitem=.FALSE.,printwtitemres=.FALSE.,printwtitemdist=.FALSE.,printwtitemdistres=.FALSE.


end module unitdef

\name{spatstat.core-internal} 
\title{Internal spatstat.core functions}
\alias{[.localpcfmatrix}
\alias{[.rat}
\alias{accumulateStatus}
\alias{active.interactions}
\alias{adaptcoef}
\alias{adjust.ratfv}
\alias{affine.msr}
\alias{ang2rad}
\alias{areadelta2}
\alias{as.data.frame.bw.optim}
\alias{as.data.frame.fv}
\alias{as.isf}
\alias{as.ppm.rppm}
\alias{assemble.plot.objects}       
\alias{augment.msr}
\alias{bandwidth.is.infinite}
\alias{BartCalc}
\alias{bermantestCalc}
\alias{bermantestEngine}
\alias{bigvaluerule}
\alias{bind.ratfv}
\alias{blankcoefnames}
\alias{bt.frame}
\alias{bw.optim}
\alias{calc.DR}
\alias{calc.NNIR}
\alias{calc.SAVE}
\alias{calc.SIR}
\alias{calc.TSE}
\alias{cannot.update}
\alias{censtimeCDFest}
\alias{check.separable}
\alias{check.testfun}
\alias{circticks}
\alias{clarkevansCalc}
\alias{coef.summary.kppm}
\alias{coef.summary.ppm}
\alias{coef.summary.slrm}
\alias{coef.vblogit} 
\alias{compatible.rat}
\alias{compileCDF}
\alias{condSimCox}
\alias{conform.ratfv}
\alias{CressieReadStatistic}
\alias{CressieReadSymbol}
\alias{CressieReadName}
\alias{cutoff2Dkernel}
\alias{CVforPCF}
\alias{damaged.ppm}
\alias{data.mppm}
\alias{deltasuffstat}
\alias{Deviation}
\alias{dfbetas.ppmInfluence}
\alias{densitycrossEngine}
\alias{densitypointsEngine}
\alias{diagnose.ppm.engine}
\alias{digestCovariates}
\alias{digital.volume}
\alias{dim.fasp}               
\alias{dim.msr}
\alias{dimnames.fasp}               
\alias{dimnames<-.fasp}
\alias{dimnames.msr}
\alias{distributecbind}
\alias{doMultiStraussHard}
\alias{dppDpcf}
\alias{dppmFixAlgorithm}
\alias{dppmFixIntensity}
\alias{envelopeEngine}
\alias{envelopeProgressData}
\alias{envelopeTest}
\alias{envelope.hasenvelope}
\alias{envelope.matrix}
\alias{equalpairs}          
\alias{evalCovar}
\alias{evalCovar.ppm}
\alias{evalCovar.slrm}
\alias{evalCovariate}
\alias{evalInteraction}
\alias{evalInterEngine}
\alias{evalPairPotential}
\alias{evaluate2Dkernel}
\alias{ExpSmoothLog}
\alias{expandDot}
\alias{extractAIC.slrm}
\alias{extractAtomicQtests}
\alias{family.glmmPQL} 
\alias{family.vblogit} 
\alias{fasp}
\alias{f3engine}
\alias{f3Cengine}
\alias{fill.coefs}
\alias{findCovariate}
\alias{fii}
\alias{findbestlegendpos}
\alias{findcbind}
\alias{flatfname}
\alias{flipxy.msr}
\alias{forbid.logi}
\alias{FormatFaspFormulae}
\alias{fvexprmap}
\alias{fvlabels}
\alias{fvlabels<-}
\alias{fvlabelmap}
\alias{fvlegend}
\alias{g3engine}
\alias{g3Cengine}
\alias{getdataname}
\alias{getglmdata}
\alias{getglmfit}
\alias{getglmsubset}
\alias{getppmdatasubset}
\alias{getppmOriginalCovariates}
\alias{getSumFun}
\alias{geyercounts}
\alias{geyerdelta2}
\alias{GLMpredict}
\alias{good.correction.K}
%\alias{gridadjacencymatrix} %DoNotExport
\alias{hackglmmPQL}
\alias{hasenvelope}
\alias{hasglmfit}
\alias{hierarchicalordering}
\alias{hiermat}
\alias{ho.engine}
\alias{illegal.iformula}
\alias{implemented.for.K}
\alias{impliedpresence}
\alias{impliedcoefficients}
\alias{influence.ppmInfluence}
\alias{instantiate.interact}
\alias{interactionfamilyname}
\alias{intermaker}
\alias{is.atomicQtest}
\alias{is.expandable.ppm}
\alias{is.interact}
\alias{is.marked.mppm}    
\alias{is.marked.msr}    
\alias{is.marked.slrm}
\alias{is.mppm}
\alias{is.multitype.mppm}    
\alias{is.multitype.msr}    
\alias{is.multitype.slrm}    
\alias{is.poisson.mppm}
\alias{is.scov}
\alias{k3engine}
\alias{Kborder.engine}
\alias{Knone.engine}
\alias{Krect.engine}
\alias{Kount}
\alias{Kwtsum}               
\alias{Kpcf.kppm}
\alias{Kmodel.slrm}
\alias{killinteraction}
\alias{km.rs.opt}
\alias{kppmComLik}
\alias{kppmMinCon}
\alias{kppmPalmLik}
\alias{kppmCLadap}
\alias{labels.ppm}
\alias{leverage.ppmInfluence}
\alias{localKengine}
\alias{localKmultiEngine}
\alias{localpcfengine}
\alias{localpcfmatrix}
\alias{logi.engine}
\alias{logLik.vblogit}
\alias{lookup2DkernelInfo}
\alias{LurkEngine}
\alias{makefvlabel}
\alias{mapInterVars}
\alias{maskLaslett}
\alias{match2DkernelName}
\alias{match.kernel}
\alias{mctestSigtraceEngine}
\alias{meanlistfv}
\alias{model.se.image}
\alias{modelFrameGam}
\alias{mpl.engine}
\alias{mpl.get.covariates}
\alias{mpl.prepare}
\alias{mpl.usable}
\alias{names<-.fv}
\alias{newformula}
\alias{newstyle.coeff.handling}
\alias{nncleanEngine}
\alias{nndcumfun}
\alias{no.trend.ppm}
\alias{objsurfEngine}
\alias{optimConverged}
\alias{optimStatus}
\alias{optimNsteps}
\alias{outdated.interact}
\alias{oversize.quad}    
\alias{pairs.listof}
\alias{pairs.solist}
\alias{PairPotentialType}
\alias{partialModelMatrix}
\alias{pcf3engine}
\alias{pcfmulti.inhom}
\alias{pcfmodel.slrm}
\alias{ploterodewin}
\alias{ploterodeimage}
\alias{plot.addvar}
\alias{plot.bw.frac}
\alias{plot.bw.optim}
\alias{plot.localpcfmatrix}
\alias{plot.lurk}
\alias{plot.minconfit}
\alias{plot.parres}
\alias{plot.plotpairsim}
\alias{plot.qqppm}
\alias{plot.spatialcdf}
\alias{PoisSaddle}
\alias{PoisSaddleArea}
\alias{PoisSaddleGeyer}
\alias{PoisSaddlePairwise}
\alias{polyLaslett}
\alias{PPMmodelmatrix}
\alias{ppm.default}
\alias{ppmCovariates}
\alias{ppmDerivatives}
\alias{ppmInfluenceEngine}
\alias{predict.profilepl}
\alias{predict.vblogit}
\alias{prefixfv}
\alias{printStatus}
\alias{printStatusList}
\alias{print.addvar}    
\alias{print.bt.frame}
\alias{print.bw.frac}
\alias{print.bw.optim}
\alias{print.diagppm}
\alias{print.densityfun}
\alias{print.detpointprocfamily}
\alias{print.detpointprocfamilyfun}
\alias{print.envelope}
\alias{print.fasp}       
\alias{print.fv}       
\alias{print.fvfun}       
\alias{print.hasenvelope}       
\alias{print.hierarchicalordering}
\alias{print.influence.ppm}       
\alias{print.interact}       
\alias{print.intermaker}       
\alias{print.isf}
\alias{print.laslett}
\alias{print.leverage.ppm}
\alias{print.localpcfmatrix}
\alias{print.lurk}
\alias{print.minconfit}
\alias{print.mppm}
\alias{print.msr}
\alias{print.parres}
\alias{print.plotpairsim}
\alias{print.plotppm}
\alias{print.profilepl}
\alias{print.quadrattest}
\alias{print.qqppm}
\alias{print.rat}
\alias{print.rppm}
\alias{print.Smoothfun}       
\alias{print.summary.mppm}
\alias{print.summary.ssf}
\alias{print.summary.slrm}
\alias{print.vblogit}
\alias{quad.mppm}
\alias{quadBlockSizes}
\alias{quadrat.testEngine}
\alias{ratfv}
\alias{reach.slrm}
\alias{rebadge.as.crossfun}
\alias{rebadge.as.dotfun}
\alias{rebadge.fv}
\alias{rebadgeLabels}
\alias{reconcile.fv}
\alias{reduceformula}
\alias{reincarnate.interact}
\alias{RelevantDeviation}
\alias{rename.fv}
\alias{rescale.msr}
\alias{resid4plot}
\alias{resid1plot}
\alias{resid1panel}
\alias{resolve.2D.kernel}
\alias{resolveEinfo}
\alias{resolve.foxall.window}
\alias{resolve.lambda}
\alias{resolve.lambda.cross}
\alias{rhohatEngine}
\alias{rhohatCalc}
\alias{rmax.Rigid}
\alias{rmax.rule}
\alias{rocData}
\alias{rocModel}
\alias{roseContinuous}
\alias{rotate.msr}
\alias{safeFiniteValue}
\alias{safePositiveValue}
\alias{scalardilate.msr}
\alias{scanmeasure}
\alias{scanmeasure.ppp}
\alias{scanmeasure.im}
\alias{scanBinomLRTS}
\alias{scanPoisLRTS}
\alias{second.moment.calc}
\alias{second.moment.engine}
\alias{sewpcf}
\alias{sewsmod}
\alias{shift.influence.ppm}              
\alias{shift.leverage.ppm}              
\alias{shift.quadrattest}              
\alias{shift.msr}              
\alias{signalStatus}
\alias{simulate.profilepl}
\alias{simulrecipe}              
\alias{slr.prepare}
\alias{slrAssemblePixelData}
\alias{slrmInfluence}
\alias{Smooth.solist}
\alias{smoothcrossEngine}              
\alias{smoothpointsEngine}
\alias{spatstatDPPModelInfo}
\alias{spatialCDFframe}
\alias{spatialCDFtest}
\alias{spatialCDFtestCalc}
\alias{sphere.volume}
\alias{splitHybridInteraction}
\alias{sp.foundclass}
\alias{sp.foundclasses}
\alias{strausscounts}
\alias{suffloc}
\alias{suffstat.generic}
\alias{suffstat.poisson}
\alias{summary.envelope}
\alias{summary.mppm}
\alias{summary.msr}
\alias{summary.profilepl}
\alias{summary.vblogit}
\alias{tweak.coefs}
\alias{tweak.fv.entry}
\alias{tweak.ratfv.entry}
\alias{twostage.test}
\alias{twostage.envelope}
\alias{unitname.msr}
\alias{unitname<-.msr}
\alias{update.ippm}
\alias{update.msr}
\alias{validate2Dkernel}
\alias{validate.angles}        
\alias{validate.weights}        
\alias{vanilla.fv}
\alias{varcountEngine}
%\alias{vblogit} %DoNotExport
%\alias{vblogit.fmla} %DoNotExport
\alias{versionstring.interact}
\alias{versionstring.ppm}
\alias{weightedclosepairs}
\alias{windows.mppm}
\alias{X2testEngine}




%%%%%%%
\description{
  Internal spatstat.core functions.
}
\usage{
\method{[}{localpcfmatrix}(x, i, \dots)
\method{[}{rat}(x, \dots)
accumulateStatus(x, stats)
active.interactions(object)
adaptcoef(new.coef, fitcoef, drop)
adjust.ratfv(f, columns, numfactor, denfactor)
\method{affine}{msr}(X, mat, vec, \dots)
ang2rad(ang, unit, start, clockwise)
areadelta2(X, r, \dots, sparseOK)
assemble.plot.objects(xlim, ylim, \dots, lines, polygon)
as.isf(object)
\method{as.ppm}{rppm}(object)
\method{as.data.frame}{bw.optim}(x, \dots)
\method{as.data.frame}{fv}(x, \dots)
augment.msr(x, \dots, sigma, recompute)
bandwidth.is.infinite(sigma)
BartCalc(fY, fK)
bermantestCalc(fram, which, alternative, \dots)
bermantestEngine(model, covariate, which, alternative, \dots,
                 modelname, covname, dataname)
bigvaluerule(objfun, objargs, startpar, \dots)
bind.ratfv(x, numerator, denominator, labl, desc, preferred,
           ratio, quotient)
blankcoefnames(x)
bt.frame(Q, trend, interaction, \dots, covariates,
         correction, rbord, use.gam, allcovar)
bw.optim(cv, h, iopt, \dots, cvname, hname, criterion,
         optimum, warnextreme, hargnames, unitname)
calc.DR(COV, z, Dim)
calc.NNIR(COV, z, pos, Dim)
calc.SAVE(COV, z, Dim)
calc.SIR(COV, z)
calc.TSE(COV, z, pos, Dim1, Dim2)
cannot.update(\dots)
censtimeCDFest(o, cc, d, breaks, \dots,
     KM, RS, HAN, RAW, han.denom, tt, pmax)
check.separable(dmat, covname, isconstant, fatal)
check.testfun(f, f1, X)
circticks(R, at, unit, start, clockwise, labels)
clarkevansCalc(X, correction, clipregion, working)
\method{coef}{summary.kppm}(object, \dots)
\method{coef}{summary.ppm}(object, \dots)
\method{coef}{summary.slrm}(object, \dots)
\method{coef}{vblogit}(object, \dots) 
\method{compatible}{rat}(A, B, \dots) 
compileCDF(D, B, r, \dots, han.denom, check)
condSimCox(object, nsim, \dots, window, n.cond, w.cond,
           giveup, maxchunk, verbose, drop)
conform.ratfv(x)
CressieReadStatistic(OBS,EXP,lambda,normalise,named)
CressieReadSymbol(lambda)
CressieReadName(lambda)
cutoff2Dkernel(kernel, sigma, varcov, \dots, scalekernel, cutoff, fatal)
CVforPCF(bw, stuff)
damaged.ppm(object)
data.mppm(x)
deltasuffstat(model, \dots,
              restrict, dataonly, sparseOK, quadsub,
              force, warn.forced, verbose, use.special)
Deviation(x, ref, leaveout, n, xi)
\method{dfbetas}{ppmInfluence}(model, \dots)
densitycrossEngine(Xdata, Xquery, sigma, \dots,
                    kernel, scalekernel,
                    weights, edge, varcov,
                    diggle, sorted, cutoff)
densitypointsEngine(x, sigma, \dots,
                    kernel, scalekernel,
                    weights, edge, varcov,
                    leaveoneout, diggle, sorted, spill, cutoff,
                    debug)
diagnose.ppm.engine(object, \dots, type, typename, opt,
                         sigma, rbord, compute.sd, compute.cts,
                         envelope, nsim, nrank,
                         rv, oldstyle, splineargs, verbose)
digestCovariates(\dots, W)
digital.volume(range, nval, vside)
\method{dim}{fasp}(x)
\method{dim}{msr}(x)
\method{dimnames}{fasp}(x)
\method{dimnames}{fasp}(x) <- value
\method{dimnames}{msr}(x)
distributecbind(x)
doMultiStraussHard(iradii, hradii, types)
dppDpcf(model, \dots)
dppmFixIntensity(DPP, lambda, po)
dppmFixAlgorithm(algorithm, changealgorithm, clusters, startpar)
envelopeEngine(X, fun, simul,
           nsim, nrank, \dots, funargs, funYargs,
           verbose, clipdata, 
           transform, global, ginterval, use.theory,
           alternative, scale, clamp,
           savefuns, savepatterns, saveresultof,
           weights,
           nsim2, VARIANCE, nSD,
           Yname,
           maxnerr, rejectNA, silent, maxerr.action,
           internal, cl,
           envir.user, expected.arg, do.pwrong,
           foreignclass, collectrubbish)
envelopeProgressData(X, fun, \dots, exponent,
                     alternative, leaveout, scale, clamp,
                     normalize, deflate, rmin, 
                     save.envelope, savefuns, savepatterns)
envelopeTest(X, \dots, exponent, alternative,
            rinterval, leaveout, scale, clamp, tie.rule,
            interpolate, save.interpolant,
            save.envelope, savefuns, savepatterns,
            Xname, badXfatal, verbose)
\method{envelope}{hasenvelope}(Y, \dots, Yname)
\method{envelope}{matrix}(Y, \dots, rvals, observed, theory, funX,
  nsim, nsim2, jsim, jsim.mean,
  type, alternative, scale, clamp, csr, use.theory, nrank, ginterval, nSD,
  savefuns, check, Yname, do.pwrong, weights, precomputed, gaveup)
equalpairs(U, X, marked=FALSE)
evalCovar(model, covariate, \dots)
\method{evalCovar}{ppm}(model, covariate, \dots, lambdatype,
          dimyx, eps, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset) 
\method{evalCovar}{slrm}(model, covariate, \dots, lambdatype,
          jitter, jitterfactor,
          modelname, covname, dataname, subset) 
evalCovariate(covariate, locations)
evalInteraction(X,P,E,interaction,correction,splitInf,\dots,
                precomputed,savecomputed)
evalInterEngine(X,P,E,interaction,correction,splitInf,\dots,
                Reach,precomputed,savecomputed)
evalPairPotential(X,P,E,pairpot,potpars,R)
evaluate2Dkernel(kernel, x, y, sigma, varcov, \dots, scalekernel)
ExpSmoothLog(X, \dots, at, weights)
expandDot(f, dotvars)
\method{extractAIC}{slrm}(fit, scale = 0, k = 2, \dots)
extractAtomicQtests(x)
\method{family}{glmmPQL}(object, \dots)
\method{family}{vblogit}(object, \dots)
fasp(fns, which, formulae, dataname, title, rowNames, colNames, checkfv)
f3engine(x, y, z, box, vside, range, nval, correction)
f3Cengine(x, y, z, box, vside, rmax, nrval)
fill.coefs(coefs, required)
findCovariate(covname, scope, scopename=NULL)
findbestlegendpos(\dots)
findcbind(root, depth, maxdepth)
fii(interaction, coefs, Vnames, IsOffset, vnameprefix)
flatfname(x)
\method{flipxy}{msr}(X)
forbid.logi(object)
FormatFaspFormulae(f, argname)
fvexprmap(x)
fvlabels(x, expand=FALSE)
fvlabels(x) <- value
fvlabelmap(x, dot=TRUE)
fvlegend(object, elang)
g3engine(x, y, z, box, rmax, nrval, correction)
g3Cengine(x, y, z, box, rmax, nrval)
getdataname(defaultvalue, \dots, dataname)
getglmdata(object, drop=FALSE)
getglmfit(object)
getglmsubset(object)
getppmdatasubset(object)
getppmOriginalCovariates(object)
getSumFun(abbreviation, classname, ismarked, fatal)
geyercounts(U,X,r,sat,Xcounts,EqualPairs)
geyerdelta2(X,r,sat,\dots,sparseOK, correction)
GLMpredict(fit, data, coefs, changecoef, type)
good.correction.K(X)
%gridadjacencymatrix(dims)
hackglmmPQL(fixed, random, family, data, correlation, weights,
            control, niter, verbose, subset, \dots, reltol)
hasenvelope(X, E)
hasglmfit(object)
hierarchicalordering(i, s)
hiermat(x, h)
ho.engine(model, \dots, nsim, nrmh, start, control, verb)
illegal.iformula(ifmla, itags, dfvarnames)
implemented.for.K(correction, windowtype, explicit)
impliedpresence(tags, formula, df, extranames=character(0))
impliedcoefficients(object, tag, new.coef)
\method{influence}{ppmInfluence}(model, \dots)
instantiate.interact(x, par)
interactionfamilyname(object)
intermaker(f, blank)
is.atomicQtest(x)
\method{is.expandable}{ppm}(x)  
is.interact(x)
\method{is.marked}{mppm}(X, \dots)
\method{is.marked}{msr}(X, \dots)
\method{is.marked}{slrm}(X, \dots)
is.mppm(x)
\method{is.multitype}{mppm}(X, \dots)
\method{is.multitype}{msr}(X, \dots)
\method{is.multitype}{slrm}(X, \dots)
\method{is.poisson}{mppm}(x)
is.scov(x)
k3engine(x, y, z, box, rmax, nrval, correction)
Kborder.engine(X, rmax, nr, correction, weights, ratio)
Knone.engine(X, rmax, nr, weights, ratio)
Krect.engine(X, rmax, nr, correction, weights,
             ratio, fname, use.integers)
Kount(dIJ, bI, b, breaks)
Kwtsum(dIJ, bI, wIJ, b, w, breaks, fatal)
Kpcf.kppm(model, what)
\method{Kmodel}{slrm}(model, \dots)
killinteraction(model)
km.rs.opt(o, cc, d, breaks, KM, RS)
kppmComLik(X, Xname, po, clusters, control, stabilize, weightfun, rmax,
           algorithm, DPP, \dots, pspace)
kppmMinCon(X, Xname, po, clusters, control, stabilize, statistic, statargs,
           algorithm, DPP, \dots, pspace)
kppmPalmLik(X, Xname, po, clusters, control, stabilize, weightfun, rmax, 
           algorithm, DPP, \dots, pspace)
kppmCLadap(X, Xname, po, clusters, control, weightfun, 
           rmax, epsilon, DPP, algorithm, \dots,
           startpar, globStrat)
\method{labels}{ppm}(object, \dots)
\method{leverage}{ppmInfluence}(model, \dots)
localKengine(X, \dots, wantL, lambda, rmax, correction, verbose, rvalue)
localKmultiEngine(X, from, to,
           lambdaFrom, lambdaTo, \dots, rmax, wantL,
           correction, verbose, rvalue, sigma, varcov,
           lambdaX, update, leaveoneout,
           Iexplain, Jexplain, Ikey, Jkey, miss.update, miss.leave)
localpcfengine(X, \dots, delta, rmax, nr, stoyan, lambda, rvalue)
localpcfmatrix(X, i, \dots, lambda, delta, rmax, nr, stoyan)
logi.engine(Q, trend, interaction, \dots, 
            covariates, subsetexpr, clipwin,
            correction, rbord, covfunargs, allcovar, 
            vnamebase, vnameprefix, justQ, savecomputed, precomputed,
            VB)
\method{logLik}{vblogit}(object, \dots)
lookup2DkernelInfo(kernel)
LurkEngine(object, type, cumulative, plot.sd, 
           quadpoints, wts, Z, subQset, 
           covvalues, resvalues, 
           clip, clipwindow, cov.is.im, covrange, 
           typename, covname,
           cl, clenv,
           oldstyle, check,
           verbose, nx, splineargs,
           envelope, nsim, nrank, Xsim,
           internal, checklength)
makefvlabel(op, accent, fname, sub, argname)
mapInterVars(object, subs, mom)
maskLaslett(X, \dots, eps, dimyx, xy, oldX, verbose, plotit)
match.kernel(kernel)
match2DkernelName(kernel)
mctestSigtraceEngine(R, devdata, devsim, \dots,
     interpolate, confint, alpha, exponent, unitname)
meanlistfv(z, \dots)
model.se.image(fit, W, \dots, what)
modelFrameGam(formula, \dots)
mpl.engine(Q, trend, interaction, \dots,
         covariates, subsetexpr, clipwin, covfunargs, correction,
	 rbord, use.gam, gcontrol,
         GLM, GLMfamily, GLMcontrol, famille,
         forcefit, nd, eps, allcovar, callstring,
         precomputed, savecomputed, preponly,
         rename.intercept, justQ, weightfactor)
mpl.get.covariates(covariates, locations, type, covfunargs, need.deriv)
mpl.prepare(Q, X, P, trend, interaction, covariates, 
            want.trend, want.inter, correction, rbord, Pname,
            callstring, \dots,
            subsetexpr,
            covfunargs, allcovar, precomputed, savecomputed,
            vnamebase, vnameprefix, warn.illegal, warn.unidentifiable,
            weightfactor, skip.border, clip.interaction, splitInf)
mpl.usable(x)
\method{names}{fv}(x) <- value
newformula(old, change, eold, enew, expandpoly, dotvars)
newstyle.coeff.handling(object)
nncleanEngine(kthNND, k, d, \dots, tol, maxit,
              plothist, lineargs, verbose, Xname)
nndcumfun(X, \dots, r)
no.trend.ppm(x)
objsurfEngine(objfun, optpar, objargs,
              \dots, dotargs, objname,
              ngrid, xlim, ylim, ratio, verbose)
optimConverged(x)
optimStatus(x, call)
optimNsteps(x)
outdated.interact(object)
oversize.quad(Q, \dots, nU, nX, p)
PairPotentialType(pairpot)
\method{pairs}{listof}(\dots, plot=TRUE)
\method{pairs}{solist}(\dots, plot=TRUE)
\method{plot}{localpcfmatrix}(x, \dots)
PoisSaddle(beta, fi)
PoisSaddleArea(beta, fi)
PoisSaddleGeyer(beta, fi)
PoisSaddlePairwise(beta, fi)
polyLaslett(X, \dots, oldX, verbose, plotit)
PPMmodelmatrix(object, data, \dots, subset, Q, keepNA, irregular,
               splitInf)
printStatus(x, errors.only)
printStatusList(stats)
partialModelMatrix(X,D,model,callstring,\dots)
pcf3engine(x, y, z, box, rmax, nrval, correction, delta)
pcfmulti.inhom(X, I, J, lambdaI = NULL, lambdaJ = NULL, \dots,
               r = NULL, breaks = NULL,
               kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
               correction = c("translate", "Ripley"),
               sigma = NULL, varcov = NULL,
               Iname = "points satisfying condition I",
               Jname = "points satisfying condition J")
\method{pcfmodel}{slrm}(model, \dots)
ploterodewin(W1, W2, col.edge, col.inside, do.plot, \dots)
ploterodeimage(W, Z, \dots, Wcol, rangeZ, colsZ, do.plot)
\method{plot}{addvar}(x, \dots, do.points=FALSE)
\method{plot}{bw.frac}(x, \dots)
\method{plot}{bw.optim}(x, \dots, showopt, optargs)
\method{plot}{lurk}(x, \dots, shade)
\method{plot}{minconfit}(x, \dots)
\method{plot}{parres}(x, \dots)
\method{plot}{plotpairsim}(x, \dots)
\method{plot}{qqppm}(x, \dots, limits=TRUE,
           monochrome=spatstat.options('monochrome'),
           limcol=if(monochrome) "black" else "red")
\method{plot}{spatialcdf}(x, \dots, xlab, ylab)
\method{ppm}{default}(Q, trend, interaction,
       \dots, covariates, data, covfunargs, subset, clipwin,
       correction, rbord, use.gam, method, forcefit, emend, project,
       prior.mean, prior.var,
       nd, eps, gcontrol, nsim, nrmh, start, control,
       verb, callstring)
ppmCovariates(model)
ppmDerivatives(fit, what, Dcovfun, loc, covfunargs)
ppmInfluenceEngine(fit, what, \dots, iScore, iHessian, iArgs,
              drop, method, fine, precomputed, sparseOK,
              fitname, multitypeOK, entrywise, matrix.action,
              dimyx, eps,
              geomsmooth)
\method{predict}{vblogit}(object, newdata, type, se.fit, dispersion,
                          terms, na.action, \dots)
\method{predict}{profilepl}(object, \dots)
prefixfv(x, tagprefix, descprefix, lablprefix, whichtags)
\method{print}{addvar}(x, \dots)
\method{print}{bt.frame}(x, \dots)
\method{print}{bw.frac}(x, \dots)
\method{print}{bw.optim}(x, \dots)
\method{print}{densityfun}(x, \dots)
\method{print}{diagppm}(x, \dots)
\method{print}{detpointprocfamily}(x, \dots)
\method{print}{detpointprocfamilyfun}(x, \dots)
\method{print}{envelope}(x, \dots)
\method{print}{fasp}(x, \dots)
\method{print}{fv}(x, \dots, tight)
\method{print}{fvfun}(x, \dots)
\method{print}{hasenvelope}(x, \dots)
\method{print}{hierarchicalordering}(x, \dots)
\method{print}{influence.ppm}(x, \dots)
\method{print}{interact}(x, \dots, family, brief, banner)       
\method{print}{intermaker}(x, \dots)
\method{print}{isf}(x, \dots)
\method{print}{laslett}(x, \dots)
\method{print}{leverage.ppm}(x, \dots)
\method{print}{localpcfmatrix}(x, \dots)
\method{print}{lurk}(x, \dots)
\method{print}{minconfit}(x, \dots)
\method{print}{mppm}(x, \dots)
\method{print}{msr}(x, \dots)
\method{print}{parres}(x, \dots)
\method{print}{plotppm}(x, \dots)
\method{print}{plotpairsim}(x, \dots)
\method{print}{profilepl}(x, \dots)
\method{print}{quadrattest}(x, \dots)
\method{print}{qqppm}(x, \dots)
\method{print}{rat}(x, \dots)
\method{print}{rppm}(x, \dots)
\method{print}{Smoothfun}(x, \dots)
\method{print}{summary.mppm}(x, \dots, brief)
\method{print}{summary.slrm}(x, \dots)
\method{print}{summary.ssf}(x, \dots)
\method{print}{vblogit}(x, \dots)
quad.mppm(x)
quadBlockSizes(nX, nD, p, nMAX, announce)
ratfv(df, numer, denom, \dots, ratio)
\method{reach}{slrm}(x, \dots)
rebadge.as.crossfun(x, main, sub, i, j)
rebadge.as.dotfun(x, main, sub, i)
rebadge.fv(x, new.ylab, new.fname, tags, new.desc, new.labl, new.yexp,
           new.dotnames, new.preferred, new.formula, new.tags)
rebadgeLabels(x, new.fname)
reconcile.fv(\dots)
reduceformula(fmla, deletevar, verbose)
RelevantDeviation(x, alternative, clamp, scaling)
rename.fv(x, fname, ylab, yexp)
\method{rescale}{msr}(X, s, unitname)
resolveEinfo(x, what, fallback, warn, atomic)
resolve.foxall.window(X, Y, W, warn.trim)
resolve.lambda(X, lambda, \dots, sigma, varcov, leaveoneout, update)
resolve.lambda.cross(X, I, J, lambdaI, lambdaJ, \dots, lambdaX,
                     sigma, varcov, leaveoneout, update, lambdaIJ,
                     Iexplain, Jexplain, miss.update, miss.leave,
                     caller)
roseContinuous(ang, rad, unit, \dots,
               start, clockwise, main, labels, at, do.plot)
\method{rotate}{msr}(X, angle, \dots, centre)
quadrat.testEngine(X, nx, ny, alternative, method, conditional, CR, 
     \dots, nsim, Xcount, xbreaks, ybreaks, tess, fit, df.est, Xname, fitname)
reincarnate.interact(object)
resid4plot(RES, plot.neg, plot.smooth,
           spacing, outer, srange, monochrome, main,
           xlab, ylab, rlab, col.neg, col.smooth, \dots)
resid1plot(RES, opt, plot.neg, plot.smooth,
              srange, monochrome, main,
              add, show.all, do.plot, col.neg, col.smooth, \dots)
resid1panel(observedX, observedV,
            theoreticalX, theoreticalV, theoreticalSD,
            xlab,ylab, \dots, do.plot)
resolve.2D.kernel(\dots,
            sigma, varcov, x, mindist, adjust, bwfun, allow.zero)
rhohatEngine(model, covariate, reference, volume, \dots,
               subset, weights, method, horvitz, smoother,
               resolution, evalCovarArgs,
               n, bw, adjust, from, to, 
               bwref, covname, covunits, confidence,
               breaks,
               modelcall, callstring)
rhohatCalc(ZX, Zvalues, lambda, denom, \dots,
           weights, lambdaX,
           method, horvitz, smoother,
           n, bw, adjust, from, to, 
           bwref, covname, confidence, breaks, positiveCI, markovCI,
           covunits, modelcall, callstring, savestuff)
rmax.Rigid(X, g)
rmax.rule(fun, W, lambda)
rocData(covariate, nullmodel, \dots, high)
rocModel(lambda, nullmodel, \dots, high)
safeFiniteValue(x, default)
safePositiveValue(x, default)
\method{scalardilate}{msr}(X, f, \dots)
scanmeasure(X, \dots)
\method{scanmeasure}{ppp}(X, r, \dots, method) 
\method{scanmeasure}{im}(X, r, \dots)
scanPoisLRTS(nZ, nG, muZ, muG, alternative)
scanBinomLRTS(nZ, nG, muZ, muG, alternative)
second.moment.calc(x, sigma, edge, what, \dots,
                   varcov, expand, obswin, npts, debug)
second.moment.engine(x, sigma, edge, what, \dots,
      kernel, scalekernel, obswin, varcov, npts, debug, fastgauss)
sewpcf(d, w, denargs, lambda2area, divisor)
sewsmod(d, ff, wt, Ef, rvals, method="smrep", \dots, nwtsteps=500)
\method{shift}{influence.ppm}(X, \dots)
\method{shift}{leverage.ppm}(X, \dots)
\method{shift}{msr}(X, \dots)
\method{shift}{quadrattest}(X, \dots)
signalStatus(x, errors.only)
\method{simulate}{profilepl}(object, \dots)
simulrecipe(type, expr, envir, csr, pois, constraints)
slr.prepare(CallInfo, envir, data, dataAtPoints, splitby, clip)
slrAssemblePixelData(Y, Yname, W, covimages, dataAtPoints, pixelarea)
slrmInfluence(model, what, \dots)
\method{Smooth}{solist}(X, \dots)
smoothcrossEngine(Xdata, Xquery, values, sigma, \dots,
                    weights, varcov,
                    kernel, scalekernel, sorted, cutoff)
smoothpointsEngine(x, values, sigma, \dots,
                   kernel, scalekernel,
                   weights, varcov, leaveoneout, sorted, cutoff, debug)
spatstatDPPModelInfo(model)
spatialCDFframe(model, covariate, \dots, jitter)
spatialCDFtest(model, covariate, test, \dots,
         dimyx, eps, interpolate, jitter,
         nsim, verbose, modelname, covname, dataname)
spatialCDFtestCalc(fra, test, \dots, details)
sphere.volume(range, nval = 10)
splitHybridInteraction(coeffs, inte)
sp.foundclass(cname, inlist, formalname, argsgiven)             
sp.foundclasses(cnames, inlist, formalname, argsgiven)
strausscounts(U,X,r,EqualPairs)
suffloc(object)
suffstat.generic(model, X, callstring)
suffstat.poisson(model, X, callstring)
\method{summary}{envelope}(object,\dots)
\method{summary}{mppm}(object, \dots, brief=FALSE)
\method{summary}{msr}(object, \dots)
\method{summary}{profilepl}(object, \dots)
\method{summary}{vblogit}(object, \dots)
tweak.coefs(model, new.coef)
tweak.fv.entry(x, current.tag, new.labl, new.desc, new.tag)
tweak.ratfv.entry(x, \dots)
twostage.test(X, \dots, exponent, nsim, nsimsub,
              alternative, reuse, leaveout, interpolate,
              savefuns, savepatterns, verbose, badXfatal, testblurb)
twostage.envelope(X, \dots, nsim, nsimsub, nrank,
                  alternative, reuse, leaveout, interpolate,
                  savefuns, savepatterns, verbose, badXfatal, testlabel)
\method{unitname}{msr}(x) 
\method{unitname}{msr}(x) <- value 
\method{update}{ippm}(object, \dots, envir) 
\method{update}{msr}(object, \dots) 
validate2Dkernel(kernel, fatal)
validate.angles(angles, unit, guess)
validate.weights(x, recip, how, allowzero, allowinf)
vanilla.fv(x)
varcountEngine(g, B, lambdaB, f, R, what)
%vblogit(y, X, offset, eps, m0, S0, S0i, xi0, verb, maxiter, \dots)
%vblogit.fmla(formula, offset, data, subset, weights, verbose, epsilon, \dots)
versionstring.interact(object)
versionstring.ppm(object)
weightedclosepairs(X, r, correction, what)
windows.mppm(x)
X2testEngine(OBS, EXP, \dots, method, CR, df, nsim,
     conditional, alternative, testname, dataname)



}
\details{
  These internal \pkg{spatstat.core} functions should not be called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat.core} to the next.
}
\value{
The return values of these functions are not documented,
and may change without warning.
}
\keyword{internal}

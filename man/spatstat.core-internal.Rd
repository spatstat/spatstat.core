\name{spatstat.core-internal} 
\title{Internal spatstat.core functions}
\alias{[.localpcfmatrix}
\alias{[.rat}
\alias{adjust.ratfv}
\alias{ang2rad}
\alias{as.data.frame.bw.optim}
\alias{as.data.frame.fv}
\alias{assemble.plot.objects}       
\alias{bandwidth.is.infinite}
\alias{BartCalc}
\alias{bermantestCalc}
\alias{bermantestEngine}
\alias{bind.ratfv}
\alias{blankcoefnames}
\alias{bt.frame}
\alias{bw.optim}
\alias{calc.DR}
\alias{calc.NNIR}
\alias{calc.SAVE}
\alias{calc.SIR}
\alias{calc.TSE}
\alias{censtimeCDFest}
\alias{check.testfun}
\alias{circticks}
\alias{clarkevansCalc}
\alias{compatible.rat}
\alias{compileCDF}
\alias{conform.ratfv}
\alias{CressieReadStatistic}
\alias{CressieReadSymbol}
\alias{CressieReadName}
\alias{cutoff2Dkernel}
\alias{CVforPCF}
\alias{Deviation}
\alias{densitycrossEngine}
\alias{densitypointsEngine}
\alias{digestCovariates}
\alias{digital.volume}
\alias{dim.fasp}               
\alias{dimnames.fasp}               
\alias{dimnames<-.fasp}
\alias{distributecbind}
\alias{envelopeEngine}
\alias{envelopeProgressData}
\alias{envelopeTest}
\alias{envelope.hasenvelope}
\alias{envelope.matrix}
\alias{evaluateCovariate}
\alias{evaluateCovariateAtPoints}
\alias{evaluateCovariateAtPixels}
\alias{evaluate2Dkernel}
\alias{ExpSmoothLog}
\alias{exactppm}
\alias{extractAtomicQtests}
\alias{fasp}
\alias{f3engine}
\alias{f3Cengine}
\alias{findbestlegendpos}
\alias{findcbind}
\alias{flatfname}
\alias{FormatFaspFormulae}
\alias{fvexprmap}
\alias{fvlabels}
\alias{fvlabels<-}
\alias{fvlabelmap}
\alias{fvlegend}
\alias{g3engine}
\alias{g3Cengine}
\alias{getSumFun}
\alias{good.correction.K}
\alias{hasenvelope}
\alias{implemented.for.K}
\alias{is.atomicQtest}
\alias{is.poisson.exactppm}
\alias{is.scov}
\alias{is.stationary.exactppm}
\alias{k3engine}
\alias{Kborder.engine}
\alias{Knone.engine}
\alias{Krect.engine}
\alias{Kount}
\alias{Kwtsum}               
\alias{km.rs.opt}
\alias{localKengine}
\alias{localKmultiEngine}
\alias{localpcfengine}
\alias{localpcfmatrix}
\alias{lookup2DkernelInfo}
\alias{makefvlabel}
\alias{maskLaslett}
\alias{match2DkernelName}
\alias{match.kernel}
\alias{mctestSigtraceEngine}
\alias{meanlistfv}
\alias{names<-.fv}
\alias{nncleanEngine}
\alias{pairs.listof}
\alias{pairs.solist}
\alias{pcf3engine}
\alias{pcfmulti.inhom}
\alias{plot.bw.frac}
\alias{plot.bw.optim}
\alias{plot.localpcfmatrix}
\alias{plot.plotpairsim}
\alias{plot.spatialcdf}
\alias{polyLaslett}
\alias{predict.exactppm}
\alias{prefixfv}
\alias{print.bw.frac}
\alias{print.bw.optim}
\alias{print.densityfun}
\alias{print.envelope}
\alias{print.exactppm}
\alias{print.fasp}       
\alias{print.fv}       
\alias{print.fvfun}       
\alias{print.hasenvelope}       
\alias{print.laslett}
\alias{print.localpcfmatrix}
\alias{print.plotpairsim}
\alias{print.quadrattest}
\alias{print.rat}
\alias{print.Smoothfun}       
\alias{print.summary.ssf}
\alias{quadrat.testEngine}
\alias{ratfv}
\alias{rebadge.as.crossfun}
\alias{rebadge.as.dotfun}
\alias{rebadge.fv}
\alias{rebadgeLabels}
\alias{reconcile.fv}
\alias{RelevantDeviation}
\alias{rename.fv}
\alias{resolve.2D.kernel}
\alias{resolve.foxall.window}
\alias{resolve.lambda}
\alias{resolve.reciplambda}
\alias{resolve.lambda.cross}
\alias{resolveEinfo}
\alias{rhohatEngine}
\alias{rhohatCalc}
\alias{rmax.Rigid}
\alias{rmax.rule}
\alias{rocData}
\alias{rocModel}
\alias{roseContinuous}
\alias{scanmeasure}
\alias{scanmeasure.ppp}
\alias{scanmeasure.im}
\alias{scanBinomLRTS}
\alias{scanPoisLRTS}
\alias{second.moment.calc}
\alias{second.moment.engine}
\alias{sewpcf}
\alias{sewsmod}
\alias{shift.quadrattest}              
\alias{Smooth.solist}
\alias{smoothcrossEngine}              
\alias{smoothpointsEngine}
\alias{spatialCDFframe}
\alias{spatialCDFtest}
\alias{spatialCDFtestCalc}
\alias{spatialCovariateEvidence}
\alias{spatialCovariateEvidence.exactppm}
\alias{spatialCovariateEvidence.ppm}
\alias{spatialCovariateEvidence.ppp}
\alias{spatialCovariateEvidence.slrm}
\alias{sphere.volume}
\alias{summary.envelope}
\alias{tweak.fv.entry}
\alias{tweak.ratfv.entry}
\alias{twostage.test}
\alias{twostage.envelope}
\alias{validate2Dkernel}
\alias{validate.angles}        
\alias{validate.weights}        
\alias{vanilla.fv}
\alias{weightedclosepairs}
\alias{X2testEngine}
\alias{accumulateStatus}
\alias{active.interactions}
\alias{adaptcoef}
\alias{affine.msr}
\alias{areadelta2}
\alias{as.isf}
\alias{as.ppm.rppm}
\alias{augment.msr}
\alias{bigvaluerule}
\alias{cannot.update}
\alias{check.separable}
\alias{coef.summary.kppm}
\alias{coef.summary.ppm}
\alias{coef.summary.slrm}
\alias{coef.vblogit} 
\alias{condSimCox}
\alias{damaged.ppm}
\alias{data.mppm}
\alias{deltasuffstat}
\alias{dfbetas.ppmInfluence}
\alias{diagnose.ppm.engine}
\alias{dim.msr}
\alias{dimnames.msr}
\alias{doMultiStraussHard}
\alias{dppDpcf}
\alias{dppmFixAlgorithm}
\alias{dppmFixIntensity}
\alias{DPPSaddle}
\alias{DPPSaddlePairwise}
\alias{equalpairs}          
\alias{evalInteraction}
\alias{evalInterEngine}
\alias{evalPairPotential}
\alias{evalPairwiseTerm}
\alias{expandDot}
\alias{extractAIC.slrm}
\alias{fakefii}
\alias{family.hackglmmPQL} 
\alias{family.vblogit} 
\alias{fill.coefs}
\alias{findCovariate}
\alias{fii}
\alias{flipxy.msr}
\alias{forbid.logi}
\alias{formula.hackglmmPQL} 
\alias{getdataname}
\alias{getglmdata}
\alias{getglmfit}
\alias{getglmsubset}
\alias{getppmdatasubset}
\alias{getppmOriginalCovariates}
\alias{geyercounts}
\alias{geyerdelta2}
\alias{GLMpredict}
\alias{hackglmmPQL}
\alias{hasglmfit}
\alias{hierarchicalordering}
\alias{hiermat}
\alias{ho.engine}
\alias{illegal.iformula}
\alias{impliedpresence}
\alias{impliedcoefficients}
\alias{influence.ppmInfluence}
\alias{instantiate.interact}
\alias{interactionfamilyname}
\alias{intermaker}
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
\alias{Kpcf.kppm}
\alias{Kmodel.slrm}
\alias{killinteraction}
\alias{kppmComLik}
\alias{kppmMinCon}
\alias{kppmPalmLik}
\alias{kppmCLadap}
\alias{labels.ppm}
\alias{leverage.ppmInfluence}
\alias{logi.engine}
\alias{logLik.vblogit}
\alias{LurkEngine}
\alias{mapInterVars}
\alias{Mayer}
\alias{model.se.image}
\alias{modelFrameGam}
\alias{mpl.engine}
\alias{mpl.get.covariates}
\alias{mpl.prepare}
\alias{mpl.usable}
\alias{newformula}
\alias{newstyle.coeff.handling}
\alias{nndcumfun}
\alias{no.trend.ppm}
\alias{objsurfEngine}
\alias{optimConverged}
\alias{optimStatus}
\alias{optimNsteps}
\alias{outdated.interact}
\alias{oversize.quad}    
\alias{PairPotentialType}
\alias{partialModelMatrix}
\alias{pcfmodel.slrm}
\alias{ploterodewin}
\alias{ploterodeimage}
\alias{plot.addvar}
\alias{plot.lurk}
\alias{plot.minconfit}
\alias{plot.parres}
\alias{plot.qqppm}
\alias{PoisSaddle}
\alias{PoisSaddleArea}
\alias{PoisSaddleGeyer}
\alias{PoisSaddlePairwise}
\alias{PPMmodelmatrix}
\alias{ppm.default}
\alias{ppmCovariates}
\alias{ppmDerivatives}
\alias{ppmInfluenceEngine}
\alias{predict.profilepl}
\alias{predict.vblogit}
\alias{printStatus}
\alias{printStatusList}
\alias{print.addvar}    
\alias{print.bt.frame}
\alias{print.diagppm}
\alias{print.detpointprocfamily}
\alias{print.detpointprocfamilyfun}
\alias{print.hierarchicalordering}
\alias{print.influence.ppm}       
\alias{print.interact}       
\alias{print.intermaker}       
\alias{print.isf}
\alias{print.leverage.ppm}
\alias{print.lurk}
\alias{print.minconfit}
\alias{print.mppm}
\alias{print.msr}
\alias{print.parres}
\alias{print.plotppm}
\alias{print.profilepl}
\alias{print.qqppm}
\alias{print.rppm}
\alias{print.summary.mppm}
\alias{print.summary.slrm}
\alias{print.vblogit}
\alias{quad.mppm}
\alias{quadBlockSizes}
\alias{reach.slrm}
\alias{reduceformula}
\alias{reincarnate.interact}
\alias{rescale.msr}
\alias{resid4plot}
\alias{resid1plot}
\alias{resid1panel}
\alias{rotate.msr}
\alias{SaddleApprox}
\alias{safeFiniteValue}
\alias{safePositiveValue}
\alias{scalardilate.msr}
\alias{shift.influence.ppm}              
\alias{shift.leverage.ppm}              
\alias{shift.msr}              
\alias{signalStatus}
\alias{simulate.profilepl}
\alias{simulrecipe}              
\alias{slr.prepare}
\alias{slrAssemblePixelData}
\alias{slrmInfluence}
\alias{spatstatDPPModelInfo}
\alias{splitHybridInteraction}
\alias{sp.foundclass}
\alias{sp.foundclasses}
\alias{strausscounts}
\alias{stripGLMM}
\alias{suffloc}
\alias{suffstat.generic}
\alias{suffstat.poisson}
\alias{summary.mppm}
\alias{summary.msr}
\alias{summary.profilepl}
\alias{summary.vblogit}
\alias{tweak.coefs}
\alias{unitname.msr}
\alias{unitname<-.msr}
\alias{update.ippm}
\alias{update.msr}
\alias{varcountEngine}
\alias{versionstring.interact}
\alias{versionstring.ppm}
\alias{windows.mppm}




%%%%%%%
\description{
  Internal spatstat.core functions.
}
\usage{

\method{[}{localpcfmatrix}(x, i, \dots)
\method{[}{rat}(x, \dots)
adjust.ratfv(f, columns, numfactor, denfactor)
ang2rad(ang, unit, start, clockwise)
\method{as.data.frame}{bw.optim}(x, \dots)
\method{as.data.frame}{fv}(x, \dots)
assemble.plot.objects(xlim, ylim, \dots, lines, polygon)
bandwidth.is.infinite(sigma)
BartCalc(fY, fK)
bermantestCalc(fram, which, alternative, \dots)
bermantestEngine(model, covariate, which, alternative, \dots,
                 modelname, covname, dataname)
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
censtimeCDFest(o, cc, d, breaks, \dots,
     KM, RS, HAN, RAW, han.denom, tt, pmax, fname, fexpr)
check.testfun(f, f1, X)
circticks(R, at, unit, start, clockwise, labels)
clarkevansCalc(X, correction, clipregion, working)
\method{compatible}{rat}(A, B, \dots)
compileCDF(D, B, r, \dots, han.denom, check)
conform.ratfv(x)
CressieReadStatistic(OBS,EXP,lambda,normalise,named)
CressieReadSymbol(lambda)
CressieReadName(lambda)
cutoff2Dkernel(kernel, sigma, varcov, \dots, scalekernel, cutoff, fatal)
CVforPCF(bw, stuff)
Deviation(x, ref, leaveout, n, xi)
densitycrossEngine(Xdata, Xquery, sigma, \dots,
                    kernel, scalekernel,
                    weights, edge, varcov,
                    diggle, sorted, cutoff)
densitypointsEngine(x, sigma, \dots,
                    kernel, scalekernel,
                    weights, edge, varcov,
                    leaveoneout, diggle, sorted, spill, cutoff,
                    debug)
digestCovariates(\dots, W)
digital.volume(range, nval, vside)
\method{dim}{fasp}(x)
\method{dimnames}{fasp}(x)
\method{dimnames}{fasp}(x) <- value
distributecbind(x)
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
evaluateCovariate(covariate, locations, \dots)
evaluateCovariateAtPixels(covariate, locations, \dots, types, eps, dimyx)
evaluateCovariateAtPoints(covariate, locations, \dots, allow.column)
evaluate2Dkernel(kernel, x, y, sigma, varcov, \dots, scalekernel)
ExpSmoothLog(X, \dots, at, weights)
exactppm(X, baseline, \dots, subset, eps, dimyx)
extractAtomicQtests(x)
fasp(fns, which, formulae, dataname, title, rowNames, colNames, checkfv)
f3engine(x, y, z, box, vside, range, nval, correction)
f3Cengine(x, y, z, box, vside, rmax, nrval)
findbestlegendpos(\dots)
findcbind(root, depth, maxdepth)
flatfname(x)
FormatFaspFormulae(f, argname)
fvexprmap(x)
fvlabels(x, expand=FALSE)
fvlabels(x) <- value
fvlabelmap(x, dot=TRUE)
fvlegend(object, elang)
g3engine(x, y, z, box, rmax, nrval, correction)
g3Cengine(x, y, z, box, rmax, nrval)
getSumFun(abbreviation, classname, ismarked, fatal)
good.correction.K(X)
hasenvelope(X, E)
implemented.for.K(correction, windowtype, explicit)
is.atomicQtest(x)
\method{is.poisson}{exactppm}(x)
is.scov(x)
\method{is.stationary}{exactppm}(x)
k3engine(x, y, z, box, rmax, nrval, correction)
Kborder.engine(X, rmax, nr, correction, weights, ratio)
Knone.engine(X, rmax, nr, weights, ratio)
Krect.engine(X, rmax, nr, correction, weights,
             ratio, fname, use.integers)
Kount(dIJ, bI, b, breaks)
Kwtsum(dIJ, bI, wIJ, b, w, breaks, fatal)
km.rs.opt(o, cc, d, breaks, KM, RS)
localKengine(X, \dots, wantL, lambda, rmax, correction, verbose, rvalue)
localKmultiEngine(X, from, to,
           lambdaFrom, lambdaTo, \dots, rmax, wantL,
           correction, verbose, rvalue, sigma, varcov,
           lambdaX, update, leaveoneout,
           Iexplain, Jexplain, Ikey, Jkey)
localpcfengine(X, \dots, delta, rmax, nr, stoyan, lambda, rvalue)
localpcfmatrix(X, i, \dots, lambda, delta, rmax, nr, stoyan)
lookup2DkernelInfo(kernel)
makefvlabel(op, accent, fname, sub, argname)
maskLaslett(X, \dots, eps, dimyx, xy, oldX, verbose, plotit)
match.kernel(kernel)
match2DkernelName(kernel)
mctestSigtraceEngine(R, devdata, devsim, \dots,
     interpolate, confint, alpha, exponent, unitname)
meanlistfv(z, \dots)
\method{names}{fv}(x) <- value
nncleanEngine(kthNND, k, d, \dots, tol, maxit,
              plothist, lineargs, verbose, Xname)
\method{pairs}{listof}(\dots, plot=TRUE)
\method{pairs}{solist}(\dots, plot=TRUE)
pcf3engine(x, y, z, box, rmax, nrval, correction, delta)
pcfmulti.inhom(X, I, J, lambdaI, lambdaJ, \dots,
               lambdaX, r, breaks, kernel, bw, stoyan, correction,
               sigma, varcov, update, leaveoneout, 
               Iname, Jname)
\method{plot}{bw.frac}(x, \dots)
\method{plot}{bw.optim}(x, \dots, showopt, optargs)
\method{plot}{localpcfmatrix}(x, \dots)
\method{plot}{plotpairsim}(x, \dots)
\method{plot}{spatialcdf}(x, \dots, xlab, ylab)
polyLaslett(X, \dots, oldX, verbose, plotit)
\method{predict}{exactppm}(object, locations, \dots, eps, dimyx)
prefixfv(x, tagprefix, descprefix, lablprefix, whichtags)
\method{print}{bw.frac}(x, \dots)
\method{print}{bw.optim}(x, \dots)
\method{print}{densityfun}(x, \dots)
\method{print}{envelope}(x, \dots)
\method{print}{exactppm}(x, \dots)
\method{print}{fasp}(x, \dots)
\method{print}{fv}(x, \dots, tight)
\method{print}{fvfun}(x, \dots)
\method{print}{hasenvelope}(x, \dots)
\method{print}{laslett}(x, \dots)
\method{print}{localpcfmatrix}(x, \dots)
\method{print}{plotpairsim}(x, \dots)
\method{print}{quadrattest}(x, \dots)
\method{print}{rat}(x, \dots)
\method{print}{Smoothfun}(x, \dots)
\method{print}{summary.ssf}(x, \dots)
quadrat.testEngine(X, nx, ny, alternative, method, conditional, CR, 
     \dots, nsim, Xcount, xbreaks, ybreaks, tess, fit, df.est, Xname, fitname)
ratfv(df, numer, denom, \dots, ratio)
rebadge.as.crossfun(x, main, sub, i, j)
rebadge.as.dotfun(x, main, sub, i)
rebadge.fv(x, new.ylab, new.fname, tags, new.desc, new.labl, new.yexp,
           new.dotnames, new.preferred, new.formula, new.tags)
rebadgeLabels(x, new.fname)
reconcile.fv(\dots)
RelevantDeviation(x, alternative, clamp, scaling)
rename.fv(x, fname, ylab, yexp)
resolve.2D.kernel(\dots,
            sigma, varcov, x, mindist, adjust, bwfun, allow.zero)
resolve.foxall.window(X, Y, W, warn.trim)
resolve.lambda(X, lambda, \dots,
               sigma, varcov, leaveoneout, update, check)
resolve.reciplambda(X, lambda, reciplambda, \dots,
                    sigma, varcov, leaveoneout, update, check)
resolve.lambda.cross(X, I, J, lambdaI, lambdaJ, \dots, lambdaX,
                     sigma, varcov, leaveoneout, update, lambdaIJ,
                     Iexplain, Jexplain)
resolveEinfo(x, what, fallback, warn, atomic)
rhohatEngine(model, covariate, reference, volume, \dots,
               subset, weights, method, horvitz, smoother,
               resolution, spatCovarArgs,
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
roseContinuous(ang, rad, unit, \dots,
               start, clockwise, main, labels, at, do.plot)
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
\method{shift}{quadrattest}(X, \dots)
\method{Smooth}{solist}(X, \dots)
smoothcrossEngine(Xdata, Xquery, values, sigma, \dots,
                    weights, varcov,
                    kernel, scalekernel, sorted, cutoff)
smoothpointsEngine(x, values, sigma, \dots,
                   kernel, scalekernel,
                   weights, varcov, leaveoneout, sorted, cutoff, debug)
spatialCDFframe(model, covariate, \dots, jitter)
spatialCDFtest(model, covariate, test, \dots,
         dimyx, eps, interpolate, jitter,
         nsim, verbose, modelname, covname, dataname)
spatialCDFtestCalc(fra, test, \dots, details)
spatialCovariateEvidence(model, covariate, \dots)
\method{spatialCovariateEvidence}{exactppm}(model, covariate, \dots, lambdatype,
          dimyx, eps, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset, clip.predict) 
\method{spatialCovariateEvidence}{ppm}(model, covariate, \dots, lambdatype,
          dimyx, eps, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset, clip.predict) 
\method{spatialCovariateEvidence}{ppp}(model, covariate, \dots, lambdatype,
          dimyx, eps, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset, clip.predict) 
\method{spatialCovariateEvidence}{slrm}(model, covariate, \dots, lambdatype,
          jitter, jitterfactor,
          modelname, covname, dataname, subset) 
sphere.volume(range, nval = 10)
\method{summary}{envelope}(object,\dots)
tweak.fv.entry(x, current.tag, new.labl, new.desc, new.tag)
tweak.ratfv.entry(x, \dots)
twostage.test(X, \dots, exponent, nsim, nsimsub,
              alternative, reuse, leaveout, interpolate,
              savefuns, savepatterns, verbose, badXfatal, testblurb)
twostage.envelope(X, \dots, nsim, nsimsub, nrank,
                  alternative, reuse, leaveout, interpolate,
                  savefuns, savepatterns, verbose, badXfatal, testlabel)
validate2Dkernel(kernel, fatal)
validate.angles(angles, unit, guess)
validate.weights(x, recip, how, allowzero, allowinf)
vanilla.fv(x)
weightedclosepairs(X, r, correction, what)
X2testEngine(OBS, EXP, \dots, method, CR, df, nsim,
     conditional, alternative, testname, dataname)

accumulateStatus(x, stats)
active.interactions(object)
adaptcoef(new.coef, fitcoef, drop)
\method{affine}{msr}(X, mat, vec, \dots)
areadelta2(X, r, \dots, sparseOK)
as.isf(object)
\method{as.ppm}{rppm}(object)
augment.msr(x, \dots, sigma, recompute)
bigvaluerule(objfun, objargs, startpar, \dots)
cannot.update(\dots)
check.separable(dmat, covname, isconstant, fatal)
\method{coef}{summary.kppm}(object, \dots)
\method{coef}{summary.ppm}(object, \dots)
\method{coef}{summary.slrm}(object, \dots)
\method{coef}{vblogit}(object, \dots) 
condSimCox(object, nsim, \dots, window, n.cond, w.cond,
           giveup, maxchunk, verbose, drop)
damaged.ppm(object)
data.mppm(x)
deltasuffstat(model, \dots,
              restrict, dataonly, sparseOK, quadsub,
              force, warn.forced, verbose, use.special)
\method{dfbetas}{ppmInfluence}(model, \dots)
diagnose.ppm.engine(object, \dots, type, typename, opt,
                         sigma, rbord, compute.sd, compute.cts,
                         envelope, nsim, nrank,
                         rv, oldstyle, splineargs, verbose)
\method{dim}{msr}(x)
\method{dimnames}{msr}(x)
doMultiStraussHard(iradii, hradii, types)
dppDpcf(model, \dots)
dppmFixIntensity(DPP, lambda, po)
dppmFixAlgorithm(algorithm, changealgorithm, clusters, startpar)
DPPSaddle(beta, fi)
DPPSaddlePairwise(beta, fi)
equalpairs(U, X, marked=FALSE)
evalInteraction(X,P,E,interaction,correction,splitInf,\dots,
                precomputed,savecomputed)
evalInterEngine(X,P,E,interaction,correction,splitInf,\dots,
                Reach,precomputed,savecomputed)
evalPairPotential(X,P,E,pairpot,potpars,R)
evalPairwiseTerm(fint, d)
expandDot(f, dotvars)
\method{extractAIC}{slrm}(fit, scale = 0, k = 2, \dots)
fakefii(model)
\method{family}{hackglmmPQL}(object, \dots)
\method{family}{vblogit}(object, \dots)
fill.coefs(coefs, required)
findCovariate(covname, scope, scopename=NULL)
fii(interaction, coefs, Vnames, IsOffset, vnameprefix)
\method{flipxy}{msr}(X)
forbid.logi(object)
\method{formula}{hackglmmPQL}(x, \dots)
getdataname(defaultvalue, \dots, dataname)
getglmdata(object, drop=FALSE)
getglmfit(object)
getglmsubset(object)
getppmdatasubset(object)
getppmOriginalCovariates(object)
geyercounts(U,X,r,sat,Xcounts,EqualPairs)
geyerdelta2(X,r,sat,\dots,sparseOK, correction)
GLMpredict(fit, data, coefs, changecoef, type)
hackglmmPQL(fixed, random, family, data, correlation, weights,
            control, niter, verbose, subset, \dots, reltol)
hasglmfit(object)
hierarchicalordering(i, s)
hiermat(x, h)
ho.engine(model, \dots, nsim, nrmh, start, control, verb)
illegal.iformula(ifmla, itags, dfvarnames)
impliedpresence(tags, formula, df, extranames=character(0))
impliedcoefficients(object, tag, new.coef)
\method{influence}{ppmInfluence}(model, \dots)
instantiate.interact(x, par)
interactionfamilyname(object)
intermaker(f, blank)
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
Kpcf.kppm(model, what)
\method{Kmodel}{slrm}(model, \dots)
killinteraction(model)
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
logi.engine(Q, trend, interaction, \dots, 
            covariates, subsetexpr, clipwin,
            correction, rbord, covfunargs, allcovar, 
            vnamebase, vnameprefix, justQ, savecomputed, precomputed,
            VB)
\method{logLik}{vblogit}(object, \dots)
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
mapInterVars(object, subs, mom)
Mayer(fi, exponent)
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
newformula(old, change, eold, enew, expandpoly, dotvars)
newstyle.coeff.handling(object)
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
partialModelMatrix(X,D,model,callstring,\dots)
\method{pcfmodel}{slrm}(model, \dots)
ploterodewin(W1, W2, col.edge, col.inside, do.plot, \dots)
ploterodeimage(W, Z, \dots, Wcol, rangeZ, colsZ, do.plot)
\method{plot}{addvar}(x, \dots, do.points=FALSE)
\method{plot}{lurk}(x, \dots, shade)
\method{plot}{minconfit}(x, \dots)
\method{plot}{parres}(x, \dots)
\method{plot}{qqppm}(x, \dots, limits=TRUE,
           monochrome=spatstat.options('monochrome'),
           limcol=if(monochrome) "black" else "red")
PoisSaddle(beta, fi)
PoisSaddleArea(beta, fi)
PoisSaddleGeyer(beta, fi)
PoisSaddlePairwise(beta, fi)
PPMmodelmatrix(object, data, \dots, subset, Q, keepNA, irregular,
               splitInf)
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
printStatus(x, errors.only)
printStatusList(stats)
\method{print}{addvar}(x, \dots)
\method{print}{bt.frame}(x, \dots)
\method{print}{diagppm}(x, \dots)
\method{print}{detpointprocfamily}(x, \dots)
\method{print}{detpointprocfamilyfun}(x, \dots)
\method{print}{hierarchicalordering}(x, \dots)
\method{print}{influence.ppm}(x, \dots)
\method{print}{interact}(x, \dots, family, brief, banner)       
\method{print}{intermaker}(x, \dots)
\method{print}{isf}(x, \dots)
\method{print}{leverage.ppm}(x, \dots)
\method{print}{lurk}(x, \dots)
\method{print}{minconfit}(x, \dots)
\method{print}{mppm}(x, \dots)
\method{print}{msr}(x, \dots)
\method{print}{parres}(x, \dots)
\method{print}{plotppm}(x, \dots)
\method{print}{profilepl}(x, \dots)
\method{print}{qqppm}(x, \dots)
\method{print}{rppm}(x, \dots)
\method{print}{summary.mppm}(x, \dots, brief)
\method{print}{summary.slrm}(x, \dots)
\method{print}{vblogit}(x, \dots)
quad.mppm(x)
quadBlockSizes(nX, nD, p, nMAX, announce)
\method{reach}{slrm}(x, \dots)
reduceformula(fmla, deletevar, verbose)
reincarnate.interact(object)
\method{rescale}{msr}(X, s, unitname)
resid4plot(RES, plot.neg, plot.smooth,
           spacing, outer, srange, monochrome, main,
           xlab, ylab, rlab, col.neg, col.smooth, \dots)
resid1plot(RES, opt, plot.neg, plot.smooth,
              srange, monochrome, main,
              add, show.all, do.plot, col.neg, col.smooth, \dots)
resid1panel(observedX, observedV,
            theoreticalX, theoreticalV, theoreticalSD,
            xlab,ylab, \dots, do.plot)
\method{rotate}{msr}(X, angle, \dots, centre)
SaddleApprox(beta, fi, approx)
safeFiniteValue(x, default)
safePositiveValue(x, default)
\method{scalardilate}{msr}(X, f, \dots)
\method{shift}{influence.ppm}(X, \dots)
\method{shift}{leverage.ppm}(X, \dots)
\method{shift}{msr}(X, \dots)
signalStatus(x, errors.only)
\method{simulate}{profilepl}(object, \dots)
simulrecipe(type, expr, envir, csr, pois, constraints)
slr.prepare(CallInfo, envir, data, dataAtPoints, splitby, clip)
slrAssemblePixelData(Y, Yname, W, covimages, dataAtPoints, pixelarea)
slrmInfluence(model, what, \dots)
spatstatDPPModelInfo(model)
splitHybridInteraction(coeffs, inte)
sp.foundclass(cname, inlist, formalname, argsgiven)             
sp.foundclasses(cnames, inlist, formalname, argsgiven)
strausscounts(U,X,r,EqualPairs)
stripGLMM(object)
suffloc(object)
suffstat.generic(model, X, callstring)
suffstat.poisson(model, X, callstring)
\method{summary}{mppm}(object, \dots, brief=FALSE)
\method{summary}{msr}(object, \dots)
\method{summary}{profilepl}(object, \dots)
\method{summary}{vblogit}(object, \dots)
tweak.coefs(model, new.coef)
\method{unitname}{msr}(x) 
\method{unitname}{msr}(x) <- value 
\method{update}{ippm}(object, \dots, envir) 
\method{update}{msr}(object, \dots) 
varcountEngine(g, B, lambdaB, f, R, what)
versionstring.interact(object)
versionstring.ppm(object)
windows.mppm(x)



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

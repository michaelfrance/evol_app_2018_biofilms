*all 3 factors significant ;

proc mixed data=ab_res.biofilm;
	class antibiotic time;
	model log_fraction = antibiotic|time/outp=ab_res.biofilm_r;
	lsmeans antibiotic*time; 
run;

*sub model for kanamycin;

proc mixed data=ab_res.biofilm;
	where antibiotic="kan";
	class antibiotic time;
	model log_fraction = time/outp=ab_res.biofilm_r;
	lsmeans time;
    estimate "0 day versus 15" time -1 1 0 0 0 0 0 0 0 0; 
	estimate "15 day versus 30" time 0 -1 0 0 0 1 0 0 0 ;
	estimate "30 day versus 45" time 0 0 0 0 0 -1 1 0 0 0 ;
	estimate "30 day versus 50" time 0 0 0 0 0 -1 0 1 0 0 ;
	estimate "30 day versus 75" time 0 0 0 0 0 -1 0 0 1 0 ;
	estimate " 75 versus plank 25" time 0 0 0 0 0 0 0 0 -1 1;
run;

*rif';
proc mixed data=ab_res.biofilm;
	where antibiotic="rif";
	class antibiotic time;
	model log_fraction = time/outp=ab_res.biofilm_r;
	lsmeans time;
    estimate "0 day versus 15" time -1 1 0 0 0 0 0 0 0 0; 
	estimate "15 day versus 30" time 0 -1 0 0 1 0 0 0 ;
	estimate "30 day versus 45" time 0 0  0 0 -1 1 0 0 0 ;
	estimate "30 day versus 50" time 0 0  0 0 -1 0 1 0 0 ;
	estimate "30 day versus 75" time 0 0  0 0 -1 0 0 1 0 ;
	estimate " 75 versus plank 25" time 0  0 0 0 0 0 0 -1 1;
run;

 * check the usual assumptions, assumptions met ;
 data ab_res.biofilm_r ; set ab_res.biofilm_r ;
  sqabsresid = sqrt(abs(Resid)) ;
  run ;
 proc plot data = ab_res.biofilm_r ;
  plot Resid*Pred sqabsresid*Pred ;

 proc capability noprint data = ab_res.biofilm_r lineprinter ;
  var Resid ;
  qqplot Resid /normal(mu = est sigma = est symbol='.') square  ;
  run ;

proc ttest data=ab_res.repl;
	class antibiotic;
	var repl_mut_rate;
run;

proc ttest data=ab_res.stat;
	class antibiotic;
	var mut_rate;
run;

proc ttest data=ab_res.fitness;
	class antibiotic;
	var rel_fit;
run;

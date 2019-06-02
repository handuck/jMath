/*
 Copyright (C) 2018 Sang-Gook Han(handuckjs@gmail.com)
 This file is part of jMath

 jMath is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 jMath is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License  
 along with this program. If not, see http://www.gnu.org/licenses/.
*/

(function(jMath){

if ( jMath.stat ) return;

jMath.stat = {

	// binomial distribution

	binopdf: function(x,n,p){
		return jMath.nchoosek(n,x) * Math.pow(p,x) * Math.pow(1-p,n-x);
	},

	binocdf: function(x,n,p){
		var s = 0;
		for ( var i = 0 ; i <= x ; i++ )
		{
			s += jMath.stat.binopdf(i,n,p);
		}
		return s;
	},

	binoinv: function(alpha, n, p)
	{
		var c = 0;
		if ( alpha > 0.5 )
		{
			alpha = 1 - alpha;
			for ( var x = n ; x >= 0 ; x-- )
			{
				c += jMath.stat.binopdf(x,n,p);
				if ( c >= alpha ) return x;
			}
			return n;
		}
		else
		{
			for ( var x = 0 ; x <= n ; x++ )
			{
				c += jMath.stat.binopdf(x,n,p);
				if ( c >= alpha ) return x;
			}		
			return 0;
		}
	},

	binostat : function(n,p) {
		return [ n*p, n*p*(1-p), Math.sqrt(n*p*(1-p)) ];
	},

	// poisson distribution

	poisspdf: function(x,lambda)
	{
		return Math.pow(lambda,x) * Math.exp(-lambda) / jMath.factorial(x);
	},

	poisscdf: function(x,lambda)
	{
		var s = 0;
		for ( var i = 0 ; i <= x ; i++ )
		{
			s += jMath.stat.poisspdf(i,lambda);
		}
		return s;
	},

	poissinv: function(p,lambda)
	{
		var x;
		if ( lambda < 0 || p < 0 || p > 1 ) return Number.NaN;
		else if ( lambda > 0 && p == 1 ) return Number.POSITIVE_INFINITY;
		if ( p > 0 && p < 1 )
		{
			if ( lambda > 10 )
			{
				x = -Math.SQRT2 * jMath.stat.erfcinv(2*p);
				var count = Math.max( Math.ceil( Math.sqrt(lambda) * x +lambda), 0);
				if ( 1-jMath.gammainc( lambda, count+1) < p )
				{
					// upward
					while( 1-jMath.gammainc(lambda, count+2) <  p) count++;
				}
				else
				{
					// downward
					while( 1-jMath.gammainc(lambda, count) >= p ) count--;
				}
				x = count;
			}
			else
			{
				for( x = 0; 1-jMath.gammainc( lambda, x+1 ) < p ; x++);
			}
		}
		return x;
	},
	
	poisstat: function(lambda)
	{
		return [lambda, lambda, Math.sqrt(lambda)];	
	},

	// hypergeometric distribution

	hygepdf: function(x,n,k,m)
	{
		return jMath.nchoosek(k,x) * jMath.nchoosek(n-k,m-x) / jMath.nchoosek(n,m);
	},

	hygecdf: function(x,n,k,m)
	{
		var s = 0;
		for ( var i = 0 ; i <= x ; i++ )
		{
			s += jMath.stat.hygepdf(i,n,k,m);
		}
		return s;
	},

	hygestat: function(n,k,m)
	{
		var mu = k*m/n;
		var sigma = mu * (n-k)*(n-m)/(n*(n-1));
		return [ k*m/n, sigma, Math.sqrt(sigma) ];
	},

	// error function erf(x)
	erf : function(x) {
		var cof = [
				-1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2,
				-9.561514786808631e-3, -9.46595344482036e-4, 3.66839497852761e-4,
				4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6,
				1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8,
				6.529054439e-9, 5.059343495e-9, -9.91364156e-10,
				-2.27365122e-10, 9.6467911e-11, 2.394038e-12,
				-6.886027e-12, 8.94487e-13, 3.13092e-13,
				-1.12708e-13, 3.81e-16, 7.106e-15,
				-1.523e-15, -9.4e-17, 1.21e-16,
				-2.8e-17
			],
			j = cof.length - 1,
			isneg = false,
			d = 0,
			dd = 0,
			t, ty, tmp, res;
		if(x < 0) {
			x = -x;
			isneg = true;
		}
		t = 2 / (2 + x);
		ty = 4 * t - 2;
		for(; j > 0; j--) {
			tmp = d;
			d = ty * d - dd + cof[j];
			dd = tmp;
		}
		res = t * Math.exp(-x*x + 0.5 * (cof[0] + ty * d) - dd);
		return (isneg) ? res - 1 : 1 - res;
	},

	// complmentary error function erfc(x)
	erfc : function(x) {
		return 1 - jMath.stat.erf(x);
	},

	// inverse of the complementary error function
	erfcinv : function(p) {
		var j = 0,
			x, err, t, pp;
		if(p >= 2) return -100;
		if(p <= 0) return 100;
		pp = (p < 1) ? p : 2 - p;
		t = Math.sqrt(-2 * Math.log(pp / 2));
		x = -0.70711 * ((2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t);
		for(; j < 2; j++) {
			err = jMath.stat.erfc(x) - pp;
			x += err / (1.12837916709551257 * Math.exp(-x * x) - x * err);
		}
		return (p < 1) ? x : -x;
	},
			

	// normal distribution
	normpdf: function(x,mu,sigma)
	{
		return Math.exp(-Math.pow((x - mu)/sigma,2)/2) / (sigma * Math.sqrt(2*Math.PI) );
	},

	normcdf: function(x,mu,sigma)
	{
		var t = (x - mu)/(Math.sqrt(2)*sigma);
		return (1 + jMath.stat.erf(t))/2;
	},

	norminv : function(p, mu, sigma) {
		return -1.41421356237309505 * sigma * jMath.stat.erfcinv(2 * p) + mu;
	},

	// exponential distribution
	exppdf: function(x,mu)
	{
		return Math.exp(-x/mu)/mu;
	},

	expcdf: function(x,mu)
	{
		return 1 - Math.exp(-x/mu);	
	},

	expinv: function(p,mu)
	{
		return -Math.log(1-p)*mu;
	},

	// uniform distribution
	unifpdf: function(x,a,b)
	{
		a = a || 0;
		b = b || 1;
		return a <= x && x <= b ? 1/(b-a) : 0;
	},

	unifcdf: function(x,a,b)
	{
		a = a || 0;
		b = b || 1;
		if ( x < a ) return 0;	
		else if ( b <= x ) return 1;
		return (x-a)/(b-a);
	},

	unifinv: function(p,a,b)
	{
		a = a || 0;
		b = b || 1;
		return a + p * (b- a);
	},

	// Student t-distribution
	tpdf: function(x,df)
	{
		var df1_2 = (df+1)/2;
		return jMath.gamma(df1_2) /jMath.gamma(df/2)
				/Math.sqrt(df*Math.PI)
				/Math.pow(1+(x*x)/df, df1_2);
	},

	tcdf: function(x,df)
	{
        var dof2 = df / 2;
        return jMath.betainc((x + Math.sqrt(x * x + df)) / (2 * Math.sqrt(x * x + df)), dof2, dof2);
	},

	tinv : function(p, df) {
		var x = jMath.betainv(2 * Math.min(p, 1 - p), 0.5 * df, 0.5);
		x = Math.sqrt(df* (1 - x) / x);
		return (p > 0.5) ? x : -x;
	},

	chi2pdf: function(x,df)
	{
		return Math.exp((df / 2 - 1) * Math.log(x) - x / 2 - (df / 2) * Math.log(2) - jMath.gammaln(df / 2));
	},

	chi2cdf: function(x,df)
	{
		return jMath.gammainc(x/2,df / 2);
	},

	chi2inv: function(p,df)
	{
		return 2 * jMath.gammaincinv(p, 0.5 * df);
	},

	// df1: level
	// df2: degree of freedom for sum square within
	stdrcdf : function(x,df1,df2)
	{
		// Accuracy can be increased by use of a finer grid.  Increase
		// jmax, kmax and 1/step proportionally.
		var jmax = 15;         // controls maximum number of steps
		var kmax = 15;         // controls maximum number of steps
		var step = 0.45;       // node spacing
		var vmax = 120;        //max d.f. for integration over chi-square

		// Handle illegal or trivial values first.
		var result = Number.NaN;
		var v = df2;
		var q = x;
		var r = df1;
		var xx = 0;
		var c, qw;
		// Compute constants, locate midpoint, adjust steps.
		var g = step / Math.pow(r,0.2);
		if (v > vmax)
		{
			c = Math.log(r * g / Math.sqrt(2*Math.PI));
		}
		else
		{
			h = step / Math.sqrt(v);
			var v2 = v * 0.5;
			c = Math.log(Math.sqrt(2/Math.PI) * r * g * h) - v2 + v2*Math.log(v2) - jMath.gammaln(v2);

			var hj= jMath(-jmax + ':' + jmax).transpose()['*'](h);
			var ehj = hj.exp();
			qw = ehj['*'](q);
			vw = hj['+']( ehj['.^'](2)['-'](1)['*'](-0.5) )['*'](v);
		}
		// Compute integral by summing the integrand over a
		// two-dimensional grid centered approximately near its maximum.
		var gk = jMath(-kmax+':'+kmax)['*'](g)['+']( 0.5 * Math.log(r) );
		var w0 = gk['.^'](2)['*'](-0.5)['+'](c); // c - 0.5 * gk.^2 
		var pz = gk.map( function(v){
		 	return jMath.stat.normcdf(-v,0,1);  
		});
		// For regular cdf, use integrand as in AS 190.
		if (v > vmax)
		{
			// don't integrate over chi-square
			x = gk.map( function(v, r, c){
			 	return jMath.stat.normcdf( q - v, 0, 1) - pz[r][c];  
			})['.^'](r-1);
			return w0.exp()['.*'](x).sum(3);
		}
		else
		{
			// integrate over chi-square
			x = qw.repeatColumn(2*kmax)['-'](gk.repeatRow(2*jmax)).map(function(v,r,c){
				return jMath.stat.normcdf(v,0,1) - pz[0][c];
			})['.^'](r-1);
			xx = w0.repeatRow(2*jmax)['+'](vw.repeatColumn(2*kmax)).exp()['.*'](x).sum(3);
		}
		return xx;
	},

	// df1: level
	// df2: degree of freedom for sum square within
	stdrinv:function(p,df1,df2)
	{
		var jmax = 20;
		var pcut = 0.00001;
		var tiny = 0.000001;
		var upper = p > 0.99;
		var p0 = upper ? 1-p : p

		//  Obtain initial values
		var q2;
		var v = df2;
		var r = df1;
		var q1 = jMath.stat.qtrng0(p, v, r);
   		var p1 = jMath.stat.stdrcdf(q1, df1, df2);
		if ( upper ) p1 = 1 - p1;
		var xx = q1;
		if (Math.abs(p1-p0) >= pcut*p0)
		{
			var p2;
	    	if (p1 > p0)
			{	
				p2 = Math.max(.75*p0, p0-.75*(p1-p0));
			}
			if (p1 < p0)
			{
				p2 = p0 + (p0 - p1) * (1 - p0) / (1 - p1) * 0.75;
			}
			if (upper)
			{
				q2 = jMath.stat.qtrng0(1-p2, v, r);
			}
      		else
			{
		        q2 = jMath.stat.qtrng0(p2, v, r);
			}
		     // Refine approximation
			for ( var j = 1 ; j < jmax ; j++ )
			{
				p2 = jMath.stat.stdrcdf(q2, df1, df2);
				if ( upper ) p2 = 1 - p2;
				var e1 = p1 - p0;
				var e2 = p2 - p0;
				var d = e2 - e1;
				xx = (q1 + q2) / 2;
				if ( Math.abs(d) > tiny*p0)
				{
					xx = (e2 * q1 - e1 * q2) / d;
				}
				if (Math.abs(e1) >= Math.abs(e2))
				{
					q1 = q2;
					p1 = p2;
				}
				if (Math.abs(p1 - p0) < pcut*p0)
				{
				   	break;
				}
				q2 = xx;
			}
		}
		return xx;
	},

	// r: level
	// v: degree of freedom for sum square within
	qtrng0 : function(p,v,r)
	{
		// Algorithm AS 190.2  Appl. Statist. (1983) Vol.32, No.2
		// Calculates an initial quantile p for a studentized range
		// distribution having v degrees of freedom and r samples
		// for probability p, p.gt.0.80 .and. p.lt.0.995.
		var t=jMath.stat.norminv(0.5 + 0.5 * p, 0,1);
		if (v < 120)
		{
			t = t + 0.25 * ( Math.pow(t,3) + t) / v;
		}
		q = 0.8843 - 0.2368 * t;
		if (v < 120)
		{
			q = q - (1.214/v) + (1.208*t/v);
		}
		return t * (q * Math.log(r-1) + 1.4142);
	},


	fpdf: function(x,df1,df2)
	{
   		return  (x >= 0) ?
	         Math.sqrt((Math.pow(df1 * x, df1) * Math.pow(df2, df2)) / (Math.pow(df1 * x + df2, df1 + df2))) / (x * jMath.beta(df1/2, df2/2)) : undefined;
	},

	fcdf: function(x,df1,df2)
	{
        return jMath.betainc((df1 * x) / (df1 * x + df2), df1 / 2, df2 / 2);
	},

	finv: function(p,df1,df2)
	{
		return df2 / (df1 * (1 / jMath.betainv(p, df1 / 2, df2 / 2) - 1));
	},

	ci_h: function(result)
	{
		var z;
		switch( result.type )
		{
			case 0:
				result.pvalue = jMath.stat.normcdf(result.mean,result.H0,result.meanStd);
				if ( result.mean > result.H0 ) result.pvalue = 1 - result.pvalue;
				result.pvalue *= 2;
				z = jMath.stat.norminv( result.alpha/2, 0, 1 );
				result.zscores = [ z, -z ];
				result.ci = [ result.H0 + z * result.meanStd, result.H0 - z * result.meanStd];
				break;
			case 1:
				result.pvalue = 1 - jMath.stat.normcdf(result.mean,result.H0,result.meanStd);
				z = jMath.stat.norminv( 1 - result.alpha, 0, 1 );
				result.zscores = [ Number.NAGATIVE_INFINITY, z ];
				result.ci = [ Number.NAGATIVE_INFINITY, result.H0 + z * result.meanStd];
				break;
			case -1:
				result.pvalue = jMath.stat.normcdf(result.mean,result.H0,result.meanStd);
				z = jMath.stat.norminv( result.alpha, 0, 1 );
				result.zscores = [ z, Number.POSITIVE_INFINITY ];
				result.ci = [  H0 + z * result.meanStd, Number.POSITIVE_INFINITY];
				break;
		}
		return result;
	},

	ci_ht: function(result)
	{
		var z;
		switch( result.type )
		{
			case 0:
				result.pvalue = jMath.stat.tcdf(result.testscore, result.df );
				if ( result.mean > result.H0 ) result.pvalue = 1 - result.pvalue;
				result.pvalue *= 2;
				z = jMath.stat.tinv( result.alpha/2, result.df );
				result.tscores = [ z, -z ];
				result.ci = [ result.H0 + z * result.meanStd, result.H0 - z * result.meanStd];
				break;
			case 1:
				result.pvalue = 1 - jMath.stat.tcdf(result.testscore, result.df);
				z = jMath.stat.tinv( 1-result.alpha, result.df);
				result.tscores = [ Number.NAGATIVE_INFINITY, z ];
				result.ci = [ Number.NAGATIVE_INFINITY, result.H0 + z * result.meanStd];
				break;
			case -1:
				result.pvalue = jMath.stat.tcdf(result.testscore, result.df);
				z = jMath.stat.tinv( result.alpha, result.df);
				result.tscores = [ z, Number.POSITIVE_INFINITY ];
				result.ci = [  result.H0 + z * result.meanStd, Number.POSITIVE_INFINITY];
				break;
		}
		return result;
	},

	// type : 0 two-tail, 1 upper tail, -1 lower tail
	test: function(alpha,type,x,H0,sigma,n,N)
	{
		alpha = alpha || 0.05;
		var testTypes = [ 'lowertail', 'twotail', 'uppertail' ];
		if ( typeof type == 'string' )
		{
			type = testTypes.indexOf(type);
			if ( type == -1 ) type = 0;
			else type--;
		}
		var result = {
			H0: H0,
			alpha: alpha,
			popSigma: sigma,
			sampleSize: n,
			popSize: N,
			type: type,
		}
		var cf = N && N * 0.05 < n ? Math.sqrt( (N-n)/(N-1) ) : 1;
		sigma = cf * sigma / Math.sqrt(n);
		result.mean = x;
		result.meanStd = sigma;
		result.testscore = (result.mean- result.H0)/sigma;
		return jMath.stat.ci_h(result);
	},

	test_t: function(alpha,type,x,H0,sigma,n,N)
	{
		alpha = alpha || 0.05;
		var testTypes = [ 'lowertail', 'twotail', 'uppertail' ];
		if ( typeof type == 'string' )
		{
			type = testTypes.indexOf(type);
			if ( type == -1 ) type = 0;
			else type--;
		}
		var result = {
			H0: H0,
			alpha: alpha,
			sampleSize: n,
			type: type,
			df: n-1,
			mean: x
		};
		var cf = N && n > 0.05 * N ?  Math.sqrt( (N-n)/(N-1) ) : 1;
		result.sigma = sigma;
		sigma = cf * sigma / Math.sqrt(n);
		result.meanStd = sigma;
		result.testscore = (result.mean- result.H0)/sigma;
		return jMath.stat.ci_ht(result);
	},

	test_p: function(alpha,type,p,H0,n,N)
	{
		var result = jMath.stat.test( alpha, type, p, H0, Math.sqrt(H0 * (1-H0)), n, N);	
		result.p = p;
		return result;
	},

	test_v: function(alpha,type,sigma2,H0,n)
	{
		var testTypes = [ 'lowertail', 'twotail', 'uppertail' ];
		var result = {
			H0: H0,
			H1: sigma2,
			alpha: alpha,
			df: n-1,
			chi2: (n-1)*sigma2/H0,
			chi2crit: jMath.stat.chi2inv( 1-alpha, n-1)
		}
		if ( typeof type == 'string' )
		{
			type = testTypes.indexOf(type);
			if ( type == -1 ) type = 0;
			else type--;
		}
		result.type = type;
		var z;
		switch( result.type )
		{
			case 0:
				result.pvalue = jMath.stat.chi2cdf(result.chi2, result.df );
				if ( result.sigma2 > result.H0 ) result.pvalue = 1 - result.pvalue;
				result.pvalue *= 2;
				result.chi2crit = [ jMath.stat.chi2inv( alpha/2, n-1), jMath.stat.chi2inv( 1-alpha/2,n-1)];
				break;
			case 1:
				result.pvalue = 1 - jMath.stat.chi2cdf(result.chi2, result.df);
				result.chi2crit = jMath.stat.chi2inv( 1-alpha, n-1);
				break;
			case -1:
				result.pvalue = jMath.stat.chi2cdf(result.chi2, result.df);
				result.chi2crit = jMath.stat.chi2inv( alpha, n-1);
				break;
		}
		return result;
	},

	compare: function(alpha,type,s1,s2,H0)
	{
		alpha = alpha || 0.05;
		var testTypes = [ 'lowertail', 'twotail', 'uppertail' ];
		if ( typeof type == 'string' )
		{
			type = testTypes.indexOf(type);
			if ( type == -1 ) type = 0;
			else type--;
		}
		var result = {
			H0: type == 0 || H0 == undefined ? 0 : H0,
			alpha: alpha,
			popSigma: [s1.sigma,s2.sigma],
			type: type
		}
		result.mean = s1.mean - s2.mean;
		result.meanStd = Math.sqrt(s1.sigma * s1.sigma/s1.n + s2.sigma*s2.sigma/s2.n);
		result.testscore = (result.mean - result.H0)/result.meanStd;
		result.samples = [s1, s2];
		jMath.stat.ci_h(result);
		return result;
	},

	compare_p: function(alpha,type,s1,s2,H0)
	{
		alpha = alpha || 0.05;
		var testTypes = [ 'lowertail', 'twotail', 'uppertail' ];
		if ( typeof type == 'string' )
		{
			type = testTypes.indexOf(type);
			if ( type == -1 ) type = 0;
			else type--;
		}
		var result = {
			H0: type == 0 || H0 == undefined ? 0 : H0,
			alpha: alpha,
			type: type,
			samples: [ s1, s2 ]
		}
		s1.p = s1.x/s1.n;
		s2.p = s2.x/s2.n;
		var p = (s1.x + s2.x)/(s1.n+s2.n);
		result.mean = s1.p - s2.p;
		result.meanStd = Math.sqrt(p * (1-p) * (1/s1.n + 1/s2.n));
		result.testscore = (result.mean - result.H0)/result.meanStd;
		jMath.stat.ci_h(result);
		var sigma = Math.sqrt( (s1.p * ( 1 - s1.p ))/s1.n + (s2.p * ( 1 - s2.p ))/s2.n);
		result.sigma = sigma;
		var z = jMath.stat.norminv( result.alpha/2, 0, 1 );
		result.ci_diff = [ result.mean + z * sigma, result.mean - z * sigma];
		return result;
	},

	compare_t: function(alpha,type,s1,s2,eqPvar,H0)
	{
		alpha = alpha || 0.05;
		var testTypes = [ 'lowertail', 'twotail', 'uppertail' ];
		if ( typeof type == 'string' )
		{
			type = testTypes.indexOf(type);
			if ( type == -1 ) type = 0;
			else type--;
		}
		var result = {
			H0: type == 0 || H0 == undefined ? 0 : H0,
			alpha: alpha,
			type: type,
			eqPvar: eqPvar,
		}
		s1.df = s1.n - 1;
		s2.df = s2.n - 1;
		result.samples = [ s1, s2 ];
		result.mean = s1.mean - s2.mean;
		var sigma;
		if ( !eqPvar )
		{
			var ms1 = s1.sigma*s1.sigma/s1.n;
			var ms2 = s2.sigma*s2.sigma/s2.n;
			var sigma = ms1 + ms2;
			result.df = Math.floor( (sigma * sigma)/( ms1*ms1/s1.df + ms2*ms2/s2.df ) );
			sigma = Math.sqrt(sigma);
		}
		else
		{
			result.df = s1.df + s2.df;
			var sigma = (s1.df*s1.sigma*s1.sigma + s2.df*s2.sigma*s2.sigma)/result.df;
			sigma = Math.sqrt(sigma*(1/s1.n + 1/s2.n));
		}

		result.meanStd = sigma;
		result.testscore = (result.mean - result.H0)/sigma;

		result.type = 0;
		jMath.stat.ci_ht(result);
		result.ci[0] += result.mean;
		result.ci[1] += result.mean;
		result.ci_diff = result.ci;

		result.type = type;
		return jMath.stat.ci_ht(result);
	},

	compare_v: function(alpha,type,s1,s2)
	{
		var testTypes = [ 'twotail', 'uppertail' ];
		if ( typeof type == 'string' )
		{
			type = testTypes.indexOf(type);
			if ( type == -1 ) type = 0;
		}
		var reverse = false;
		// the variance of s1 must be greater than that of s2
		if ( s1.sigma < s2.sigma )
		{
			var t = s1;
			s1 = s2;
			s2 = t;
			reverse = true;
		}
		s1.df = s1.n - 1;
		s2.df = s2.n - 1;
		var result = {
			alpha: alpha,
			samples: [ s1, s2 ],	
			F: s1.sigma*s1.sigma/s2.sigma/s2.sigma,
			type: type,
			reverse: reverse 
		};
		switch( result.type )
		{
			case 0:
				result.pvalue = 2 *(1-jMath.stat.fcdf(result.F, s1.df, s2.df));
				result.Fcrit = jMath.stat.finv(1-alpha/2, s1.df, s2.df );
				break;
			case 1:
				result.pvalue = 1 - jMath.stat.fcdf(result.F, s1.df, s2.df);
				result.Fcrit = jMath.stat.finv(1-alpha, s1.df, s2.df );
				break;
		}

		return result;		
	},

	ci: function(alpha,mu,sigma,n,N)
	{
		alpha = alpha || 0.05;
		var cf = N && N * 0.05 < n ? Math.sqrt( (N-n)/(N-1) ) : 1;
		sigma = cf * sigma / Math.sqrt(n);
		var z = jMath.stat.norminv(alpha/2, 0, 1);
		return [ mu + z * sigma, mu - z * sigma];
	},

	ci_t: function(alpha,mu,sigma,n,N)
	{
		alpha = alpha || 0.05;
		var cf = N && N * 0.05 < n ? Math.sqrt( (N-n)/(N-1) ) : 1;
		sigma = cf * sigma / Math.sqrt(n);
		var z = jMath.stat.tinv( alpha/2, n-1);
		return [ z * sigma, -z * sigma];
	},

	ci_p: function(alpha,p,n,N)
	{
		alpha = alpha || 0.05;
		var cf = N && N * 0.05 < n ? Math.sqrt( (N-n)/(N-1) ) : 1;
		var z = jMath.stat.norminv(alpha/2, 0, 1);
		sigma = cf * Math.sqrt( p * (1-p) / n );
		return [ p + z * sigma, p - z * sigma];
	},

	ci_var : function(alpha,sigma,n)
	{
		alpha = alpha || 0.05;
		var df = n - 1;
		var v = sigma * sigma;
		return [df*v/jMath.stat.chi2inv(1-alpha/2,df), 
			    df*v/jMath.stat.chi2inv(alpha/2, df) ];
	},

	ci_std : function(alpha,sigma,n)
	{
		var ret = jMath.stat.ci_var(alpha,sigma,n);
		return ret.map( function(v){
			return Math.sqrt(v);
		});
	},

	numSamples: function(alpha,me,sigma,isProp)
	{
		alpha = alpha || 0.05;	
		if ( isProp )
		{
			sigma = Math.sqrt(sigma * ( 1 - sigma ));
		}
		var z = jMath.stat.norminv(alpha/2, 0, 1);
		return Math.ceil(Math.pow(z*sigma/me,2));
	}
};

jMath.fn.ci = function(alpha,sigma,N)
{
	return jMath.stat.ci(alpha,this.mean()[0][0],sigma,this.cols,N);
}

jMath.fn.ci_t = function(alpha,N)
{
	alpha = alpha || 0.05;
	var cf = N && this.cols > 0.05 * N ?  Math.sqrt( (N-this.cols)/(N-1) ) : 1;
	var mean = this.mean()[0][0];
	var std = cf * this.std(false)[0][0] / Math.sqrt(this.cols);
	var t = jMath.stat.tinv( alpha/2, this.cols -1 );
	return [ mean + t*std, mean - t*std];
}

jMath.fn.ci_p = function(alpha,func,N)
{
	var p = this[0].map(func).reduce(function(acc,v){
		return acc+v;
	})/this.cols;		
	return jMath.stat.ci_p(alpha,p,this.cols, N);
}

jMath.fn.ci_var = function(alpha)
{
	alpha = alpha || 0.05;
	var v = this.var(false,3);
	var df = this.cols-1;
	return [df*v/jMath.stat.chi2inv(1-alpha/2,df), 
		    df*v/jMath.stat.chi2inv(alpha/2, df) ];
}

jMath.fn.ci_std = function(alpha)
{
	return this.ci_var(alpha).map(function(v){
		return Math.sqrt(v);	
	});
}

jMath.fn.test_p = function(alpha,type,H0,func,N)
{
	var p = this[0].map(func).reduce(function(acc,v){
		return acc+v;
	})/this.cols;		
	return jMath.stat.ci_p(alpha,type,p,H0, this.cols, N);
}

jMath.fn.test_v = function(alpha,type,H0)
{
	var result = jMath.stat.test_v(alpha,type, this.var()[0][0],H0,this.cols);
	result.mean = this.mean()[0][0];
	return result;
}

jMath.fn.test = function(alpha,type,H0,N)
{
	var x = this.mean(3);
	var sigma = this.std(false,3);
	return jMath.stat.test_t(alpha,type,x,H0,sigma,this.rows*this.cols);
}

// row: level
// col: frequencies
jMath.fn.anova_p = function(alpha)
{
	var sumR = this.sum(1);
	var sumC = this.sum(2);

	// the elements on expect should be 5 or more
	var expect = sumR['./'](sumR.sum(3));

	var expF = sumC['*'](expect);
	var chi2 = this['-'](expF)['.^'](2)['./'](expF);

	var result = {
		alpha: alpha,
		chi2: chi2.sum(3),
		chi2crit: jMath.stat.chi2inv(1-alpha, this.rows-1),
		df: this.rows - 1,
		prop: this['./'](sumC.repeatColumn(this.rows-1)),
		expect: expect
	};
	result.pvalue = 1-jMath.stat.chi2cdf( result.chi2, result.df );
	return result;
}

jMath.fn.fitness = function(alpha,df,expectfn)
{
	var sumR = this.sum(3);	
	var result = {
		alpha : alpha,
		df: df
	};
	var expect;
	if ( Array.isArray(expectfn) )
	{
		result.expect = jMath(expectfn);
	}
	else if ( expectfn instanceof m$.fn.init )
	{
		result.expect = expectfn;
	}
	else 
	{
		expect = [];
		for ( var i = 0 ; i < this.cols; i++ )
		{
			expect[i] = expectfn(i);
		}
		result.expect = jMath(expect);
	}
	expect = result.expect['*'](sumR); 
	result.expectFreq = expect;
	result.removed = [];
	// the elements on expect should be 5 or more
	// otherwise, it will be merged into an adjacent element.

	for( var i = expect.cols - 1 ; i >= 0 ; i-- )
	{
		if ( expect[0][i] < 5 )
		{
			if ( i == 0 )
			{
				expect[0][i+1] += expect[0][i];
				this[0][i+1] += this[0][i];
			}
			else
			{
				expect[0][i-1] += expect[0][i];
				this[0][i-1] += this[0][i];
			}			
			expect[0][i] = 0;
			this[0][i] = 0;
			result.df--;
			result.removed[i] = true;
		}	
		else
		{
			result.removed[i] = false;
		}
	}

	var chi2 = this['-'](expect)['.^'](2)['./'](expect);
	for ( var i = 0 ; i < expect.cols ; i++ )
	{
		if ( expect[0][i] == 0 )
		{
			chi2[0][i] = 0;
		}
	}
//	console.log(this.toString());
//	console.log(expect.toString())
//	console.log(chi2.toString());
	result.chi2 = chi2.sum(3);
	result.chi2crit = jMath.stat.chi2inv(1-alpha, result.df);
	result.chi2norm = chi2['./']( result.chi2 == 0 ? 1 : result.chi2 );
	result.prop = this['./'](sumR);
	result.pvalue = 1-jMath.stat.chi2cdf( result.chi2, result.df );
	return result;		
}		

jMath.fn.fitness_poiss = function(alpha,lambda)
{
	alpha = alpha || 0.05;
	var df = this.cols - 1;
	var sum = 0;
	var total = 0;
	for ( var i = 0 ; i < this.cols ; i++ )
	{
		total += this[0][i];
		sum += i * this[0][i];	
	}
	var slambda = sum/total;
	if ( arguments.length < 2 )
	{
		lambda = slambda;
		df--;
	}
	var len = this.cols;
	var result = this.fitness( alpha, df, function(i){
		if ( i == len - 1 )
		{
			return 1 - jMath.stat.poisscdf( i-1, lambda );
		}
		return jMath.stat.poisspdf(i,lambda);
	});
	result.lambda = lambda;
	result.sample = {
		lambda : slambda
	};
	return result;
}		

jMath.fn.fitness_bino = function(alpha,p)
{
	alpha = alpha || 0.05;
	var n = this.cols-1;
	var df = n;
	var sum = 0;
	var total = 0;
	for ( var i = 0 ; i < n ; i++ )
	{
		total += this[0][i];
		sum += i * this[0][i];	
	}
	var sp = sum/total/n;
	if ( arguments.length < 2 )
	{
		p = sp;
	}
	var result = this.fitness( alpha, df, function(x){
		return jMath.stat.binopdf(x,n,p);
	});
	result.sample = {
		p: sp,
		mean: n * sp
	};
	result.p = p;
	result.mean = n * p;
	return result;
}

jMath.fn.fitness_norm = function(alpha,bins,mu,sigma)
{
	alpha = alpha || 0.05;

	var smu = this.mean()[0][0];
	var ssigma = this.std()[0][0];
	var dfComp = 0;

	if ( mu === undefined )
	{
		mu = smu;	
		dfComp++;
	}
	if ( sigma === undefined )
	{
		sigma = ssigma;
		dfComp++;
	}

	var zscores = this.map( function(v){
		return (v - mu)/sigma;
	});
	if ( bins instanceof m$.fn.init )
	{
		bins = bins.toArray()[0];
	}
	else if ( typeof bins == 'number' )
	{
		var s = -bins/2;
		var interval = 1;
		if ( s < -3 )
		{
			interval = 6/bins;
			s = -3;
		}
		var len = bins;
		bins = [s];
		for ( var i = 1 ; i < len ; i++ )
		{
			bins[i] = bins[i-1] + interval;
		}
	}
	var hist = jMath(zscores.freqdist(bins));

	var interval = bins[1] - bins[0];

	var probs = bins.map( function(v){
		return jMath.stat.normcdf(v+interval,0,1);
	});
	probs[probs.length-1] = 1;
	var df = probs.length - 1 - dfComp;

	for ( var i = bins.length-1 ; i > 0; i-- )
	{
		probs[i] -= probs[i-1];
	}
	var samples = hist.slice(':', 1).transpose();
	var result = samples.fitness( alpha, df, function(i){
		return probs[i];
	});
	result.samples = {
		mu: smu,
		sigma: ssigma
	};
	result.interval = interval;
	result.histogram = hist;
	return result;
}

jMath.fn.test_indep = function(alpha)
{
	alpha = alpha || 0.05;
	result = {
		alpha : alpha,
		df : (this.cols-1)*(this.rows-1)
	};
	var rowS = this.sum(1);
	var colS = this.sum(2);
	var sumR = this.sum(3);
	var exp = colS['*'](rowS)['./'](sumR);
	var chi2 = this['-'](exp)['.^'](2)['./'](exp);
	result.chi2 = chi2.sum(3);
	result.chi2crit = jMath.stat.chi2inv(1-alpha, result.df);
	result.expect = exp;
	result.pvalue = 1-jMath.stat.chi2cdf( result.chi2, result.df );
	return result;		
}

// style 
// 0: diff pop var
// 1: equal pop var
// 2: dependent sample
jMath.fn.compare = function(alpha,type,sample,style,H0)
{
	style = style || 0;
	if ( style == 2 )
	{
		var diff = this['-'](sample);
		return jMath.stat.test_t(alpha,type,diff.mean()[0][0],
			   0, diff.std()[0][0], diff.cols);
	}
	return jMath.stat.compare_t(alpha,type,{
		mean: this.mean()[0][0],
		sigma: this.std()[0][0],
		    n: this.cols
	},{
		mean: sample.mean()[0][0],
		sigma: sample.std()[0][0],
		    n: sample.cols
	},style,H0);
}

jMath.fn.compare_v = function(alpha,type,sample)
{
	return jMath.stat.compare_v( alpha, type, {
		sigma: this.std()[0][0],
		n: this.cols
	},{
		sigma: sample.std()[0][0],
		n: sample.cols
	});
}

jMath.sampling = function(numsamples,numpops,type,opts)
{
	type = type || 'random';
	var idxs =  [];
	switch( type )
	{
		case 'random':
			for( var i = 0; i < numsamples ; i++ )
			{
				var v = Math.min(numpops-1, Math.round(Math.random() * numpops));
				if ( idxs.indexOf(v) == -1 )
				{
					idxs.push(v);
				}
			}
			break;
		case 'systematic':
			// opts : frequency
			var v = Math.min( numpops-1, Math.round(Math.random() * numpops));
			for ( var i = 0 ; i < numsamples ; i++ )
			{
				idxs.push( v );
				v = (v + opts)%numpops;
			}
			break;
		case 'stratified':
			// opts: group ratios	
			for ( var i = 0 ; i < opts.length; i++ )
			{
				idxs.push( jMath.sampling( numsamples * opts[i], numpops[i] * opts[i], 'random' ) );	
			}
			break;
	}
	return idxs;
}

jMath.stat.anova1 = function(alpha)
{
	var result = {
		means: [],
		mean: [],
		sigma: [],
		alpha: alpha,
		numSamples: 0,
		between: { ss: 0, df: 0, ms: 0 },
		within: { ss: 0, df: 0, ms: 0 },
		total: { ss: 0, df: 0 },
	};
	var sum = 0;
	var samples = Array.prototype.slice.call( arguments, 1);
	result.numGroups = samples.length;
	for ( var i = 0 ; i < samples.length; i++)
	{
		result.numSamples += samples[i].cols;
		sum += samples[i].sum(3);
		result.means.push( samples[i].mean()[0][0] );
		result.sigma.push( samples[i].std()[0][0] );
	}
	result.total.df = result.numSamples - 1;
	result.between.df = samples.length - 1;
	result.within.df = result.numSamples - samples.length;
	result.mean = sum/result.numSamples;
	for ( var i = 0 ; i < samples.length; i++)
	{
		result.total.ss += samples[i]['-'](result.mean)['.^'](2).sum(3);
		result.between.ss += Math.pow(result.means[i] - result.mean,2) * samples[i].cols;
	}
	result.within.ss= result.total.ss - result.between.ss;	

	result.total.ms = result.total.ss/result.total.df;
	result.between.ms = result.between.ss/result.between.df;
	result.within.ms = result.within.ss/result.within.df;
	result.F = result.between.ms/result.within.ms;
	result.Fcrit = jMath.stat.finv(1-alpha,result.between.df, result.within.df)
	result.pvalue = 1-jMath.stat.fcdf( result.F, result.between.df, result.within.df);

	result.multicompare = { };
	var Q = jMath.stat.stdrinv(1-alpha, samples.length, result.within.df);
	for ( var i = 0 ; i < result.means.length ; i++ )
	{
		var s1 = 1/samples[i].cols;
		for ( var j = i+1 ; j < result.means.length; j++ )
		{
			var s2 = 1/samples[j].cols;
			result.multicompare[i+'-'+j] = {
				diff: Math.abs(result.means[i] - result.means[j]),
				CR: Q * Math.sqrt( result.within.ms/2 * (s1 + s2))
			};
		}
	}
	return result;
}

jMath.stat.anovarbl = function(alpha)
{
	var result = {
		means: {
			within: [],
			block: [],
			total: 0,
		},
		sigmas: {
			within: [],
			block: []
		},
		alpha: alpha,
		numSamples: 0,
		between: { ss: 0, df: 0, ms: 0 },
		block: { ss: 0, df: 0, ms: 0 },
		error: { ss: 0, df: 0, ms: 0 },
		total: { ss: 0, df: 0 }
	};
	var sum = 0;
	var samples = jMath();
	for ( var i = 1 ; i < arguments.length ; i++ )
	{
		samples[';='](arguments[i]);
	}
	result.numSamples = samples.rows * samples.cols;
	sum = samples.sum(3);

	result.means.within = samples.mean(2);
	result.means.block = samples.mean(1);
	result.sigmas.within = samples.std(2);
	result.sigmas.block = samples.std(1);

	result.total.df = result.numSamples - 1;
	result.between.df = samples.rows - 1;
	result.block.df = samples.cols - 1;
	result.error.df = result.between.df * result.block.df;
	result.means.total = sum/result.numSamples;

	result.total.ss = samples['-']( result.means.total )['.^'](2).sum(3);
	result.between.ss = result.means.within['-'](result.means.total)['.^'](2)['*'](samples.cols).sum(3);
	result.block.ss = result.means.block['-'](result.means.total)['.^'](2)['*'](samples.rows).sum(3);
	result.error.ss= result.total.ss - result.between.ss - result.block.ss;	

	result.total.ms = result.total.ss/result.total.df;
	result.between.ms = result.between.ss/result.between.df;
	result.block.ms = result.block.ss/result.block.df;
	result.error.ms = result.error.ss/result.error.df;

	result.F = {
		between:{
			value: result.between.ms/result.error.ms,
			C: jMath.stat.finv(1-alpha,result.between.df, result.error.df),
		},
		block: {
			value: result.block.ms/result.error.ms,
			C: jMath.stat.finv(1-alpha,result.block.df, result.error.df),
		}
	};
	result.F.between.pvalue = 1-jMath.stat.fcdf( result.F.between.value, 
						result.between.df, 
						result.error.df);
	result.F.block.pvalue = 1-jMath.stat.fcdf( result.F.block.value, 
						result.block.df, result.error.df);

	result.multicompare = { };
	var Q = jMath.stat.stdrinv(1-alpha, samples.rows, result.error.df);
	var CR = Q * Math.sqrt( result.error.ms/samples.cols );
	for ( var i = 0 ; i < result.means.within.length ; i++ )
	{
		for ( var j = i+1 ; j < result.means.within.length; j++ )
		{
			result.multicompare[i+'-'+j] = {
				diff: Math.abs(result.means.within[i] - result.means.within[j]),
				CR: CR
			};
		}
	}
	return result;
}

jMath.stat.anova2 = function(alpha, samples)
{
	var result = {
		means: {
			factorA: [],
			factorB: [],
			cell: [],
			total: 0
		},
		sigmas: {
			factorA: [],
			factorB: [],
			cell: []
		},
		size: [ samples.length, samples[0].length ],
		alpha: alpha,
		numSamples: 0,
		factorA: { ss: 0, df: 0, ms: 0 },
		factorB: { ss: 0, df: 0, ms: 0 },
		interact: { ss: 0, df: 0, ms: 0 },
		within: { ss: 0, df: 0, ms: 0 },
		total: { ss: 0, df: 0 }
	};
	var sum = 0;
	var totalSize = 0;
	var numColumns = [];
	for ( var j = 0 ; j < samples[0].length; j++ )
	{
		result.means.factorB[j] = 0;
		numColumns[j] = 0;
	}
	for ( var i = 0 ; i < samples.length; i++ )
	{
		result.means.cell[i] = [];
		result.sigmas.cell[i] = [];
		var sum = 0;
		var num = 0;
		for ( var j = 0 ; j < samples[i].length; j++ )
		{
			var s = samples[i][j].sum(3);
			result.means.factorB[j] += s;
			numColumns[j] += samples[i][j].cols;
			sum += s;
			num += samples[i][j].cols;
			result.within.df += samples[i][j].cols - 1;
			result.means.cell[i][j] = samples[i][j].mean()[0][0];
			result.sigmas.cell[i][j] = samples[i][j].std()[0][0];
		}
		result.means.total += sum;
		totalSize += num;
		result.means.factorA[i] = sum/num;
	}
	result.means.total /= totalSize;

	for ( var j = 0 ; j < samples[0].length; j++ )
	{
		result.means.factorB[j] /= numColumns[j];
	}

	result.total.df = totalSize - 1;
	result.factorA.df = result.size[0] - 1;
	result.factorB.df = result.size[1] - 1;
	result.interact.df = result.factorA.df * result.factorB.df;

	for ( var i = 0 ; i < samples.length; i++ )
	{
		for ( var j = 0 ; j < samples[i].length; j++ )
		{
			result.total.ss += samples[i][j]['-'](result.means.total)['.^'](2).sum(3);
			result.factorA.ss += Math.pow(result.means.factorA[i] - result.means.total,2) * samples[i][j].cols;
			result.interact.ss += Math.pow( result.means.cell[i][j] - result.means.factorA[i] - result.means.factorB[j] + result.means.total, 2 )
									* samples[i][j].cols;
		}
	}
	for ( var j = 0 ; j < samples[0].length ; j++ )
	{
		for ( var i = 0 ; i < samples.length; i++ )
		{
			result.factorB.ss += Math.pow(result.means.factorB[j] - result.means.total,2) * samples[i][j].cols;
		}
	}
	result.within.ss = result.total.ss - result.factorA.ss - result.factorB.ss - result.interact.ss;

	result.total.ms = result.total.ss/result.total.df;
	result.factorA.ms = result.factorA.ss/result.factorA.df;
	result.factorB.ms = result.factorB.ss/result.factorB.df;
	result.interact.ms = result.interact.ss/result.interact.df;
	result.within.ms = result.within.ss/result.within.df;

	result.F = {
		factorA:{
			value: result.factorA.ms/result.within.ms,
			C: jMath.stat.finv(1-alpha,result.factorA.df, result.within.df),
		},
		factorB:{
			value: result.factorB.ms/result.within.ms,
			C: jMath.stat.finv(1-alpha,result.factorB.df, result.within.df),
		},
		interact:{
			value: result.interact.ms/result.within.ms,
			C: jMath.stat.finv(1-alpha,result.interact.df, result.within.df),
		},
	};
	result.F.factorA.pvalue = 1-jMath.stat.fcdf( result.F.factorA.value, 
						result.factorA.df, result.within.df);
	result.F.factorB.pvalue = 1-jMath.stat.fcdf( result.F.factorB.value, 
						result.factorB.df, result.within.df);
	result.F.interact.pvalue = 1-jMath.stat.fcdf( result.F.interact.value, 
						result.interact.df, result.within.df);

	result.multicompare = {
  		factorA: {},
		factorB: {},
	};
	var Q = jMath.stat.stdrinv(1-alpha, result.size[0], result.within.df);
	result.multicompare.factorA.CR = Q * Math.sqrt( result.within.ms/(result.size[1]*samples[0][0].cols));
	for ( var i = 0 ; i < result.size[0] ; i++ )
	{
		for ( var j = i+1 ; j < result.size[0] ; j++ )
		{
			result.multicompare.factorA[i+'-'+j] =
				Math.abs(result.means.factorA[i] - result.means.factorA[j]);
		}
	}

	Q = jMath.stat.stdrinv(1-alpha, result.size[1], result.within.df);
	result.multicompare.factorB.CR = Q * Math.sqrt( result.within.ms/(result.size[0]*samples[0][0].cols));

	for ( var i = 0 ; i < result.size[1] ; i++ )
	{
		for ( var j = i+1 ; j < result.size[1] ; j++ )
		{
			result.multicompare.factorB[i+'-'+j] =
				Math.abs(result.means.factorB[i] - result.means.factorB[j]);
		}
	}

	return result;
}

jMath.fn.corrcoef = function(alpha)
{
	alpha = alpha || 0.05;
	var sigma = this.std(true,1);
	var mean = this.mean(1);
	var cov = this['-'](mean.repeatRow( this.rows - 1 ));
	var s = 0;
	for ( var r = 0 ; r < cov.length; r++)
	{
		s += cov[r][0] * cov[r][1];
	}
	var result = {
		alpha : alpha,
		r: s/(sigma[0][0]*sigma[0][1]*cov.length),
		df: this.rows - 2,
		mean: mean[0],
		sigma: sigma[0]
	}	
	var denum = Math.sqrt( (1-result.r*result.r)/result.df );
	result.testscore = result.r/denum;
	result.pvalue = jMath.stat.tcdf( result.testscore, result.df );
	if ( result.r> 0 ) result.pvalue = 1 - result.pvalue;
	result.pvalue *= 2;
	var t = jMath.stat.tinv( result.alpha/2, result.df );
	result.tscores = [ t, -t ];

	// standard error
	denum = 1/Math.sqrt( this.rows - 3 );
	var z = jMath.stat.norminv( result.alpha/2, 0, 1 );
	var fisherz = Math.atanh(result.r); 
	result.ci = [ fisherz +  z * denum, fisherz - z*denum ].map( Math.tanh );
	return result;
}

jMath.fn.regressline = function(alpha, extra)
{
	alpha = alpha || 0.05;
	var result = this.corrcoef(alpha);
	var X = this.slice(':', 0);
	var Y = this.slice(':', 1);
	result.n = this.rows;
	result.r = {
		value: result.r,
		ci: result.ci,
		tscores : result.tscores,
		tvalue:  result.testscore
	};
	delete result.ci;
	delete result.pvalue;
	delete result.tscores;
	delete result.testscore;

	result.line = [0,0];
	result.line[1] = {
		value: result.r.value * result.sigma[1]/result.sigma[0]
	};
	result.line[0] = {
		value: result.mean[1] - result.line[1].value * result.mean[0]
	};
	result.predict = {
		value: X['*'](result.line[1].value)['+'](result.line[0].value)
	};
	result.residues = Y['-'](result.predict.value);
	result.anova = {
		sse: result.residues['.^'](2).sum(3), 
		ssr: result.predict.value['-'](result.mean[1])['.^'](2).sum(3)
	};
	result.anova.sst = result.anova.sse + result.anova.ssr;
	result.r2 ={
		value: result.anova.ssr/result.anova.sst,
		f : result.anova.ssr/(result.anova.sse/result.df),
		fcrit : jMath.stat.finv( 1-alpha, 1, result.df )
	}
	result.r2.pvalue = 1 - jMath.stat.fcdf( result.r2.f, 1, result.df );

	result.se = Math.sqrt(result.anova.sse / result.df);
	var talpha2 = jMath.stat.tinv( 1-alpha/2, result.df );
	result.tscores = [-talpha2, talpha2];
	result.predict.ci = [];
	result.predict.pi = [];
	for ( var i = 0 ; i < this.rows ; i++ )
	{
		var yhat = result.predict.value[i][0];
		var sx = (1+ Math.pow((this[i][0]-result.mean[0])/result.sigma[0],2))/result.n;
		var s = result.se * talpha2 * Math.sqrt(sx);
		result.predict.ci.push( [ yhat - s, yhat +s ] );
		s = result.se * talpha2 * Math.sqrt(1+sx);
		result.predict.pi.push( [ yhat - s, yhat +s ] );
	}
	if ( extra )
	{
		result.extra = [];
		for ( var i = 0 ; i < extra.length ; i++ )
		{
			var ev = {};
			var yhat = result.line[0].value + result.line[1].value * extra[i];
			ev.yhat = yhat;	
			var sx = (1+ Math.pow((extra[i]-result.mean[0])/result.sigma[0],2))/result.n;
			var s = result.se * talpha2 * Math.sqrt(sx);
			ev.ci = [ yhat - s, yhat +s ];
			s = result.se * talpha2 * Math.sqrt(1+sx);
			ev.pi = [ yhat - s, yhat +s ];
			result.extra.push( ev );
		}
	}
	result.line[1].se = result.se/(Math.sqrt(result.n)*result.sigma[0]);
	result.line[0].se = result.line[1].se*Math.sqrt(X['.^'](2).mean(3));
	for( var i = 0 ; i < result.line.length; i++ )
	{
		result.line[i].t = result.line[i].value / result.line[i].se;
		result.line[i].ci = [
			result.line[i].value + result.tscores[0]* result.line[i].se,
			result.line[i].value + result.tscores[1] * result.line[i].se ];
		if ( result.line[i].t > 0 )
		{
			result.line[i].pvalue = 1 - jMath.stat.tcdf( result.line[i].t, result.df );
		}
		else
		{
			result.line[i].pvalue = jMath.stat.tcdf( result.line[i].t, result.df );
		}
		result.line[i].pvalue *= 2;
	}
	return result;
}

jMath.regress = function( X, Y, alpha, extra )
{
	X = jMath.ones(Y.rows, 1)[':='](X);
	var XT = X.transpose();
	var XXinv = XT['*'](X).inv();
	var B = XXinv['*'](XT)['*'](Y);
	var mY = Y.mean(3);
	var Yhat = X['*'](B);
	var SSE = Y['-'](Yhat)['.^'](2).sum(3);
	var SST = Y['-'](mY).pow(2).sum(3);
	var SSR = SST - SSE;
	var rdf = X.cols - 1;
	var df = Y.rows - rdf - 1;

	var MSE = SSE/df;
	var MST = SST/(Y.rows-1);
	var MSR = SSR/rdf;

	var r2 = 
	{
		value: SSR/SST,
		adjust: 1 - MSE/MST,
			F: MSR/MSE,
		Fcrit: jMath.stat.finv(0.95,rdf,df)
	};
	r2.pvalue = 1 - jMath.stat.fcdf( r2.F, rdf, df );
	var se = Math.sqrt(MSE);
	var C = XXinv['*'](MSE);

	alpha = alpha || 0.05;

	var result = {
		B: {
			value: B,
			se: C.diag().sqrt().transpose(),
			df: df,
		},
	    C: C,
		anova: {
			SST: SST,
			SSR: SSR,
			SSE: SSE,
			MST: MST,
			MSR: MSR,
			MSE: MSE // Standard Error
		},
		df: {
			SSR: rdf,
			SSE: df,
			SST: Y.rows - 1
		},
		r2 : r2,
		r: Math.sqrt( r2.value ),
	};
	result.B.t = result.B.value['./']( result.B.se );
	result.B.tcrit = [ jMath.stat.tinv( alpha/2, result.B.df ), 
						 jMath.stat.tinv( 1-alpha/2, result.B.df) ];
	result.B.ci = [];
	result.B.pvalue = [];
	for ( var i = 0 ; i  < result.B.value.rows ; i++ )
	{
		result.B.ci.push([
			result.B.value[i][0] + result.B.tcrit[0]*result.B.se[i][0],
			result.B.value[i][0] + result.B.tcrit[1]*result.B.se[i][0]
		]);
		var pv = jMath.stat.tcdf( result.B.t[i][0], result.B.df );
		if ( result.B.t[i][0] > 0 )
		{
			pv = 1-pv;
		}
		result.B.pvalue.push(2*pv);
	}

	if ( extra )
	{
		var ex;
		if ( typeof extra == 'boolean' )
		{
			eX = X;
		}
		else
		{
			eX = jMath.ones(extra.length, 1)[':='](jMath(extra));
		}
		var yHat = eX['*'](result.B.value);
		result.yHat = X.removeColumns(0)[':='](yHat);
		result.residues = Y['-'](yHat);
		// autocorrelation: 0:+, 2: no, 4: -
		result.durbinWatson = result.residues.diff().pow(2).sum(3) / result.residues.pow(2).sum(3);

		var b= eX['*'](result.C)['*'](eX.transpose()).diag().transpose();
		var eSE = b.sqrt();
		var eSEP = b['+'](MSE).sqrt();
		var ME = jMath(result.B.tcrit).transpose()['*'](eSE.transpose());
		var MEP = jMath(result.B.tcrit).transpose()['*'](eSEP.transpose());
		result.extra = [
		];
		for ( var i = 0 ; i < yHat.rows; i++ )
		{
			var v = {
				yhat: yHat[i][0],
				ci: [ yHat[i][0] + ME[0][i], yHat[i][0] + ME[1][i] ],
				pi: [ yHat[i][0] + MEP[0][i], yHat[i][0] + MEP[1][i] ]
			};
			result.extra.push(v);
		}
	}
	return result;
}

jMath.fn.regress = function(alpha, extra)
{
	var Y = this.slice(':', this.cols - 1);
	var X = this.slice(':', '0:' + (this.cols-2));
	return jMath.regress(X,Y,alpha,extra);
}

jMath.fn.vif = function()
{
	var result = [];
	var Xs = this.slice(':', '0:' + (this.cols-2));
	for ( var i = 0 ; i < this.cols - 1 ; i++ )
	{
		var list = [];
		var Y = Xs.slice(':', i);		
		var X = jMath.ones(Y.rows, 1)[':=']( Xs.removeColumns(i) );
		var XT = X.transpose();
		var XXinv = XT['*'](X).inv();
		var B = XXinv['*'](XT)['*'](Y);
		var Yhat = X['*'](B);
		var mY = Y.mean(3);
		var Yhat = X['*'](B);
		var SSE = Y['-'](Yhat)['.^'](2).sum(3);
		var SST = Y['-'](mY).pow(2).sum(3);
		var SSR = SST - SSE;
		result.push( 1/(1-SSR/SST));
	}
	return result;
}

jMath.fn.modelbuilding = function(method)
{
	var alpha =0.05;
	method = method || 'stepwise';
	var result;
	var Y = this.slice(':', this.cols-1);		
	// General Stepwise Regression
	if ( method == 'stepwise' || method == "forward" )
	{
		var corrcoeff = [];
		var X;
		var skipCoeffCheck = method == 'forward'; 
		for ( var i = 0 ; i < this.cols - 1; i++ )
		{
			X = this.slice(':', i);
			result = jMath.regress(X,Y,alpha);
			corrcoeff.push([result.r, i, result.r2.pvalue]);
		}
		corrcoeff = corrcoeff.sort(function(v1,v2){
			if ( v1[0] > v2[0] ) return -1;
			else if ( v1[0] < v2[0] ) return 1;
			return 0;
		});
		if ( corrcoeff[0][2] > alpha ) return null;
		X = this.slice(':', corrcoeff[0][1] );
		for ( var i = 1 ; i < corrcoeff.length; i++ )
		{
			X = X[':='](this.slice(':', corrcoeff[i][1]));	
			result = jMath.regress(X,Y,alpha);
			var ok = true;
			if ( result.r2.pvalue > alpha )
			{
				ok = false;
			}
			else if ( !skipCoeffCheck )
			{
				for ( var j = 1 ; ok && j < result.B.pvalue.length; j++ )
				{
					if ( result.B.pvalue[j] > alpha )
					{
						ok = false;
						break;
					}
				}
			}
			if ( !ok )
			{
				X = X.removeColumns( X.cols - 1 );
				corrcoeff[i][1] = -1;
			}
		}
		result = jMath.regress(X,Y,alpha);
		result.indexs = corrcoeff.filter(function(v){
			return v[1] != -1;	
		}).map( function(v){
			return v[1];		
		});
	}
	else if ( method == 'bestsubset' )
	{
		var all = this.regress();
		var numCols = this.cols-1;
		var num = Math.pow(2, numCols);
		var r2as = [];
		var bits;
		var X;
		for ( var i = 1 ; i < num ; i++ )
		{
			bits = i;		
			X = null;
			var cnt = 0;
			for ( var j = 0 ; j < numCols; j++ )
			{
				if ( (bits & 1) == 1 )
				{
					if ( X == null ) X = this.slice(':', j);
					else X[':=']( this.slice(':',j) );
					cnt++;
				}
				bits >>= 1;
			}	
			var result = jMath.regress(X,Y,alpha);
			var Cp = result.anova.SSE/all.anova.MSE - ( this.rows - 2 * (cnt+1) );
			r2as.push( [result.r2.adjust, i, cnt, Cp] );
		}
		r2as = r2as.filter(function(v){
			return v[3] <= v[2] + 1;
		}).sort(function(v1,v2){
			if ( v1[0] > v2[0] ) return -1;
			else if ( v1[0] < v2[0] ) return 1;
			return 0;
		});
		num = r2as.length - 1;
		for( var i = 0 ; i < num; i++ )
		{
			var diff = r2as[i][0] - r2as[i+1][0];
			if ( diff > 0.01 )
			{
				num = i+1;
				break;
			}
		}
		if ( num )
		{
			r2as = r2as.slice(0,num).sort(function(v1,v2){
				if ( v1[2] < v2[2] ) return -1;
				else if ( v1[2] > v2[2] ) return 1;
				return 0;
			});
		}
		bits = r2as[0][1];
		var indexs = [];
		X = null;
		for ( var j = 0 ; j < numCols ; j++ )
		{
			if ( (bits & 1) == 1 )
			{
				if ( X == null ) X = this.slice(':', j);
				else X[':=']( this.slice(':',j) );
				indexs.push(j);
			}
			bits >>= 1;
		}
		result = jMath.regress(X,Y,alpha,true);
		result.indexs = indexs;
	}
	if ( X != null )
	{
		if ( result.residues )
		{
			result.residues = {
			  	value: result.residues,
				corrcoef : []
			}
		}
		else
		{
			var X = jMath.ones(Y.rows, 1)[':='](X);
			result.residues = {
			  	value: Y['-'](X['*'](result.B.value)),
				corrcoef : []
			};
		}
		for ( var i = 0 ; i < result.indexs.length; i++ )
		{
			X = this.slice( ':', result.indexs[i] )[':='](result.residues.value.abs());
			result.residues.corrcoef.push( X.corrcoef(alpha) );
		}
	}
	return result;
}

jMath.fn.forecast = function(type,p)
{
	var list = [];
	var mad = 0;
	if ( type == 'sma' )
	{
		p = p || 3;	
		for ( var i = p ; i < this.rows ; i++ )
		{
			var sum = 0;
			for ( var j = i - p ; j < i ; j++ )
			{
				sum += this[j][1];
			}
			sum /= p;
			list.push( [ this[i][0], sum] );
			mad += Math.abs(this[i][1] - sum);
		}
		mad /= this.rows - p;
	}
	else if ( type == 'wma' )
	{
		var ps  = p.reduce(function(a,v){ return a+v; });
		for ( var i = p.length ; i < this.rows ; i++ )
		{
			var sum = 0;
			for ( var j = i - p.length, k=p.length-1 ; j < i ; j++, k-- )
			{
				sum += this[j][1] * p[k];
			}
			sum /= ps
			list.push( [ this[i][0], sum] );
			mad += Math.abs(this[i][1] - sum);
		}
		mad /= this.rows - p.length;
	}
	else if ( type == 'exp' )
	{
		if ( typeof p == 'number' ) p = [ p ];
		for ( var i = 0 ; i < p.length; i++ )
		{
			list.push( [ this[i][0], this[i][1] ] );
		}
		for ( var i = p.length  ; i < this.rows ; i++  )
		{
			list[i] = [ this[i][0], list[i-1][1] ];
			for ( var j = i - p.length, k=p.length-1 ; j < i ; j++, k-- )
			{
				list[i][1] += p[k] * ( this[j][1] - list[j][1] );
			}
			mad += Math.abs(this[i][1] - list[i][1]);
		}
		mad /= this.rows - p.length;
		for ( var i = 0 ; i < p.length; i++ )
		{
			list.splice(0,p.length);
		}
	}
	else if ( type == 'expta' )
	{
		var alpha = p;
		var beta = arguments[2];
		list.push( [ this[0][0], this[0][1] ] );
		var T = 0;
		var FIT = this[0][1];
		for ( var i = 1  ; i < this.rows ; i++  )
		{
			FIT = list[i-1][1] + T;
			list[i] = [ this[i][0], FIT + alpha * ( this[i-1][1] - FIT )];
			T = beta * ( list[i][1] - list[i-1][1] ) +  ( 1 - beta ) * T;
			mad += Math.abs(this[i][1] - list[i][1]);
		}
		mad /= this.rows - 1;
		list.splice(0,1);
	}
	else if ( type == 'season' )
	{
		var sIdx = p>>1;
		var numRepeat = Math.floor(this.rows/p);
		var nsf = jMath.zeros(numRepeat, p);
		if ( p%2 == 0 )
		{
			var len = p - 1;
			for ( var i = sIdx; i < this.rows - sIdx; i++ )
			{
				var j = i-sIdx;
				var s = this[j++][1];
				for ( var k = 0 ; k < len ; j++, k++ )
				{
					s += 2*this[j][1];
				}
				s += this[j][1];
				var r = Math.floor(i/p);
				var c = i%p;
				nsf[r][c] = this[i][1]/(s/(2*p));
			}
		}
		else
		{
			for ( var i = sIdx; i < this.rows - sIdx ; i++ )
			{
				var s = 0;
				for ( var j = i-sIdx, k = 0 ; k < p ; j++, k++ )
				{
					s += this[j][1];
				}
				var r = Math.floor(i/p);
				var c = i%p;
				nsf[r][c] = this[i][1]/(s/p);
			}
		}
		var Y = this.slice(':',1);
		var X = this.slice(':',0);
		numRepeat--;
		nsf = nsf.sum()['./'](numRepeat);
	   	nsf = nsf['./'](nsf.sum(3))['*'](p).transpose();
		var repnsf = nsf.repeatRow(numRepeat);
		var result = jMath.regress( X, Y['./'](repnsf), 0.05, true);
		result.interval = p;
		result.nsf = nsf;
		for ( var i = 0 ; i < this.rows; i++ )
		{
			result.yHat[i][1] *= repnsf[i][0];
			result.residues[i][0] = this[i][1] - result.yHat[i][1];
		}
		console.log(result.yHat.toString());
		result.mad = result.residues.abs().mean(3);
		return result;
	}
	else if ( type == 'dummy')
	{
		var X = this.slice(':', 0);
		for ( var i = 0 ; i < this.rows ; i++ )
		{
			var v = i%p;				
			for ( var k = 1 ; k < p ; k++ )
			{
				X[i][k] = 0;
			}
			if ( v )
			{
				X[i][v] = 1;
			}
		}
		X.cols = p;
		var XY = X[':=']( this.slice(':',1));
		var result = XY.modelbuilding('bestsubset');
		result.mad = result.residues.value.abs().mean(3);
		result.XY = XY;
		return result;
	}
	var result = jMath(list);
	result.mad = mad;
	return result;
}

})(jMath);

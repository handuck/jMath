(function(jStat,Math){

jStat.extend( jStat.fn, {
	stdev_s : function(whole) {
		var r = this.rows();
		var c = this.cols();
		if ( r == 1 ) return this.stdev(true);
		else
		{
			if ( whole ) 
			{
				var n = r * c ;
				var coeff = Math.sqrt( n / ( n-1) );
				return this.stdev(true) * coeff;
			}
			else
			{
				var n = this.rows();
				var coeff = Math.sqrt( n / ( n-1) );
				return this.stdev().map(function(v){
					return v * coeff;
				});
			}
		}
	},
	coeffvar_s : function() {
		if ( this.rows() == 1 )
		{
			return this.stdev(true) / this.mean();
		}
		else
		{
			var std = this.stdev_s();
	    	return this.mean().map(function(v,idx){
				return std[idx]/v;
			});
		}		
	},

	sort: function() { 
		return jStat( this[0].sort(function(arg1,arg2){
				if ( arg1 == arg2 ) return 0;
				else if ( arg1 > arg2 ) return 1;
				return -1;
		}) );
	},

	percentile: function(p) {
		var list = this.sort();
		var v = (this[0].length-1) * p;
		var idx = Math.floor(v);
		var decimal = v - idx;
		return list[0][idx] + decimal * ( list[0][idx+1] - list[0][idx] );
	},

	percentrank: function(v){
		var list = this.sort();
		if ( list[0][0] > v || list[0][list[0].length-1] < v ) return Number.NaN;
		for( var i = 1 ; i < list[0].length ; i++ )
		{
			if ( list[0][i]	>= v )
			{
				break;
			}
		}
		i--;
		var idx = i + (v-list[0][i])/(list[0][i+1] - list[0][i]);
		return idx / (list[0].length-1);
	}
             
});

})(jStat,Math);

(function(global){

global.jMath = global.m$ = function()
{
	return new global.jMath.fn.init(arguments);
}

m$.epsilon = 2.220446049250313e-16;

m$.fn = m$.prototype = {
	length: 0,
	constructor: global.m$,
	splice: Array.prototype.splice,
	init : function(args)
	{
		if ( args.length == 0 )
		{
			this.length = 0;
			this.rows = 0;
			this.cols = 0;
		}
		else if ( typeof(args[0]) == 'string' )
		{
			this.parse(args[0]);
		}
		else
		{
			var list;
			if ( typeof(args[0]) == 'number' )
			{
				list = Array.prototype.slice.call(args,0);						
			}
			else
			{
				list = args[0];
			}
			this.rows = list.length;
			this.cols = list[0].length;
			if ( !this.cols )
			{
				this.rows = 1;		
				this.cols = list.length;
				this[0] = list;
			}
			else
			{
				for( var i = 0 ; i < list.length; i++ )
				{
					this[i] = list[i];
					if ( this.cols < list[i].length )
					{
						this.cols = list[i].length;
					}
				}
			}
		}
		this.length = this.rows;
	},

	parseElement: function(value,dir)
	{
		if ( value.indexOf(':') != -1 )
		{
			var list = [];
			var v = value.split(/:/);
			var init, step, end;
			if ( v.length == 2 )
			{
				init = v[0].length > 0 ? parseFloat(v[0]) : 0;
				end = v[0].length > 0 ? parseFloat(v[1]) : NaN;
				if ( dir && isNaN(end) )
				{
					end = dir == 1 ? this.rows-1 : this.cols-1;
				}
				step = 1;
			}	
			else
			{
				init = parseFloat(v[0]);
				step = parseFloat(v[1]);
				end = parseFloat(v[2]);
			}
			if( step > 0 )
			{
				for ( var v = init ; v <= end ; v+=step )
				{
					list.push( v );
				}
			}
			else if ( step < 0 )
			{
				for ( var v = init ; v >= end ; v+=step )
				{
					list.push( v );
				}
			}
			return list;
		}
		else
		{
			return [parseFloat(value)];
		}
	},

	parse: function(value)
	{
		var list = null;
		var rows = value.split(/;/);	
		this.rows = rows.length;	
		this.cols = 0;
		for ( var i = 0 ; i < rows.length; i++ )
		{
			var cols = rows[i].trim().split(/\s+/);		
			this[i] = [];
			for ( var j = 0 ; j < cols.length; j++ )
			{
				var list = this.parseElement(cols[j]);
				for ( var c = 0 ; c < list.length; c++ )
				{
					this[i].push( list[c] );
				}
				if ( this.cols < this[i].length )
				{
					this.cols = this[i].length;
				}
			}
		}
	},

	toString: function(v)
	{
		var str = "";
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			if ( v )
			{
				str += parseFloat(this[r][0]).toFixed(v);
			}
			else
			{
				str += this[r][0];
			}
			for( var c = 1 ; c < this.cols; c++ )
			{
				if ( v )
				{
					str += "\t" + parseFloat(this[r][c]).toFixed(v);
				}
				else
				{
					str += "\t" + this[r][c];
				}
			}
			str += "\n";
		}
		return str.trimRight();
	},

	toArray: function()
	{
		return Array.prototype.slice.call(this,0); 
	},

	valueOf: function()
	{
		if ( this.rows == 0 || this.cos == 0 ) return 0;
		else if ( this.rows == 1 && this.cols == 1 ) return this[0][0];	
		return this;
	}
}

m$.fn.init.prototype = m$.fn;

m$.fn.removeColumns = function(c)
{
	var list = [];
	var cols = arguments.length == 1 ? [c] : Array.prototype.slice.call(arguments, 0);
	for ( var j = 0 ; j < this.rows ; j++ )
	{
		var sub = [];
		for ( var k = 0 ; k < this.cols; k++ )
		{
			if ( cols.indexOf(k) != -1 ) continue;	
			sub.push( this[j][k] );
		}
		list.push( sub );
	}
	return jMath(list);
}

m$.fn.slice = function( rIdxs, cIdxs )
{
	var rRange, cRange; 
	if ( typeof rIdxs == 'number' )
	{
		rRange = [rIdxs];
	}
	else if ( Array.isArray(rIdxs) )
	{
		rRange = rIdxs;
	}
	else
	{
		var range = rIdxs.split(/:/).map( function(v){
			return isNaN(parseFloat(v)) ? v : parseFloat(v);
		});
		if ( range.length == 1 )
		{
			range[0] = range[1] = parseInt(rIdxs);
		}
		else
		{
			if ( range[0] === '' ) range[0] = 0;
			if ( range[1] === '' || range[1] === 'end' ) 
			{
				range[1] = this.rows - 1;
			}
		}
		rRange = [];
		for ( var i = range[0] ; i <= range[1] ; i++ )
		{
			rRange.push( i );
		}
	}
	if ( typeof cIdxs == 'number' )
	{
		cRange = [cIdxs];
	}
	else if ( Array.isArray(cIdxs) )
	{
		cRange = cIdxs;
	}
	else
	{
		var range = cIdxs.split(/:/).map( function(v){
			return isNaN(parseFloat(v)) ? v : parseFloat(v);
		});
		if ( range.length == 1 ) 
		{
			range[0] = range[1] = parseInt(cIdxs);
		}
		if ( range[0] === '' ) range[0] = 0;
		if ( range[1] === '' || range[1] === 'end' ) 
		{
			range[1] = this.cols - 1;
		}
		cRange = [];
		for ( var i = range[0] ; i <= range[1] ; i++ )
		{
			cRange.push( i );
		}
	}
	var list = [];
	for ( var i = 0 ; i < rRange.length ; i++ )
	{
		var item = [];
		for ( var j = 0 ; j < cRange.length; j++ )
		{
			item.push(this[ rRange[i] ][ cRange[j] ]);
		}
		list.push(item);
	}
	return m$(list);
}

m$.fill = function(rows,cols, cbfunc)
{
	var list = [];
	for ( var r = 0 ; r < rows ; r++)
	{
		list.push([]);
		for ( var c = 0 ; c < cols ; c++ )
		{
			list[r][c] = cbfunc(r,c);
		}
	}
	return m$(list);
}

m$.eye = function(rows,cols)
{
	var list = [];
	rows = rows || 1;
	cols = cols || rows;
	for ( var r = 0 ; r < rows ; r++)
	{
		list.push([]);
		for ( var c = 0 ; c < cols ; c++ )
		{
			list[r][c] = r == c ? 1 : 0;
		}
	}
	return m$(list);
}

m$.zeros = function(rows,cols)
{
	return m$.fill(rows,cols,function(r,c){
		return 0;	
	});
}

m$.ones = function(rows,cols)
{
	return m$.fill(rows,cols,function(r,c){
		return 1;	
	});
}

m$.rand = function(rows,cols,from,end)
{
	from = from || 0;
	end = end || 1;
	var range = end - from;
	return m$.fill(rows,cols,function(r,c){
		return from + Math.random() * range;
	});
}

m$.randn = function(rows, cols, mu, sigma) 
{
	rows = rows || 1;
	cols = cols || 1;
	mu = mu || 0;
	sigma = sigma || 1;
	var u, v, x, y, q, mat;
	return m$.fill(rows,cols, function(r,c){
		var u, v, x, y, q;
		do {
			u = Math.random();
			v = 1.7156 * (Math.random() - 0.5);
			x = u - 0.449871;
			y = Math.abs(v) + 0.386595;
			q = x*x + y * (0.19600 * y - 0.25472 * x);
		} while(q > 0.27597 && (q > 0.27846 || v*v > -4 * Math.log(u) * u*u));
		return mu + (v/u)*sigma;
	});
}

m$.fn.map = function(cbfunc)
{
	var list = [];
	for ( var r = 0 ; r < this.rows ; r++)
	{
		list.push([]);
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			list[r][c] = cbfunc(this[r][c],r,c);
		}
	}
	return m$(list);
}

m$.fn.forEach = function(cbfunc)
{
	for ( var r = 0 ; r < this.rows ; r++)
	{
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			cbfunc(this[r][c],r,c);
		}
	}
}

m$.fn.diag = function()
{
	var list = [];
	var size = Math.min( this.rows, this.cols );
	for ( var i = 0 ; i < size ; i++ )
	{
		list.push( this[i][i] );
	}
	return jMath(list);
}

m$.fn['.^'] = function(x)
{
	return this.map( function(v,r,c){
		return Math.pow(v,x);
	});
}

m$.fn['^'] = function(x)
{
	if ( this.rows != this.cols ) 
	{
		throw new RangeError('must be square matrix');
	}
	var ret = this;
	for ( var i = 1 ; i < x ; i++ )
	{
		ret = ret['*'](this);	
	}
	return ret;
}

m$.fn['.*'] = function(x)
{
	if ( Array.isArray(x) )
	{
		x = m$(x);
	}
	if ( x instanceof m$.fn.init )
	{
		return this.map( function(v,r,c){
			return v * x[r][c];
		});
	}
	else
	{
		return this.map( function(v,r,c){
			return v * x;
		});
	}
}

m$.fn['*'] = m$.fn.mul = function(x)
{
	if ( Array.isArray(x) )
	{
		x = m$(x);
	}
	if ( x instanceof m$.fn.init )
	{
		if ( this.cols != x.rows )
		{
			throw new RangeError("the number of column must be the same as the number of rows");
		}
		var list = [];
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			list.push([]);
			for ( var c = 0; c < x.cols ; c++ )
			{
				var sum = 0;
				for ( var k = 0 ; k < x.rows ; k++ )
				{
					sum += this[r][k] * x[k][c];
				}
				list[r].push(sum);
			}
		}
		return m$(list);
	}
	else
	{
		return this.map( function(v,r,c){
			return v * x;
		});
	}
}

m$.fn['+'] =  m$.fn.add = function(x)
{
	if ( Array.isArray(x) )
	{
		x = m$(x);
	}
	if ( x instanceof m$.fn.init )
	{
		return this.map( function(v,r,c){
			return v + x[r][c];
		});
	}
	else
	{
		return this.map( function(v,r,c){
			return v + x;
		});
	}
}

m$.fn['-'] = m$.fn.sub = function(x)
{
	if ( Array.isArray(x) )
	{
		x = m$(x);
	}
	if ( x instanceof m$.fn.init )
	{
		return this.map( function(v,r,c){
			return v - x[r][c];
		});
	}
	else
	{
		return this.map( function(v,r,c){
			return v - x;
		});
	}
}

m$.fn['./'] = function(x)
{
	if ( Array.isArray(x) )
	{
		x = m$(x);
	}
	if ( x instanceof m$.fn.init )
	{
		return this.map( function(v,r,c){
			return v / x[r][c];
		});
	}
	else
	{
		return this.map( function(v,r,c){
			return v / x;
		});
	}
}

m$.fn['/'] = function(x)
{
	if ( Array.isArray(x) )
	{
		x = m$(x);
	}
	if ( x instanceof m$.fn.init )
	{
		return this['*'](x.inv());
	}
	else
	{
		return this.map( function(v,r,c){
			return v / x;
		});
	}
}

m$.fn.inv = function()
{
	var b = jMath.eye(this.rows, this.cols);
	var c = jMath.gauss_jordan(this, b);
	var obj = [];
	for ( var i = 0 ; i < this.rows; i++) 
	{
		obj[i] = [];
		for ( var j = this.cols - 1; j < c[0].length; j++)
		{
			obj[i][j - this.cols] = c[i][j];
		}
	}
	return jMath(obj);
}

m$.fn.det = function()
{
 	var alen = this.rows;
	var alend = alen * 2;
	var vals = new Array(alend);
	var rowshift = alen - 1;
	var colshift = alend - 1;
	var mrow = rowshift - alen + 1;
	var mcol = colshift;
	var result = 0;
    // check for special 2x2 case
	if (alen === 2) {
		return this[0][0] * this[1][1] - this[0][1] * this[1][0];
	}
	for (var i = 0; i < alend; i++) {
		vals[i] = 1;
	}
	for (var i = 0; i < alen; i++) {
		for (var j = 0; j < alen; j++) {
			vals[(mrow < 0) ? mrow + alen : mrow ] *= this[i][j];
			vals[(mcol < alen) ? mcol + alen : mcol ] *= this[i][j];
			mrow++;
			mcol--;
		}
		mrow = --rowshift - alen + 1;
		mcol = --colshift;
	}
	for (i = 0; i < alen; i++) {
		result += vals[i];
	}
	for (; i < alend; i++) {
		result -= vals[i];
	}
	return result;		
}

// TODO: not working properly.
m$.fn.qr = function() {
	var Q = jMath.zeros(this.rows, this.cols);
	var R = jMath.zeros(this.cols, this.cols);
	var m = Math.min( this.rows, this.cols );
	var s = 0;
	for ( var c = 0 ; c < this.cols ; c++)
	{
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			Q[r][c] = this[r][c];
		}
		for ( var j = 0 ; j < c; j++ )
		{
			s = 0;
			for ( var k = 0 ; k < this.rows; k++ )
			{
				s += this[k][c] * Q[k][j]
			}
			for ( var k = 0 ; k < this.rows; k++ )
			{
				Q[k][c] -= s * Q[k][j];
			}
		}
		s = 0;
		for ( r = 0 ; r < this.rows; r++ )
		{
			s += Q[r][c] * Q[r][c];
		}
		s = Math.sqrt(s);
		for ( r = 0 ; r < this.rows; r++ )
		{
			Q[r][c] /= s;
		}
	}
	for ( var c = 0 ; c < R.cols ; c++ )
	{
		for ( var r = 0 ; r <= c ; r++ )
		{
			for ( var i = 0 ; i < R.rows ; i++ )
			{
				R[r][c] += this[i][c] * Q[i][r];			
			}
		}
	}
	return { Q: Q, R: R}
}


m$.fn.svd = function()
{
    var temp;
//Compute the thin SVD from G. H. Golub and C. Reinsch, Numer. Math. 14, 403-420 (1970)
	var prec= m$.epsilon; //Math.pow(2,-52) // assumes double prec
	var tolerance= 1.e-64/prec;
	var itmax= 50;
	var c=0;
	var i=0;
	var j=0;
	var k=0;
	var l=0;
	
	var u= this.clone();
	var m= u.rows;
	
	var n= u.cols;
	
	if (m < n) throw "Need more rows than columns"
	
	var e = new Array(n);
	var q = new Array(n);
	for (i=0; i<n; i++) e[i] = q[i] = 0.0;
	var v = jMath.zeros(n,n);
//	v.zero();
	
 	function pythag(a,b)
 	{
		a = Math.abs(a)
		b = Math.abs(b)
		if (a > b)
			return a*Math.sqrt(1.0+(b*b/a/a))
		else if (b == 0.0) 
			return a
		return b*Math.sqrt(1.0+(a*a/b/b))
	}

	//Householder's reduction to bidiagonal form

	var f= 0.0;
	var g= 0.0;
	var h= 0.0;
	var x= 0.0;
	var y= 0.0;
	var z= 0.0;
	var s= 0.0;
	
	for (i=0; i < n; i++)
	{	
		e[i]= g;
		s= 0.0;
		l= i+1;
		for (j=i; j < m; j++) 
			s += (u[j][i]*u[j][i]);
		if (s <= tolerance)
			g= 0.0;
		else
		{	
			f= u[i][i];
			g= Math.sqrt(s);
			if (f >= 0.0) g= -g;
			h= f*g-s
			u[i][i]=f-g;
			for (j=l; j < n; j++)
			{
				s= 0.0
				for (k=i; k < m; k++) 
					s += u[k][i]*u[k][j]
				f= s/h
				for (k=i; k < m; k++) 
					u[k][j]+=f*u[k][i]
			}
		}
		q[i]= g
		s= 0.0
		for (j=l; j < n; j++) 
			s= s + u[i][j]*u[i][j]
		if (s <= tolerance)
			g= 0.0
		else
		{	
			f= u[i][i+1]
			g= Math.sqrt(s)
			if (f >= 0.0) g= -g
			h= f*g - s
			u[i][i+1] = f-g;
			for (j=l; j < n; j++) e[j]= u[i][j]/h
			for (j=l; j < m; j++)
			{	
				s=0.0
				for (k=l; k < n; k++) 
					s += (u[j][k]*u[i][k])
				for (k=l; k < n; k++) 
					u[j][k]+=s*e[k]
			}	
		}
		y= Math.abs(q[i])+Math.abs(e[i])
		if (y>x) 
			x=y
	}
	
	// accumulation of right hand gtransformations
	for (i=n-1; i != -1; i+= -1)
	{	
		if (g != 0.0)
		{
		 	h= g*u[i][i+1]
			for (j=l; j < n; j++) 
				v[j][i]=u[i][j]/h
			for (j=l; j < n; j++)
			{	
				s=0.0
				for (k=l; k < n; k++) 
					s += u[i][k]*v[k][j]
				for (k=l; k < n; k++) 
					v[k][j]+=(s*v[k][i])
			}	
		}
		for (j=l; j < n; j++)
		{
			v[i][j] = 0;
			v[j][i] = 0;
		}
		v[i][i] = 1;
		g= e[i]
		l= i
	}
	
	// accumulation of left hand transformations
	for (i=n-1; i != -1; i+= -1)
	{	
		l= i+1
		g= q[i]
		for (j=l; j < n; j++) 
			u[i][j] = 0;
		if (g != 0.0)
		{
			h= u[i][i]*g
			for (j=l; j < n; j++)
			{
				s=0.0
				for (k=l; k < m; k++) s += u[k][i]*u[k][j];
				f= s/h
				for (k=i; k < m; k++) u[k][j]+=f*u[k][i];
			}
			for (j=i; j < m; j++) u[j][i] = u[j][i]/g;
		}
		else
			for (j=i; j < m; j++) u[j][i] = 0;
		u[i][i] += 1;
	}
	
	// diagonalization of the bidiagonal form
	prec= prec*x
	for (k=n-1; k != -1; k+= -1)
	{
		for (var iteration=0; iteration < itmax; iteration++)
		{	// test f splitting
			var test_convergence = false
			for (l=k; l != -1; l+= -1)
			{	
				if (Math.abs(e[l]) <= prec)
				{	test_convergence= true
					break 
				}
				if (Math.abs(q[l-1]) <= prec)
					break 
			}
			if (!test_convergence)
			{	// cancellation of e[l] if l>0
				c= 0.0
				s= 1.0
				var l1= l-1
				for (i =l; i<k+1; i++)
				{	
					f= s*e[i]
					e[i]= c*e[i]
					if (Math.abs(f) <= prec)
						break
					g= q[i]
					h= pythag(f,g)
					q[i]= h
					c= g/h
					s= -f/h
					for (j=0; j < m; j++)
					{	
						y= u[j][l1]
						z= u[j][i]
						u[j][l1] =  y*c+(z*s)
						u[j][i] = -y*s+(z*c)
					} 
				}	
			}
			// test f convergence
			z= q[k]
			if (l== k)
			{	//convergence
				if (z<0.0)
				{	//q[k] is made non-negative
					q[k]= -z
					for (j=0; j < n; j++)
						v[j][k] = -v[j][k]
				}
				break  //break out of iteration loop and move on to next k value
			}
			if (iteration >= itmax-1)
				throw 'Error: no convergence.'
			// shift from bottom 2x2 minor
			x= q[l]
			y= q[k-1]
			g= e[k-1]
			h= e[k]
			f= ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
			g= pythag(f,1.0)
			if (f < 0.0)
				f= ((x-z)*(x+z)+h*(y/(f-g)-h))/x
			else
				f= ((x-z)*(x+z)+h*(y/(f+g)-h))/x
			// next QR transformation
			c= 1.0
			s= 1.0
			for (i=l+1; i< k+1; i++)
			{	
				g= e[i]
				y= q[i]
				h= s*g
				g= c*g
				z= pythag(f,h)
				e[i-1]= z
				c= f/z
				s= h/z
				f= x*c+g*s
				g= -x*s+g*c
				h= y*s
				y= y*c
				for (j=0; j < n; j++)
				{	
					x= v[j][i-1]
					z= v[j][i]
					v[j][i-1] = x*c+z*s
					v[j][i] = -x*s+z*c
				}
				z= pythag(f,h)
				q[i-1]= z
				c= f/z
				s= h/z
				f= c*g+s*y
				x= -s*g+c*y
				for (j=0; j < m; j++)
				{
					y= u[j][i-1]
					z= u[j][i]
					u[j][i-1] = y*c+z*s
					u[j][i] = -y*s+z*c
				}
			}
			e[l]= 0.0
			e[k]= f
			q[k]= x
		} 
	}
		
	//vt= transpose(v)
	//return (u,q,vt)
	for (i=0;i<q.length; i++) 
	  if (q[i] < prec) q[i] = 0
	  
	//sort eigenvalues	
	for (i=0; i< n; i++)
	{	 
	//writeln(q)
	 for (j=i-1; j >= 0; j--)
	 {
	  if (q[j] < q[i])
	  {
	//  writeln(i,'-',j)
	   c = q[j]
	   q[j] = q[i]
	   q[i] = c
	   for(k=0;k<u.length;k++) { temp = u[k][i]; u[k][i] = u[k][j]; u[k][j] = temp; }
	   for(k=0;k<v.length;k++) { temp = v[k][i]; v[k][i] = v[k][j]; v[k][j] = temp; }
//	   u.swapCols(i,j)
//	   v.swapCols(i,j)
	   i = j	   
	  }
	 }	
	}

	var qm = jMath.zeros( u.cols, v.rows );
	for ( var i = 0 ; i < q.length; i++ )
	{
		qm[i][i] = q[i];
	}
	
	return {U:u,S:qm,V:v}
}


m$.fn["'"] = m$.fn.transpose = function()
{
	var list = [];
	for ( var c = 0 ; c < this.cols ; c++ )
	{
		list.push([]);
		for ( var r = 0 ; r < this.rows; r++ )
		{
			list[c][r] = this[r][c];
		}
	}
	return m$(list);
}

// append column
m$.fn[':='] = m$.fn.appendToColumn = function()
{
	if ( this.length = 0 )
	{
		this[0] = [[]];
	}
	for ( var i = 0 ; i < arguments.length; i++ )
	{
		var mat = arguments[i];
		this.rows = Math.max( this.rows, mat.rows );
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			if ( !this[r] )
			{
				this[r] = [];
				for ( var c = 0 ; c < this.cols ; c++ )
				{
					this[r].push(undefined);
				}
			}
			if ( !mat[r] )
			{
				mat[r] = [];
				for ( var c = 0 ; c < this.cols ; c++ )
				{
					mat[r].push(undefined);
				}
			}
			this[r] = this[r].concat(mat[r]);
		}
		this.cols += mat.cols;
	}
	this.length = this.rows;
	return this;
}

// append row
m$.fn[';='] = m$.fn.appendToRow = function()
{
	if ( this.length == 0 )
	{
		this[0] = [[]];
	}
	for ( var i = 0 ; i < arguments.length; i++ )
	{
		var mat = arguments[i];
		this.cols = Math.max( this.cols , mat.cols );
		for ( var r = 0 ; r < mat.rows ; r++ )
		{
			var rIdx = r + this.rows;
			this[rIdx] = [];
			for ( var c = 0 ; c < mat.cols ; c++ )
			{
				this[rIdx][c] = mat[r][c];
			}
		}
		this.rows += mat.rows;
	}
	this.length = this.rows;
	return this;
}

m$.joinByRow = function()
{
	var first;
	var list;
	if ( arguments.length == 1 )
	{
		first = arguments[0][0].clone();
		list = arguments[0].slice(1);
	}
	else
	{
		first = arguments[0].clone();
		list = Array.prototype.slice.call( arguments, 1 );
	}
	return first.appendToRow.apply( first, list );
}

m$.joinByColumn = function()
{
	var first;
	var list;
	if ( arguments.length == 1 )
	{
		first = arguments[0][0].clone();
		list = arguments[0].slice(1);
	}
	else
	{
		first = arguments[0].clone();
		list = Array.prototype.slice.call( arguments, 1 );
	}
	return first.appendToColumn.apply( first, list );
}

m$.fn.clone = function()
{
	var list = [];
	for ( var i = 0 ; i < this.rows ; i++ )
	{
		var sub = [];
		for( var j = 0 ; j < this.cols ; j++ )
		{
			sub.push(this[i][j]);
		}
		list.push(sub);
	}
	return jMath( list );
}

m$.fn.repeatColumn = function(len)
{
	var list = this.toArray();
	for( var i = 0 ; i < len ; i++ )
	{
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			for ( var r = 0 ; r < this.rows; r++ )
			{
				list[r].push( list[r][c] );
			}
		}
	}
	return jMath(list);
}

m$.fn.repeatRow = function(len)
{
	var list = this.toArray();
	for( var i = 0 ; i < len ; i++ )
	{
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			list.push( Array.prototype.slice.call( list[r], 0) );
		}
	}
	return jMath(list);
}

// 삼각함수


Math.atanh = function(v)
{
	return 0.5 * Math.log((1+v)/(1-v));
}

Math.tanh = function(v)
{
	var e = Math.exp(2*v);
	return (e-1)/(e+1);
}

m$.fn.sin = function() { return this.map( Math.sin ); }
m$.fn.cos = function() { return this.map( Math.cos ); }
m$.fn.tan = function() { return this.map( Math.tan ); }
m$.fn.acos = function() { return this.map( Math.acos ); }
m$.fn.asin = function() { return this.map( Math.asin ); }
m$.fn.atan = function() { return this.map( Math.atan ); }
m$.fn.atanh = function() { return this.map( Math.atanh ); }
m$.fn.tanh = function() { return this.map( Math.tanh ); }

// real to int
m$.fn.ceil = function() { return this.map( Math.ceil ); }
m$.fn.floor = function() { return this.map( Math.floor ); }
m$.fn.abs = function() { return this.map( Math.abs ); }

// power
m$.fn.sqrt = function() { return this.map( Math.sqrt ); }
m$.fn.log = function() { return this.map( Math.log ); }
m$.fn.exp = function() { return this.map( Math.exp ); }
m$.fn.pow = function(y) { return this.map( function(x){
 	return Math.pow(x,y);
}); }


m$.fn.sum = function(dir)
{
	var ret;
	if ( dir != 3 )
	{
		if ( this.rows == 1 ) dir = 2;
		else if( this.cols == 1 ) dir = 1;
		else dir = dir || 1;	
	}

	if ( dir == 1 ) 	// column
	{ 
		ret = jMath.zeros( 1, this.cols );
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			for ( var r = 0 ; r < this.rows ; r++ )
			{
				ret[0][c] += this[r][c];
			}
		}
	}
	else if ( dir == 2 ) // row
	{
		ret = jMath.zeros( this.rows, 1);
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			for ( var c = 0 ; c < this.cols ; c++ )
			{
				ret[r][0] += this[r][c];
			}
		}
	}
	else if ( dir == 3 ) // all
	{
		ret = 0;
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			for ( var c = 0 ; c < this.cols ; c++ )
			{
				ret += this[r][c];
			}
		}
	}
	return ret;
}

m$.fn.__dir = function(dir)
{
	if ( this.rows == 1 ) dir = 2;
	else if( this.cols == 1 ) dir = 1;
	else dir = dir || 1;	
	return dir;
}

m$.fn.prod = function(dir)
{
	dir = this.__dir(dir);
	var ret;
	if ( dir == 1 ) 	// column
	{ 
		ret = jMath.ones( 1, this.cols );
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			for ( var r = 0 ; r < this.rows ; r++ )
			{
				ret[0][c] *= this[r][c];
			}
		}
	}
	else if ( dir == 2 ) // row
	{
		ret = jMath.ones( this.rows, 1);
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			for ( var c = 0 ; c < this.cols ; c++ )
			{
				ret[r][0] *= this[r][c];
			}
		}
	}
	else if ( dir == 3 )
	{
		ret = 1;
		for ( var r = 0; r < this.rows;  r++)
		{
			for( var c = 0 ; c < this.cols ; c++ )
			{
				ret *= this[r][c];
			}
		}
	}
	return ret;
}

m$.fn.diff = function(n,dir)
{
	n = n || 1;
	dir = this.__dir(dir);
	var ret;
	if ( dir == 1 ) 	// column
	{ 
		ret = jMath.zeros( this.rows-n, this.cols );
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			for ( var r = 0 ; r < ret.rows ; r++ )
			{
				ret[r][c] = this[r][c] - this[r+n][c];
			}
		}
	}
	else if ( dir == 2 ) // row
	{
		ret = jMath.zeros( this.rows, this.cols-n);
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			for ( var c = 0 ; c < ret.cols ; c++ )
			{
				ret[r][c] = this[r][c] - this[r][c+n];
			}
		}
	}
	return ret;
}

m$.fn.trace = function()
{
	var m = Math.min( this.rows, this.cols );
	var s = 0;
	for ( var i = 0 ;  i < m ; i++ )
	{
		s += this[i][i];
	}
	return s;
}

m$.fn.mean = function(dir)
{
	var len;
	if ( dir != 3 )
	{
		dir = this.__dir(dir);
		len = dir == 1 ? this.rows : this.cols;
		return this.sum(dir)['./'](len);
	}
	else
	{
		return this.sum(dir)/(this.rows*this.cols);
	}
}

m$.fn.geomean = function(dir)
{
	var len;
	if ( dir == 3 )
	{
		var mul = 1;
		for ( var r=0 ; r < this.rows ; r++ )
		{
			for( var c=0; c < this.cols; c++ )
			{
				mul *= this[r][c];
			}
		}
		return Math.pow( mul, 1/(this.rows*this.cols));
	}
	else
	{
		dir = this.__dir(dir);
		var ret;
		var size = 0;
		if ( dir == 1 )
		{
			size = this.rows;
			ret = jMath.ones( 1, this.cols );
			for ( var c = 0 ; c < this.cols ; c++ )
			{
				for ( var r=0 ; r < this.rows ; r++ )
				{
					ret[0][c] *= this[r][c];
				}
			}
		}
		else
		{
			size = this.cols;
			ret = jMath.ones( this.rows, 1);
			for ( var r = 0 ; r < this.rows ; r++ )
			{
				for ( var c=0 ; c < this.cols; c++ )
				{
					ret[r][0] *= this[r][c];
				}
			}
		}
		size = 1/size;
		return ret.map(function(v){
			return Math.pow(v,size);
		});
	}
}

m$.fn.var = function(pop,dir)
{
	if ( dir != 3 )
	{
		dir = this.__dir(dir);
	}
	var m = this.mean(dir);
	var ret;
	var n;
	if ( dir == 1 )
	{
		ret = this.map( function(v,r,c){
			return Math.pow(v - m[0][c],2);
		}).sum(dir);
		n = this.rows;
	}
	else if ( dir == 2 )
	{
		ret = this.map( function(v,r,c){
			return Math.pow(v - m[r][0],2);
		}).sum(dir);
		n = this.cols;
	}
	else if ( dir == 3 )
	{
		ret = this.map( function(v,r,c){
			return Math.pow(v - m,2);
		}).sum(3);
		n = this.rows * this.cols;
		return ret/(pop ? n : n - 1);
	}
	return ret['./'](pop ? n : n - 1);
}

m$.fn.std = function(pop,dir)
{
	if ( dir == 3 )
	{
		return Math.sqrt(this.var(pop,dir));
	}
	return this.var(pop,dir).sqrt();
}

m$.fn.sort = function(dir)
{
	dir = this.__dir(dir);
	if ( dir == 1 ) 
	{
		ret = this.transpose();	
	}
	else if ( dir == 2 )
	{
		ret = this;
	}
	for ( var i = 0 ; i < ret.length; i++ )
	{
		ret[i] = ret[i].sort( function(arg1,arg2){ return arg1 - arg2; });	
	}
	return dir == 1 ? ret.transpose() : ret;
}

m$.fn.median = function(dir)
{
	if ( dir != 3 )
	{
		dir = this.__dir(dir);
	}
	if ( dir == 1 )
	{
		return this.transpose().median(2).transpose();
	}
	else if ( dir == 2 )
	{
		var sret = this.sort(dir);		
		var ret;
   		ret = jMath.zeros( this.rows, 1);
		var p = Math.floor(this.cols/2);
		if ( this.cols % 2 == 0 ) 
		{
			for( var r = 0; r < this.rows ; r++ )
			{
				ret[r][0] = (sret[r][p] + sret[r][p-1])/2;
			}
		}
		else
		{
			for( var r = 0; r < sret.rows ; r++ )
			{
				ret[r][0] = sret[r][p];
			}
		}
	}
	else if ( dir == 3 )
	{
		var list = this[0];
		for ( var i = 1 ; i < this.rows ; i++ )
		{
			list = list.concat( this[i] );
		}
		return jMath(list).median(2)[0][0];
	}
	return ret;
}

m$.fn.max = function(dir)
{
	if ( dir != 3 )
	{
		dir = this.__dir(dir);
	}
	var ret;
	if ( dir == 1 )
	{
		ret = jMath.zeros(2, this.cols);
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			ret[0][c] = this[0][c];			
			ret[1][c] = 0;
		}
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			for ( var r = 1 ; r < this.rows ; r++ )
			{
				if ( ret[0][c] < this[r][c] )
				{
					ret[0][c] = this[r][c];
					ret[1][c] = r;
				}
			}
		}
	}
	else if ( dir == 2 )
	{
		ret = jMath.zeros(this.rows, 1);
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			ret[r][0] = this[r][0];			
			ret[r][1] = 0;
		}
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			for ( var c = 1 ; c < this.cols ; c++ )
			{
				if ( ret[r][0] < this[r][c] )
				{
					ret[r][0] = this[r][c];
					ret[r][1] = c;
				}
			}
		}
	}
	else if ( dir == 3 )
	{
		ret = this[0][0];
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			for ( var c = 0 ; c < this.cols ; c++ )
			{
				if ( ret < this[r][c] )
				{
					ret = this[r][c];
				}
			}
		}
	}
	return ret;
}

m$.fn.min = function(dir)
{
	if ( dir != 3 )
	{
		dir = this.__dir(dir);
	}
	var ret;
	if ( dir == 1 )
	{
		ret = jMath.zeros(2, this.cols);
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			ret[0][c] = this[0][c];			
			ret[1][c] = 0;
		}
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			for ( var r = 1 ; r < this.rows ; r++ )
			{
				if ( ret[0][c] > this[r][c] )
				{
					ret[0][c] = this[r][c];
					ret[1][c] = r;
				}
			}
		}
	}
	else if ( dir == 2 )
	{
		ret = jMath.zeros(this.rows, 1);
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			ret[r][0] = this[r][0];			
			ret[r][1] = 0;
		}
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			for ( var c = 1 ; c < this.cols ; c++ )
			{
				if ( ret[r][0] > this[r][c] )
				{
					ret[r][0] = this[r][c];
					ret[r][1] = c;
				}
			}
		}
	}
	else if ( dir == 3 )
	{
		ret = this[0][0];
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			for ( var c = 0 ; c < this.cols ; c++ )
			{
				if ( ret > this[r][c] )
				{
					ret = this[r][c];
				}
			}
		}
	}
	return ret;
}

m$.fn.range = function(dir)
{
	if ( dir != 3 )
	{
		dir = this.__dir(dir);
	}
	var max = this.max(dir);
	var min = this.min(dir);

	if ( dir == 1 )
	{
		return m$( [ min[0], max[0] ] );
	}
	else if ( dir == 2 )
	{
		var list = [];
		for ( var i = 0 ; i < this.rows ; i++ )
		{
			list.push( [  min[i][0], max[i][0] ] );
		}
		return m$(list);
	}
	else if ( dir == 3 )
	{
		return [ min, max ];	
	}
}


m$.fn.freqdist = function(idxfunc)
{
	var list = {};
	var args = Array.prototype.slice.call( arguments, 1);
	args.unshift(0);
	this.forEach( function(v,r,c){
		if (idxfunc)
		{
			if ( Array.isArray(idxfunc) )
			{
				var i;
				for ( i = 1 ; i < idxfunc.length; i++ )
				{
					if ( v < idxfunc[i] )
					{
						v = idxfunc[i-1];
						break;
					}
				}		
				if ( i == idxfunc.length )
				{
					v = idxfunc[idxfunc.length-1];
				}
			}
			else
			{
				args[0] = v;
				v = idxfunc.apply(this,args);
			}
		}
		if ( !list[v] ) 
		{
			list[v] = 1;
		}
		else
		{
			list[v]++;
		}
	}, this);
	var keys = Object.keys(list).sort(function(a1,a2){ 
		if ( parseFloat(a1) != NaN )
		{
			a1 = parseFloat(a1);
			a2 = parseFloat(a2);
		}
		if ( a1 == a2 ) return 0;
		else if ( a1 > a2 ) return 1;
		return -1;
	});
	var dist = [];
	for ( var i = 0 ; i < keys.length; i++ )
	{
		var k = keys[i];
		dist.push( [ k, list[k] ] );
	}
	return jMath(dist);
}

m$.fn.relfreqdist = function()
{
	var list = [];
	var sum = 0;
	for ( var i = 0 ; i < this.rows ; i++ )
	{
		sum += this[i][1];
	}
	for ( var i = 0 ; i < this.rows ; i++ )
	{
		list.push( [this[i][0], this[i][1]/sum] );
	}		
	return jMath(list);
}

m$.fn.mode = function()
{
	var fdist = this.freqdist();
	var max = fdist.max();
	var mlist = [];
	for ( var i = 0 ; i < fdist.rows; i++ )
	{
		if ( fdist[i][1] == max[0][1] )
		{
			mlist.push( fdist[i][0] );
		}
	}
	return mlist.length == fdist.rows ? null: m$(mlist);
}

m$.fn.wsum = function()
{
	var s = 0;
	var freq = 0;
	for ( var r = 0 ; r < this.rows ; r++ )
	{
		s +=  this[r][0] * this[r][1];
		freq += this[r][1];
	}	
	return m$([s, freq]);
}

m$.fn.wmean = function()
{
	var avg = this.wsum();
	avg[0][0] /= avg[0][1];
	return avg;
}

m$.fn.wvar = function(pop)
{
	var s = 0;
	var avg = this.wmean();
	for ( var r = 0 ;  r < this.rows ; r++ )
	{
		s += Math.pow( this[r][0] - avg[0][0], 2 ) * this[r][1];
	}
	return m$( [s /(pop ? avg[0][1] : (avg[0][1]-1)) , avg[0][1]] );
}

m$.fn.wstd = function(pop)
{
	var v = this.wvar(pop);
	v[0][0] = Math.sqrt(v[0][0]);
	return m$(v);
}

m$.fn.dpdmean = function()
{
	var m = 0;
	for ( var i = 0 ; i < this.rows ; i++ )
	{
		m += this[i][0] * this[i][1];
	}	
	return m;
}

m$.fn.dpdvar = function()
{
	var m = this.dpdmean();	
	var v = 0;
	for ( var i = 0 ;i  < this.rows ; i++ )
	{
		v += Math.pow( this[i][0] - m, 2 ) * this[i][1];
	}
	return v;
}

m$.fn.dpdstd = function()
{
	return Math.sqrt(this.dpdvar());
}

m$.fn.percentile = function(ps,dir)
{
	dir = this.__dir(dir);					
	var sobj = this.sort(dir);
	var list = [];
	if ( typeof ps == 'number' )
	{
		ps = [ps];
	}
	if ( dir == 1 )
	{
		for ( var k = 0 ; k < ps.length; k++ )
		{
			var pos = ( this.rows - 1 ) * ps[k];
			var i = Math.floor(pos);
			var real;
			switch( ps[k] )
			{
				case 0: 
					real = 0;
					break;
				case 1:
					i = this.rows - 2;
					real = 1;
					break;
				default:
					real = pos - i;
					break;
			}
			list[k] = [];
			for ( var c = 0 ; c < this.cols ; c++ )
			{
				list[k][c] = sobj[i][c] + real *  ( sobj[i+1][c] - sobj[i][c] );
			}
		}
	}
	else if ( dir == 2 )
	{
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			list[r] = [];
		}
		for ( var k = 0 ; k < ps.length; k++ )
		{
			var pos = ( this.cols - 1 ) * ps[k];
			var i = Math.floor(pos);
			var real;
			switch( ps[k] )
			{
				case 0: 
					real = 0;
					break;
				case 1:
					i = this.cols- 2;
					real = 1;
					break;
				default:
					real = pos - i;
					break;
			}
			for ( var r = 0 ; r < this.rows ; r++ )
			{
				var v = sobj[r][i] + real * ( sobj[r][i+1] - sobj[r][i] );
				list[r].push(v);
			}
		}
	}
	return m$(list);
}

m$.fn.percentrank = function(v,dir)
{
	dir = this.__dir(dir);					
	var sobj = this.sort(dir);
	var list = [];
	if ( dir == 1 )
	{
		for ( var c = 0 ; c < this.cols ; c++ )
		{
			for ( var r = 0 ; r < this.rows; r++)
			{
				if(  sobj[r][c] > v )				
				{
					r--;
					list[c] = r;
					break;
				}
			}
			if( r == -1 || r == this.rows )
			{
				list[c] = NaN;
			}
			else
			{
				list[c] = r + (v-sobj[r][c])/( sobj[r+1][c] - sobj[r][c]);
				list[c] /= (sobj.rows-1);
			}
		}
	}
	else if ( dir == 2 )
	{
		for ( var r = 0 ; r < this.rows ; r++ )
		{
			for ( var c = 0 ; c < this.cols; c++)
			{
				if( sobj[r][c] > v )				
				{
					c--;
					list[r] = c;
					break;
				}
			}
			if( c == -1 || c == this.cols )
			{
				list[r] = [NaN];
			}
			else
			{
				list[r] = c + (v-sobj[r][c])/( sobj[r][c+1] - sobj[r][c]);
				list[r] = [ list[r] / (this.cols-1) ];
			}
		}
	}
	return m$(list);
}

m$.fn.quartiles = function(dir)
{
	return this.percentile([0,0.25,0.5,0.75,1],dir);
}

m$.zscore = function(v,mu,sigma)
{
	return (v - mu)/sigma;
}

m$.fn.zscore = function(v,dir)
{
	if ( arguments.length <= 1 )
	{
		dir = v;
		if ( dir == 3 )
		{
			var mu = this.mean(dir);
			var sigma = this.std(false,dir);
			return this.map(function(v) {
				return (v - mu)/sigma;
			} )[0];
		}		
		else
		{
			dir = this.__dir(dir);	
			var mu = this.mean(dir);
			var sigma = this.std(false,dir);
			if ( dir == 1 )
			{
				return this.map(function(v,r,c){
					return (v - mu[0][c])/sigma[0][c];
				});
			}
			else
			{
				return this.map(function(v,r,c){
					return (v - mu[r][0])/sigma[r][0];
				});
			}
		}
	}
	else
	{
		if ( dir == 3 )
		{
			var mu = this.mean(dir);
			var sigma = this.std(false,dir);
			return (v - mu)/sigma;
		}
		else
		{
			dir = this.__dir(dir);	
			var mu = this.mean(dir).map(function(m){
				return v - m;
			});
			var sigma = this.std(false,dir);
			return mu['./'](sigma);
		}
	}
}

jMath.factorial = function(x)
{
	return jMath.gamma(x+1);
}

jMath.nchoosek = function(n,x)
{
	return jMath.gamma(n+1)/(jMath.gamma(n-x+1) * jMath.gamma(x+1));
}

// error function

jMath.erf = function(x)
{
	var a1 =  0.254829592;
	var a2 = -0.284496736;
	var a3 =  1.421413741;
	var a4 = -1.453152027;
	var a5 =  1.061405429;
	var p  =  0.3275911;

	var sign = 1;
	if (x < 0) 
	{
		sign = -1;
	}
	x = Math.abs(x);

	var t = 1.0/(1.0 + p*x);
	var y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.exp(-x*x);

	return sign*y;
}

jMath.erfinv = function(x)
{
	var z;
	var a  = 0.147;                                                   
	var the_sign_of_x;
	if(x==0) 
	{
		the_sign_of_x = 0;
	} 
	else if(x>0)
	{
		the_sign_of_x = 1;
	}
	else 
	{
		the_sign_of_x = -1;
	}

	if(0 != x) 
	{
		var ln_1minus_x_sqrd = Math.log(1-x*x);
		var ln_1minusxx_by_a = ln_1minus_x_sqrd / a;
		var ln_1minusxx_by_2 = ln_1minus_x_sqrd / 2;
		var ln_etc_by2_plus2 = ln_1minusxx_by_2 + (2/(Math.PI * a));
		var first_sqrt = Math.sqrt((ln_etc_by2_plus2*ln_etc_by2_plus2)-ln_1minusxx_by_a);
		var second_sqrt = Math.sqrt(first_sqrt - ln_etc_by2_plus2);
		z = second_sqrt * the_sign_of_x;
	} 
	else 
	{ // x is zero
		z = 0;
	}
	return z;
 }

// gamma function
jMath.gamma = function(x) {
	var p = [
			-1.716185138865495, 24.76565080557592, -379.80425647094563,
			629.3311553128184, 866.9662027904133, -31451.272968848367,
			-36144.413418691176, 66456.14382024054
		],
		q = [
			-30.8402300119739, 315.35062697960416, -1015.1563674902192,
			-3107.771671572311, 22538.118420980151, 4755.8462775278811,
			-134659.9598649693, -115132.2596755535
		],
		fact = false,
		n = 0,
		xden = 0,
		xnum = 0,
		y = x,
		i, z, yi, res, sum, ysq;
	if(y <= 0) {
		res = y % 1 + 3.6e-16;
		if (res) {
			fact = (!(y & 1) ? 1 : -1) * Math.PI / Math.sin(Math.PI * res);
			y = 1 - y;
		} else {
			return Infinity;
		}
	}
	yi = y;
	if (y < 1) {
		z = y++;
	} else {
		z = (y -= n = (y | 0) - 1) - 1;
	}
	for (i = 0; i < 8; ++i) {
		xnum = (xnum + p[i]) * z;
		xden = xden * z + q[i];
	}
	res = xnum / xden + 1;
	if (yi < y) {
		res /= yi;
	} else if (yi > y) {
		for (i = 0; i < n; ++i) {
			res *= y;
			y++;
		}
	}
	if (fact) {
		res = fact / res;
	}
	return res;
}

jMath.gammaln = function(x) {
	var j = 0,
	cof = [ 76.18009172947146, -86.50532032941677, 24.01409824083091,
		-1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 ],
	ser = 1.000000000190015,
	xx, y, tmp;
	tmp = (y = xx = x) + 5.5;
	tmp -= (xx + 0.5) * Math.log(tmp);
	for (; j < 6; j++) ser += cof[j] / ++y;
	return Math.log(2.5066282746310005 * ser / xx) - tmp;
}


// lower incomplete gamma function P(a,x)
jMath.gammainc = function(x,a)
{
	var aln = jMath.gammaln(a);
	var ap = a, sum = 1 / a, del = sum, b = x + 1 - a, c = 1 / 1.0e-30;
	var d = 1 / b, h = d, i = 1;
	// calculate maximum number of itterations required for a
	var ITMAX = -~(Math.log((a >= 1) ? a : 1 / a) * 8.5 + a * 0.4 + 17);
	var an, endval;

	if (x < 0 || a <= 0) {
		return NaN;
	} else if (x < a + 1) {
		for (; i <= ITMAX; i++) {
			sum += del *= x / ++ap;
		}
		return sum * Math.exp(-x + a * Math.log(x) - (aln));
	}
	for (; i <= ITMAX; i++) {
		an = -i * (i - a);
		b += 2;
		d = an * d + b;
		c = b + an / c;
		d = 1 / d;
		h *= d * c;
	}
	return 1 - h * Math.exp(-x + a * Math.log(x) - (aln));
}

jMath.gammaincinv = function(p, a) 
{
	var j = 0, a1 = a - 1, EPS = 1e-8, gln = jMath.gammaln(a);
	var x, err, t, u, pp, lna1, afac;
	if(p >= 1) return Math.max(100, a + 100 * Math.sqrt(a));
	if(p <= 0) return 0;
	if(a > 1) {
		lna1 = Math.log(a1);
		afac = Math.exp(a1 * (lna1 - 1) - gln);
		pp = (p < 0.5) ? p : 1 - p;
		t = Math.sqrt(-2 * Math.log(pp));
		x = (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t;
		if(p < 0.5) x = -x;
		x = Math.max(1e-3, a * Math.pow(1 - 1 / (9 * a) - x / (3 * Math.sqrt(a)), 3));
	} 
	else 
	{
		t = 1 - a * (0.253 + a * 0.12);
		if(p < t) x = Math.pow(p / t, 1 / a);
		else x = 1 - Math.log(1 - (p - t) / (1 - t));
	}
	for(; j < 12; j++) 
	{
		if(x <= 0) return 0;
		err = jMath.gammainc(x,a) - p;
		if(a > 1) t = afac * Math.exp(-(x - a1) + a1 * (Math.log(x) - lna1));
		else t = Math.exp(-x + a1 * Math.log(x) - gln);
		u = err / t;
		x -= (t = u / (1 - 0.5 * Math.min(1, u * ((a - 1) / x - 1))));
		if(x <= 0) x = 0.5 * (x + t);
		if(Math.abs(t) < EPS * x) break;
	}
	return x;
}

// beta function
jMath.beta = function(x, y) {
	// ensure arguments are positive
	if (x <= 0 || y <= 0) return undefined;
	// make sure x + y doesn't exceed the upper limit of usable values
	return (x + y > 170) ?  Math.exp(jMath.betaln(x, y)) : jMath.gamma(x) * jMath.gamma(y) / jMath.gamma(x + y);
}
	
// natural logarithm of beta function
jMath.betaln = function(x, y) {
	return jMath.gammaln(x) + jMath.gammaln(y) - jMath.gammaln(x + y);
}

jMath.betacdf = function(x,a,b) 
{
	var fpmin = 1e-30,
		m = 1,
		m2, aa, c, d, del, h, qab, qam, qap;
	// These q's will be used in factors that occur in the coefficients
	qab = a + b;
	qap = a + 1;
	qam = a - 1;
	c = 1;
	d = 1 - qab * x / qap;
	if(Math.abs(d) < fpmin) d = fpmin;
	d = 1 / d;
	h = d;
	for (; m <= 100; m++) {
		m2 = 2 * m;
		aa = m * (b - m) * x / ((qam + m2) * (a + m2));
		// One step (the even one) of the recurrence
		d = 1 + aa * d;
		if(Math.abs(d) < fpmin) d = fpmin;
		c = 1 + aa / c;
		if(Math.abs(c) < fpmin) c = fpmin;
		d = 1 / d;
		h *= d * c;
		aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
		// Next step of the recurrence (the odd one)
		d = 1 + aa * d;
		if(Math.abs(d) < fpmin) d = fpmin;
		c = 1 + aa / c;
		if(Math.abs(c) < fpmin) c = fpmin;
		d = 1 / d;
		del = d * c;
		h *= del;
		if(Math.abs(del - 1.0) < 3e-7) break;
	}
	return h;
}

jMath.betainc = function(x,a,b)
{
	// Factors in front of the continued fraction.
	var bt = (x === 0 || x === 1) ?  0 :
		Math.exp(jMath.gammaln(a + b) - jMath.gammaln(a) -
		jMath.gammaln(b) + a * Math.log(x) + b *
		Math.log(1 - x));
	if(x < 0 || x > 1) return false;
	if(x < (a + 1) / (a + b + 2))
	{
		// Use continued fraction directly.
		return bt * jMath.betacdf(x, a, b) / a;
	}
	// else use continued fraction after making the symmetry transformation.
	return 1 - bt * jMath.betacdf(1 - x, b, a) / b;
}

jMath.betainv = function(p,a,b)
{
	var EPS = 1e-8,
		a1 = a - 1,
		b1 = b - 1,
		j = 0,
		lna, lnb, pp, t, u, err, x, al, h, w, afac;
	if(p <= 0) return 0;
	if(p >= 1) return 1;
	if(a >= 1 && b >= 1) {
		pp = (p < 0.5) ? p : 1 - p;
		t = Math.sqrt(-2 * Math.log(pp));
		x = (2.30753 + t * 0.27061) / (1 + t* (0.99229 + t * 0.04481)) - t;
		if(p < 0.5) x = -x;
		al = (x * x - 3) / 6;
		h = 2 / (1 / (2 * a - 1)  + 1 / (2 * b - 1));
		w = (x * Math.sqrt(al + h) / h) - (1 / (2 * b - 1) - 1 / (2 * a - 1)) * (al + 5 / 6 - 2 / (3 * h));
		x = a / (a + b * Math.exp(2 * w));
	} else {
		lna = Math.log(a / (a + b));
		lnb = Math.log(b / (a + b));
		t = Math.exp(a * lna) / a;
		u = Math.exp(b * lnb) / b;
		w = t + u;
		if(p < t / w) x = Math.pow(a * w * p, 1 / a);
		else x = 1 - Math.pow(b * w * (1 - p), 1 / b);
	}
	afac = -jMath.gammaln(a) - jMath.gammaln(b) + jMath.gammaln(a + b);
	for(; j < 10; j++) {
		if(x === 0 || x === 1) return x;
		err = jMath.betainc(x, a, b) - p;
		t = Math.exp(a1 * Math.log(x) + b1 * Math.log(1 - x) + afac);
		u = err / t;
		x -= (t = u / (1 - 0.5 * Math.min(1, u * (a1 / x - b1 / (1 - x)))));
		if(x <= 0) x = 0.5 * (x + t);
		if(x >= 1) x = 0.5 * (x + t + 1);
		if(Math.abs(t) < EPS * x && j > 0) break;
	}
	return x;
}

jMath.gauss_jordan = function(a, b) {
	var m = a.clone()[':='](b);
	var h = m.rows;
	var w = m.cols;
	// find max pivot
	for (var y = 0; y < h; y++) {
		var maxrow = y;
		for (var y2 = y+1; y2 < h; y2++) {
			if (Math.abs(m[y2][y]) > Math.abs(m[maxrow][y]))
			maxrow = y2;
		}
		var tmp = m[y];
		m[y] = m[maxrow];
		m[maxrow] = tmp
		for (var y2 = y+1; y2 < h; y2++) {
			c = m[y2][y] / m[y][y];
			for (var x = y; x < w; x++) {
				m[y2][x] -= m[y][x] * c;
			}
		}
	}
	// backsubstitute
	for (var y = h-1; y >= 0; y--) {
		c = m[y][y];
		for (var y2 = 0; y2 < y; y2++) {
			for (var x = w-1; x > y-1; x--) {
				m[y2][x] -= m[y][x] * m[y2][y] / c;
			}
		}
		m[y][y] /= c;
		for (var x = h; x < w; x++) {
			m[y][x] /= c;
		}
	}
	return m;
}

})(window);

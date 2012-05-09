
Plogist : Pattern {
	var <>r=1.5, <>xi=0.1, <>length=inf;

	*new {| r=1.5, xi=0.1, length=inf |
		^super.newCopyArgs(r, xi, length)
	}
	
	embedInStream {|inval|
		var xn = xi;
		var rStrm = r.asStream;
		var rVal;
		length.do {
			rVal = rStrm.next(inval);
			if (rVal.isNil, { ^inval });
			xn = 1.0 - (rVal * xn.squared);
			inval = xn.yield;
		};
		^inval;
	}
}

Pcml : ListPattern {
	var <>r=1.5, <>g=0.1;
	
	*new {| list, r=1.5, g=0.1, repeats=inf |
		^super.new(list, repeats).r_(r).g_(g);
	}
	
	f {|r,x| ^1.0 - (r * x.squared) }//logistic

	evolve {| prev, r, g |
		var next = Array.newClear(prev.size), halfG = g * 0.5;
		prev.size.do {|i|
			next[i] = ((1.0 - g) * this.f(r, prev[i])) 
						+ (halfG * (this.f(r, prev.wrapAt(i+1)) + this.f(r, prev.wrapAt(i-1))));
		};
		^next;
	}

	embedInStream {|inval|
		var items = list.copy;
		var rStrm = r.asStream, gStrm = g.asStream;
		var rVal, gVal;
		repeats.do {
			rVal = rStrm.next(inval);
			gVal = gStrm.next(inval);
			if (rVal.isNil || gVal.isNil, { ^inval });
			items = this.evolve(items, rVal, gVal);
			inval = items.yield;
		}
		^inval;
	}
}

Pgcm : Pcml {
	
	evolve {| prev, r, g |
		var next = Array.newClear(prev.size), global = g / prev.size * prev.sum;
		prev.do {|item,i|
			next[i] = ((1.0 - g) * this.f(r, item)) + global;
		};
		^next;
	}
}

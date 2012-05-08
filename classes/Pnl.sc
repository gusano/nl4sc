
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

	evolve {| copy, r, g |
		copy.size.do {|i|
			copy[i] = ((1.0 - g) * this.f(r, copy[i])) 
						+ (0.5 * g * (this.f(r, copy.wrapAt(i+1)) + this.f(r, copy.wrapAt(i-1))));
		};
		^copy;
	}

	embedInStream {|inval|
		var items = list.copy;
		var rStrm = r.asStream, gStrm = g.asStream;
		var rVal, gVal;
		repeats.do {
			rVal = rStrm.next(inval);
			gVal = gStrm.next(inval);
			if (rVal.isNil || gVal.isNil, { ^inval });
			inval = yield(this.evolve(items, rVal, gVal));
		}
		^inval;
	}
}

Pgcm : Pcml {
	
	evolve {| copy, r, g |
		var sum = copy.sum, reciprocal = 1 / copy.size;
		copy.do {|item,i|
			copy[i] = ((1.0 - g) * this.f(r,item)) + ((g * reciprocal) * sum);
		};
		^copy;
	}
}

//TODO 3D attractors Plorenz, Plangford, Prossler, Prabinovich?

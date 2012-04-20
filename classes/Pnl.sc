
Plogist : Pattern {
	var <>r=1.5, <>xi=0.1, <>length=inf;

	*new {| r=1.5, xi=0.1, length=inf |
		^super.newCopyArgs(r, xi, length)
	}
	storeArgs {^[r, xi, length] }
	
	embedInStream {|inval|
		var rStrm = r.asStream, rVal;
		var xn = xi;
		length.do {
			rVal = rStrm.next(inval);
			if (rVal.isNil, { ^inval });
			xn = 1.0 - (rVal * xn.squared);
			inval = xn.yield;
		};
		^inval;
	}
}

//no good
//Pcml : ListPattern {
//	var <>r=1.5, <>g=0.1, <>length=inf;
//	
//	*new {| list, r=1.4, g=0.1, length=inf |
//		^super.newCopyArgs(list, r, g, length)
//	}
//	storeArgs {^[list, r, g, length] }
//	
//	f {|i,r| ^1.0 - ( r * list[i].squared ) }
//
//	evolve {|r,g|
//		list.size.do {|i| 
//			list[i] = ((1.0 - g) * this.f(i, r)) + ((g / list.size) * list.sum);
//		};
//		^list;
//	}
//		
//	embedInStream {|inval|
//		var item;
//		var rStrm = r.asStream;
//		var gStrm = g.asStream;
//		length.do {
//			item = this.evolve(rStrm.next(inval), gStrm.next(inval));
//			inval = item.embedInStream(inval);
//		}
//		^inval;
//	}
	
//}
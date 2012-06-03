
Plogist : Pattern {
	var <>r, <>xi, <>length;

	*new {| r=1.5, xi=0.1, length=inf |
		^super.newCopyArgs(r, xi, length);
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
	var <>r, <>g;
	
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
		var next = Array.newClear(prev.size), nG = g / prev.size, sum = 0;
		prev.do {|item,i| sum = sum + this.f(r, item) };
		prev.do {|item,i|
			next[i] = ((1.0 - g) * this.f(r, item)) + (nG * sum);
		};
		^next;
	}
}

Plorenz : Pattern {
	var <>h, <>s, <>r, <>b, <>xi, <>yi, <>zi, <>length;
	//s: the fluid viscosity of a substance to its thermal conductivity
	//r: the difference in temperature between the top and bottom of the gaseous system.
	//b: the width to height ratio of the box which is being used to hold the gas

	*new {| h=0.01, s=10, r=28, b=2.667, xi=1, yi=0, zi=0, length=inf |
		^super.newCopyArgs(h, s, r, b, xi, yi, zi, length);
	}
	
	embedInStream {|inval|
		var xn = xi, yn = yi, zn = zi;
		var sStrm = s.asStream, rStrm = r.asStream, bStrm = b.asStream;
		var sVal, rVal, bVal;
		length.do {
			sVal = sStrm.next(inval);
			rVal = rStrm.next(inval);
			bVal = bStrm.next(inval);
			if (sVal.isNil || rVal.isNil || bVal.isNil, { ^inval });
			xn = xn + (h * sVal * (yn - xn));
			yn = yn + (h * (xn * (rVal - zn) - yn));
			zn = zn + (h * ((xn * yn) - (bVal * zn)));
			inval = [xn, yn, zn].yield;
		};
		^inval;
	}
}
/*
p = Plorenz(0.05).asStream;
{ p.next }.plot;
Array.fill(10, { p.next }).plot;
Plotter
Pbind(\freq, p.next).play;
(
Ê Ê var width = 500, height = 400, rate = 0.005;
Ê Ê var w, u;

Ê Ê w = Window("3d canvas demo", Rect(128, 64, width, height), false)
Ê Ê Ê Ê .front;

Ê Ê u = Canvas3D(w, Rect(0, 0, width, height))
Ê Ê Ê Ê .background_(Color.black)
Ê Ê Ê Ê .scale_(200)
Ê Ê Ê Ê .perspective_(0.4)
Ê Ê Ê Ê .distance_(2);

Ê Ê // add triangular spiral
Ê Ê u.add(Canvas3DItem()
Ê Ê Ê Ê .color_(Color.green)
Ê Ê Ê Ê .width_(1.5)
Ê Ê Ê Ê .paths_([900.collect { p.next * 0.1 }])
Ê Ê );

Ê Ê // add cube
Ê Ê u.add(Canvas3DItem.cube
Ê Ê Ê Ê .color_(Color.white)
Ê Ê Ê Ê .width_(2)
Ê Ê );

Ê Ê // spin canvas on mouse move
Ê Ê u.mouseMoveAction = {|v,x,y|
Ê Ê Ê Ê u.transforms = [
Ê Ê Ê Ê Ê Ê Canvas3D.mRotateY(x / -200 % 2pi),
Ê Ê Ê Ê Ê Ê Canvas3D.mRotateX(y / 200 % 2pi)
Ê Ê Ê Ê ];
Ê Ê Ê Ê u.refresh;
Ê Ê };

Ê Ê u.mouseMoveAction.value(nil, 50, 50); // initial rotation
)
*/
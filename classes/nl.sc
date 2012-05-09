
NLUGen : UGen {
	*categories {^#["UGens>nl4sc"] }
}

Logist0 : NLUGen {
	*ar {| freq=22050, r=1.5, xi=0.1, mul=1, add=0 |
		^this.multiNew(\audio, freq, r, xi).madd(mul, add);
	}
}
Logist1 : Logist0 {}
Logist2 : Logist0 {}

CML0 : NLUGen {
	*ar {| freq=22050, r=1.5, g=0.5, xi=0.1, mul=1.0, add=0.0 |
		^this.multiNew(\audio, freq, r, g, xi).madd(mul, add);
	}
}
CML1 : CML0 {}
CML2 : CML0 {}

GCM0 : NLUGen {
	*ar {| freq=22050, r=1.5, g=0.5, mul=1.0, add=0.0 |
		^this.multiNew(\audio, freq, r, g).madd(mul, add);
	}
}
GCM1 : GCM0 {}
GCM2 : GCM0 {}

HCM0 : NLUGen {
	*ar {| freq=22050, r=1.5, g=0.5, mul=1.0, add=0.0 |
		^this.multiNew(\audio, freq, r, g).madd(mul, add);
	}
}
HCM1 : HCM0 {}
HCM2 : HCM0 {}

Nagumo : NLUGen {
	*ar {| uh=0.01, vh=0.01, pulse=0 mul=1, add=0 |
		^this.multiNew(\audio, uh, vh, pulse).madd(mul, add);
	}
}
FIS : NLUGen {
	*ar {| r=3.5, xi=0.1, n=3, mul=1, add=0 |
		^this.multiNew(\audio, r, xi, n).madd(mul, add);
	}
}
TLogist : NLUGen {
	*kr {| r=1.5, xi=0.1, trg=0, mul=1, add=0 |
		^this.multiNew(\control, r, xi, trg).madd(mul, add)
	}
}

NLMultiOutUGen : MultiOutUGen {
	*categories {^#["UGens>nl4sc"] }
}

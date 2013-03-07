s.boot;
{ CML0.ar(SampleRate.ir/1.3, MouseY.kr(1.5, 1.9), MouseX.kr(0, 1)) * 0.1 }.scope;
{ CML1.ar(SampleRate.ir/10, MouseY.kr(1.5, 1.9), MouseX.kr(0, 1)) * 0.1 }.scope;
{ CML3.ar(SampleRate.ir/10, MouseY.kr(1.5, 1.9), MouseX.kr(0, 1)) * 0.1 }.scope;

{ GCM0.ar(SampleRate.ir/1.3, MouseY.kr(1.5, 1.9), MouseX.kr(0, 1)) * 0.1 }.scope;
{ GCM1.ar(SampleRate.ir/10, MouseY.kr(1.5, 1.9), MouseX.kr(0, 1)) * 0.1 }.scope;
{ GCM3.ar(SampleRate.ir/10, MouseY.kr(1.5, 1.9), MouseX.kr(0, 1)) * 0.1 }.scope;

{ HCM0.ar(SampleRate.ir/1.3, MouseY.kr(1.5, 1.9), MouseX.kr(0, 1)) * 0.1 }.scope;
{ HCM1.ar(SampleRate.ir/10, MouseY.kr(1.5, 1.9), MouseX.kr(0, 1)) * 0.1 }.scope;
{ HCM3.ar(SampleRate.ir/10, MouseY.kr(1.5, 1.9), MouseX.kr(0, 1)) * 0.1 }.scope;

{Out.ar(1, SinOsc.ar(200)) }.play;
//plug&play with pulse
{ Nagumo.ar(0.1, 0.001, LFPulse.ar(20, 0, MouseX.kr(0.0, 0.5), MouseY.kr(0.01, 0.5)), 0.1) }.play;
{ var sig = Nagumo.ar(0.1, 0.001, LFPulse.ar(1000, 0, MouseX.kr(0.0, 0.5), MouseY.kr(0.01, 0.75)), 0.1); [sig, Delay1.ar(sig)] }.play;
{ Nagumo.ar(0.1, 0.1, LFPulse.ar(4410, 0, MouseX.kr(0.0, 0.5), MouseY.kr(0.01, 0.5)), 0.1) }.play;
{ Nagumo.ar(0.1, 0.1, LFPulse.ar(4410, 0, MouseX.kr(0.0, 0.5), MouseY.kr(0.01, 0.5)), 0.1) }.play;

//non-pulse, watch-out. noisy!
{ Nagumo.ar(0.1, 0.00001, Impulse.ar(10, 0, MouseY.kr(0.01, 0.5))) }.play;
{ Nagumo.ar(0.1, 0.001, Impulse.ar(10, 0, MouseY.kr(0.01, 0.5))) }.play;
{ var sig = Nagumo.ar(0.1, 0.001, SinOsc.ar(20, LFNoise2.kr(10), MouseY.kr(0.01, 0.5)));[sig, Delay1.ar(sig)]}.play;
{ Nagumo.ar(0.1, 0.1, SinOsc.ar(50, 0, MouseY.kr(0.01, 0.5))) }.play;
{ Nagumo.ar(0.1, 0.1, BrownNoise.ar(MouseY.kr(0.01, 0.5))) }.play;
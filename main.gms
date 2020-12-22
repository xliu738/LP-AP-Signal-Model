option limcol=0, limrow=0;
option decimals=8;

set A 'set of attributes' /'latitude', 'longitude'/;
set N 'set of APs' /ap1*ap13/;

Table coord(A,N) 'table of decimal degree coordinates to work from'
$ONDELIM
$INCLUDE data1.csv
$OFFDELIM
;

alias(M,N);

parameter arcs(N,N);
arcs(M,N)$(ord(M) ne ord(N)) = yes;


* Haversine formula to convert to euclidean distance
scalar
  E_rad 'radius of Earth in km' /6373/;

parameters
  dlat(M,N) 'distance in latitude values for AP pairs in radians',
  dlong(M,N) 'distance in longitude values for AP pairs in radians',
  angle(M,N) 'squared angular distance of the decimal degrees',
  euc_dist(M,N) 'euclidean distance in meters';

dlat(M,N)$arcs(M,N) = pi/180 * (coord('latitude',N) - coord('latitude',M));
dlong(M,N)$arcs(M,N) = pi/180 * (coord('longitude',N) - coord('longitude',M));
angle(M,N)$arcs(M,N) = power(sin(dlat(M,N)/2), 2)+
                        cos(coord('latitude',M))*
                        cos(coord('latitude',N))*
                        power(sin(dlong(M,N)/2), 2);
euc_dist(M,N)$arcs(M,N) = E_rad * 2 * arctan2( sqrt(angle(M,N)), sqrt(1-angle(M,N))) * 1000;


free variables
  signal_strength;
  
positive variables
  x(M) 'power to set for AP M in mW',
  y(M,N) 'power received at AP N from AP M in mW';
  
scalar
  noise_thresh 'power under which signal will be treated as noise in mW' /1e-9/,
  wavelen 'wavelength of signal waves in m' /0.125/;
  
equations
  obj,
  pathloss(M,N),
  noise(N);

obj..
  signal_strength =e= sum((M,N)$arcs(M,N), y(M,N));
  
pathloss(M,N)$arcs(M,N)..
  y(M,N) =e= x(M) * power(wavelen,2) / (power((4*pi),2) * power(euc_dist(M,N),2));
  
noise(N)..
  sum(M, y(M,N))*1e9 =l= noise_thresh *1e9; 

x.up(M) = 1e-3;
x.lo(M) = 1e-9;
  
model main /all/;
solve main using lp maximizing signal_strength;

parameter upbound(M,N);
upbound(M,N) = noise_thresh * power((4*pi),2) * power(euc_dist(M,N),2) / power(wavelen,2);

parameter signal_dBm(M);
signal_dBm(M) = 10*log10(x.l(M));

option signal_dBm:1:0:1;
display signal_dBm;
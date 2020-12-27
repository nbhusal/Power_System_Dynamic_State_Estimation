function mpc=case9_new
%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data 
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
% I, IDE, PL, QL, GL, BL, area, VM,  VA,  BASEKVA, ZONE, Vmax	Vmin
mpc.bus = [ 
        1  3 0 0  0   0   1  1.1     0          16.5    1    1.1     0.9
        2  2 0 0  0   0   1  1.09736 0          18      1    1.1     0.9
        3  2 0 0  0   0   1  1.08662 0          13.8    1    1.1     0.9
        4  1 0 0  0   0   1  1       0          230     1    1.1     0.9
        5  1 125  50 0 0  1  1       0          230     1    1.1     0.9                                  
        6  1 90   30 0 0  1  1       0          230     1    1.1     0.9
        7  1 0 0  0   0   1  1       0          230     1    1.1     0.9
        8  1 100  35 0 0  1  1       0          230     1    1.1     0.9              
        9  1 0 0  0   0   1  1       0          230     1    1.1     0.9         
];
%% generator data
%I - Bus number
%ID - Machine identifier (0-9, A-Z)
%PG - MW output
%QG - MVAR output
%QT - Max MVAR
%QB - Min MVAR
%VS - Voltage setpoint
%IREG - Remote controlled bus index (must be type 1), zero to control own
%voltage, and must be zero for gen at swing bus
%MBASE - Total MVA base of this machine (or machines), defaults to system
%MVA base.
%ZR,ZX - Machine impedance, pu on MBASE
%RT,XT - Step up transformer impedance, p.u. on MBASE
%GTAP - Step up transformer off nominal turns ratio
%STAT - Machine status, 1 in service, 0 out of service
%RMPCT - Percent of total VARS required to hold voltage at bus IREG
%to come from bus I - for remote buses controlled by several generators
%PT - Max MW
%PB - Min MW

%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
%     I, PG,    QG,QT,  QB,   VS,     MBASE, status,Pmax, Pmin Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen=[ 
      1  160    0 300 -300  1.1      100   1      800  10  0  0	0	0	0	0	0	0	0	0	0	
      2  134.32 0 300 -300  1.09736  100   1      163  10  0  0	0	0	0	0	0	0	0	0	0	
      3  94.19  0 300 -300  1.08662  100   1      85   10  0  0	0	0	0	0	0	0	0	0	0	];


%% Branch Data 
%Branch records, ending with a record with from bus of zero
%I,J,CKT,R,X,B,RATEA,RATEB,RATEC,RATIO,ANGLE,GI,BI,GJ,BJ,ST
%I - From bus number
%J - To bus number
%CKT - Circuit identifier (two character) not clear if integer or alpha
%R - Resistance, per unit
%X - Reactance, per unit
%B - Total line charging, per unit
%RATEA - MVA rating A
%RATEB, RATEC - Higher MVA ratings
%RATIO - Transformer off nominal turns ratio
%ANGLE - Transformer phase shift angle
%GI,BI - Line shunt complex admittance for shunt at from end (I) bus, pu.
%GJ,BJ - Line shunt complex admittance for shunt at to end (J) bus, pu.
%ST - Initial branch status, 1 - in service, 0 - out of service
        %I, J,  R,     X,       B,   RATEA,RATEB,RATEC,RATIO,  ANGLE,  Status Angmin Angmax 
mpc.branch=[ 1  4   0      0.0576   0      400 500 600     16.5/230 0      1	-360	360;
         2  7   0      0.0625   0      400 500 600     18.0/230 0      1	-360	360;
         3  9   0      0.0586   0      400 500 600     13.8/230 0      1	-360	360;
         4  6   0.017  0.092  2*0.079  400 500 600     1        0      1	-360	360; 
         4  5   0.01   0.085  2*0.088  400 500 600     1        0      1	-360	360; 
         5  7   0.032  0.161  2*0.153  400 500 600     1        0      1	-360	360; 
         7  8   0.0085 0.072  2*0.0745 400 500 600     1        0      1	-360	360;
         8  9   0.0119 0.1008 2*0.1045 400 500 600     1        0      1	-360	360;
         6  9   0.039  0.17   2*0.179  400 500 600     1        0      1	-360	360;
];
  
end 

clc;
clear;
close all;  
F = 4;
m = 500;
a = F/m;
inV = a*60;
deltaV = [0; 0; 0];
mu = 3.986004418*10^14;
%% Initial simulate
tleFile = "ISS_OLD.TLE";
startTime = datetime(2023,09,04);
stopTime = datetime(2023,10,04);
sampleTime =3600;
sc1 = satelliteScenario(startTime,stopTime,sampleTime);

satSGP4 = satellite(sc1,tleFile, ...    
    "Name","ISSSGP4", ...
    "OrbitPropagator","sgp4");


[positionSGP4,velocitySGP4] = states(satSGP4,"CoordinateFrame", "ecef");
r_ecef = positionSGP4(:,size(positionSGP4, 2))';
v_ecef = velocitySGP4(:,size(positionSGP4, 2))';
magnitude = norm(v_ecef);
% Calculate the unit vector
unit_vector = v_ecef / magnitude;
[r_eci, v_eci] = ecef2eci(stopTime, r_ecef, v_ecef);
%% give the ISS a delta-V 
r_ijk = r_eci;
r = 6790000;
v = (norm(v_eci));
[a, ecc1, incl, RAAN, argp, nu1, truelon, arglat, lonper] = ijk2keplerian(r_ijk, v_eci);
ecc1 = sprintf('%07d', round(ecc1 * 1e6));
deltaV_f = (sqrt(mu*((2/r)-(1/6.796586702118547e+06))) - v);
deltaV = [0; deltaV_f; 0];
deltaV = [0; 1; 0];
deltaV_magnitude = norm(deltaV);
v_ijk = v_eci+(deltaV_magnitude * unit_vector)';
[positionSGP4_n1, velocitySGP4_n1] = convert_positions(positionSGP4, velocitySGP4);
positionSGP4_n_mean = mean(positionSGP4_n1);
[a, ecc, incl, RAAN, argp, nu, truelon, arglat, lonper] = ijk2keplerian(r_ijk, v_ijk);
T = (((2*pi)/(sqrt(mu)))*a^(3/2))/60;
%mean motion 
n = 1/T; %revolutions per min
n = n*60*24; %revolutions per day
%nu (rad)
nu = nu*(pi/180);
%Mean anomaly
M_e = [2*atan(sqrt((1 - ecc)/(1 + ecc))*tan(nu/2))-((ecc*sqrt(1 - ecc^2)*sin(nu))/(1 + ecc*cos(nu)))];
%Mean anomaly (deg)
M_e = M_e* (180/pi) ;   
%% update TLE file 
satellite_name = 'ISS (ZARYA)';
tle_filename = 'ISS.tle';

% Convert to formatted string
ecc = sprintf('%07d', round(ecc * 1e6));

fid = fopen(tle_filename, 'w');
fprintf(fid, '%s\n', satellite_name);
fprintf(fid, '1 25544U 98067A   23246.98284133  .00013126  00000-0  23590-3 0  9990\n');    
fprintf(fid, '2 25544 %8.4f %8.4f %s %8.4f %8.4f %11.8f1847', incl, RAAN, ecc, argp, M_e, n);
fclose(fid);
fprintf('新的TLE檔案已生成:%s\n',tle_filename);
%% input new TLE to SGP4

tleFile = "ISS.tle";
startTime = datetime(2023,10,04);
stopTime = datetime(2023,11,04);
sampleTime =3600;
sc2 = satelliteScenario(startTime,stopTime,sampleTime);

satSGP42 = satellite(sc2,tleFile, ...    
    "Name","ISSSGP4", ...
    "OrbitPropagator","sgp4");


[positionSGP4_n,velocitySGP4_n] = states(satSGP42,"CoordinateFrame", "ecef");
[positionSGP4_n2, velocitySGP4_n2] = convert_positions(positionSGP4_n, velocitySGP4_n);
positionSGP4_n_mean2 = mean(positionSGP4_n2);
new = [positionSGP4_n1,positionSGP4_n2];

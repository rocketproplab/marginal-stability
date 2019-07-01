clear all; close all; clc;


% This program is adapted from Huzel and Huang's "Modern Engineering for
% Design of Liquid-Propellant Rocket Engines" and uses the calculations for
% tank dimensions found in Chapter 8

%% Lists all required constants and values required for math-ing
P = input('Propellant Tank Pressure (psi): '); %input tank pressure in psig

D = 12; %tank outer diameter in inches
R = D./2; %tank radius, in
S = 40000; %yield strength of 6061-T6 alloy in psi
Sc = 56000; %compressive yield strength of 6061-T6 in psi
FoS = 1.5; %factor of safety
v = 0.33; %Poisson's ratio of 6061-T6 alloy, unitless
ew = 0.85; %approximate efficiency of welded joints
E = 1e7; %6061-T6 modulus of elasticity, psi
rho_alum = 0.0975; %density of aluminum, lbm./in.^3
rho_lox = 0.04122124; %density of LOx, lbm./in.^3
rho_rp1 = 0.0292631; %density of RP-1, lbm./in.^3
k = 2; %ratio of a./b, which should be 2:1 due to the semi-elliptical tank head design
K = 1.2; %obtained from graph on page 292 of H&H book
Fc = 3305; %estimated critical loading on rocket during flight, lbf
g = 32.2; %gravitational acceleration, ft./s.^3

mix_ratio = 2.23; %fuel mix ratio, for prop mass calcs and volume
m_prop = 1180.799; %propellant mass in lbm--given from prop team
m_rp1 = m_prop ./ (1+mix_ratio); %computes needed amount of RP-1 based on mixture ratio
m_lox = m_prop - m_rp1; %remaining propellant must be LOx
vol_lox = m_lox ./ rho_lox;%use density of propellants to compute volumes
vol_lox = vol_lox .* 1.02; %adds an extra 2% of volume for ullage space
vol_rp1 = m_rp1 ./ rho_rp1;
vol_rp1 = vol_rp1 .* 1.02; %also adds an extra 2% of volume for ullage

%not required unless volume is given instead of mass
% vol_lox = 16748.707; %volume of oxidizer in in.^3
% vol_rp1 = 10579.792; %volume of fuel in in.^3
% d_pressline = 0.5; %outer diameter of pressurant lines through tanks, in
% d_annulus = 2; %outer diameter of annular lining, in
% d_rp1fill = 1; %outer diameter of fill line to rp-1, in
d_pressline = 0;
d_annulus = 0;
d_rp1fill = 0;
%as of new design, all plumbing is external :)
t_barlows = (P.*D)./(2.*S./FoS);

%% Calculates needed tank thickness based on loading characteristics
% Algebraic manipulation of eqn. 8-32 on p. 293 required to caclculate
% required cylindrical wall thickness based on loading criteria

tc_vec = 0:0.001:0.25; %vector of wall thicknesses for testing
LHS = 4.*(Fc.*FoS).*(1-v.^2) ./ (E.*pi); %the left-hand side of the eqn. on p.293
RHS = (36.*tc_vec.^3 - tc_vec.^5) ./ (6-tc_vec).^3;

indexer = find(RHS <= LHS); %finds intersection of two curves
tc = tc_vec(max(indexer)); %saves tc--cylindrical wall thickness

%% test to see whether barlows is a bit more accurate than the other
if t_barlows >= tc
    tc = t_barlows;
end
%tc = t_barlows;
%% Uses required cylindrical wall thickness to estimate bulkhead thickness

a = R - tc; %a is the internal radius, or major axis length of ellipse
b = a./2; %minor axis length, or head height--property of 2:1 shape
Eprime = 2.*k + (1./sqrt(k.^2-1)).*log((k+sqrt(k.^2-1))./(k-sqrt(k.^2-1)));
%eprime is the "design factor" described on p. 292, for weight calcs later
%on

tk = (K.*P.*a) ./ (S.*ew ./ FoS); %knuckle thickness, in
tcr = (P.*R) ./ (2.*S.*ew ./ FoS); %crown thickness, in
tequiv = (tk + tcr) ./ 2; %average bulkhead thickness for weight calcs

%% Redimensions tanks according to volume and wetted area eqns. given in H&H

v_endcap = (2.*pi.*a.^2.*b)./3; %endcap volume capacity, in.^3
leftovervol_lox = vol_lox - 2.*v_endcap; %volume left for cylindrical section
leftovervol_rp1 = vol_rp1 - 2*v_endcap; %vol left for rp-1 cylindrical section
%this has been adjusted due to the fact that the RP-1 tank is now its own
%entity, and not a common bulkhead.  The addition of a convex bulkhead adds
%volume to the tank itself and removes the need to exclude the bulkhead
%volume
vcyl_lox = leftovervol_lox; %lox cylinder volume needed to hold remaining propellant
vcyl_rp1 = leftovervol_rp1; %rp-1 cylinder volume needed to hold remaining propellant

cylheight_lox = vcyl_lox ./ (pi.*a.^2 - pi.*d_rp1fill.^2 - 2.*pi.*(d_pressline ./ 2).^2); %computes height of cylinders for each propellant
cylheight_rp1 = vcyl_rp1 ./ (pi.*a.^2 - pi.*(d_pressline ./ 2).^2 - pi.*(d_annulus ./ 2).^2);

intertank_airframe_height = 12; %height of intertank airframe in inches
overall_height = cylheight_lox + cylheight_rp1 + 4.*(b+tcr) + intertank_airframe_height; %head height, with cylindrical sections and two end caps

%% Computes weights of tank segments also based on given equations
%pi.*(outer_diameter.^2 - inner_diameter.^2) .* height .* rho_alum ./ 4
% weight of bulkheads in lbm for calculating--from SolidWorks
loxupper_weight = 5.83;
loxlower_weight = 5.17;
rp1upper_weight = 5.29;
rp1lower_weight = 5.65;
w_endcaps = loxupper_weight + loxlower_weight + rp1upper_weight + rp1lower_weight;

%weight of respective cylindrical sections
w_cyl_rp1 = pi.*(R.^2 - a.^2) .* cylheight_rp1 .* rho_alum;
w_cyl_lox = pi.*(R.^2 - a.^2) .* cylheight_lox .* rho_alum;

% all plumbing weights are no longer included
dry_weight = w_endcaps + w_cyl_rp1 + w_cyl_lox;
wet_weight = dry_weight + m_prop;



%% Pressurant Tank Sizing Code
prop_vol = vol_lox + vol_rp1;
prop_press = P;
gamma = 5/3; %ratio of specific heats for Helium
press_temp = 536.67; %pressurant temp in Rankine
R = 10.73159*12^3; %gas constant in imperial units, in3-psi / R-lb-mol
molar_mass = 4.003;
initial_press = 3600; %initial pressure in bottle, psi
final_press = P; %final pressure at burnout, psi
press_tank_diam = 11; %outer diameter of pressurant tank, in, based on Ray's SW model


press_mass = ((prop_press * prop_vol) / (R * press_temp)) * (gamma / (1-(final_press/initial_press))) * molar_mass;
press_vol = (press_mass * (1/molar_mass) * R * press_temp) / initial_press;
% Pressurant Tank sizing calcs also excluded since they suck

%% Print the results for ease of design and troubleshooting if needed
fprintf('Propellant Volume: %f in.^3\n', vol_lox+vol_rp1)
fprintf('Tank Pressure: %f psi\n', P)
fprintf('Cylindrical Wall Thickness: %f in\n', tc)
fprintf('Cylindrical Thickness--Barlow''s Eqn.(comparison): %f in\n', (P.*D)./(2.*S./FoS))
%I threw Barlow's estimate in the printed response to make sure that these
%calculations, which take into account the axial loading of max Q and
%thrust, are larger than just plain internal pressure estimates.  This is
%just a handy troubleshooting feature./sanity check to make sure everything
%is at least in the ballpark
fprintf('Equivalent Bulkhead Thickness: %f in\n', tequiv)
fprintf('LOx Cylindrical Wall Height: %f in\n', cylheight_lox)
fprintf('RP-1 Cylindrical Wall Height: %f in\n', cylheight_rp1)
fprintf('Approximate Overall Tank Height: %f in\n', overall_height)
fprintf('Approximate Dry Mass: %f lbm\n', dry_weight)
fprintf('Approximate Wet Mass: %f lbm\n', wet_weight)
fprintf('Pressurant Volume and Mass Required: %f in^3 and %f lbm\n', press_vol, press_mass)
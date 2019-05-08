% tylka.m, Matthias Brennwald
%
% This program calculates the approximated power response of a loudspeaker using the method described by JG Tylka and EY Choueiri in "On the Calculation of Full and Partial Directivity Indices", available for download at https://www.princeton.edu/3D3A/Publications/Tylka_3D3A_DICalculation.pdf

% Revisions:
% - 2 May 2019: first version
% - 3 May 2019: changed computation of omega(n): just calculate first quarter of the whole series, then use periodicity in omega to complete the sequence as described in Tylka
% - 4 May 2019: fixed an issue where the H_m,n SPL data for a given orbit/angle were loaded from data file that belonged to a different orbit/angle. Also made the code more flexible by automagically figuring out the number of points per orbit, and the number of points per SPL curve.

% Notes:
%
% (1) While Tylka describes the calculation method for the directivity index (DI), the power response is just the denominator term of equation (1) in the Tylka paper. This term and it's equivalents throughout the Tylka paper are used here to determine the power response.
%
% (2) The summation indices n and m in the Tylka paper are defined as n = 0...N-1 and m=0,1 (zero based). For use with GNU Octave (or Matlab) in this program, this convention is changed to n = 1...N and m = 1,2.
%
% (3) This program loads the H_m,n data from ASCII data files reflecting the SPL responses measured on two orbits around a loudspeaker (one horizontal and hone vertical orbit), as described by Tylka. The points are distributed equally along the two orbits.
%
% (4) The angular positions of the SPL data in the ASCII files is specified in the file name, which follows the 'SomeFileNameXYZ <orb> <ang>.txt' format. The file name defines if the data belongs the horizontal orbit (<orb>="hor") or vertical orbit (<orb>="ver"), and provides the angle (<ang>) of the point on the orbit where the SPL data was measured. The angle is measured in degrees relative to the "listening axis", and is positive upwards (vertical orbit) or to the right (horizontal orbit, viewed from the front of the speaker). Negative angles are used for points on the left or bottom of the speaker.
% Examples:
%   - 'somefile hor 25.txt' contains the SPL data measured on the horizontal orbit at 25° to the right of the main listening axis (as viewed from the front of the speaker).
%   - 'somefile hor -70.txt' contains the SPL data measured on the horizontal orbit at 70 to the left of the main listening axis.
%   - 'somefile hor -105.txt' contains the SPL data measured on the horizontal orbit at 105 to the left of the main listening axis (and slightly behind the speaker).
%   - 'somefile ver 65.txt' contains the SPL data on the vertical orbit at 65° above the main listening axis.
%   - 'somefile ver -155.txt' contains the SPL data on the vertical orbit at 155° below the main listening axis (and behind the speaker);
%
% (5) For purposes of testing and demonstration, a set of ASCII data files is provided with this program in the directory "testdata". The data in these files were generated using the diffraction calculator in Vituix CAD.



% Ask user for data directory:
disp ('Select directory containing the data files...')
D = uigetdir (); % directory containing the ASCII data files
F = readdir(D); % list of files in directory D



% Determine base name of data files:
B = []; i = 1;
while isempty(B)
	if i > length(F)
		error (sprintf('Could not find any *.txt file in directory %s.',D))
		break
	end
	[dd,nn,ee] = fileparts (F{i});
	if strcmp (ee,'.txt')
		if any ( j = findstr (nn,' hor ') )
			B = nn(1:j-1);
		elseif any ( j = findstr (nn,' ver ') )
			B = nn(1:j-1);
		end
	end
	i = i+1;
end
disp (sprintf('Base name of data files: %s',B))

% filter matching files:
k = strfind (F,B);
l = repmat (0,size(k));
for i = 1:length(l)
	l(i) = ~isempty(k{i});
end
F = F(find(l)); % only keep file names that match the base name




% Determine data size (length of SPL curves, angular spacing):

disp ('Determine size and spacing of data...')

% determine angles from file names
deg_hor = deg_ver = [];
for i = 1:length(F)
	if findstr(F{i},' hor ')
		u = strsplit (F{i},' hor '){2};
		deg_hor = [ deg_hor str2num(strsplit(u,'.txt'){1}) ];
	elseif findstr(F{i},' ver ')
		u = strsplit (F{i},' ver '){2};
		deg_ver = [ deg_hor str2num(strsplit(u,'.txt'){1}) ];
	end
end
deg_hor = unique (deg_hor);
deg_ver = unique (deg_ver);
n_hor = length (deg_hor);
n_ver = length (deg_ver);
deg_hor_delta = diff ( deg_hor );
deg_ver_delta = diff ( deg_ver );
if min(deg_hor_delta) < max(deg_hor_delta)
	error ('Horizontal angles are not equally spaced!')
elseif min(deg_ver_delta) < max(deg_ver_delta)
	error ('Vertical angles are not equally spaced!');
end
deg_hor_delta = deg_hor_delta(1);
deg_ver_delta = deg_ver_delta(1);
if deg_hor_delta ~= deg_ver_delta
	error ('Angular spacing on horizontal and vertical orbits are not the same!')
end
deg_delta = deg_hor_delta;
disp (sprintf('Angular spacing: %g degrees',deg_delta))
N = 360 / deg_delta; % number of points per orbit


% Determine length of SPL curves:
[x1,x2,x3] = textread ( [D filesep F{i}] , '%f%f%f' , 'headerlines',1 );
L = length (x1); % number of data points in the SPL response curves



% Read data files:
disp ('Reading data files...')
deg = [ [0:N/2-1]*deg_delta  180-[N/2+1:N]*deg_delta ]; % theta(n) values corresponding to H(m,n,:) as used by Tylka equation (2)

% Init array of SPL data:
H = repmat(NaN,2,N,L); % H(m,n) is the SPL response at points m,n, length of each SPL response curve is L

% load files for all theta(m,n) / phi(m,n):
for m = 1:2
	for n = 1:N
		if m == 1
			orbit = 'hor';
		elseif m == 2
			orbit = 'ver';
		end
		fn = sprintf('%s %s %i.txt',B,orbit,deg(n)); % file name
		if deg(n) == -180
			fn = strrep (fn,' -180',' 180');
		end
		if deg(n) == 0 || deg(n) == -180
			if strcmp(orbit,'ver')
				fn = strrep (fn,' ver ',' hor ');
			end
		end
		disp (sprintf('   m=%i, n=%i: loading file for orbit = %s, angle = %i deg. (%s)',m,n,orbit,deg(n),fn));
		[x1,x2,x3] = textread ( [D filesep fn] , '%f%f%f' , 'headerlines',1 ); % x1: frequency, x2: SPL in dB, x3: phase in degrees
		H(m,n,:) = 10.^(x2/20); % convert SPL from dB to linear scale
	end % n
end % m
f = x1; % frequency



% for debugging / code testing:
if exist('kk')
	% Override H(m,n) from data files with H(m,n)=1 for n = kk and H(m,n)=0 elsewhere.
	% The resulting power response should be flat and be equal to 2*omega(kk).

	H = repmat (0,size(H)); % replace H by uniform polar response of a perfect point source with flat SPL 
	H(:,kk,:) = 1;
end

% calculate omega_n coefficients (following equations 9...12 in the Tylka paper):
% note: the code below was verified to give the same numbers as those listed in Tab. 1 in the Tylka paper for N = 72 / delta_theta = 5°.

function om = __omega (theta) % equation (9) in Tylka paper
	om = 2*pi*(1-cos(theta));
end % function

delta_theta = 2*pi/N; % angular spacing between measurement in each orbit

% determine summation coefficients for use in equation (3); see also Fig. 2
omega = repmat (NaN,1,N/4+1); % first N/4+1 values
omega(1) = 1/4/pi * __omega(delta_theta/2) / 2; % "cap" centered at theta = phi = 0, equation (10)
for n = 2:N/4+1
	omega(n) = 1/4/pi * ( __omega((n-1)*delta_theta+delta_theta/2) - __omega((n-1)*delta_theta-delta_theta/2) ) / 4  ;
end

% remaining values for n = N/4+2...N:
omega = [ omega fliplr(omega(1:end-1)) ];
omega = [ omega omega(2:end-1) ];




%%%%%% Numbers from Vituix CAD:
%%%%%% https://www.diyaudio.com/forums/software-tools/307910-vituixcad-post5779704.html
v = [ 0.000237559492495 , 0.000950402661548 , 0.001893572184973 , 0.0028223304807 , 0.003729609137299 , 0.004608503216363 , 0.005452323803261 , 0.006254648913813 , 0.007009372369471 , 0.007710750269022 , 0.008353444703147 , 0.008932564379134 , 0.009443701846565 , 0.009882967040682 , 0.01024701688812 , 0.010533080749721 , 0.010738981506769 , 0.010863152130176 , 0.010904647606522 , 0.010863152130176 , 0.010738981506769 , 0.010533080749721 , 0.01024701688812 , 0.009882967040682 , 0.009443701846565 , 0.008932564379134 , 0.008353444703147 , 0.007710750269022 , 0.007009372369471 , 0.006254648913813 , 0.005452323803261 , 0.004608503216363 , 0.003729609137299 , 0.0028223304807 , 0.001893572184973 , 0.000950402661548 , 0.000118779746248 ];

%% omega = [ v(1:end-1) v(1:end-1) ]; % use the Vituix coefficients instead



% calculate power response using the denominator of equation (3)
Sn  = sum ( omega .* H.^2 , 2 ); % summation over n=1:N
Smn = Sn(1,:) + Sn(2,:); % summation over m=1,2
PR  = 10*log10(4*pi*Smn); % power response in dB; 4pi term is the required to normalise to the surface of the unit sphere


% debugging info:
if exist('kk')
	disp(Smn(1))
	disp(2*omega(kk))
end


% plot polar response diagram (horizontal orbit) of input data:

% convert H to dB:
SPL = 20*log10( reshape( H(1,:),N,length(f) ) ); % matrix of SPL data (horizontal orbit), in dB

% rearrange for monotonous horizontal angles:
[u,k] = sort(deg); SPL = SPL(k,:);

% contour plot:
figure(1)
contourf(f,[-180 deg_hor],[SPL(end,:);SPL],linspace(min(min(SPL)),max(max(SPL)),30))
set(gca,'xscale','log')
set(gca,'ytick',[-180:90:180]);
axis([20 20e3 -180 180])
colormap ('jet');
xlabel ('Frequency (Hz)')
ylabel ('Horizontal angle (degrees)')
title ('Polar SPL response (horizontal, dB-SPL)')



% plot power response:
figure(2)
semilogx (f,PR);
xlim([20 20e3])
xlabel ('Frequency (Hz)')
ylabel ('Power response (dB)')
title ('Power Response according to Tylka eqn. (3)')

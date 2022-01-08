% Proracun mosta Metodom Konacnih Elemenata
% Primer proracuna Metodom Konacnih Elemenata sa Ojler Bernulijevim Gredama 
% implementiran u jednom GNU Octave fajlu.
% Autor: Milos D. Petrasinovic <mpetrasinovic@mas.bg.ac.rs>
% Proracun strukture letelica
% Masinski fakultet, Univerzitet u Beogradu
% Katedra za vazduhoplovstvo, Struktura letelica
% https://vazmfb.com
% Beograd, 2022
%
% ---------------
%
% Copyright (C) 2022 Milos Petrasinovic <info@vazmfb.com>
%  
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as 
% published by the Free Software Foundation, either version 3 of the 
% License, or (at your option) any later version.
%   
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%   
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% ---------------
clear('all'), clc, close('all'), tic
disp([' --- ' mfilename ' --- ']);

% - Ulazni podaci
% Dimenzije mosta
H = 3000; % [mm] visina mosta
B = 2500; % [mm] sirina mosta
L = 12000; % [mm] duzina mosta

% Karaktetristike poprecnog preseka
Ds = 60; % [mm] spoljni precnik cevi
t = 5; % [mm] debljina zida cevi
% Karakteristike materijala
E = 210e9; % [Pa] Modul elasticnosti
ni = 0.3; % [-] Poasonov koeficijent

% Koordinate cvorova
% nc(cvor, 1:3) = [x, y, z];
% Leva strana konstrukcije
nc = [0, 0, 0; L/4, 0, 0; L/2, 0, 0; L*3/4, 0, 0; L, 0, 0; 
  L/4, 0, H; L/2, 0, H;L*3/4, 0, H;];
% Desna strana konstrukcije
nc = [nc; nc+[0, B, 0];]; 

% Cvorovi elemenata 
% en(element, :) = [cvor1, cvor2];
% Leva strana konstrukcije
en = [1, 2; 1, 6; 6, 2; 2, 7; 6, 7; 2, 3; 
  3, 7; 7, 8; 3, 4; 7, 4; 4, 8; 4, 5; 8, 5;]; 
% Desna strana konstrukcije
en = [en; en+[8, 8]]; 
% Rasponke
en = [en; (1:8).', (9:16).']; 

% Povezivanje cvorova
% ec(i, 1:8) = [cvor1, cvor2, Tx, Ty, Tz, Rx, Ry, Rz];
% 0 - slobodan
% 1 - povezan
ec = [];

% Konturni uslovi
% bc(cvor, 1:6) = [Tx, Ty, Tz, Rx, Ry, Rz];
% 0 - slobodan
% 1 - ogranicen
bc([1, 9], :) = repmat([1, 1, 1, 0, 0, 0], 2, 1);
bc([5, 13], :) = repmat([0, 1, 1, 0, 0, 0], 2, 1);

% Spoljasnja opterecenja
% F(cvor, 1:6) = [Fx[N], Fy[N], Fz[N], Mx[Nm], My[Nm], Mz[Nm]];
F(3, :) = [0, 0, -20000, 0, 0, 0];
F(11, :) = [2000, 0, -20000, 0, 0, 0];

% Dodatne promenljive
s = [1000, 200, 300, 3*10^-2, 1, 1000, 0.2, 100]; % Koeficijenti za skaliranje  

% - Definisanje matrice krutosti
d2r = 180/pi; % stepeni u radijane
ne = size(en, 1); % broj elemenata
nn = size(nc, 1); % broj cvorova
nec = size(ec, 1); % broj veza cvorova
dof = 6*nn; % broj stepeni slobode
D = zeros(dof, 1); % pomeranja
bc(end+1:nn, :) = zeros(length(size(bc, 1)+1:nn), 6); % konturni uslovi
F(end+1:nn, :) = zeros(length(size(F, 1)+1:nn), 6); % spoljasnja opterecenja
F = reshape(F.', [], 1); 
K = zeros(dof, dof); % matrica krutosti

% Stepeni slobode cvorova
rdof = find(bc')'; % ograniceni stepeni slobode
fdof = find(~bc')'; % slobodni stepeni slobode

% Povezivanje cvorova
Cp = eye(6);
Np = 0; % broj dodatnih stepeni slobode
if(nec)
    Np = size(ec(ec(:, 3:8)>0), 1);
end
lmn = (1:Np)+dof; % dodatni stepeni slobode

Q = zeros(Np, 1); % vektor za prosirenje cvornih opterecenja
C = zeros(Np, dof); % matrica za prosirenje matrice krutosti
j = 0;
for i = 1:nec
  Npi = size(ec(ec(i, 3:8)>0), 2);
  C(j+1:j+Npi, (ec(i, 1)-1)*6+(1:6)) = Cp(ec(i, 3:8)>0, :);
  C(j+1:j+Npi, (ec(i, 2)-1)*6+(1:6)) = -Cp(ec(i, 3:8)>0, :);
  j = j+Npi;
end

% Karakteristike materijala i poprecnog preseka
G = E/(2*(1+ni)); % Modul klizanja

% Karakteristike poprecnog preseka (kruzna cev)
Ds = Ds/s(1); % [m]
t = t/s(1); % [m]
A = (Ds^2-(Ds-2*t)^2)*pi/4; % [m^2] Povrsina poprecnog preseka
Iy = (Ds^4-(Ds-2*t)^4)*pi/64; % [m^4] Aksijalni moment inercije za osu y
Iz = Iy; % [m^4] Aksijalni moment inercije za osu z
J = (Ds^4-(Ds-2*t)^4)*pi/32; % [m^4] Polarni moment inercije

% Odredjivanje matrice krutosti
nc = nc/s(1); % [m]
for i=1:ne 
    j = en(i, :);       
    edof = [6*j(1)-5 6*j(1)-4 6*j(1)-3 ... % stepeni slobode elementa
           6*j(1)-2 6*j(1)-1 6*j(1)...
           6*j(2)-5 6*j(2)-4 6*j(2)-3 ...
           6*j(2)-2 6*j(2)-1 6*j(2)]; 
    L_e = sqrt((nc(j(2), 1)-nc(j(1), 1))*(nc(j(2), 1)-... % duzina elementa
        nc(j(1), 1))+(nc(j(2), 2)-nc(j(1), 2))*(nc(j(2), 2)-...
        nc(j(1), 2))+(nc(j(2), 3)-nc(j(1), 3))*(nc(j(2), 3)-nc(j(1), 3)));
    ki = [E*A/L_e, 12*E*Iz/(L_e^3), 6*E*Iz/(L_e^2), 4*E*Iz/L_e, ...
        2*E*Iz/L_e, 12*E*Iy/(L_e^3), 6*E*Iy/(L_e^2), ...
        4*E*Iy/L_e, 2*E*Iy/L_e, G*J/L_e];
    a = diag([ki(1), ki(2), ki(6)]);
    b(2, 3) = ki(3); b(3, 2) = -ki(7);
    c = diag([ki(10), ki(8), ki(4)]);
    d = diag([-ki(10), ki(9), ki(5)]);
    
    % Matrica krutosti elementa u lokalnom ks elementa
    k_e = [a, b, -a, b; b.', c, b, d; -a.', b.', a, -b; b.', d.', -b.', c];
  
    if nc(j(1), 1) == nc(j(2), 1) && nc(j(1), 2) == nc(j(2), 2)
        if nc(j(2), 3) > nc(j(1), 3)
            r = [0, 0, 1; 0, 1, 0; -1, 0, 0];
        else
            r = [0, 0, -1; 0, 1, 0; 1, 0, 0];
        end
    else
        CXx = (nc(j(2), 1)-nc(j(1), 1))/L_e;
        CYx = (nc(j(2), 2)-nc(j(1), 2))/L_e;
        CZx = (nc(j(2), 3)-nc(j(1), 3))/L_e;
        i_eXY = sqrt(CXx^2+CYx^2);
        CXy = -CYx/i_eXY;
        CYy = CXx/i_eXY;
        CZy = 0;
        CXz = -CXx*CZx/i_eXY;
        CYz = -CYx*CZx/i_eXY;
        CZz = i_eXY;
        r = [CXx, CYx, CZx; CXy, CYy, CZy; CXz, CYz, CZz];
    end
    
    T = [r, zeros(3, 9);  % Matrica transformacije
          zeros(3), r, zeros(3, 6);
          zeros(3, 6), r, zeros(3); 
          zeros(3, 9), r];
        
    K(edof, edof) = K(edof, edof)+T.'*k_e*T; % Matrica krutosti elementa
end  

% - Resavanje MKE
Kp = [K, C.'; C, zeros(Np, Np)]; % prosirena matrica krutosti
Fp = [F; Q]; % prosiren vektor spoljasnjih opterecenja

% Resavanje redukovane jednacine konacnih elemenata
D1 = Kp([fdof, lmn], [fdof, lmn])\Fp([fdof, lmn]); 

D([fdof, lmn]) = D1;
Rp = Kp*D; % prosiren vektor reakcija veza
lam = D((end+1-Np):end); % Lagranzevi mnozioci
R = Rp(rdof); % reakcije veza
D = reshape(D(1:(end-Np)), 3, []); % pomeranja cvorova
DT = D(:, 1:2:end); % translatorna pomeranja
DR = D(:, 2:2:end); % ugaona pomeranja

% - Prikaz rezultata
% Prikaz vrednosti pomeranja i reakcija veza
disp(' -------------------- ');
disp(' Pomeranja')
i = 1:dof/2;
DTv = [reshape(repmat(1:nn, 3, 1), 1, []); ...
  reshape(repmat(1:3, nn, 1).', [], 1).'; DT(i)*s(1)];
disp(' Cvor | Komponenta | Pomeranje [mm]');
fprintf(' %3d | %3d | %14.10f\n', DTv);

DRv = [reshape(repmat(1:nn, 3, 1), 1, []); ...
  reshape(repmat(4:6, nn, 1).', [], 1).'; DR(i)/d2r];
disp(' Cvor | Komponenta | Pomeranje [deg]');
fprintf(' %3d | %3d | %14.10f\n', DRv);

% Prikaz vrednosti reakcija
disp(' -------------------- ');
disp(' Reakcije veza')
nrdof = mod(rdof, 6); % ograniceni stepen slobode u cvoru
nrdof(nrdof == 0) = ones(length(find(nrdof == 0)), 1)*6;
nr = (rdof-nrdof)/6+1; % cvor sa ogranicenjem
FTv = [nr(nrdof < 4).', nrdof(nrdof < 4).', R(nrdof < 4)].';
disp(' Cvor | Komponenta | Reakcije [N]');
fprintf(' %3d | %3d | %14.10f\n', FTv);

FRv = [nr(nrdof > 3).', nrdof(nrdof > 3).', R(nrdof > 3)].';
disp(' Cvor | Komponenta | Reakcije [Nm]');
fprintf(' %3d | %3d | %14.10f\n', FRv);

% - Prikaz proracunskog modela
% Prikaz polaznog modela sa opterecenjima
disp(' -------------------- ');
disp(' Prikaz proracunskog modela... ');
drawArrow = @(x, y, z, varargin) quiver3(x(1), y(1), z(1), ...
    x(2)+10^-5, y(2)+10^-5, z(2)+10^-5, 0, varargin{:});   
drawArrowMarker = @(x, y, z, m, varargin) [plot3([x(1); x(1)+x(2)], ...
  [y(1); y(1)+y(2)], [z(1); z(1)+z(2)], '-', varargin{:}), ...
  plot3(x(1)+x(2), y(1)+y(2), z(1)+z(2), m, varargin{:})];    

nc = nc*s(1);
bms = [[nc(en(:, 1), 1), nc(en(:, 2), 1)], ...
    [nc(en(:, 1), 2), nc(en(:, 2), 2)], ...
    [nc(en(:, 1), 3), nc(en(:, 2), 3)]]; % grede
  
figure(1);
box on, grid on, hold on
c = get(gca, 'colororder');
plot3(bms(:, 1:2).', bms(:, 3:4).', bms(:, 5:6).', ...
    'LineWidth', 2, 'Color', c(1, :));
plot3(nc(:, 1), nc(:, 2), nc(:, 3), 'o', 'LineWidth', 2, ...
    'MarkerFaceColor', c(2, :), 'Color', c(2, :));
text(nc(:, 1)+s(8), nc(:, 2)+s(8), nc(:, 3), ...
   num2str((1:nn).'), 'FontSize', 18, 'Color', c(2, :));
rotate3d, view(45, 15), axis equal, set(gca, 'FontSize', 18);
xlim([-s(6), L+s(6)]), ylim([-s(6), B+s(6)]), zlim([-s(6), H+s(6)])

% Prikaz konturnih uslova
for i=1:length(rdof)
  p = zeros(1, 3);
  if(nrdof(i) < 4) 
    p(nrdof(i)) = 1*s(3);
    drawArrowMarker([nc(nr(i), 1), p(1)], [nc(nr(i), 2), p(2)], ...
      [nc(nr(i), 3), p(3)], '+', 'Color', 'b', 'LineWidth', 2, ...
      'MarkerFaceColor', 'b');
  else
    p(nrdof(i)-3) = -1*s(3);
    drawArrowMarker([nc(nr(i), 1), p(1)], [nc(nr(i), 2), p(2)], ...
      [nc(nr(i), 3), p(3)], 'o', 'Color', 'r', 'LineWidth', 2, ...
      'MarkerFaceColor', 'r');
  end
end

% Prikaz spoljasnjih opterecenja
for i=1:nn
  p = zeros(1, 3);
  if(norm(F((i-1)*6+1:(i-1)*6+3)) > 0)
    p = (F((i-1)*6+1:(i-1)*6+3))*s(4);
    drawArrow([nc(i, 1), p(1)], [nc(i, 2), p(2)], ...
      [nc(i, 3), p(3)], 'Color', 'b', 'LineWidth', 2, ...
      'MarkerFaceColor', 'b', 'MaxHeadSize', s(7));
  end
  if(norm(F((i-1)*6+4:(i-1)*6+6)) > 0)
    p = (F((i-1)*6+4:(i-1)*6+6))*s(5);
    drawArrow([nc(i, 1), p(1)], [nc(i, 2), p(2)], ...
      [nc(i, 3), p(3)], 'Color', 'r', 'LineWidth', 2, ...
      'MarkerFaceColor', 'r', 'MaxHeadSize', s(7));
  end
end

% Prikaz deformisanog modela sa opterecenjima
disp(' Prikaz deformisanog modela... ');
ncd = nc+(DT.'.*s(1)).*s(2);
bmsd = [[ncd(en(:, 1), 1), ncd(en(:, 2), 1)], ...
    [ncd(en(:, 1), 2), ncd(en(:, 2), 2)], ...
    [ncd(en(:, 1), 3), ncd(en(:, 2), 3)]]; % grede

figure(2);
box on, grid on, hold on
plot3(bms(:, 1:2).', bms(:, 3:4).', bms(:, 5:6).', '--', ...
    'LineWidth', 1, 'Color', c(1, :));
plot3(bmsd(:, 1:2).', bmsd(:, 3:4).', bmsd(:, 5:6).', ...
    'LineWidth', 2, 'Color', c(2, :));
plot3(ncd(:, 1), ncd(:, 2), ncd(:, 3), 'o', 'LineWidth', 2, ...
    'MarkerFaceColor', c(2, :), 'Color', c(2, :));
rotate3d, view(45, 15), axis equal, set(gca, 'FontSize', 18);
xlim([-s(6), L+s(6)]), ylim([-s(6), B+s(6)]), zlim([-s(6), H+s(6)])

% - Kraj programa
disp(' -------------------- ');
disp(' Svaciji i prema svakom jednaki, korisni, podignuti');
disp(' uvek smisleno, na mestu na kom se ukrstava najveci broj');
disp(' ljudskih potreba, istrajniji su od drugih gradjevina');
disp(' i ne sluze nicem sto je tajno ili zlo.');
disp('  - Mostovi, Ivo Andric');
disp(' -------------------- ');
disp(' Program je uspesno izvrsen... ');
disp([' Potrebno vreme: ' num2str(toc, '%.2f') ' sekundi']);
disp(' -------------------- ');

function mat = Rz(alfa)
% --------------------------------------------------
% sestavi rotacni matici 3x3 kolem osy z o uhel alfa
% IN:  alfa ... uhel rotace v [rad]
% OUT: mat  ... rotacni matice Rz 3x3
% --------------------------------------------------

if nargin ~= 1
	error('Incorrect number of input arguments.');
end

mat = [cos(alfa) sin(alfa) 0;
      -sin(alfa) cos(alfa) 0;
          0         0      1];

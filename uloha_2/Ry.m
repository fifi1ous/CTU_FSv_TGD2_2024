function mat = Ry(alfa)
% --------------------------------------------------
% sestavi rotacni matici 3x3 kolem osy y o uhel alfa
% IN:  alfa ... uhel rotace v [rad]
% OUT: mat  ... rotacni matice Ry 3x3
% --------------------------------------------------

if nargin ~= 1
	error('Incorrect number of input arguments.');
end

mat = [cos(alfa) 0 -sin(alfa);
          0      1     0;
       sin(alfa) 0  cos(alfa)];

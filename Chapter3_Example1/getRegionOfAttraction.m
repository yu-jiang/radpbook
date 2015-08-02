function D = getRegionOfAttraction(p)
%%

% The set for approximation is set as
%\Omega = {x1,x2,x3,x4 | x1\in(-0.5, 0.5),
%                        x2\in(-5,5),
%                        x3\in(-0.2,0.2),
%                        x4\in(-10,10)}

% We need to find D, by solving
%
% max D
% st. V(x)\le D => x\in\Omega
%
% One necessary condition is that if x in boundary of \Omega
% => V(x) >= D


% Boundary checking
D = 0;
xs = 1;
dx1 = 0.1;
dx2 = 1;
dx3 = 0.05;
dx4 = 2;

while all(xs) == 1
	D = D + 0.01;
	% 1. Check the boundary of x1
	for x1 = [-0.5 0.5]
		for x2 = -5:dx2:5
			for x3 = -0.2:dx3:0.2
				for x4 = -10:dx4:10
					% xs = [xs;[x1 x2 x3 x4]];
					xs = [xs; p'*Phi_fun([x1 x2 x3 x4])'>= D];
				end
			end
		end
	end
	
	% 2. Check the boundary of x2
	for x1 = -0.5:dx1:0.5
		for x2 = [-5 5]
			for x3 = -0.2:dx3:0.2
				for x4 = -10:dx4:10
					% xs = [xs;[x1 x2 x3 x4]];
					xs = [xs; p'*Phi_fun([x1 x2 x3 x4])'>= D];
				end
			end
		end
	end
	
	% 3. Check the boundary of x3
	for x1 = -0.5:dx1:0.5
		for x2 = -5:dx2:5
			for x3 = [-0.2 0.2]
				for x4 = -10:dx4:10
					% xs = [xs;[x1 x2 x3 x4]];
					xs = [xs; p'*Phi_fun([x1 x2 x3 x4])'>= D];
				end
			end
		end
	end
	
	% 4. Check the boundary of x4
	for x1 = -0.5:dx1:0.5
		for x2 = -5:dx2:5
			for x3 = -0.2:dx3:0.2
				for x4 = [-10 10]
					% xs = [xs;[x1 x2 x3 x4]];
					xs = [xs; p'*Phi_fun([x1 x2 x3 x4])'>= D];
				end
			end
		end
	end
end

D = D - 0.01;
end
%%


classdef (Sealed) paramMgr < handle
	% Singleton class to store parameters
	
	properties
		D1 = 1;
		T1 = 5;
 		B12 = 0.06;
 		B11 = 0.05;
 		B22 = 0.01;
 		w0 = 100*pi;
 		H1 = 3;
 		angle10 = 1;
 		angle20 = 1.2;
		
 		D2 = 2;
 		T2 = 5;
 		H2 = 3;
 		Ef1 = 2;
 		Ef2 = 3;
		
 		K1 = [-1 -2 0];
		KM = zeros(1,3);
    	Q = 300*eye(3);
 		x0 = [0.1;0;0];
		z0 = [0.5;-1;0];
		
		K2 = [-10 -10 -10];
		
		A1
		B1
		A2
		B2				
		
		% Initial gain
		K0 = [2.2361    0.2891];
		
		% For exploration noise
		w = [-1.8902    0.2456    0.4173 -0.6776   -0.6973   -0.6715  0.1009   -0.4142   -1.2529 -1.2968];
		oscTstart = 2;
		connectTstart = 5;
		
		% Learning
		T = 0.1;
	end
	
	methods (Access = private)
		function this = paramMgr()
			% disp('Creating a new parameter manager');
			Initialize(this)
	    end
	end
	
	methods (Static)
		function this = getInstance()
			persistent localObj
			if isempty(localObj) || isvalid(localObj)
				localObj = paramMgr;
			end
			this = localObj;
		end
	end	
end

function Initialize(this)
this.A1 = [0     1         0;
	       0  -this.D1/2/this.H1  this.w0/2/this.H1;
	       0    0         -1/this.T1];
this.B1 = [0;	0;	1/this.T1];
% this.K1 = lqr(this.A1,this.B1,this.Q,1);
this.A2 = [0     1         0;
		   0  -this.D2/2/this.H2     this.w0/2/this.H2;
		   0    0         -1/this.T2];
	   
this.B2 = [0;   0;	1/this.T2];
%this.K2 = lqr(this.A2,this.B2,0.5*[1  0 0;0 1 0;0 0 1],10);

end
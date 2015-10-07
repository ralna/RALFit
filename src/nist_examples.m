% Find derivates of the NIST test examples
% (requires the symbolic toolbox)

syms v1 v2 v3 v4 v5 v6 v7 real
syms x x1 x2 pi 

v = [v1 v2 v3 v4 v5 v6 v7];

% generate the jacobians and hessians....
e{1} = (v1 * x).^v2; var{1} = v(1:2);
e{2} = v1 * exp(-v2*x); var{2} = v(1:2);
e{3} = v1/(1-(1+v2*x/2)^2); var{3} = v(1:2);
e{4} = v2/(1-(1+2*v2*x)^(.5));  var{4} = v(1:2);
e{5} = v1*v2*x*((1+v2*x)^(-1)); var{5} = v(1:2);
e{6} = v1*x1*exp(-v2*x2); var{6} = v(1:2); 
e{7} = atan(v1/(x-v2))/pi; var{7} = v(1:2);
e{8} = v2*cos( 2*pi*x/v1 ); var{8} = v(1:2);
e{9} = v2*sin( 2*pi*x/v1 ); var{9} = v(1:2);
e{10} = v1*(x^2+x*v2) / (x^2+x*v3+v4); var{10} = v(1:4);
e{11} = v1/(1+exp(v2-v3*x));  var{11} = v(1:3);
e{12} = v1 * exp(v2/(x+v3)); var{12} = v(1:3);
e{13} = (v1/v2) * exp(-0.5*((x-v3)/v2)^2);  var{13} = v(1:3);
e{14} = v1/((1+exp(v2-v3*x))^(1/v4));  var{14} = v(1:4);
e{15} = v1*(v2+x)^(-1/v3); var{15} = v(1:3);
e{16} = exp(-v1*x)/(v2+v3*x); var{16} = v(1:3);
e{17} = v1*exp( -(x-v2)^2 / v3^2 ); var{17} = v(1:3);
e{18} = (v1 + v2*x + v3*x^2) / (1 + v4*x + v5*x^2); var{18} = v(1:5);
e{19} = (v1 + v2*x + v3*x^2 + v4*x^3) / (1 + v5*x + v6*x^2 + v7*x^3);
var{19} = v(1:7);

for i = 1:length(e)
    grad{i} = jacobian(e{i},var{i}).';
    hess{i} = jacobian(grad{i},var{i});
end

% output in fortran
fileID = fopen('nist_derivates.txt','w');
for i = 1:length(e) % loop over the problems
    e_string = char(e{i});
    e_string = strrep(e_string,'^','**');
    fprintf(fileID,'e%i = %s \n',i, e_string);   
    for j = 1:length(grad{i})
        jac_string = char(grad{i}(j));
        jac_string = strrep(jac_string,'^','**');
        fprintf(fileID,'j%i(%i) = %s \n',i,j,jac_string); 
    end
    
    for j = 1:length(grad{i})
        for k = 1:length(grad{i})
            hess_string = char(hess{i}(j,k));
            hess_string = strrep(hess_string,'^','**');
            fprintf(fileID,'h%i(%i,%i) = %s \n',i,j,k,hess_string); 
        end
    end
    
    fprintf(fileID,'\n');   
    
end
fclose(fileID);
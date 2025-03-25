%SYSTEM PARAMETERS
rho=15/2.1; g_u=0.96; g_v=0.45; bU=1*2.1;% the product of b and u^U   
kappa=1.0003; delta=6.5*10^-6; L=1.001; mu=1.01; eps=0.1;
%INEQUALITIES INITIALIZATION AND SOLUTION
Gamma=[0 0 0 1-mu 1; 0 1-g_u g_u L*g_v*rho L*g_v*rho;
-(1-mu*kappa*g_v)*rho/kappa 1-g_u g_u L*g_v*rho L*g_v*rho;
 (1-g_u)*bU 0 0 0 0];
Theta=[0; 1-(1+delta)*g_u;-1+(1-delta)*g_u; eps];
e1=[1 0 0 0]; e2=[0 1 0 0]; e3=[0 0 1 0]; e4=[0 0 0 1];
setlmis([]); Psi=lmivar(2,[5 1]);
lmiterm([-1 1 1 Psi],.5*[1,0,0,0,0],1,'s');
lmiterm([-2 1 1 Psi],.5*[0,1,0,0,0],1,'s');
lmiterm([-3 1 1 Psi],.5*[0,0,1,0,0],1,'s');
lmiterm([-4 1 1 Psi],.5*[0,0,0,1,0],1,'s');
lmiterm([-5 1 1 Psi],.5*[0,0,0,0,1],1,'s');
lmiterm([6 1 1 Psi],.5*e1*Gamma,1,'s'); lmiterm([7 1 1 Psi],.5*e2*Gamma,1,'s');
lmiterm([8 1 1 Psi],.5*e3*Gamma,1,'s'); lmiterm([9 1 1 Psi],.5*e4*Gamma,1,'s');
lmiterm([-6 1 1 0],e1*Theta); lmiterm([-7 1 1 0],e2*Theta);
lmiterm([-8 1 1 0],e3*Theta); lmiterm([-9 1 1 0],e4*Theta);
LM1=getlmis;[tmin,xfeas] = feasp(LM1, [0 1000 1e22 0 0]); 
%if tmin < 0, then the inequalities (4.9) is feasible
PS1 = [1 0 0 0 0]*dec2mat(LM1,xfeas,Psi); 
PS2 = [0 1 0 0 0]*dec2mat(LM1,xfeas,Psi);
PS3 = [0 0 1 0 0]*dec2mat(LM1,xfeas,Psi);
PS4 = [0 0 0 1 0]*dec2mat(LM1,xfeas,Psi);
PS5 = [0 0 0 0 1]*dec2mat(LM1,xfeas,Psi);
K=PS1^-1; lambda1=PS2; lambda2=PS3; l=PS4^-1; mu=mu; tau=PS5;
clear classes
clear;
clc;
close all;
N=4;
n=2;
d=0.2;
yita=0.002;
alpha=0.4;
sita=0.5;
w=0.8;
l3=1.2;
l4=1.8;
A=[-2,0.5;1,-0.7];
D=[-1;0];
C=[1,-0.7];
L=[1,-1,0,0;-1,2,-1,0;0,-1,2,-1;0,0,-1,1];
IN=eye(4);
In=eye(3);
G=sdpvar(2,1,'full');
H=sdpvar(2,2);
H3=kron(IN,H);
H4=kron(IN,H);
P=sdpvar(n,n,'symmetric');
P1=kron(IN,P);

X11=alpha*P1+kron(IN,H'*A-G*C)+(kron(IN,H'*A-G*C))';
X12=P1-H3'-kron(IN,H'*A-G*C)'*l4;
X13=-kron(IN,G);
X14=l3*H3'*kron(IN,D);

X22=-l4*H4-l4*H4';
X23=-l4*kron(IN,G);
X24=H4'*kron(IN,D)*l4;

X33=(-alpha/(2*sita^2))*eye(N);
X34=zeros(N,N);

X44=(-alpha/(2*w^2))*eye(N);

X=[X11,X12,X13,X14;
    X12',X22,X23,X24;
    X13',X23',X33,X34;
    X14',X24',X34',X44];
% Y=[-0.1*eye(2),C*H-HB*C;
%     H'*C'-C'*HB',-eye(2)];
Lmi=[P>=0,X<=0];%,Y<=0
solvesdp(Lmi);
pres = checkset(Lmi)
if sum(pres>=0)==size(pres,1)
     disp('Positive solution.')
else
    disp('No solution.');
end
G1=inv(double(H'))*double(G)
P1=double(P);

   
  

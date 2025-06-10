clear classes
clear;
clc;
close all;
N=4;
n=2;
h=0.2;
yita=0.02;
alpha=0.4;
sita=0.4;
w=0.8;
l1=1;
l2=0.2;
A=[-1.5,0.5;0.2,-1.2];
B=[-0.5;-1];
D=[-1;0];
C=[1,0;
    0,1];
G1=[    -0.0265,    0.0194;
    0.0057 ,  -0.0052];
L=[1,-1,0,0;-1,2,-1,0;0,-1,2,-1;0,0,-1,1];

E=[1,-1,0,0;
    1,0,-1,0;
    1,0,0,-1];

F=[0,0,0;
    -1,0,0;
    0,-1,0;
    0,0,-1];
IN=eye(4);
In=eye(3);
K=sdpvar(1,2,'full');
H=sdpvar(2,2);
H1=kron(In,H);
H2=H1;
Q=sdpvar(n,n,'symmetric');
Q1=kron(In,Q);
R=sdpvar(n,n,'symmetric');
R1=kron(In,R);
O=sdpvar(2,2,'symmetric');
O1=kron(In,O);
z=sdpvar(n,n,'symmetric');
Z=kron(In,z);
x=exp(alpha*h)-1;

X11=alpha*Q1+R1+kron(In,A)*H1*l1+l1*H1'*kron(In,A')-(alpha/x)*Z+l1^2*kron(E*E',G1*(C*C')*G1');
X12=Q1-l1*H1+l2*kron(In,H'*A');
X13=(alpha/x)*Z+l1'*kron(E*L*F,B*K);
X14=zeros((N-1)*n,(N-1)*n);
X15=l1'*kron(E*L*F,B*K);
X16=kron(E,G1*C)*l1;

X22=-l2*H2-l2*H2'+h*Z;
X23=l2'*kron(E*L*F,B*K);
X24=zeros((N-1)*n,(N-1)*n);
X25=l2'*kron(E*L*F,B*K);
X26=l2'*kron(E,G1*C);

X33=yita*kron(F'*(L'*L)*F,O)-2*(alpha/x)*Z;
X34=(alpha/x)*Z;
X35=-yita*kron(F'*(L'*L)*F,O);
X36=zeros((N-1)*n,(N)*n);

X44=-exp(-alpha*h)*R1-(alpha/x)*Z;
X45=zeros((N-1)*n,(N-1)*n);
X46=zeros((N-1)*n,(N)*n);

X55=yita*kron(F'*L'*IN*L*F,O)-kron(In,O);
X56=zeros((N-1)*n,(N)*n);

X66=(-alpha/(sita^2))*kron(IN,eye(n));

X=[X11,X12,X13,X14,X15,X16;
    X12',X22,X23,X24,X25,X26;
    X13',X23',X33,X34,X35,X36;
    X14',X24',X34',X44,X45,X46;
    X15',X25',X35',X45',X55,X56;
    X16',X26',X36',X46',X56',X66;];
Lmi=[Q>=0,R>=0,O>=0,z>=0,X<=0];
solvesdp(Lmi);
pres = checkset(Lmi)
if sum(pres>=0)==size(pres,1)
     disp('Positive solution.')
else
    disp('No solution.');
end
K1=double(K)*inv(double(H));
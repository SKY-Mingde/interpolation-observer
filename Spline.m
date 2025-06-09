%三次样条插值的matlab实现
function [x2,y2]=Spline_3(x,y,String,dx)  
n=length(x);
if n<=1
    disp('数组长度过短!')
    return
end
if n~=len(y)
    disp('两个数组长度必须相同！')
    return
end
 
if nargin==2
    dx=zeros(1,2);
    dx(1)=(y(2)-y(1))/(x(2)-x(1));
    dx(2)=(y(n)-y(n-1))/(x(n)-x(n-1));
else
    if nargin ==4
        if String == 'dx'
           
        end
    end
end
%========预分配内存空间===============================
h=zeros(n-1,1);  %步长
delta_y=zeros(n-1,1);   %增量
%解对应的系数
A=zeros(n-1,1);
B=zeros(n-1,1);
C=zeros(n-1,1);
D=zeros(n-1,1);
%G和f为需要解的n*n矩阵及其数值
G=zeros(n,n);  %是 n*n 矩阵,右端有n个二阶导数值
f=zeros(n,1);  % f是右边的矩阵
%=========================
 
for i=1:1:n-1
    h(i)=x(i+1)-x(i);
    delta_y(i)=y(i+1)-y(i); %每次y的增量
end
%   A是y(i)而B是导数
 

    G(1,1)=2*h(1);
    G(1,2)=h(1);
    G(n,n-1)=h(n-1); 
    G(n,n)=2*h(n-1);
    a=dx(1);
    b=dx(2);
%=============构造矩阵=====================
    %填写G矩阵中其余的参数()
    for i=2:1:n-1  %从2行到n-1行的数据
        G(i,i-1)=h(i-1);
        G(i,i)=2*(h(i)+h(i-1));
        G(i,i+1)=h(i);
    end
    
    f(1)=6*( delta_y(1)/h(1) - a  ) ;    
    f(n)=6*(  b-delta_y(n-1)/h(n-1)  );
    for i = 2:1:n-1
        f(i)=6*(delta_y(i)/h(i)-delta_y(i-1)/h(i));
    end
    %m为二阶导数值
    m=G\f;   
    %是G\f（实际上是inv(G)*f  ）解出m 的值）
    m=m';
    
%=====================================================
    %计算每段样条插值的系数：
    for i=1:1:n-1
        A(i)=y(i);
        B(i)=  (y(i+1)-y(i))/h(i) -h(i)*m(i)/2 - h(i)/6*(m(i+1)-m(i));
        C(i)=m(i)/2;
        D(i)=(m(i+1)-m(i))/(6*h(i));
    end
    
%*******注意:从这里设置插值后绘制函数的精度*****可以有不同的要求,

    x_2=min(x):0.05:max(x);  
    x_2=x_2';   
    y2=zeros(len(x_2),1);
    j=1;
    hold on
        
    for i=1:1:len(x_2)
        while x_2(i)>x(j+1)  
            j=j+1;
        end
        y2(i)=A(j)   +B(j)*(x_2(i)-x(j))  +C(j)*(x_2(i)-x(j))^2  +  D(j)*(x_2(i)-x(j))^3;
    end
 
x2=x_2;
 
if nargout==0
     plot(x2,y2)
end
 
end
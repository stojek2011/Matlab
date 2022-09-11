function n = ADTK51
global  t eps w
eps=0.08;
a=1;
w=sqrt(2^2-eps^2)
k=100;  
h=pi/(k*w);
i=1;
l=1;
dk=2*pi*eps/w;
T=100;
%время
t=[0:h:T*pi/w];
n=length(t)
A=0.5*a*(1+exp(-dk/2))/(-1+exp(-dk/2))
x10=0.001;
x20=A+a;
figure(1)
i=2;

while i<n
hold on; grid on
xlabel('x');
ylabel('dx');

title('Фазовый портрет системы');
  d=xx1(i);
   if(x10>0)
   if(((i)==l*k))
     if(((mod(i/k,2)==0)))
      x10=(x20*x2(i+1)+a)/(x1(i+1));
      l=l+1;
    end
    if(((mod(i/k,2)~=0)))
       x20=(x10*x1(i+1)-a)/x2(i+1)  ;
       l=l+1;
       end 
   end 
   end
  if(x10<0)
  if(((i)==l*k))
     if(((mod(i/k,2)~=0)))
     x20=(x10*x1(i+1)-a)/x2(i+1)  ;
       l=l+1;
     
    end
    if(((mod(i/k,2)==0)))
      x10=(x20*x2(i-1)-a)/(x1(i+1));
      l=l+1;
      end 
  end 
  end
  if d <=0
 line([x10*x1(i) x10*x1(i+1)],[x10*xx1(i) x10*xx1(i+1)])
end
 if d >0
  line([ (x20*x2(i)+a) (x20*x2(i+1)+a)],[x20*xx2(i) x20*xx2(i+1)]);
 end
drawnow;
i=i+1;
end
end
function xi = x1(i)
global t eps w


x1=exp(-eps*t(i))*(cos(w*t(i))+(eps/w)*sin(w*t(i)));
xi=x1;

end
function xi = x2(i)
global t eps w


x2=exp(-eps*t(i))*(cos((w*t(i)))+(eps/w)*sin(w*t(i)));
xi=x2;

end
function xi = xx1(i)

global t eps w


x1=exp(-eps*t(i))*(-w*sin(w*t(i))+(eps)*cos(w*t(i)))-eps*exp(-eps*t(i))*(cos(w*t(i))+(eps/w)*sin(w*t(i)));
xi=x1;

end
function xi = xx2(i)
global  t eps w
xx2=exp(-eps*t(i))*(-w*sin(w*t(i))+(eps)*cos(w*t(i)))-eps*exp(-eps*t(i))*(cos(w*t(i))+(eps/w)*sin(w*t(i)));
xi=xx2;

end

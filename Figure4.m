%This code provides a simple prototypical 
%implmenetation of robustness test with
%integral quadratic constraints.
%
%The theory is exposed in the paper
%C.W. Scherer
%   Dissipativity and integral quadratic constraints: 
%   Tailored computational robustness tests for 
%   complex interconnections
%IEEE Control Systems Magazine 42 (3), 115-139
%
%This paper is also available on arXiv under https://doi.org/10.48550/arXiv.2105.07401
%All references in the code are related to this paper.
%
%It calls robinv.m and requires the following toolboxes to run:
%
%Control System Toolbox
%Robust Control Toolbox
%Yalmip
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specifically, this file reproduces the results of
%the example in Section 
%
%"The Benefit of Dynamic Integral Quadratic Constraints: An Example"
%
%as depicted in Figure 4.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Given system
G0=ss([-3 -2;1 0],[1;0],[0 -1],0);

%Poles a for filters in Zames-Falb multipliers
pole1=10;
pole2=100;

%grid for parameter alpha
alv=linspace(5,50,200);

ov=[];ind=0;
for al=alv;
    ind=ind+1;
    p.sys=[1;1]*G0*al*[1 1];    
    p.P0=[0 1;1 -2];
    
    %Nominal invariance
    p.type='nom';
    p.psi=[];
    s=robinv(p);
    ov(1,ind)=s.ov;

    %Static multiplier only
    p.type='zf';
    s=robinv(p);
    ov(2,ind)=s.ov;
    
    %Dynamic multiplier with pole1
    p.psi=ss(-pole1,pole1,-1,1);
    s=robinv(p);
    ov(3,ind)=s.ov;
    
    %Dynamic multiplier with pole2
    p.psi=ss(-pole2,pole2,-1,1);
    s=robinv(p);
    ov(4,ind)=s.ov
end;

figure(1);clf
ymax=20;
plot(alv,ov(1,:),'k',alv,ov(2,:),'b',alv,ov(3,:),'r',alv,ov(4,:),'g');grid on;
p0=alv(min(find(ov(2,:)==Inf)));
if ~isempty(p0);    
    h=line([p0;p0],[0;ymax],'Color','b','LineStyle',':');
end;
p1=alv(min(find(ov(3,:)==Inf)));
if ~isempty(p1);
    h=line([p1;p1],[0;ymax],'Color','r','LineStyle',':');
end;
p2=alv(min(find(ov(4,:)==Inf)));
if ~isempty(p2);
    h=line([p2;p2],[0;ymax],'Color','g','LineStyle',':');
end;

xlabel('Parameter $\alpha$','interpreter','latex');
ylabel('$\sqrt{{\rm trace}(Y)}$','interpreter','latex')
a=axis;a(4)=ymax;axis(a)

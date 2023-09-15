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
%Example 23 as depicted in Figure 5 with the following correction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Case 1 as described in Example 23 of IEEE Control Systems Magazine version
%
%The system in Fig. 3(a) with performance output e=Gz reads
%z=Gw+d, w=delta z, e=Gz.
%The uncertainty delta is contained the set 
%as in Theorem 21 with P0=[1 0;0 -1].
%
%The interconnection translates into 
%
%z=Gw+d, w=delta z
%e=G*(Gw+d)=G^2w+Gd
%
%The computed results resemble those in Figure 5 and can be 
%generated with this code as seen below. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Case 2 as described in Example 23 of latest paper version on arXiv.
%
%The exact results in Figure 5 are generated for the
%interconnction in Figure 1 with K=1 and the performance
%output Ge instead of e. 
%
%This results in the description 
%z=G(d+w), w=delta z, e=Gz 
%
%This reads as
% 
%z=Gw+Gd, w=delta z
%e=G*(Gw+Gd)
%
%The computed results are eactly those in Figure 5.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%System in Example 23 for al=1
G=ss([-1 1;0 -1],[0;1],[-2 1],0);

%Grid for parameter al of system in Example 23
alv=linspace(.9,1,200);
ov=[];

v=[];ind=0;
for al=alv;
    ind=ind+1;
    Gal=G*al;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Case 1: Configuration in IEEE Control Systems Magazine version
    %p.sys=[1;Gal]*(Gal*[1 0]+[0 1]);

    %Case 2: Configuration in arXiv version 
    p.sys=[1;Gal]*Gal*[1 1];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Matrix in (58) to describe uncertainty constraint 
    p.P0=[1 0;0 -1];

    %Nominal invariance
    p.type='nom';
    p.psi=ss([],[],[],[]);   
    s=robinv(p);
    ov(1,ind)=s.ov;

    %Static multiplier
    p.type='D';
    p.psi=ss([],[],[],1);
    s=robinv(p);
    ov(2,ind)=s.ov;
    
    %Dynamic multiplier 
    s=zpk('s');
    p.psi=ss([1;(s+.9)^2/(s+1)^2]);
    s=robinv(p);
    ov(3,ind)=s.ov;        
    
    %Dynamic multiplier implemented with conservatism as in [26]   
    s=robinv_Ref26(p);
    ov(4,ind)=s.ov      
end;
%%
figure(1);clf
ymax=10;
plot(alv,ov(2,:),'b',alv,ov(3,:),'r',alv,ov(4,:),'k');grid on;
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

function s=robinv_Ref26(p);
%function s=robinv_Ref26(p);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code is a slight modification of robinv.m.
%
%It modifies the implementation of dynamic D-scalings as described
%in Example 23 of the paper
% 
%C.W. Scherer
%Dissipativity and integral quadratic constraints: Tailored 
%computational robustness tests for complex interconnections
%IEEE Control Systems Magazine 42 (3), 115-139
%
%This conservative implementation is proposed in [26]. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For a description about how to apply the code we refer to robinv.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set default type to p.type='nom'
if ~isfield(p,'type')
    p.type='nom';
else
    if ~strcmp(p.type,'zf') & ~strcmp(p.type,'D')
        p.type='nom';
    end;
end;

%Dimensions of channels
[k,m]=size(p.sys);
nz=1;ne=k-nz;
nw=1;nd=m-nw;
%kp is number of filters
if ~strcmp(p.type,'nom');
    kp=size(p.psi,1);
end;

%Build system from [w;d] to [z;w] and to [e;d] as in (40)
%Output dimensions: nzO=1, nw=1, ne, nd
sys=[p.sys;eye(m)];  %system from [w;d] to [z;e;w;d]

%Permutation gives system from inuts [w;d] to combined outputs [z;w;e;d] 
sys=sys([1:nz k+(1:nz) nz+(1:ne) end-nd+1:end],:); 
n=size(sys.a,1);

%Set up LMIs
lmi=[];
%%% Nominal Case - No Uncertainties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p.type,'nom');
    Psi=ss([],[],[],[]); %identitiy filter in (41)
    P=[];                %empty filter in (41)
    Ze=zeros(n);         %no terminal cost matrix
    %extract system from d to [e;d] (ignore input w, outputs z,w)
    sys=sys(nz+nw+(1:(nd+ne)),nw+1:end); 
end;
%%% Build Zames Falb multiplier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p.type,'zf');

    %Define positive real matrix as in (53)
    Ppr=[0 1;1 0];   
    
    %Define conic combination coefficients as in (53)
    %Allow for more than one filter if compared to (53)
    la=sdpvar(1,kp+1);

    %Compose Psi and P as in (53) (with different ordering)
    %
    %Sector bound multiplier (static) is always included
    Psi=ss([],[],[],[blkdiag(eye(nz),eye(nz))]);
    P=la(1)*p.P0;    
    lmi=lmi+[la(1)>=0];
    
    %Adjoin Zames-Falb multiliers according to Lemma 18
    for j=1:kp;
        Psi=[Psi;blkdiag(p.psi(j,:)*eye(nz),eye(nz))];
        P=blkdiag(P,la(j+1)*Ppr);
        lmi=lmi+[la(j+1)>=0];        
    end
    %Terminal cost is zero
    Ze=zeros(size(Psi.a,1)+n);
end;

%%% Build D multiplier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p.type,'D');

    %Define Psi=diag(psi,psi) as required in Theorem 21
    Psi=blkdiag(p.psi,p.psi);

%%%%%%%%%%%%%%%% MODIFICATION of robinv.m starts here %%%%%%%%%%%%%%%%%%%%%    
    %
    %Impose constraint on M
    %
    %Extract state-space data
    [Ap,Bp,Cp,Dp]=ssdata(p.psi);
    np=size(Ap,1);
    %mp=size(Dp,2);
    %Define certificate for FDI
    %Xp=sdpvar(np,np);
    
    %Define general middle matrix M
    M=sdpvar(kp,kp);

    %Adjoin KYP LMI (59)
    %lmi=lmi+[[Ap'*Xp+Xp*Ap Xp*Bp;Bp'*Xp zeros(mp)]-[Cp Dp]'*M*[Cp Dp]<=0];

    %Multiplier parameter P as in Corollary 22
    P=kron(p.P0,M);
    
    %Terminal cost matrix Z as in Corollary 22
    Z=kron(p.P0,zeros(np));
    
    %Fill up with zeros to match with dimension of X in (44)
    Ze=blkdiag(Z,zeros(n));
%%%%%%%%%%%%%%%% MODIFICATION of robinv.m stops here %%%%%%%%%%%%%%%%%%%%%%   

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set up common ingredients of test for all types p.type

%Define Psi-filtered system as in (40) and required for (41)
sysf=blkdiag(Psi,eye(ne+nd))*sys;

%System sysf has output v (see (35)) and e, d with dimensions nv, ne, nd
nv=size(Psi,1);

%Extract information about filtered system
[Af,Bf,Cf,Df]=ssdata(sysf);
nf=size(Af,1);
[kf,mf]=size(sysf);

%Output matrix related to output e
Cfe=Cf(nv+(1:ne),:);

%Cannot allow a direct feedthrough term for invariance
%(D_e and D_ed in (38) must vanish)
if norm(Df(nv+(1:ne),:))>0;error('[De Ded] must vanish.');end;

%Define special performance index matrix such that
%(40) reads as (43)
Pp=[zeros(ne) zeros(ne,nd);zeros(nd,ne) -eye(nd)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up main LMI (41) 
%
%KYP certificate denoted as Xf
Xf=sdpvar(nf,nf);

%Introduce variable Y for robust invariance as in (44)
Y=sdpvar(ne,ne);

%Combine IQC matrix P with performance index Pp to deifne (41)
Pe=blkdiag(P,Pp);

%Adjoin main LMI (41) to exisitign LMI system 
lmi=lmi+[[Af'*Xf+Xf*Af Xf*Bf;Bf'*Xf zeros(mf)]+[Cf Df]'*Pe*[Cf Df]<=0];

%Adjoin the invariance LMI (44) 
lmi=lmi+[[Y Cfe;Cfe' Xf+Ze]>=0];

%Minimie the trace of Y under the LMI constraints
%Call the LMI solver in Robust Control Toolbox from Yalmip
%A direct implementation is more efficient
h=solvesdp(lmi,trace(Y),sdpsettings('solver','lmilab','verbose',0));
if h.problem>0;
    s.ov=inf;
else
    %Return information about validity of LMIs
    s.lmi=checkset(lmi);
    %Return optimal value
    s.ov=sqrt(trace(double(Y)));    
end
end
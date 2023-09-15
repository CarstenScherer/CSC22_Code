function s=robinv(p);
%function s=robinv(p);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code provides a simple prototypical 
%implmenetation of the test for robust invariance with
%integral quadratic constraints.
%
%The theory is exposed in the paper
%
%C.W. Scherer
%Dissipativity and integral quadratic constraints: Tailored 
%computational robustness tests for complex interconnections
%IEEE Control Systems Magazine 42 (3), 115-139
%
%This paper is also available on arXiv under 
%https://doi.org/10.48550/arXiv.2105.07401
%All references in the code are related to this paper.
%
%Specifically, the code is tailored towards 
%generating the results as depicted in Figs. 4 and 5.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Required inputs are collected in structure p with following fields: 
%
%p.sys:  Stable ss system with input [w;d], output [z;e] as in (38).
%        In this implementation, w and z are assumed to be of
%        dimension 1.
%
%Two types of uncertainties are covered and encoded as follows:
%
%p.type: 'zf' (Zames-Falb) or 'D' (Dynamic D-scalings)
%p.psi: ss column vector with stable SISO filters as follows:
%       p.type='nom': Nominal problem - no uncertainties 
%       p.type='zf':  Zames Falb filters h as in Section 
%                     "The Benefit of Dynamic Integral 
%                      Quadratic Constraints: An Example"
%       p.type='D':   General filters psi as in Theorem 21
%
%p.P0:  2x2 index matrix as follows:
%       p.type='zf':  P_{m,L} for sector bound matrix as in (8) 
%       p.type='D':   P_0 for region of uncertainty frequency 
%                     response in (58)
%
%Output information is encded in structure s with following fields:
%
%s.ov:   Square root of minimial trace of Y
%s.lmi:  Checkset information from Yalmip about overall LMI system 

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

    %Impose constraint on M as in Theorem 21 with KYP
    %
    %Extract state-space data
    [Ap,Bp,Cp,Dp]=ssdata(p.psi);
    np=size(Ap,1);
    mp=size(Dp,2);
    %Define certificate for FDI
    Xp=sdpvar(np,np);
    
    %Define general middle matrix M
    M=sdpvar(kp,kp);
    
    %Adjoin KYP LMI (59)
    lmi=lmi+[[Ap'*Xp+Xp*Ap Xp*Bp;Bp'*Xp zeros(mp)]-[Cp Dp]'*M*[Cp Dp]<=0];

    %Multiplier parameter P as in Corollary 22
    P=kron(p.P0,M);
    
    %Terminal cost matrix Z as in Corollary 22
    Z=kron(p.P0,Xp);
    
    %Fill up with zeros to match with dimension of X in (44)
    Ze=blkdiag(Z,zeros(n));
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
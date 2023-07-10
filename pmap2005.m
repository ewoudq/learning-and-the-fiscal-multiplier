function [pA,pB,Omega,Pi] = pmap2005(a,m,c,rho,GvarG,loc_z,loc_v,A,B,D,F,G,nz,nv)
%PMAP Maps the restricted Perceived Law of Motion (PLM) to the projected
%Actual Law of Motion (ALM).
%
%   The 'p-map', or 'projected T-map', maps the restricted Perceived Law of
%   Motion (PLM) to the projected Actual Law of Motion (ALM). See Guse,
%   2008, JECD.
%
%   Inputs a and c are the coefficient matrices of the PLM, _including the
%   coeffients that are restricted to be equal to zero_.

% preparation
n=nv+nz;
%loc_z=loc_z(1:nz,:);
%loc_v=loc_z(nz+1:end,:);
l=sum(loc_z(:)~=0)+sum(loc_v(:)~=0);
% Unprojected ALM

% EQ: JEDC
%Ta=A+B*a+D*(eye(size(c))+c)*a;
%Tc=B*c+D*c^2+F;

Ta=A+B*a+D*(eye(size(m))+m)*a;
Tm=B*m+D*m^2+F;
Tc=B*c+D*(m*c+c*rho)+G*rho;

%% The moment matrices
% Omega =variance matrix of the SUR information matrix
% Pi = covariance matrix cov(W,X)

% EQ: EY=(eye(size(Tc,1))-Tc)^(-1)*Ta;
%  EQ: Pay attention to the notation in WP. $$Y_t$$ on p. 8 is not the same
% as on p. 20. Here I use the notation from p. 8.

TA=[Ta; zeros(nv,1)];
TY= [Tm Tc;
    zeros(nv,nz) rho];
EY=(eye(size(TY,1))-TY)^(-1)*TA;
% method 1
% EQ: vecvarY=(eye(n*n)-kron(Tc,Tc))^(-1)*GvarG(:);
vecvarY=(eye(n*n)-kron(TY,TY))^(-1)*GvarG(:);
varY=reshape(vecvarY,[n,n]);

% method 2: numerically
% #EQ: Requires MATLAB 'Control system toolbox' (only available via Athena)
%varY2=dlyap(Tc,GvarG);

EZ=EY(1:nz);
EV=EY(nz+1:end);

%i=1;
%for j=1:2    
    %for k=1:2
        %omegajki=[];
        %pijki=[];
        %EX=[];
        %for z=1:nz
            %loci=loc_z(z,:);
            %loci=loci(loci~=0);
            %rows=(j-1)*nz+1:(j-1)*nz+nv;
            %columns=(k-1)*nz+1:(k-1)*nz+nv;
            %varYjk=varY(rows,columns);
            %pivarYjk=varYjk(loci,:);
            %varYjk=varYjk(loci,loci);
            %omegajki=blkdiag(omegajki,varYjk);
            %pijki=blkdiag(pijki,pivarYjk);
            %Ez=EZ(loci);
            %Ev=EV(loci);
            %Ex=[Ez;Ev]';
            %EX=blkdiag(EX,Ex);
        %end
        %omega(:,:,i)=omegajki;
        %pi(:,:,i)=pijki;        
        %i=i+1;
    %end
%end
%Omega=[omega(:,:,1) omega(:,:,2);
%    omega(:,:,3) omega(:,:,4)];
%Pi=[pi(:,:,1) pi(:,:,2);
%    pi(:,:,3) pi(:,:,4)];

% New

% EX
% Ew
Ew=[];
for i=1:nz
    Ewi=EZ(logical(loc_z(i,:)))';
    Ew=blkdiag(Ew,Ewi);
end

% Ev
Ev=[];
for i=1:nz
    Evi=EV(logical(loc_v(i,:)))';
    Ev=blkdiag(Ev,Evi);
end

% EX
EX=[Ew Ev];

% Omega
% $$Omega_11
varYjk=varY(1:nz,1:nz);
Omega11=[];
for i=1:nz
    % Omega_11i
    Omega11i=varYjk(logical(loc_z(i,:)),logical(loc_z(i,:)));
    Omega11=blkdiag(Omega11,Omega11i);
end
% $$Omega_12
varYjk=varY(1:nz,nz+1:end);
Omega12=[];
for i=1:nz
    % Omega_12i
    Omega12i=varYjk(logical(loc_z(i,:)),logical(loc_v(i,:)));
    Omega12=blkdiag(Omega12,Omega12i);
end
% $$Omega_21
varYjk=varY(nz+1:end,1:nz);
Omega21=[];
for i=1:nz
    % Omega_21i
    Omega21i=varYjk(logical(loc_v(i,:)),logical(loc_z(i,:)));
    Omega21=blkdiag(Omega21,Omega21i);
end
% $$Omega_22
varYjk=varY(nz+1:end,nz+1:end);
Omega22=[];
for i=1:nz
    % Omega_22i
    Omega22i=varYjk(logical(loc_v(i,:)),logical(loc_v(i,:)));
    Omega22=blkdiag(Omega22,Omega22i);
end
Omega=[Omega11 Omega12;
    Omega21 Omega22];

% Pi
% $$Pi_11
varYjk=varY(1:nz,1:nz);
Pi11=[];
for i=1:nz
    % Pi_11i
    Pi11i=varYjk(logical(loc_z(i,:)),:);
    Pi11=blkdiag(Pi11,Pi11i);
end
% $$Pi_12
varYjk=varY(1:nz,nz+1:end);
Pi12=[];
for i=1:nz
    % Pi_12i
    Pi12i=varYjk(logical(loc_z(i,:)),:);
    Pi12=blkdiag(Pi12,Pi12i);
end
% $$Pi_21
varYjk=varY(nz+1:end,1:nz);
Pi21=[];
for i=1:nz
    % Pi_21i
    Pi21i=varYjk(logical(loc_v(i,:)),:);
    Pi21=blkdiag(Pi21,Pi21i);
end
% $$Pi_22
varYjk=varY(nz+1:end,nz+1:end);
Pi22=[];
for i=1:nz
    % Pi_22i
    Pi22i=varYjk(logical(loc_v(i,:)),:);
    Pi22=blkdiag(Pi22,Pi22i);
end
Pi=[Pi11 Pi12;
    Pi21 Pi22];

%% p-map
% p(B)
%EQ JECD Tx=Tc(:);
Tmvec=Tm';
Tcvec=Tc';
Tmvec=Tmvec(:);
Tcvec=Tcvec(:);
Tx=[Tmvec;Tcvec];

pB=Omega^(-1)*Pi*Tx;

%Y=[x;pi;g;u];
% P(A)

pA=EZ-EX*pB; %pA=EY-EX'*pB;
end


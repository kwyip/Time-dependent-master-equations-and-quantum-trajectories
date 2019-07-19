%kawayip@usc.edu
%8-qubit_chain
function ame_4qubit
%warning('off','MATLAB:eigs:TooManyRequestedEigsForRealSym')
% % Define Pauli matrices
sX = [0 1; 1 0];
sY = [0 -1i; 1i 0];
sZ = [1 0; 0 -1];
unit = speye(2);
natom = 4;

neval = 2^natom;
nevaltruc = 16;

% Define tensor product of Pauli matrices
sX_1 = kron(sX, speye(2^(natom-1)));
sX_2 = kron(kron(speye(2),sX),speye(2^(natom-2)));
sX_3 = kron(kron(speye(2^2),sX),speye(2^(natom-3)));
sX_4 = kron(speye(2^(natom-1)),sX);

sZ_1 = kron(sZ, speye(2^(natom-1)));
sZ_2 = kron(kron(speye(2),sZ),speye(2^(natom-2)));
sZ_3 = kron(kron(speye(2^2),sZ),speye(2^(natom-3)));
sZ_4 = kron(speye(2^(natom-1)),sZ);

A_1 = sZ_1;
A_2 = sZ_2;
A_3 = sZ_3;
A_4 = sZ_4;

sZsZII = kron(kron(kron(sZ,sZ),unit),unit);
IsZsZI = kron(kron(kron(unit,sZ),sZ),unit);
IIsZsZ = kron(kron(kron(unit,unit),sZ),sZ);


plus = [1/sqrt(2); 1/sqrt(2)];

dlm = dlmread('DW2000_parameters.txt');
slist = dlm(:,1).';
A_s = dlm(:,2).';
B_s = dlm(:,3).';
A_sp1 = @(s)interp1(slist,A_s,s);
B_sp1 = @(s)interp1(slist,B_s,s);


%Temperature and coupling parameters
beta = (1/2.6)*1e-9;
wc = 8*pi*1e9;
gsq2pi = (1.27)*1e-4*2*pi; %g^2*pi
betainv = 2.6*1e9; %1/beta



%Pure state initialization
psi0 = sparse(kron(kron(kron(plus,plus),plus),plus));
rho0 = psi0*psi0';
%rho = rho0(:);
rho = rho0;


tf = 5e-6;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%timestep between each ode executions
dt_me = tf/1000;
tstep_me = 0:dt_me:tf;
fidelitylist_me = zeros(1, numel(tstep_me));

for index = 1:numel(tstep_me)
    index    
    Hs = -1e9.*(1-tstep_me(index).*tf.^(-1)).*(sX_1+sX_2+sX_3+sX_4) + 1e9.*tstep_me(index).*tf.^(-1).*(((-1).*sZ_1+(0.95).*sZ_2+(0.95).*sZ_3+(-1).*sZ_4) + ...
        ((-1).*sZsZII + (-1).*IsZsZI + (-1).*IIsZsZ));

    [V,D] = eig(full(Hs));
    %make sure eigenmatrix is sorted
    if ~issorted(diag(D))
        [V,D] = eig(full(Hs));
        [D,I] = sort(diag(D));
        D = diag(D);
        V = V(:, I);
        fprintf('sorted');
    end
    v = sparse(neval,nevaltruc);
    e = sparse(1,nevaltruc);
    
    for ii = 1:nevaltruc
        v(:,ii) = sparse(V(:,ii));
        e(ii) = sparse(D(ii,ii));
    end
    Hsd = v'*Hs*v;
    v0 = v(:,1);
    rhom = reshape(rho,[neval,neval]);
    rhomcb = sparse(v'*rhom*v);
    

    fidelity = rhomcb(1,1);
    fidelitylist_me(1, index) = fidelity;
    
    alindblad  = @(t, rho)lindblad(t, rho, Hsd, natom, gsq2pi, beta, betainv, wc, A_1,A_2,A_3,A_4,v,e);
    t    = [tstep_me(index), tstep_me(index) + dt_me];
    options = odeset('RelTol',1e-3,'AbsTol',1e-6);
    %[~, RHO] = ode45(alindblad, t, rhomcb(:));
    [~, RHO] = ode45(alindblad, t, rhomcb(:),options);
    % Final value of rho is initial value for next step:
    rho = RHO(end, :); 
    rho = v*reshape(rho,[nevaltruc,nevaltruc])*v';   %in computational basis
    rho = rho(:);
end


eptime = toc

figure(1)
plot(tstep_me, fidelitylist_me,'-b','LineWidth',2);
xlabel('$t$','Interpreter','latex')
xlim([0 tf])
ylabel('Ground State Population')
title(['tf: ' num2str(tf)])
print -dpdf me4qubitexample1_linear

figure(2)
plot(tstep_me./tf, fidelitylist_me,'-b','LineWidth',2); 
xlabel('$s$','Interpreter','latex')
ylabel('Ground State Population')
title(['tf: ' num2str(tf)])
print -dpdf me4qubitexample2_linear

txt1 = sprintf('me4qubitexample1_linear%d.txt',tf);
fid1 = fopen(txt1,'w');
fprintf(fid1,'%13d %8d\n',[tstep_me;fidelitylist_me]);
fclose(fid1);


txt10 = sprintf('eptime.txt');
fid10 = fopen(txt10,'w');
fprintf(fid10,'%d\n',eptime);
fclose(fid10);


%function drhodt = lindblad(~, rho, Hsd, natom, g, beta, wc,A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,v,e)
function drhodt = lindblad(~, rho, Hsd, natom, gsq2pi, beta, betainv, wc,A_1,A_2,A_3,A_4,v,e)
neval = 2^natom;
nevaltruc = 16;
rhomcb = sparse(reshape(rho, [nevaltruc,nevaltruc]));

drhodt = -1i*(Hsd*rhomcb - rhomcb*Hsd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = bsxfun(@minus, e.', e);

[sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1);%AbsTol
length(sortedOutput(sortedOutput>0));
count = 0;
for w = sortedOutput(sortedOutput>0)
    gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));
    if isnan(gamma) || isinf(gamma)
        gamma = gsq2pi*betainv;
    end
    [b, a] = ind2sub(size(X), find(abs(X - w)<= 0.1));% AbsTol
    count = count+length(b);

    Lpcomponents1 = sparse(nevaltruc,nevaltruc); 
    Lpcomponents2 = sparse(nevaltruc,nevaltruc);
    Lpcomponents3 = sparse(nevaltruc,nevaltruc);
    Lpcomponents4 = sparse(nevaltruc,nevaltruc);

    for s = 1:length(b)
      %matrixelement = v(:,i(s))'*A*v(:,j(s));
      matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));     %j<->a, i<->b in paper
      matrixelement2 = v(:,a(s))'*A_2*v(:,b(s));
      matrixelement3 = v(:,a(s))'*A_3*v(:,b(s));
      matrixelement4 = v(:,a(s))'*A_4*v(:,b(s));


      Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent2 = matrixelement2*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent3 = matrixelement3*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent4 = matrixelement4*sparse(a(s),b(s),1,nevaltruc,nevaltruc);


      Lpcomponents1 = Lpcomponents1 + Lpcomponent1;
      Lpcomponents2 = Lpcomponents2 + Lpcomponent2;
      Lpcomponents3 = Lpcomponents3 + Lpcomponent3;
      Lpcomponents4 = Lpcomponents4 + Lpcomponent4;

    end

    Lncomponents1 = Lpcomponents1';
    Lncomponents2 = Lpcomponents2';
    Lncomponents3 = Lpcomponents3';
    Lncomponents4 = Lpcomponents4';

    
    Lp1 = sqrt(gamma)*Lpcomponents1;
    Lp2 = sqrt(gamma)*Lpcomponents2;
    Lp3 = sqrt(gamma)*Lpcomponents3;
    Lp4 = sqrt(gamma)*Lpcomponents4;


    Ln1 = sqrt(gamma*exp(-beta*w))*Lncomponents1;
    Ln2 = sqrt(gamma*exp(-beta*w))*Lncomponents2;
    Ln3 = sqrt(gamma*exp(-beta*w))*Lncomponents3;
    Ln4 = sqrt(gamma*exp(-beta*w))*Lncomponents4;
  
    
    


%%%%%%%%%%Superoperator2%%%%%%%%%  
drhodt = drhodt + (Lp1*rhomcb*Lp1'-0.5*(Lp1)'*Lp1*rhomcb - 0.5*rhomcb*(Lp1)'*Lp1) + (Ln1*rhomcb*Ln1'-0.5*(Ln1)'*Ln1*rhomcb - 0.5*rhomcb*(Ln1)'*Ln1)...
    + (Lp2*rhomcb*Lp2'-0.5*(Lp2)'*Lp2*rhomcb - 0.5*rhomcb*(Lp2)'*Lp2) + (Ln2*rhomcb*Ln2'-0.5*(Ln2)'*Ln2*rhomcb - 0.5*rhomcb*(Ln2)'*Ln2) ...
    + (Lp3*rhomcb*Lp3'-0.5*(Lp3)'*Lp3*rhomcb - 0.5*rhomcb*(Lp3)'*Lp3) + (Ln3*rhomcb*Ln3'-0.5*(Ln3)'*Ln3*rhomcb - 0.5*rhomcb*(Ln3)'*Ln3)...
    + (Lp4*rhomcb*Lp4'-0.5*(Lp4)'*Lp4*rhomcb - 0.5*rhomcb*(Lp4)'*Lp4) + (Ln4*rhomcb*Ln4'-0.5*(Ln4)'*Ln4*rhomcb - 0.5*rhomcb*(Ln4)'*Ln4);
end

gamma0 = gsq2pi*betainv;
%[b0, a0] = ind2sub(size(X), find(X == 0));
[b0, a0] = ind2sub(size(X), find(abs(X - 0)<=0.1));%AbsTol

L01 = sparse(nevaltruc,nevaltruc);
L02 = sparse(nevaltruc,nevaltruc);
L03 = sparse(nevaltruc,nevaltruc);
L04 = sparse(nevaltruc,nevaltruc);

for s = 1:length(b0)
    matrixelement1 = v(:,a0(s))'*A_1*v(:,b0(s));
    matrixelement2 = v(:,a0(s))'*A_2*v(:,b0(s));
    matrixelement3 = v(:,a0(s))'*A_3*v(:,b0(s));
    matrixelement4 = v(:,a0(s))'*A_4*v(:,b0(s));

    L0component1 = matrixelement1*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component2 = matrixelement2*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component3 = matrixelement3*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component4 = matrixelement4*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L01 = L01 + L0component1;
    L02 = L02 + L0component2;
    L03 = L03 + L0component3;
    L04 = L04 + L0component4;
end

L01 = sqrt(gamma0)*L01;
L02 = sqrt(gamma0)*L02;
L03 = sqrt(gamma0)*L03;
L04 = sqrt(gamma0)*L04;



%%%%%%%%%Superoperator2%%%%%%%%%    
drhodt = drhodt + (L01*rhomcb*(L01)'-0.5*(L01)'*L01*rhomcb - 0.5*rhomcb*(L01)'*L01)...
    + (L02*rhomcb*(L02)'-0.5*(L02)'*L02*rhomcb - 0.5*rhomcb*(L02)'*L02)...
    + (L03*rhomcb*(L03)'-0.5*(L03)'*L03*rhomcb - 0.5*rhomcb*(L03)'*L03)...
    + (L04*rhomcb*(L04)'-0.5*(L04)'*L04*rhomcb - 0.5*rhomcb*(L04)'*L04);


%%%%%%%%%Superoperator2%%%%%%%%%
% drhodt = L*rhovcb;
%%%%%%%%%Superoperator2%%%%%%%%%
%whos
drhodt = drhodt(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
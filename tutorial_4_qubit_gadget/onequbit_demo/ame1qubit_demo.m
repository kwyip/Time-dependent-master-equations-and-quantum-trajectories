%kawayip@usc.edu
%1-qubit_chain
function ame1qubit_demo
%warning('off','MATLAB:eigs:TooManyRequestedEigsForRealSym')
% % Define Pauli matrices
sX = [0 1; 1 0];
sY = [0 -1i; 1i 0];
sZ = [1 0; 0 -1];
unit = speye(2);
natom = 1;

neval = 2^natom;
nevaltruc = 2;

% Define tensor product of Pauli matrices
sX_1 = sX;

sZ_1 = sZ;

A_1 = sZ_1;



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



psi0 = plus;
rho0 = psi0*psi0';
%rho = rho0(:);
rho = rho0;

e1list = [];
e2list = [];

tf = 100e-6;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%timestep between each ode executions
dt_me = tf/1000;
tstep_me = 0:dt_me:tf;
fidelitylist_me = zeros(1, numel(tstep_me));

for index = 1:numel(tstep_me)
    index    
    Hs = -1e9.*A_sp1(tstep_me(index)./tf).*(sX_1) + 1e9.*B_sp1(tstep_me(index)./tf).*((-1).*((1/4).*sZ_1));

    [V,D] = eig(full(Hs));
    %make sure eigenmatrix is sorted
    if ~issorted(diag(D))
        [V,D] = eig(full(Hs));
        [D,I] = sort(diag(D));
        D = diag(D);
        V = V(:, I);
        fprintf('sorted');
    end
    e1list = [e1list D(1,1)];
    e2list = [e2list D(2,2)];
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
    
    
    alindblad  = @(t, rho)lindblad(t, rho, Hsd, natom, gsq2pi, beta, betainv, wc, A_1,v,e);
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
plot(tstep_me, fidelitylist_me,'LineWidth',2);
xlabel('$t$','Interpreter','latex','FontSize',25)
xlim([0 tf])
ylabel('Ground State Population','FontSize',25)
title(['tf(s): ' num2str(tf)],'FontSize',9)
set(gca,'FontSize',20)
print -dpdf ame1qubit_demo_1

figure(2)
plot(tstep_me./tf, fidelitylist_me,'LineWidth',2); 
xlabel('$s$','Interpreter','latex','FontSize',25)
ylabel('Ground State Population','FontSize',25)
title(['tf(s): ' num2str(tf)],'FontSize',9)
set(gca,'FontSize',20)
print -dpdf ame1qubit_demo_2

txt1 = sprintf('me1qubitexample1_linear%d.txt',tf);
fid1 = fopen(txt1,'w');
fprintf(fid1,'%13d %8d\n',[tstep_me;fidelitylist_me]);
fclose(fid1);

figure(1)
plot(tstep_me./tf, e1list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e2list./(2*pi)/1e9, 'LineWidth',2);  
xlabel('$s$','Interpreter','latex','FontSize',25)
ylabel('GHz','FontSize',25)
set(gca,'FontSize',20)
legend('0th','1st', 'location', 'best')
title('Spectrum','Interpreter','latex')
print -dpdf 1qubit_spectrum

txt10 = sprintf('eptime.txt');
fid10 = fopen(txt10,'w');
fprintf(fid10,'%d\n',eptime);
fclose(fid10);


%function drhodt = lindblad(~, rho, Hsd, natom, g, beta, wc,A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,v,e)
function drhodt = lindblad(~, rho, Hsd, natom, gsq2pi, beta, betainv, wc,A_1,v,e)
neval = 2^natom;
nevaltruc = 2;
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

    for s = 1:length(b)
      %matrixelement = v(:,i(s))'*A*v(:,j(s));
      matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));     %j<->a, i<->b in paper



      Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);



      Lpcomponents1 = Lpcomponents1 + Lpcomponent1;


    end

    Lncomponents1 = Lpcomponents1';


    
    Lp1 = sqrt(gamma)*Lpcomponents1;



    Ln1 = sqrt(gamma*exp(-beta*w))*Lncomponents1;

  
    
    


%%%%%%%%%%Superoperator2%%%%%%%%%  
drhodt = drhodt + (Lp1*rhomcb*Lp1'-0.5*(Lp1)'*Lp1*rhomcb - 0.5*rhomcb*(Lp1)'*Lp1) + (Ln1*rhomcb*Ln1'-0.5*(Ln1)'*Ln1*rhomcb - 0.5*rhomcb*(Ln1)'*Ln1);
end

gamma0 = gsq2pi*betainv;
%[b0, a0] = ind2sub(size(X), find(X == 0));
[b0, a0] = ind2sub(size(X), find(abs(X - 0)<=0.1));%AbsTol

L01 = sparse(nevaltruc,nevaltruc);


for s = 1:length(b0)
    matrixelement1 = v(:,a0(s))'*A_1*v(:,b0(s));


    L0component1 = matrixelement1*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);

    L01 = L01 + L0component1;

end

L01 = sqrt(gamma0)*L01;




%%%%%%%%%Superoperator2%%%%%%%%%    
drhodt = drhodt + (L01*rhomcb*(L01)'-0.5*(L01)'*L01*rhomcb - 0.5*rhomcb*(L01)'*L01);


%%%%%%%%%Superoperator2%%%%%%%%%
% drhodt = L*rhovcb;
%%%%%%%%%Superoperator2%%%%%%%%%
%whos
drhodt = drhodt(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
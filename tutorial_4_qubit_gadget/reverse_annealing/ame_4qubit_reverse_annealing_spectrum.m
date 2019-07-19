%kawayip@usc.edu
%8-qubit_chain
function ame_4qubit_reverse_annealing_spectrum
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




dlm = dlmread('DW2000_parameters.txt');
% slist1 = linspace(1,0.75,400);
% slist2 = linspace(0.75,0.75,200);
% slist3 = linspace(0.75,1,400);
% slist = [slist1, slist2, slist3];

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


% %Gibbs state initialization
% Hs = -1e9.*A_sp1(0).*(sX_1+sX_2+sX_3+sX_4) + 1e9.*B_sp1(0).*((-1).*((1/4).*sZ_1) + ((-1).*sZsZII + (-1).*IsZsZI + (-1).*IIsZsZ));
% [V,D] = eig(full(Hs));
% if ~issorted(diag(D))
%     [V,D] = eig(full(Hs));
%     [D,I] = sort(diag(D));
%     D = diag(D);
%     V = V(:, I);
%     fprintf('sorted');
% end
% diag(D)
% e = zeros(1,nevaltruc);
% v = sparse(V(:,1:nevaltruc));
% for i = 1:nevaltruc
%     e(i) = sparse(D(i,i));
% end
% Z = 0;
% for i = 1:nevaltruc
%     Z = Z + exp(-beta*e(i));
% end
% 
% gibbs = 0;
% for i = 1:nevaltruc
%     gibbs = gibbs + (exp(-beta*e(i))/Z)*v(:,i)*v(:,i)';
% end
%rho = gibbs(:);
%Pure state initialization
plus = [1/sqrt(2); 1/sqrt(2)];
psi0 = sparse(kron(kron(kron(plus,plus),plus),plus));
rho0 = psi0*psi0';
%rho = rho0(:);
rho = rho0;

gaplist = [];
e1list = [];
e2list = [];
e3list = [];
e4list = [];
e5list = [];
e6list = [];
e7list = [];
e8list = [];
e9list = [];
e10list = [];
e11list = [];
e12list = [];
e13list = [];
e14list = [];
e15list = [];
e16list = [];
tf = 5e-6;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%timestep between each ode executions

step = 1000;
sstar = 0.523
inslist1 = linspace(1,sstar,step*2/5);
inslist2 = linspace(sstar,sstar,step*1/5+1);
inslist3 = linspace(sstar,1,step*2/5);
inslist = [inslist1, inslist2, inslist3];


dt_me = tf/1000;
tstep_me = 0:dt_me:tf;

fidelitylist_me = zeros(1, numel(tstep_me));

for index = 1:numel(tstep_me)
    %index
%     Hs = -1e9.*A_sp1(tstep_me(index)./tf)/2.*(sX_1+sX_2+sX_3+sX_4) + 1e9.*B_sp1(tstep_me(index)./tf)/2.*(((-1).*sZ_1+(0.95).*sZ_2+(0.95).*sZ_3+(-1).*sZ_4) + ...
%         ((-1).*sZsZII + (-1).*IsZsZI + (-1).*IIsZsZ));

    Hs = -1e9.*A_sp1(tstep_me(index)./tf)/2.*(sX_1+sX_2+sX_3+sX_4) + 1e9.*B_sp1(tstep_me(index)./tf)/2.*(((-1).*sZ_1+(0.95).*sZ_2+(0.95).*sZ_3+(-1).*sZ_4) + ...
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
    
    e1list = [e1list D(1,1)];
    e2list = [e2list D(2,2)];
    gaplist = [gaplist D(2,2)-D(1,1)];
    e3list = [e3list D(3,3)];
    e4list = [e4list D(4,4)];
    e5list = [e5list D(5,5)];
    e6list = [e6list D(6,6)];
    e7list = [e7list D(7,7)];
    e8list = [e8list D(8,8)];
    e9list = [e9list D(9,9)];
    e10list = [e10list D(10,10)];
    e11list = [e11list D(11,11)];
    e12list = [e12list D(12,12)];
    e13list = [e13list D(13,13)];
    e14list = [e14list D(14,14)]; 
    e15list = [e15list D(15,15)];
    e16list = [e16list D(16,16)];    
    
    v = sparse(neval,nevaltruc);
    e = sparse(1,nevaltruc);
    
    for ii = 1:nevaltruc
        v(:,ii) = sparse(V(:,ii));
        e(ii) = sparse(D(ii,ii));
    end
    Hsd = v'*Hs*v;
    v0 = v(:,1);
%     rhom = reshape(rho,[neval,neval]);
%     rhomcb = sparse(v'*rhom*v);
%     
% 
%     fidelity = rhomcb(1,1);
%     fidelitylist_me(1, index) = fidelity;
%     
%     alindblad  = @(t, rho)lindblad(t, rho, Hsd, natom, gsq2pi, beta, betainv, wc, A_1,A_2,A_3,A_4,v,e);
%     t    = [tstep_me(index), tstep_me(index) + dt_me];
%     options = odeset('RelTol',1e-3,'AbsTol',1e-6);
%     %[~, RHO] = ode45(alindblad, t, rhomcb(:));
%     [~, RHO] = ode45(alindblad, t, rhomcb(:),options);
%     % Final value of rho is initial value for next step:
%     rho = RHO(end, :); 
%     rho = v*reshape(rho,[nevaltruc,nevaltruc])*v';   %in computational basis
%     rho = rho(:);
end


eptime = toc

% figure(1)
% plot(tstep_me, fidelitylist_me,'-b','LineWidth',2);
% xlabel('$t$','Interpreter','latex')
% xlim([0 tf])
% ylabel('Ground State Population')
% title(['tf: ' num2str(tf)])
% 
% figure(2)
% plot(tstep_me./tf, fidelitylist_me,'-b','LineWidth',2); 
% xlabel('$s$','Interpreter','latex')
% ylabel('Ground State Population')
% title(['tf: ' num2str(tf)])
% 
% txt1 = sprintf('ame%d.txt',tf);
% fid1 = fopen(txt1,'w');
% fprintf(fid1,'%13d %8d\n',[tstep_me;fidelitylist_me]);
% fclose(fid1);

z = 1
txt1 = sprintf('e1_%f.txt', z);
fid1 = fopen(txt1,'w');
fprintf(fid1, '%.13f %.20f\n', [(tstep_me./tf).' (e1list./(2*pi)).'].');
fclose(fid1);

txt2 = sprintf('e2_%f.txt', z);
fid2 = fopen(txt2,'w');
fprintf(fid2, '%.13f %.20f\n', [(tstep_me./tf).' (e2list./(2*pi)).'].');
fclose(fid2);

txt3 = sprintf('e3_%f.txt', z);
fid3 = fopen(txt3,'w');
fprintf(fid3, '%.13f %.20f\n', [(tstep_me./tf).' (e3list./(2*pi)).'].');
fclose(fid3);

txt4 = sprintf('e4_%f.txt', z);
fid4 = fopen(txt4,'w');
fprintf(fid4, '%.13f %.20f\n', [(tstep_me./tf).' (e4list./(2*pi)).'].');
fclose(fid4);

txt5 = sprintf('e5_%f.txt', z);
fid5 = fopen(txt5,'w');
fprintf(fid5, '%.13f %.20f\n', [(tstep_me./tf).' (e5list./(2*pi)).'].');
fclose(fid5);

txt9 = sprintf('gap_%f.txt', z);
fid9 = fopen(txt9,'w');
fprintf(fid9, '%.13f %.20f\n', [(tstep_me./tf).' (gaplist./(2*pi)).'].');
fclose(fid9);

txt10 = sprintf('eptime.txt');
fid10 = fopen(txt10,'w');
fprintf(fid10,'%d\n',eptime);
fclose(fid10);

figure(1)
plot(tstep_me./tf, e1list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e2list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e3list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e4list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e5list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e6list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e7list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e8list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e9list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e10list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e11list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e12list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e13list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e14list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e15list./(2*pi)/1e9, 'LineWidth',2); 
hold on
plot(tstep_me./tf, e16list./(2*pi)/1e9, 'LineWidth',2); 
xlabel('$s$','Interpreter','latex','FontSize',25)
ylabel('GHz','FontSize',25)
set(gca,'FontSize',20)
%legend('0th','1st','2nd','3rd','4th', 'location', 'best')
title('Spectrum','Interpreter','latex')
print -dpdf ame_4qubit_reverse_annealing_spectrum

figure(2)
plot(tstep_me./tf, e1list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e2list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e3list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e4list./(2*pi)/1e9, 'LineWidth',2);  
hold on
plot(tstep_me./tf, e5list./(2*pi)/1e9, 'LineWidth',2);  
xlim([0.5 0.55])
xlabel('$s$','Interpreter','latex','FontSize',25)
ylabel('GHz','FontSize',25)
set(gca,'FontSize',20)
%legend('0th','1st','2nd','3rd','4th', 'location', 'best')
title('Spectrum','Interpreter','latex')
print -dpdf ame_4qubit_reverse_annealing_spectrum_2


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

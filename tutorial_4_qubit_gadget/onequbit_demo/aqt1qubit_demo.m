%kawayip@usc.edu
%8-qubit_chain
function aqt1qubit_demo
% cluster = parallel.cluster.Generic;
% set(cluster,'JobStorageLocation', '/home/rcf-proj2/ky/kawayip/research/fourqubitgadget_linear/qt');
% set(cluster,'HasSharedFilesystem', true);
% set(cluster,'IntegrationScriptsLocation','/usr/usc/matlab/R2018b/SlurmIntegrationScripts');
% cluster.AdditionalProperties.SlurmArgs='--time=23:59:59';

% cluster = get_LOCAL_cluster('/home/rcf-proj2/ky/kawayip/research/fourqubitgadget_linear/qt');
% 
% pool=parpool(cluster,15)

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

% 
% sZsZII = kron(kron(kron(sZ,sZ),unit),unit);
% IsZsZI = kron(kron(kron(unit,sZ),sZ),unit);
% IIsZsZ = kron(kron(kron(unit,unit),sZ),sZ);


plus = [1/sqrt(2); 1/sqrt(2)];

dlm = dlmread('DW2000_parameters.txt');
slist = dlm(:,1).';
A_s = dlm(:,2).';
B_s = dlm(:,3).';
A_sp1 = @(s)interp1(slist,A_s,s);
B_sp1 = @(s)interp1(slist,B_s,s);
% c0 = [1; 0];
% c1 = [0; 1];

%initial state 
psi0 = plus;
rho0 = psi0*psi0';
%rho = rho0(:);
rho = rho0;


%Temperature and coupling parameters
beta = (1/2.6)*1e-9;
wc = 8*pi*1e9;
gsq2pi = (1.27)*1e-4*2*pi; %g^2*pi
betainv = 2.6*1e9; %1/beta

%fidelityendlist_me = [];
dlm = dlmread('DW2000_parameters.txt');
slist = dlm(:,1).';
A_s = dlm(:,2).'; 
B_s = dlm(:,3).';
A_sp1 = @(s)interp1(slist,A_s,s);
B_sp1 = @(s)interp1(slist,B_s,s);


%tf = 2e-6;
tf = 100e-6;
hbar = 1;
tic
% Quantum trajectory by using waiting time distribution
ntraj = 5000;
dt_qt = tf/1000;
tstep_qt = 0:dt_qt:tf;
%matrix to store fidelity: traj in row, time in col
fidelitylistmat_qt = zeros(ntraj, numel(tstep_qt));



parfor n = 1:ntraj
    psi = psi0;
    unpsi = psi0;
    r = rand([1 2]);
    fidelitylist_qt = zeros(1, numel(tstep_qt));
    v = zeros(neval,nevaltruc);
    e = zeros(1,nevaltruc);
%     jlist_qt = zeros(1, 1000);
%     jorder = 1;
%     unpsinormlist_qt = zeros(1, numel(tstep_qt));
    for index = 1:numel(tstep_qt)
        %index
        %clear dp
        Hs = -1e9.*A_sp1(tstep_qt(index)./tf).*(sX_1) + 1e9.*B_sp1(tstep_qt(index)./tf).*((-1).*((1/4).*sZ_1));

        
        [V,D] = eig(full(Hs));
        if ~issorted(diag(D))
            %[V,D] = eigs(Hs,18,'sa');
            [V,D] = eig(full(Hs));
            [D,I] = sort(diag(D));
            D = diag(D);
            V = V(:, I);
            sorted = 1
        end
        for ii = 1:nevaltruc
            v(:,ii) = sparse(V(:,ii));
            e(ii) = sparse(D(ii,ii));
        end
        Hsd = v'*Hs*v;
        psicb = v'*psi;
        unpsicb = v'*unpsi;
        norm(psicb);
        %Fidelity
        fidelity = psicb(1)*psicb(1)';
        %fidelity = (v(:,1)'* psi) * (psi' * v(:,1));
        %fidelitylist = cat(2, fidelitylist, fidelity);
        fidelitylist_qt(1, index) = fidelity; 
        
        X = bsxfun(@minus, e.', e);   %matrix of \omega_{ba}
        [sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1); %AbsTol
        length(sortedOutput(sortedOutput>0));
        w_unique = length(sortedOutput);
        
        dp = zeros(1, w_unique*4);
                
        %gamma0 = 2.*g.^2.*pi.*beta.^(-1);
        gamma0 = gsq2pi*betainv;
        %[b0, a0] = ind2sub(size(X), find(X == 0));
        [b0, a0] = ind2sub(size(X), find(abs(X - 0)<=0.1));%AbsTol
        L01 = sparse(nevaltruc,nevaltruc);


        for s = 1:length(b0)
            matrixelement1 = v(:,a0(s))'*A_1*v(:,b0(s));      %j<->a, i<->b in paper


            L0component1 = matrixelement1*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);


            L01 = L01 + L0component1;

        end
        L01 = sqrt(gamma0)*L01;


        dp(1) = (psicb'*L01')*(L01*psicb)*dt_qt;



        H_eff = Hsd - (1i*hbar/2)*(L01'*L01);
        pdx = 1+1;
        
        count = 0;
        for w = sortedOutput(sortedOutput>0) 
            gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));   
            if isnan(gamma) || isinf(gamma)
                gamma = gsq2pi*betainv;
            end
            %[b, a] = ind2sub(size(X), find(X == w));
            [b, a] = ind2sub(size(X), find(abs(X - w)<= 0.1));% AbsTol
            count = count+length(b);
            Lpcomponents1 = sparse(nevaltruc,nevaltruc); 


            for s = 1:length(b)
                matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));


                Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);


                
                Lpcomponents1 = Lpcomponents1 + Lpcomponent1;


            end
            Lncomponents1 = Lpcomponents1';


                  
            Lp1 = sqrt(gamma)*Lpcomponents1;


            
            Ln1 = sqrt(gamma*exp(-beta*w))*Lncomponents1;


            
            dp(pdx) = (psicb'*Lp1')*(Lp1*psicb)*dt_qt;


            
            dp(pdx+1) = (psicb'*Ln1')*(Ln1*psicb)*dt_qt;


            
            pdx = pdx + 2;
            H_eff = H_eff - (1i*hbar/2)*(Lp1'*Lp1)...
                - (1i*hbar/2)*(Ln1'*Ln1);
        end
        length(dp);
             
        %nonzeros(dp)
        dpj = sum(dp);
        %H_eff = sparse(H_eff);
        U_eff = expm(-1i*dt_qt*H_eff/hbar);
        
        dp0 = 1 - dpj;
        
        %unpsi_prev = unpsi;
        unpsi_prev = unpsicb;
        norm2_prev = norm(unpsi_prev)^2;
        %unpsi = U_eff*unpsi;     %Evolve until the norm of unpsi becomes r1
        unpsicb = U_eff*unpsicb;
        %%%%%%%

        unpsinorm = norm(unpsicb);
        
        
        norm2_unpsi = unpsinorm^2;
        r1 = r(1);
        if unpsicb'*unpsicb > r1
            psicb = unpsicb/unpsinorm;
            %Change back to comp.basis
            unpsi = v*unpsicb;
            psi = v*psicb;
            %Change back to comp.basis
        else % Rigorously should have implemented backtrack:
             % collapse has occured:
             % find collapse time to within specified tolerance
             % ------------------------------------------------
             % Rigorously should have implemented backtrack:
             % collapse has occured:
             % find collapse time to within specified tolerance
             % ------------------------------------------------
            t_prev = tstep_qt(index);
            t_final = tstep_qt(index) + dt_qt;
            %r1;
            ii = 0;
            while ii < 5
                ii = ii + 1;
                t_guess = t_prev + (log(norm2_prev/r1)/log(norm2_prev/norm2_unpsi)) * (t_final - t_prev);
                %t_guess - t_prev
                unpsi_guess = expm(-1i*(t_guess - t_prev)*H_eff/hbar)*unpsi_prev;
                %norm2_guess = norm(unpsi_prev)^2
                norm2_guess = norm(unpsi_guess)^2;
                if abs(r1 - norm2_guess) < 0.001*r1  %error tolerance
                    break
                elseif (norm2_guess < r1)
                    t_final = t_guess;
                    norm2_unpsi = norm2_guess;
                else
                    t_prev = t_guess;
                    unpsi_prev = unpsi_guess;
                    norm2_prev = norm2_guess;
                end
            end
            %r2 = rand([1 1]);
            r2 = r(2);
            %condition = zeros(1, (nevaltruc^2 - nevaltruc + 1));
            condition = zeros(1, length(dp));
            cumsumlist = cumsum(dp);
            for m = 1:length(dp)
                condition(m) = (r2 < cumsumlist(m)/dpj);
            end
            %Find Lindblad operators corresponding to "jump"
            k = find(condition,1,'first');
            Lk = lindbladsearch(k,natom,v,e,nevaltruc);
            psicb = Lk*psicb/norm(Lk*psicb);
            %jlist = [jlist k];

            %more timestep control
            if ~(t_guess >= tstep_qt(index) && t_guess <= tstep_qt(index) + dt_qt)
                t_guess = tstep_qt(index) + dt_qt;
            end

            if t_guess > tf
                t_guess = tf;
            end
            %Change back to comp.basis
            %unpsi = v*unpsicb;
            psi = v*psicb;
            %Change back to comp.basis
                       
            Hs = -1e9.*A_sp1(t_guess./tf).*(sX_1) + 1e9.*B_sp1(t_guess./tf).*((-1).*((1/4).*sZ_1));
            [V,D] = eig(full(Hs));
            %[V,D] = eigs(Hs,18,'sa');
            if ~issorted(diag(D))
                %[V,D] = eigs(Hs,18,'sa');
                [V,D] = eig(full(Hs));
                [D,I] = sort(diag(D));
                D = diag(D);
                V = V(:, I);
            end
            for ii = 1:nevaltruc
                v(:,ii) = sparse(V(:,ii));
                e(ii) = sparse(D(ii,ii));
            end
            
            Hsd = v'*Hs*v;
            psicb = v'*psi;
            
            X = bsxfun(@minus, e', e);
            %[sortedOutput,~,~] = unique(reshape(X,[1 nevaltruc^2]));
            [sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1);%AbsTol
            w_unique = length(sortedOutput);
            dp = zeros(1, w_unique*4);
            
            %gamma0 = 2.*g.^2.*pi.*beta.^(-1);
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


            dp(1) = (psicb'*L01')*(L01*psicb)*dt_qt;


            
            H_eff = Hsd - (1i*hbar/2)*(L01'*L01);
            pdx = 4+1;
            for w = sortedOutput(sortedOutput>0)
                %pdx = pdx + 1;
                %gamma = (2*pi*g^2*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));   
                gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));    
                if isnan(gamma) || isinf(gamma)
                    %gamma = 2.*g.^2.*pi.*beta.^(-1);
                    gamma = gsq2pi*betainv;
                end
                %[b, a] = ind2sub(size(X), find(X == w));
                [b, a] = ind2sub(size(X), find(abs(X - w)<= 0.1));% AbsTol
                Lpcomponents1 = sparse(nevaltruc,nevaltruc); 


                for s = 1:length(b)
                    matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));


                    
                    Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);


                  %L = L + sqrt(natom*gamma)*nfactor*Lcomponent;
                    Lpcomponents1 = Lpcomponents1 + Lpcomponent1;


                end
                Lncomponents1 = Lpcomponents1';



                Lp1 = sqrt(gamma)*Lpcomponents1;



                Ln1 = sqrt(gamma*exp(-beta*w))*Lncomponents1;



                dp(pdx) = (psicb'*Lp1')*(Lp1*psicb)*dt_qt;



                dp(pdx+1) = (psicb'*Ln1')*(Ln1*psicb)*dt_qt;


                              
                pdx = pdx + 2;
                %dp(pdx) = (psi'*Lp')*(Lp*psi)*dt_qt;
                H_eff = H_eff - (1i*hbar/2)*(Lp1'*Lp1)...
                    - (1i*hbar/2)*(Ln1'*Ln1);
            end


            U_eff = expm(-1i*((tstep_qt(index) + dt_qt)-t_guess)*H_eff/hbar);
            %unpsi = U_eff*psi;
            unpsicb = U_eff*psicb;
            psicb = unpsicb/norm(unpsicb);
            
            
            %Change back to comp.basis
            unpsi = v*unpsicb;
            psi = v*psicb;
            %Change back to comp.basis   
            
            %r1 = rand([1 1]); % Draw a new random number    
            r = rand([1 2]); % Draw new random numbers   
        end
    end
    %fidelitylistall = cat(1, fidelitylistall, fidelitylist);
    fidelitylistmat_qt(n, :) = fidelitylist_qt;
%     jlistmat_qt(n,:) = jlist_qt;
%     unpsinormlistmat_qt(n,:) = unpsinormlist_qt;
end

if ntraj == 1
    fidelitylistavg_qt = fidelitylistmat_qt;
else
    fidelitylistavg_qt = mean(fidelitylistmat_qt);
end   
eptime = toc

figure(1)
h = plot(tstep_qt, fidelitylistavg_qt,'-','LineWidth',2); 
xlabel('$t$','Interpreter','latex','FontSize',25)
xlim([0 tf])
ylabel('Ground State Population','FontSize',25)
title(['tf: ' num2str(tf) 's, Number of trajectories: ' num2str(ntraj)],'FontSize',8)
set(gca,'FontSize',20)
print -dpdf qt1qubitexample1_linear

figure(2)
h = plot(tstep_qt./tf, fidelitylistavg_qt,'-','LineWidth',2); 
xlabel('$s$','Interpreter','latex','FontSize',25)
ylabel('Ground State Population','FontSize',25)
title(['tf: ' num2str(tf) 's, Number of trajectories: ' num2str(ntraj)],'FontSize',8)
set(gca,'FontSize',20)
print -dpdf qt1qubitexample2_linear


txt1 = sprintf('qt1qubitexample_linear_fidelity_tf%d_ntraj%d.txt',tf,ntraj);
fid1 = fopen(txt1,'w');
fprintf(fid1,'%13d %8d\n',[tstep_qt;fidelitylistavg_qt]);
fclose(fid1);

fid2 = fopen('qt4qubitexample_eptime.txt','w');
fprintf(fid2,'%d\n',eptime);
fclose(fid2);




delete(pool)
%delete(gcp('nocreate'))


% Subfunction for single-trajectory of deterministic evolution ODE 
% Fast: Constant Hamiltonian
function unpsidot = dunpsifast(unpsi,hbar,H_eff)
% No jump
unpsidot = -1i*hbar*H_eff*unpsi;



function L = lindbladsearch(k,natom,v,e,nevaltruc)
sZ = [1 0; 0 -1];
neval = 2^natom;
%nevaltruc = 2^natom;
unit = speye(2);

z = ceil(k/1);

beta = (1/2.6)*1e-9;
wc = 8*pi*1e9;
gsq2pi = (1.27)*1e-4*2*pi;
betainv = 2.6*1e9;

sZ_1 = sZ;


A_1 = sZ_1;

X = bsxfun(@minus, e.', e);
%[sortedOutput,~,~] = unique(reshape(X,[1 nevaltruc^2]));
[sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1);%AbsTol


if z == 1
   column = ceil(length(sortedOutput)/2);   
elseif logical(mod(z,2))
   column = ceil(length(sortedOutput)/2) - (z-1)/2;
else
   column = ceil(length(sortedOutput)/2) + z/2;
end



w = sortedOutput(column);
%gamma = (2*pi*g^2*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));  
gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));   
if isnan(gamma) || isinf(gamma)
    %gamma = 2.*g.^2.*pi.*beta.^(-1);
    gamma = gsq2pi*betainv;
end
%[b, a] = ind2sub(size(X), find(X == w));
[b, a] = ind2sub(size(X), find(abs(X - w)< 0.1));% AbsTol
L = zeros(nevaltruc,nevaltruc); 
A = A_1;
for s = 1:length(b)
  matrixelement = v(:,a(s))'*A*v(:,b(s));
  Lcomponent = matrixelement*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
  L = L + Lcomponent;
end
L = sqrt(gamma)*L;

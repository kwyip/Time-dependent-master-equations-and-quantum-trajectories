%kawayip@usc.edu
%8-qubit_chain
function aqt_4qubit_reverse_annealing
% cluster = parallel.cluster.Generic;
% set(cluster,'JobStorageLocation', '/home/rcf-proj2/ky/kawayip/research/fourqubitgadget_linear/qt');
% set(cluster,'HasSharedFilesystem', true);
% set(cluster,'IntegrationScriptsLocation','/usr/usc/matlab/R2018b/SlurmIntegrationScripts');
% cluster.AdditionalProperties.SlurmArgs='--time=23:59:59';

cluster = get_LOCAL_cluster('/home/rcf-proj2/ky/kawayip/research/fourqubitgadget_reverse_annealing/ssol/qt');

pool=parpool(cluster,15)

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
% c0 = [1; 0];
% c1 = [0; 1];

plus = [1/sqrt(2); 1/sqrt(2)];
zero0 = [1; 0];
one1 = [0; 1];
%psi0 = sparse(kron(kron(kron(plus,plus),plus),plus));
%psi0 = sparse(kron(kron(kron(one1,zero0),zero0),one1));
psi0 = sparse(kron(kron(kron(zero0,one1),one1),zero0));
%psi0 = sparse(kron(kron(kron(zero0,zero0),zero0),zero0));
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
tf = 5e-6;
hbar = 1;
tic
% Quantum trajectory by using waiting time distribution
ntraj = 1000;

step = 1000;
sstar = 0.525;
inslist1 = linspace(1,sstar,step*2/5);
inslist2 = linspace(sstar,sstar,step*1/5+1);
inslist3 = linspace(sstar,1,step*2/5);
inslist = [inslist1, inslist2, inslist3];

dt_qt = tf/step;
tstep_qt = 0:dt_qt:tf;

%matrix to store fidelity: traj in row, time in col
fidelitylistmat_qt = zeros(ntraj, numel(tstep_qt));
unpsinormlistmat_qt = zeros(ntraj, numel(tstep_qt));
jlistmat_qt = zeros(ntraj, 1000);


parfor n = 1:ntraj
    psi = psi0;
    unpsi = psi0;
    r = rand([1 2]);
    fidelitylist_qt = zeros(1, numel(tstep_qt));
    v = zeros(neval,nevaltruc);
    e = zeros(1,nevaltruc);
    jlist_qt = zeros(1, 1000);
    jorder = 1;
    unpsinormlist_qt = zeros(1, numel(tstep_qt));
    for index = 1:numel(tstep_qt)
        %index
        %clear dp
%         Hs = -1e9.*(1-tstep_qt(index).*tf.^(-1)).*(sX_1+sX_2+sX_3+sX_4) + 1e9.*tstep_qt(index).*tf.^(-1).*(((-1).*sZ_1+(0.95).*sZ_2+(0.95).*sZ_3+(-1).*sZ_4) + ...
%             ((-1).*sZsZII + (-1).*IsZsZI + (-1).*IIsZsZ));
        Hs = -1e9.*A_sp1(inslist(index))/2.*(sX_1+sX_2+sX_3+sX_4) + 1e9.*B_sp1(inslist(index))/2.*(((-1).*sZ_1+(0.95).*sZ_2+(0.95).*sZ_3+(-1).*sZ_4) + ...
            ((-1).*sZsZII + (-1).*IsZsZI + (-1).*IIsZsZ));
        
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
        L02 = sparse(nevaltruc,nevaltruc);
        L03 = sparse(nevaltruc,nevaltruc);
        L04 = sparse(nevaltruc,nevaltruc);

        for s = 1:length(b0)
            matrixelement1 = v(:,a0(s))'*A_1*v(:,b0(s));      %j<->a, i<->b in paper
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

        dp(1) = (psicb'*L01')*(L01*psicb)*dt_qt;
        dp(2) = (psicb'*L02')*(L02*psicb)*dt_qt;
        dp(3) = (psicb'*L03')*(L03*psicb)*dt_qt;
        dp(4) = (psicb'*L04')*(L04*psicb)*dt_qt;


        H_eff = Hsd - (1i*hbar/2)*(L01'*L01+L02'*L02+L03'*L03+L04'*L04);
        pdx = 4+1;
        
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
            Lpcomponents2 = sparse(nevaltruc,nevaltruc);
            Lpcomponents3 = sparse(nevaltruc,nevaltruc);
            Lpcomponents4 = sparse(nevaltruc,nevaltruc);

            for s = 1:length(b)
                matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));
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

            
            dp(pdx) = (psicb'*Lp1')*(Lp1*psicb)*dt_qt;
            dp(pdx+1) = (psicb'*Lp2')*(Lp2*psicb)*dt_qt;
            dp(pdx+2) = (psicb'*Lp3')*(Lp3*psicb)*dt_qt;
            dp(pdx+3) = (psicb'*Lp4')*(Lp4*psicb)*dt_qt;

            
            dp(pdx+4) = (psicb'*Ln1')*(Ln1*psicb)*dt_qt;
            dp(pdx+5) = (psicb'*Ln2')*(Ln2*psicb)*dt_qt;
            dp(pdx+6) = (psicb'*Ln3')*(Ln3*psicb)*dt_qt;
            dp(pdx+7) = (psicb'*Ln4')*(Ln4*psicb)*dt_qt;

            
            pdx = pdx + 8;
            H_eff = H_eff - (1i*hbar/2)*(Lp1'*Lp1 + Lp2'*Lp2 + Lp3'*Lp3 + Lp4'*Lp4)...
                - (1i*hbar/2)*(Ln1'*Ln1 + Ln2'*Ln2 + Ln3'*Ln3 + Ln4'*Ln4);
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
        
        unpsinormlist_qt(index) = unpsinorm;
        
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
            jlist_qt(jorder) = k;
            jorder = jorder + 1;

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
                       
%             Hs = -1e9.*(1-tstep_qt(index).*tf.^(-1)).*(sX_1+sX_2+sX_3+sX_4) + 1e9.*tstep_qt(index).*tf.^(-1).*(((-1).*sZ_1+(0.95).*sZ_2+(0.95).*sZ_3+(-1).*sZ_4) + ...
%                 ((-1).*sZsZII + (-1).*IsZsZI + (-1).*IIsZsZ));
            Hs = -1e9.*A_sp1(inslist(index))/2.*(sX_1+sX_2+sX_3+sX_4) + 1e9.*B_sp1(inslist(index))/2.*(((-1).*sZ_1+(0.95).*sZ_2+(0.95).*sZ_3+(-1).*sZ_4) + ...
                ((-1).*sZsZII + (-1).*IsZsZI + (-1).*IIsZsZ));
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

            dp(1) = (psicb'*L01')*(L01*psicb)*dt_qt;
            dp(2) = (psicb'*L02')*(L02*psicb)*dt_qt;
            dp(3) = (psicb'*L03')*(L03*psicb)*dt_qt;
            dp(4) = (psicb'*L04')*(L04*psicb)*dt_qt;

            
            H_eff = Hsd - (1i*hbar/2)*(L01'*L01+L02'*L02+L03'*L03+L04'*L04);
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
                Lpcomponents2 = sparse(nevaltruc,nevaltruc);
                Lpcomponents3 = sparse(nevaltruc,nevaltruc);
                Lpcomponents4 = sparse(nevaltruc,nevaltruc);

                for s = 1:length(b)
                    matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));
                    matrixelement2 = v(:,a(s))'*A_2*v(:,b(s));
                    matrixelement3 = v(:,a(s))'*A_3*v(:,b(s));
                    matrixelement4 = v(:,a(s))'*A_4*v(:,b(s));

                    
                    Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
                    Lpcomponent2 = matrixelement2*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
                    Lpcomponent3 = matrixelement3*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
                    Lpcomponent4 = matrixelement4*sparse(a(s),b(s),1,nevaltruc,nevaltruc);

                  %L = L + sqrt(natom*gamma)*nfactor*Lcomponent;
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


                dp(pdx) = (psicb'*Lp1')*(Lp1*psicb)*dt_qt;
                dp(pdx+1) = (psicb'*Lp2')*(Lp2*psicb)*dt_qt;
                dp(pdx+2) = (psicb'*Lp3')*(Lp3*psicb)*dt_qt;
                dp(pdx+3) = (psicb'*Lp4')*(Lp4*psicb)*dt_qt;


                dp(pdx+4) = (psicb'*Ln1')*(Ln1*psicb)*dt_qt;
                dp(pdx+5) = (psicb'*Ln2')*(Ln2*psicb)*dt_qt;
                dp(pdx+6) = (psicb'*Ln3')*(Ln3*psicb)*dt_qt;
                dp(pdx+7) = (psicb'*Ln4')*(Ln4*psicb)*dt_qt;

                              
                pdx = pdx + 8;
                %dp(pdx) = (psi'*Lp')*(Lp*psi)*dt_qt;
                H_eff = H_eff - (1i*hbar/2)*(Lp1'*Lp1 + Lp2'*Lp2 + Lp3'*Lp3 + Lp4'*Lp4)...
                    - (1i*hbar/2)*(Ln1'*Ln1 + Ln2'*Ln2 + Ln3'*Ln3 + Ln4'*Ln4);
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
    jlistmat_qt(n,:) = jlist_qt;
    unpsinormlistmat_qt(n,:) = unpsinormlist_qt;
end

if ntraj == 1
    fidelitylistavg_qt = fidelitylistmat_qt;
else
    fidelitylistavg_qt = mean(fidelitylistmat_qt);
end   
eptime = toc

figure(1)
h = plot(tstep_qt, fidelitylistavg_qt,'-');
xlabel('$t$','Interpreter','latex')
xlim([0 tf])
ylabel('Ground State Population')
title(['tf: ' num2str(tf) ', Number of trajectories: ' num2str(ntraj)])
print -dpdf qt4qubitexample1_reverse

figure(2)
h = plot(tstep_qt./tf, fidelitylistavg_qt,'-');
xlabel('$s$','Interpreter','latex')
ylabel('Ground State Population')
title(['tf: ' num2str(tf) ', Number of trajectories: ' num2str(ntraj)])
print -dpdf qt4qubitexample2_reverse

% figure(3)
% h = plot(tstep_qt./tf, unpsinormlistmat_qt,'-');
% xlabel('$s$','Interpreter','latex')
% ylabel('unpsinorm','Interpreter','latex')
% title(['tf: ' num2str(tf) ', Number of trajectories: ' num2str(ntraj)])
% print -dpdf qt4qubitexample3
txt1 = sprintf('qt4qubitexample_reverse_fidelity_tf%d_ntraj%d.txt',tf,ntraj);
fid1 = fopen(txt1,'w');
fprintf(fid1,'%13d %8d\n',[tstep_qt;fidelitylistavg_qt]);
fclose(fid1);

fid2 = fopen('qt4qubitexample_reverse_eptime.txt','w');
fprintf(fid2,'%d\n',eptime);
fclose(fid2);

fid3 = fopen('qt4qubitexample_reverse_jlistmat_cb.txt','w');
% fprintf(fid3,'%d\n',jlistmat);
% fclose(fid3);
[~, cols] = size(jlistmat_qt');
x = repmat('%d\t',1,(cols-1));
fprintf(fid3,[x,'%d\n'],jlistmat_qt);
fclose(fid3);

% fid5 = fopen('qt4qubitexample_unpsinormlistmat_cb.txt','w');
% [~, cols] = size(unpsinormlistmat_qt');
% x = repmat('%d\t',1,(cols-1));
% fprintf(fid5,[x,'%d\n'],unpsinormlistmat_qt);
% fclose(fid5);


delete(pool)
%delete(gcp('nocreate'))


% Subfunction for single-trajectory of deterministic evolution ODE 
% Fast: Constant Hamiltonian
function unpsidot = dunpsifast(unpsi,hbar,H_eff)
% No jump
unpsidot = -1i*hbar*H_eff*unpsi;



%subfunction to search for lindblad
function L = lindbladsearch(k,natom,v,e,nevaltruc)
sZ = [1 0; 0 -1];
z = ceil(k/4);

beta = (1/2.6)*1e-9;
wc = 8*pi*1e9;
gsq2pi = (1.27)*1e-4*2*pi;
betainv = 2.6*1e9;

sZ_1 = kron(sZ, speye(2^(natom-1)));
sZ_2 = kron(kron(speye(2),sZ),speye(2^(natom-2)));
sZ_3 = kron(kron(speye(2^2),sZ),speye(2^(natom-3)));
sZ_4 = kron(speye(2^(natom-1)),sZ);

A_1 = sZ_1;
A_2 = sZ_2;
A_3 = sZ_3;
A_4 = sZ_4;



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
gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));
if isnan(gamma) || isinf(gamma)
    %gamma = 2.*g.^2.*pi.*beta.^(-1);
    gamma = gsq2pi*betainv;
end
[b, a] = ind2sub(size(X), find(abs(X - w)< 0.1));% AbsTol
L = zeros(nevaltruc,nevaltruc); 
if mod(k,4) == 1
    A = A_1;
elseif mod(k,4) == 2
    A = A_2;
elseif mod(k,4) == 3
    A = A_3;
else
    A = A_4;
end

for s = 1:length(b)
  matrixelement = v(:,a(s))'*A*v(:,b(s));
  Lcomponent = matrixelement*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
  L = L + Lcomponent;
end
L = sqrt(gamma)*L;
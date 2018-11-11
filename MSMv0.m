function [optG membership rateResults resource] = MSMv0(Q, SNRdBth, NRB, Rv)

maxSEgrouping = 1; % 1: MaxSE, 2: MQS, 3: MUS
ConvexOptimization = true; % true: Convex, false: equal

T = 16; % number of tiles
V = 2; % number of videos
I = T*V; % number of tiles x number of vides
Qu = 4; % number of representations
U = zeros(I,Qu); % marginal utility
Sa_aerial = [0.0504    0.0888    0.1147    0.0766    0.2932    0.4214    0.5408    0.3543    0.3604    0.5566     0.7570    0.3085    0.1186    0.1748    0.1949    0.0745];
Sa_pole = [0.0620    0.0736    0.0951    0.0669    0.3177    0.5244    0.4228    0.2449    0.3727    0.6355    0.4120    0.1967    0.0793    0.1177    0.1391    0.0453];
Sa_gas = [0.0950    0.1487    0.1060    0.0596    0.3756    0.6837    0.5878    0.2934    0.3829    0.5482    0.3120    0.1751    0.0543    0.1034    0.0826    0.0513];
Sa = [Sa_aerial Sa_gas];

[Vn Nv] = size(Rv);
tempA = max(Rv')'; A = 1e2; 
alpha1 = 1./log(tempA/A); beta1 = tempA/A;
% utility calculation
for uu=1:length(alpha1)
    oriu(uu,:) = alpha1(uu).*log(beta1(uu).*Rv(uu,:)./tempA(uu)); 
end
% marginal utility and marginal cost
for ll=1:Nv
    if ll==1, U(:,ll) = oriu(:,ll); C(:,ll) = Rv(:,ll);
    else U(:,ll) = oriu(:,ll)-oriu(:,ll-1); C(:,ll) = Rv(:,ll) - Rv(:,ll-1); end 
end % marginal utility

SNRth = 10.^(SNRdBth/10);
c = 12*7*[0.1523 0.2344 0.3770 0.6016 0.8770 1.1758 1.4766 1.9141 2.4063 2.7305 3.3223 3.9023 4.5234 5.1152 5.5547]/(0.5e-3);


count = 0; % number of possible minimum SNR
fmar = 0.91; % FEC error margin
sigma = 10;
alpha = 0.0001;
testG = 5;

Nue = length(Q);

%% Spectral Efficiency
for ii=1:Nue
    temp(ii,:) = 1-erf(-(Q(ii)-SNRdBth)/(sigma*sqrt(2)));
    tempSE(ii,:) = fmar*c.*temp(ii,:)/2; 
    [SE(ii) MCS(ii)] = max(tempSE(ii,:));
end


%% Grouping
switch maxSEgrouping
    case 1
        % Grouping Algorithm - MaxSE
        maxG = 5; indexSE(1) = 1;
        for g = 1:maxG
            if g==1
                minSE(g) = SE(1);
                tempN = Nue;
                N{g} = Nue;
            else
                tempSE2 = minSE(g-1)*(1:tempN)+SE(indexSE(g-1):(indexSE(g-1)+tempN-1)).*(tempN:-1:1);
                %plot(tempSE2)
                [maxSE id] = max(tempSE2);
                indexSE(g) = indexSE(g-1) + id;
                minSE(g) = SE(indexSE(g));
                for gg=1:g-1
                    N{g}(gg) = indexSE(gg+1)-indexSE(gg);
                end
                N{g}(g) = Nue - indexSE(g);
                tempN = tempN-id;
            end
        end
    case 2
        % Grouping Algorithm - MQS
        maxG = 5; indexSE(1) = 1;
        for g = 1:maxG
            if g==1
                minSE(g) = SE(1);
                tempN = Nue;
                N{g} = Nue;
            else
                %tempSE2 = minSE(g-1)*(1:tempN)+SE(indexSE(g-1):(indexSE(g-1)+tempN-1)).*(tempN:-1:1);
                tempSE2 = mean(SE(indexSE(g-1):(indexSE(g-1)+tempN-1)));
                %plot(tempSE2)
                [maxSE id] = min(abs(SE(indexSE(g-1):(indexSE(g-1)+tempN-1))-tempSE2));
                indexSE(g) = indexSE(g-1) + id;
                minSE(g) = SE(indexSE(g));
                for gg=1:g-1
                    N{g}(gg) = indexSE(gg+1)-indexSE(gg);
                end
                N{g}(g) = Nue - indexSE(g);
                tempN = tempN-id;
            end
        end
    case 3
        % Grouping Algorithm - MUS
        maxG = 5; indexSE(1) = 1;
        for g = 1:maxG
            if g==1
                minSE(g) = SE(1);
                tempN = Nue;
                N{g} = Nue;
            else
                tempSE2 = mean(SE(indexSE(g-1):(indexSE(g-1)+tempN-1)));
                [maxSE id] = min(abs(SE(indexSE(g-1):(indexSE(g-1)+tempN-1))-tempSE2));
                indexSE(g) = indexSE(g-1) + id;
                minSE(g) = SE(indexSE(g));
                for gg=1:g-1
                    N{g}(gg) = indexSE(gg+1)-indexSE(gg);
                end
                N{g}(g) = Nue - indexSE(g);
                tempN = tempN-id;
            end
        end  
end

%% Iteration - 1
% Rate selection algorithm
tempU = 0; G = 0; a = 0; b = 0; groupU = 0; totalU = 0; convexU = 0; iter = 0; dU = 0; equalRB = 0; equalU = 0; minNrb = 0; 
%temp = 0; 
tempA = 0; tempDU = 0; tempMin = 0; 
while G<testG % multicasting
%while G==0 % broadcasting
    G = G + 1; dU = zeros(G,NRB); convexNrb = zeros(1,G);
    state = zeros(1,Vn);
    stateBest = 0;
    groupU = 0;
    for g = 1:G
        % initial resource allocation for TBRS
        nRB = floor(NRB/G);
        if nRB<1, break; end
        % grouping
        [state a(g) b(g) initc(g)] = TBRS5(nRB, minSE(g), U, Sa, Rv, state);
        stateG{G}(g,:) = state;
%             fity = a(g)*log(b(g)*((1:NRB)));
%             figure(g);
%             plot(1:NRB, fity); hold on;
        groupU(g) = 0;
        for kk=1:I, 
            if state(kk)>0, groupU(g) = groupU(g)+oriu(kk,state(kk)).*Sa(kk); end
        end
        %plot(nRB, groupU(g),'*'); axis([0 NRB 0 10]);
        groupU(g) = groupU(g)*N{G}(g)/Nue;
    end
    equalU = groupU;
    equalState = state;
    equalRB = ones(1,G)*NRB/G;
    bestNrb{G}(1:G) = equalRB;
    %equal
    totalU(G) = sum(groupU); % utility achieve by Equal resource allocation
    % resource allocation algorithm (convex)
    iter(G) = 0; 
    
    %% Iteration - 2
    if G>1
        convexU = groupU; % initial utilty
        stateBest = stateG{G};
        if ConvexOptimization % performs Convex optimization
            state = zeros(1,I);
            stateConvex = zeros(G,I);
            % dU/dNRB
            for g = 1:G
                dU(g,:) = a(g)*N{G}(g)./(1:NRB);
            end
            %plot(dU')
            tempMin = 10000; tempMin2 = 100000;
            % optimal condition search algorithm
            for nn = 1:NRB
                state = zeros(1,I);
                lambda(1) = dU(1,nn);
                minNrb(1) = nn;
                for g = 2:G
                    tempDU = abs(dU(g,:)-lambda(1));
                    tempDU = tempDU + 1e10*(tempDU<0);
                    [lambda(g) minNrb(g)] = min(tempDU);
                end
                if sum(minNrb)>NRB
                    convexNrb(1) = minNrb(1);
                    convexNrb(2:end) = floor((NRB-minNrb(1))*(minNrb(2:end))/sum(minNrb(2:end)));
                    if tempMin2>=abs(sum(lambda(1)-lambda(2:end)))
                        tempMin2 = abs(sum(lambda(1)-lambda(2:end)));
                        %optNrb(end) = minNrb(end)-sum(minNrb)+NRB;
                    else
                        if ~exhaustiveSearch
                            break;
                        end
                    end                        
                    %break;
                    %optNrb
                else
                    if tempMin>=abs(sum(lambda(1)-lambda(2:end))) || tempMin>10
                        tempMin = abs(sum(lambda(1)-lambda(2:end)));
                        %optNrb(end) = minNrb(end)-sum(minNrb)+NRB;
                    else
                        if ~exhaustiveSearch
                            break;
                        end
                    end
                    convexNrb = minNrb;
                end
                exhU = 0;
                for g = 1:G
                    [state a(g) b(g) initc(g)] = TBRS5(convexNrb(g), minSE(g), U, Sa, Rv, state);
                    stateConvex(g,:) = state;
                    exhU(g) = 0;
                    for kk=1:I, 
                        if state(kk)>0, exhU(g) = exhU(g)+oriu(kk,state(kk)).*Sa(kk); end
                    end
                    %plot(nRB, groupU(g),'*'); axis([0 NRB 0 10]);
                    exhU(g) = exhU(g)*N{G}(g)/Nue;
                end
                if sum(exhU)>sum(groupU)
                    groupU = exhU;
                    bestNrb{G}(1:G) = convexNrb;
                    stateG{G} = stateConvex;
                end
            end
            if sum(convexU)<sum(groupU)
                convexU = groupU;
                stateBest = stateG{G};
                iter(G) = iter(G) + 1;
            else
                groupU = convexU;
                stateG{G} = stateBest;
                break;
            end
        end
    end
    if sum(groupU)>tempU
        tempU = sum(groupU);
    else
        %break;
        %optG(NRB) = G-1;
    end
    if totalU(G)<sum(groupU)
        totalU(G) = sum(groupU);
    else
        totalU(G) = totalU(G);
        iter(G) = 0;
    end
end

[outputU optG] = max(totalU);
membership = indexSE;
rateResults = stateG{optG};
resource = bestNrb{optG};

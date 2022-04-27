close all;
load('s_box.mat');

no_bits = 8;
fb=2;%no.faulty bits
sigma = 1.5;
no_traces = 10000;
range = 2^no_bits-1;

K=0:range;
A=0:no_bits;
B=0:no_bits;

cardK = length(K);
cardA = length(A);
cardB = length(B);
pr_uni = 1/cardK;

% PART1: JOINT DISTRIBUTION ATTACK ON UNPROTECTED WITH ML DISTINGUISHER

% PART1A: Compute the theoretical distributions


Distribution_Ineff = zeros(cardB, cardK);
Distribution_Faulty = zeros(cardB, cardK);
Distribution_Eff = zeros(cardB, cardK);
Distribution_joint = zeros(cardA, cardB, cardK);


for k=K
    for a=0:255
        for f=0:2^(fb)-1
            % compute the discrete part of the leakage function 
            faultvalue=bitxor((256-(2^fb)),f);
            b=bitand(a,faultvalue);
            hc=Hamming(bitxor(s_box(a+1),k));
            hf =Hamming(bitxor(s_box(b+1),k));
            % compute the joint distribution for every k
            if (b)==(a)
                Distribution_Ineff(hf+1, k+1) = Distribution_Ineff(hf+1, k+1) + pr_uni;
            else
                Distribution_Faulty(hf+1, k+1) = Distribution_Faulty(hf+1, k+1) + pr_uni;
                Distribution_Eff(hc+1, k+1) = Distribution_Eff(hc+1, k+1) + pr_uni;
                Distribution_joint( hc+1, hf+1, k+1) = Distribution_joint( hc+1, hf+1, k+1) + (pr_uni);
            end
        end
    end
end

% PART1B: ML computation

% simulate the leakage using HW plus noise model
% we assume that we know the POIs involved in the leakage

m = randi(range+1,no_traces,1)-1;
Key = 10;
fault_vector= bitxor((256-(2^fb)),randi(((2^fb)),no_traces,1)-1);
m_faulty=bitand(m,fault_vector);
c = reshape((bitxor(s_box(m+1), Key)),no_traces,1);
c_faulty=reshape((bitxor(s_box(m_faulty+1), Key)),no_traces,1);


% here linking the leakage to the hamming weight is straightforward
L = Hamming(c) + normrnd(0,sigma,no_traces,1);
L_F = Hamming(c_faulty) + normrnd(0,sigma,no_traces,1);
%-------------------------------------------------------------------------

% estimation of sigma 
sigma_c = std(L);
sigma_f = std(L_F);

%-------------------------------------------------------------------------

% technique to compute ML distinguisher

score_Ineff = zeros(cardK,1);
score_Eff = zeros(cardK,1);
score_Faulty = zeros(cardK,1);
score_joint = zeros(cardK,1);
Key_rank_ineff=zeros(no_traces,1);
Key_rank_Eff=zeros(no_traces,1);
Key_rank_Faulty=zeros(no_traces,1);
Key_rank_joint=zeros(no_traces,1);
ni=0;
ne=0;
est=zeros(256,1);
estj=zeros(256,1);

for i=1:no_traces
    for k=K
        summation_Ineff = 0;
        summation_Eff = 0;
        summation_Faulty = 0;
        summation_joint = 0;
        if c(i)==c_faulty(i)
            ni=ni+1;
            for h = A
                 term1i = 1/sigma_c*sqrt(2*pi) * exp(-0.5 * ((L(i) - h)/sigma_c)^2 );
                 term2i = Distribution_Ineff(h+1,k+1);
                 full_termi = term1i * term2i;
                 summation_Ineff = summation_Ineff + full_termi;
            end
            est(k+1)=Distribution_Ineff(Hamming(c(i))+1,k+1)+est(k+1);
            if ni>256
                score_Ineff(k+1) = score_Ineff(k+1) + log(summation_Ineff);
            else
                score_Ineff(k+1) =  log(summation_Ineff);
            end
                
        else
            ne=ne+1;
            for h = A
                 term1e = 1/sigma_c*sqrt(2*pi) * exp(-0.5 * ((L(i) - h)/sigma_c)^2 );
                 term2e = Distribution_Eff(h+1,k+1);
                 full_terme = term1e * term2e;
                 summation_Eff = summation_Eff + full_terme;
                for h_f = B
                    if h==0
                        term1f = 1/sigma_f*sqrt(2*pi) * exp(-0.5 * ((L_F(i) - h_f)/sigma_f)^2 );
                        term2f = Distribution_Faulty(h_f+1,k+1);
                        full_termf = term1f * term2f;
                        summation_Faulty = summation_Faulty + full_termf;
                    end
                    term1j = 1/sigma_c*sqrt(2*pi) * exp(-0.5 * ((L(i) - h)/sigma_c)^2 ) * 1/sigma_f*sqrt(2*pi) * exp(-0.5 * ((L_F(i) - h_f)/sigma_f)^2 );
                    term2j = Distribution_joint(h+1,h_f+1,k+1);
                    full_termj = term1j * term2j;
                    summation_joint = summation_joint + full_termj;
                end
            end
            if ne>256
                score_Eff(k+1) = score_Eff(k+1) +log( summation_Eff);
                score_Faulty(k+1) = score_Faulty(k+1) + log(summation_Faulty);
                score_joint(k+1) = score_joint(k+1) +log( summation_joint);
            else
                score_Eff(k+1) =log(summation_Eff);
                score_Faulty(k+1) =log(summation_Faulty);
                score_joint(k+1) =log(summation_joint);
            end
            estj(k+1)=Distribution_joint(Hamming(c(i))+1, Hamming(c_faulty(i))+1,k+1)+estj(k+1);
                
        end
    end
    Key_rank_ineff(i)=sum(score_Ineff(:)>= score_Ineff(Key+1));
    Key_rank_Eff(i)=sum(score_Eff(:)>= score_Eff(Key+1));
    Key_rank_Faulty(i)=sum(score_Faulty(:)>= score_Faulty(Key+1));
    Key_rank_joint(i)=sum(score_joint(:)>= score_joint(Key+1));
end
ni=ni/256;
ne=ne/256;
score_Ineff
% max(score_Ineff)-min(score_Ineff)
% max(score_Eff)-min(score_Eff)
% max(score_joint)-min(score_joint)

hold on
plot( Key_rank_ineff(1:10000),'LineWidth',1)
% plot( Key_rank_Eff)
plot( Key_rank_Faulty(1:10000),'LineWidth',1)
plot( Key_rank_joint(1:10000),'LineWidth',1)

hold off
legend('sifa','sefFault','joint')
















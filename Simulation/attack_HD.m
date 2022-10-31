close all;
load('s_box.mat');

no_bits =8;
fb=4;%no.faulty bits
sigma = 0.7;
no_traces = 110000;
range = 2^no_bits-1;

K=0:range;
A=0:no_bits;
B=0:no_bits;

cardK = length(K);
cardA = length(A);
cardB = length(B);
pr_uni = 1/cardK;



Distribution_SIFA = zeros(cardB, cardK);
Distribution_SEFA1 = zeros(cardB, cardK);
Distribution_SEFA2 = zeros(cardA, cardB, cardK);

for k=K
    for a=0:255
        for f=0:2^(fb)-1
            % compute the discrete part of the leakage function 
            faultvalue=bitxor((256-(2^fb)),f);
            b=bitand(a,faultvalue);
            hc=Hamming(bitxor(s_box(a+1),k));
            hf =Hamming(bitxor(s_box(b+1),k));
            hd=Hamming(bitxor(s_box(b+1),s_box(a+1)));
            % compute the joint distribution for every k
            if (b)==(a)
                Distribution_SIFA(hf+1, k+1) = Distribution_SIFA(hf+1, k+1) + pr_uni;
            else
                Distribution_SEFA1(hf+1, k+1) = Distribution_SEFA1(hf+1, k+1) + pr_uni;
                Distribution_SEFA2( hc+1, hd+1, k+1) = Distribution_SEFA2( hc+1, hd+1, k+1) + (pr_uni);
            end
        end
    end
end

%PART1B: ML computation

%simulate the leakage using HW plus noise model
%we assume that we know the POIs involved in the leakage
no_key=20;
    Key_rank_SIFA_n=zeros(no_traces,no_key);%n denotes as a noisy circumstances
    Key_rank_SEFA1_n=zeros(no_traces,no_key);
    Key_rank_SEFA2_n=zeros(no_traces,no_key);
    Key_rank_sifa=zeros(no_traces,no_key);
    Key_rank_sefa1=zeros(no_traces,no_key);
    Key_rank_sefa2=zeros(no_traces,no_key); 
for no_k=1:no_key
    m = randi(range+1,no_traces,1)-1;
    Key = randi(256)-1;
    fault_vector= bitxor((256-(2^fb)),randi(((2^fb)),no_traces,1)-1);
    m_faulty=bitand(m,fault_vector);
    c = reshape((bitxor(s_box(m+1), Key)),no_traces,1);
    c_faulty=reshape((bitxor(s_box(m_faulty+1), Key)),no_traces,1);
    HD=(bitxor(c,c_faulty));
    
    % here linking the leakage to the hamming weight is straightforward
    L = Hamming(c) + normrnd(0,sigma,no_traces,1);
    L_F = Hamming(c_faulty) + normrnd(0,sigma,no_traces,1);
    L_HD=Hamming(HD) + normrnd(0,sigma,no_traces,1);
    %-------------------------------------------------------------------------
    
    % est_sifamation of sigma
    sigma_c =sigma; %std(L);
    sigma_f =sigma; %std(L_F);
    
    %-------------------------------------------------------------------------
    
    % technique to compute ML distinguisher
    
    score_SIFA = zeros(cardK,1);
    score_SEFA1 = zeros(cardK,1);
    score_SEFA2 = zeros(cardK,1);
    ni=0;
    ne=0;
    est_sifa=ones(256,1);%esti,ation of SIFA for noiseless secnario
    est_sefa1=ones(256,1);
    est_sefa2=ones(256,1);
    for i=1:no_traces
        for k=K
            summation_SIFA = 0;
            summation_SEFA1 = 0;
            summation_SEFA2 = 0;
            if c(i)==c_faulty(i)
                ni=ni+1;
                    est_sifa(k+1)=log(Distribution_SIFA(Hamming(c(i))+1,k+1))+est_sifa(k+1);
                for hc = A
                    term1i = (1/sigma_c)*sqrt(2*pi) * exp(-0.5 * ((L(i) - hc)/sigma_c)^2 );
                    term2i = Distribution_SIFA(hc+1,k+1);
                    full_termi = term1i * term2i;
                    summation_SIFA = summation_SIFA + full_termi;
                end
                if ni>256
                    score_SIFA(k+1) = score_SIFA(k+1) + log(summation_SIFA);
                else
                    score_SIFA(k+1) =  log(summation_SIFA);
                end
                
            else
                ne=ne+1;
                    est_sefa1(k+1)=log(Distribution_SEFA1(Hamming(c_faulty(i))+1,k+1))+est_sefa1(k+1);
                    est_sefa2(k+1)=log(Distribution_SEFA2(Hamming(c(i))+1, Hamming(HD(i))+1,k+1))+est_sefa2(k+1);
                for hc = A
                    for hf = B 
                        if hc==0
                            term1f = (1/sigma_f)*sqrt(2*pi) * exp(-0.5 * ((L_F(i) - hf)/sigma_f)^2 );
                            term2f = Distribution_SEFA1(hf+1,k+1);
                            full_termf = term1f * term2f;
                            summation_SEFA1 = summation_SEFA1 + full_termf;
                        end
                        term1jHD = (1/sigma_c)*sqrt(2*pi) * exp(-0.5 * ((L(i) - hc)/sigma_c)^2 ) * (1/sigma_f)*sqrt(2*pi) * exp(-0.5 * ((L_HD(i) - hf)/sigma_f)^2 );
                        term2jHD = Distribution_SEFA2(hc+1,hf+1,k+1);
                        full_termjHD = term1jHD * term2jHD;
                        summation_SEFA2 = summation_SEFA2 + full_termjHD;
                        
                    end
                end
                if ne>256
                    score_SEFA1(k+1) = score_SEFA1(k+1) + log(summation_SEFA1);
                    score_SEFA2(k+1) = score_SEFA2(k+1) +log( summation_SEFA2);
                else
                    score_SEFA1(k+1) =log(summation_SEFA1);
                    score_SEFA2(k+1) =log(summation_SEFA2);
                end
                
            end
        end
        Key_rank_SIFA_n(i,no_k)=sum(score_SIFA(:)>= score_SIFA(Key+1));
        Key_rank_SEFA1_n(i,no_k)=sum(score_SEFA1(:)>= score_SEFA1(Key+1));
        Key_rank_SEFA2_n(i,no_k)=sum(score_SEFA2(:)>= score_SEFA2(Key+1));
%%%########################################################################        
        Key_rank_sifa(i,no_k)=sum(est_sifa(:)>= est_sifa(Key+1));
        Key_rank_sefa1(i,no_k)=sum(est_sefa1(:)>= est_sefa1(Key+1));
        Key_rank_sefa2(i,no_k)=sum(est_sefa2(:)>= est_sefa2(Key+1));
    end
    ni=ni/256;
    ne=ne/256;

end
hold on
plot(mean( Key_rank_SIFA_n(1:no_traces,1:no_key),2),'red','LineWidth',1)
plot( mean (Key_rank_SEFA1_n(1:no_traces,1:no_key),2),'green','LineWidth',1)
plot( mean(Key_rank_SEFA2_n(1:no_traces,1:no_key),2),'blue','LineWidth',1)

hold off
% legend('sifa','sefFault','sefa')

hold on
plot(mean( Key_rank_sifa(1:no_traces,:),2),'red','LineWidth',1)
plot( mean (Key_rank_sefa1(1:no_traces,:),2),'green','LineWidth',1)
plot( mean(Key_rank_sefa2(1:no_traces,:),2),'blue','LineWidth',1)














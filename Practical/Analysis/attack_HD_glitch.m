load('s_box.mat');
no_bits = 8;
% fb=2;%no.faulty bits
sigma = 0;
no_traces = 40000;%95000
range = 2^no_bits-1;

K=0:range;
A=0:no_bits;
B=0:no_bits;

cardK = length(K);
cardA = length(A);
cardB = length(B);
pr_uni = 1/cardK;
no_key=1;
    Key_rank_ineff=zeros(no_traces,no_key);
    Key_rank_Eff=zeros(no_traces,no_key);
    Key_rank_Faulty=zeros(no_traces,no_key);
    Key_rank_joint=zeros(no_traces,no_key);
    Key_rank_joint_HD=zeros(no_traces,no_key);
    Key_rank_ni=zeros(no_traces,no_key);
    Key_rank_ne=zeros(no_traces,no_key);
    Key_rank_nj=zeros(no_traces,no_key);
    Key_rank_nj_HD=zeros(no_traces,no_key); 
for no_k=1:no_key
    m = data_n(1:no_traces,2);%randi(range+1,no_traces,1)-1;
    m_faulty = data_n(1:no_traces,1);
    Key = randi(256)-1;
%     fault_vector= bitxor((256-(2^fb)),randi(((2^fb)),no_traces,1)-1);
%     m_faulty=bitand(m,fault_vector);
    c = reshape((bitxor(s_box(m+1), Key)),no_traces,1);
    c_faulty=reshape((bitxor(s_box(m_faulty+1), Key)),no_traces,1);
    HD=bitxor(c,c_faulty);
    
    % here linking the leakage to the hamming weight is straightforward
    L = Hamming(c) + normrnd(0,sigma,no_traces,1);
    L_F = Hamming(c_faulty) + normrnd(0,sigma,no_traces,1);
    L_HD=Hamming(HD) + normrnd(0,sigma,no_traces,1);
    %-------------------------------------------------------------------------
    sigma=0.25;
    % estimation of sigma
    sigma_c =sigma; %std(L);
    sigma_f =sigma; %std(L_F);
    
    %-------------------------------------------------------------------------
    
    % technique to compute ML distinguisher
    
    score_Ineff = zeros(cardK,1);
    score_Eff = zeros(cardK,1);
    score_Faulty = zeros(cardK,1);
    score_joint = zeros(cardK,1);
    score_joint_HD = zeros(cardK,1);
    ni=0;
    ne=0;
    esti=ones(256,1);
    este=ones(256,1);
    estj=ones(256,1);
    estjHD=ones(256,1);
    for i=1:no_traces
        for k=K
            summation_Ineff = 0;
            summation_Eff = 0;
            summation_Faulty = 0;
            summation_joint = 0;
            summation_joint_HD = 0;
            if c(i)==c_faulty(i)
                ni=ni+1;
                    esti(k+1)=log(Distribution_Ineff(Hamming(c(i))+1,k+1))+esti(k+1);
                for h = A
                    term1i = (1/sigma_c)*sqrt(2*pi) * exp(-0.5 * ((L(i) - h)/sigma_c)^2 );
                    term2i = Distribution_Ineff(h+1,k+1);
                    full_termi = term1i * term2i;
                    summation_Ineff = summation_Ineff + full_termi;
                end
                if ni>256
                    score_Ineff(k+1) = score_Ineff(k+1) + log(summation_Ineff);
                else
                    score_Ineff(k+1) =  log(summation_Ineff);
                end
                
            else
                ne=ne+1;
                    este(k+1)=log(Distribution_Faulty(Hamming(c_faulty(i))+1,k+1))+este(k+1);
                    estj(k+1)=log(Distribution_joint(Hamming(c(i))+1, Hamming(c_faulty(i))+1,k+1))+estj(k+1);
                    estjHD(k+1)=log(Distribution_joint_HD(Hamming(c(i))+1, Hamming(HD(i))+1,k+1))+estjHD(k+1);
                for h = A
                    term1e = (1/sigma_c)*sqrt(2*pi) * exp(-0.5 * ((L(i) - h)/sigma_c)^2 );
                    term2e = Distribution_Eff(h+1,k+1);
                    full_terme = term1e * term2e;
                    summation_Eff = summation_Eff + full_terme;
                    for h_f = B
                        if h==0
                            term1f = (1/sigma_f)*sqrt(2*pi) * exp(-0.5 * ((L_F(i) - h_f)/sigma_f)^2 );
                            term2f = Distribution_Faulty(h_f+1,k+1);
                            full_termf = term1f * term2f;
                            summation_Faulty = summation_Faulty + full_termf;
                        end
                        term1j = (1/sigma_c)*sqrt(2*pi) * exp(-0.5 * ((L(i) - h)/sigma_c)^2 ) * (1/sigma_f)*sqrt(2*pi) * exp(-0.5 * ((L_F(i) - h_f)/sigma_f)^2 );
                        term2j = Distribution_joint(h+1,h_f+1,k+1);
                        full_termj = term1j * term2j;
                        summation_joint = summation_joint + full_termj;
                        term1jHD = (1/sigma_c)*sqrt(2*pi) * exp(-0.5 * ((L(i) - h)/sigma_c)^2 ) * (1/sigma_f)*sqrt(2*pi) * exp(-0.5 * ((L_HD(i) - h_f)/sigma_f)^2 );
                        term2jHD = Distribution_joint_HD(h+1,h_f+1,k+1);
                        full_termjHD = term1jHD * term2jHD;
                        summation_joint_HD = summation_joint_HD + full_termjHD;
                        
                    end
                end
                if ne>256
                    score_Eff(k+1) = score_Eff(k+1) +log( summation_Eff);
                    score_Faulty(k+1) = score_Faulty(k+1) + log(summation_Faulty);
                    score_joint(k+1) = score_joint(k+1) +log( summation_joint);
                    score_joint_HD(k+1) = score_joint_HD(k+1) +log( summation_joint);
                else
                    score_Eff(k+1) =log(summation_Eff);
                    score_Faulty(k+1) =log(summation_Faulty);
                    score_joint(k+1) =log(summation_joint);
                    score_joint_HD(k+1) =log(summation_joint_HD);
                end
                
            end
        end
        Key_rank_ineff(i,no_k)=sum(score_Ineff(:)>= score_Ineff(Key+1));
        Key_rank_Eff(i,no_k)=sum(score_Eff(:)>= score_Eff(Key+1));
        Key_rank_Faulty(i,no_k)=sum(score_Faulty(:)>= score_Faulty(Key+1));
        Key_rank_joint(i,no_k)=sum(score_joint(:)>= score_joint(Key+1));
        Key_rank_joint_HD(i,no_k)=sum(score_joint_HD(:)>= score_joint_HD(Key+1));
%%%########################################################################        
        Key_rank_ni(i,no_k)=sum(esti(:)>= esti(Key+1));
        Key_rank_ne(i,no_k)=sum(este(:)>= este(Key+1));
        Key_rank_nj(i,no_k)=sum(estj(:)>= estj(Key+1));
        Key_rank_nj_HD(i,no_k)=sum(estjHD(:)>= estjHD(Key+1));
    end
    ni=ni/256;
    ne=ne/256;
% score_Ineff
% max(score_Ineff)-min(score_Ineff)
% max(score_Eff)-min(score_Eff)
% max(score_joint)-min(score_joint)
end
hold on
plot(mean( Key_rank_ineff(1:no_traces,:),2),'red','LineWidth',1)
% plot(mean (Key_rank_Eff(1:no_traces,:),2),'LineWidth',1)
% % plot( mean (Key_rank_Faulty(1:no_traces,:),2),'green','LineWidth',1)
% % plot( mean(Key_rank_joint(1:no_traces,:),2),'blue','LineWidth',1)
plot( mean(Key_rank_joint_HD(1:no_traces,:),2),'blue','LineWidth',1)
% 
% % 
hold off
% legend('sifa','sefFault','sefa')
% % savings=zeros(2,no_traces);
% % savings(1,:)=mean( Key_rank_ineff(1:no_traces,:),2);
% % savings(2,:)=mean (Key_rank_Eff(1:no_traces,:),2);
% % savings(3,:)=mean (Key_rank_Faulty(1:no_traces,:),2);
% % savings(4,:)=mean(Key_rank_joint(1:no_traces,:),2);
% % savings(4,:)=mean(Key_rank_joint_HD(1:no_traces,:),2);
% % 
% % csvwrite('file415.csv',savings)
% 
% % hold on
% % plot(mean( Key_rank_ni(1:no_traces,:),2),'red','LineWidth',1)
% % plot(mean (Key_rank_Eff(1:no_traces,:),2),'LineWidth',1)
% % plot( mean (Key_rank_ne(1:no_traces,:),2),'green','LineWidth',1)
% % % plot( mean(Key_rank_nj(1:no_traces,:),2),'blue','LineWidth',1)
% % plot( mean(Key_rank_nj_HD(1:no_traces,:),2),'blue','LineWidth',1)
% % 
% % hold off
legend('sifa','sefa')
% savings=zeros(2,no_traces);
% savings(1,:)=mean( Key_rank_ni(1:no_traces,:),2);
% % savings(2,:)=mean (Key_rank_Eff(1:no_traces,:),2);
% % savings(3,:)=mean (Key_rank_ne(1:no_traces,:),2);
% savings(2,:)=mean(Key_rank_nj_HD(1:no_traces,:),2);
% csvwrite('file4ni7.csv',savings)





% data =  csvread('data.csv');
% data1 =  csvread('data1.csv');
% ne=79384;
% ni=20616;
% input_text=csvread('input_text4.csv');
% correct_c=csvread('correct_c4.csv');
% faulty_c=csvread('faulty_c4.csv');
% no_tempi=26435;
% no_bits = 8;
% % fb=2;%no.faulty bits
% sigma = 3;
% no_traces = 60000;
% range = 2^no_bits-1;
% K=0:range;
% A=0:no_bits;
% B=0:no_bits;
% 
% input_ineff4=csvread('input_ineff4.csv');
% Distribution_Ineff =zeros(9,256);
% Distribution_Faulty = zeros(9,256);
% Distribution_Eff = zeros(9,256);
% Distribution_joint=zeros(9,9,256);
% Distribution_joint_HD=zeros(9,9,256);
% 
% for key=1:256
%     for s=1:no_temp
%         Distribution_Eff(Hamming(correct_c(s,key))+1,key)=Distribution_Eff(Hamming(correct_c(s,key))+1,key)+1;
%         Distribution_Faulty(Hamming(faulty_c(s,key))+1,key)=Distribution_Faulty(Hamming(faulty_c(s,key))+1,key)+1;
%         Distribution_joint(Hamming(correct_c(s,key))+1,Hamming(faulty_c(s,key))+1,key)=Distribution_joint(Hamming(correct_c(s,key))+1,Hamming(faulty_c(s,key))+1,key)+1;
%         Xor=Hamming(bitxor(correct_c(s,key),faulty_c(s,key)));
%         Distribution_joint_HD(Hamming(correct_c(s,key))+1,Hamming(Xor)+1,key)=Distribution_joint_HD(Hamming(correct_c(s,key))+1,Hamming(Xor)+1,key)+1;
%     end
% end
% correct_ci=zeros(no_tempi,256);
% for key=0:255
%     for si=1:no_tempi
%         correct_ci(si,key+1)=bitxor(input_ineff4(si,1),key);
%         Distribution_Ineff(Hamming(correct_ci(si,key+1))+1,key+1)=Distribution_Ineff(Hamming(correct_ci(si,key+1))+1,key+1)+1;
%     end
% end
data_n=zeros(100000,2);
s=1;
i=1;
n=1;
while n<100000
    data_n(n:n+3,:)=data(s:s+3,:);
    s=s+4;
    data_n(n+4,:)=data1(i,:);
    i=i+1;
    n=n+5;
end

































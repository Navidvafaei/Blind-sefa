# Blind-sefa
First of all we should create cipher with a faulty value:
```
 faultvalue=bitxor((256-(2^fb)),round(rand*((2^fb)-1)));
```
This fault value will be bitand in the last round of encryption with the following code:
```
    if i_round==9
        if (faultee>0)
            state(1,1)=bitand(state(1,1),faultvalue);
        end;
    end
```
Since we should compute encryption two times we use **faultee** parameter; when a **faultee** is more than 1 encryption is **faultee**.

 We use **regularfault** function to create faulty and non faulty ciphertext.
 
 **Key_col** is denoted as an array of different keys.
```
function [key_col,cipherc,cipherf]=regularfault(sample_t,key_t,fb)
key_col=cell(1,1000);
cipherc=zeros(16,sample_t,,key_t);
cipherf=zeros(16,sample_t,,key_t);
key_n=0;
    while key_n<key_t
        key(1,1)=round(rand*255)  ;
        key(2,1)=round(rand*255)  ;
        key(3,1)=round(rand*255)  ;
        key(4,1)=round(rand*255)  ;
        key(5,1)=round(rand*255)  ;
        key(6,1)=round(rand*255)  ;
        key(7,1)=round(rand*255)  ;
        key(8,1)=round(rand*255)  ;
        key(9,1)=round(rand*255)  ;
        key(10,1)=round(rand*255) ;
        key(11,1)=round(rand*255) ;
        key(12,1)=round(rand*255) ;
        key(13,1)=round(rand*255) ;
        key(14,1)=round(rand*255) ;
        key(15,1)=round(rand*255) ;
        key(16,1)=round(rand*255) ;
        key_col{key_n} = {dec2hex(key)};
        SS=reshape(key,1,16);
        key_hex=dec2hex(SS);
        [s_box, inv_s_box, w, poly_mat, inv_poly_mat] = aes_init(key_hex);
        if w(41,1)==key_n;
            key_n=key_n+1;
            for i=1:sample_t
                faultee=1;
                plaintext(1,1)=round(rand*255);
                plaintext(2,1)=round(rand*255);
                plaintext(3,1)=round(rand*255);
                plaintext(4,1)=round(rand*255);
                plaintext(5,1)=round(rand*255);
                plaintext(6,1)=round(rand*255);
                plaintext(7,1)=round(rand*255);
                plaintext(8,1)=round(rand*255);
                plaintext(9,1)=round(rand*255);
                plaintext(10,1)=round(rand*255);
                plaintext(11,1)=round(rand*255);
                plaintext(12,1)=round(rand*255);
                plaintext(13,1)=round(rand*255);
                plaintext(14,1)=round(rand*255);
                plaintext(15,1)=round(rand*255);
                plaintext(16,1)=round(rand*255);
                faultvalue=bitxor((256-(2^fb)),round(rand*((2^fb)-1)));
                [ciphertext,] = cipher (plaintext, w, s_box, poly_mat,0,0,faultvalue);
                cipherc(:,i,key_n)=ciphertext;
                [ciphertextf,] = cipher (plaintext, w, s_box, poly_mat,0,faultee,faultvalue);
                cipherf(:,i,key_n)=ciphertextf;
            end
        end
    end
end
```
Then create a template by the following function for the attack

**temp_i** is denoted for ineffective template

**temp_e** is denoted for effective template when output ciphertext is correct

**temp_e_f** is denoted for effective template when output ciphertext is faulty

**temp_joint** is used for a joint of effective fault for both correct and faulty ciphertext

```
function[temp_i,temp_e,temp_e_f,temp_joint]=template(fb)
[s_box,~] = s_box_gen (1);
temp_i=zeros(2^(fb),256,256,2);
temp_e=zeros(2^(fb),256,256,2);
temp_h_e=zeros(9,256,2);
temp_h_e_f=zeros(9,256,2);
temp_h_i=zeros(9,256,2);
temp_h_i_f=zeros(9,256,2);
temp_h_ef=zeros(256,9,9);
for key=1:256
    for inp=1:256
        for i=1:2^(fb)
            faultvalue=bitxor((256-(2^fb)),((i-1)));
            s_out_f=sub_bytes(bitand(faultvalue,inp-1),s_box);
            s_out=sub_bytes(inp-1,s_box);
            if  s_out_f==s_out
                temp_i(i,key,inp,1)=bitxor(s_out,key-1);
                temp_h_i(sum(de2bi(temp_i(i,key,inp,1))>0)+1,key,1)=temp_h_i(sum(de2bi(temp_i(i,key,inp,1))>0)+1,key,1)+1;
                temp_h_i(sum(de2bi(inp))+1,key,2)=temp_h_i(sum(de2bi(inp-1))+1,key,2)+1;
                temp_i(i,key,inp,2)=bitxor(s_out_f,key-1);
                temp_h_i_f(sum(de2bi(temp_i(i,key,inp,2))>0)+1,key,1)=temp_h_i(sum(de2bi(temp_i(i,key,inp,2))>0)+1,key,1)+1;
                temp_h_i_f(sum(de2bi(bitxor(faultvalue,inp-1)))+1,key,2)=temp_h_i_f(sum(de2bi(bitxor(faultvalue,inp-1)))+1,key,2)+1;
            else
                temp_e(i,key,inp,1)=bitxor(s_out,key-1);
                temp_h_e(sum(de2bi(temp_e(i,key,inp,1))>0)+1,key,1)=temp_h_e(sum(de2bi(temp_e(i,key,inp,1))>0)+1,key,1)+1;
                temp_h_e(sum(de2bi(inp))+1,key,2)=temp_h_e(sum(de2bi(inp-1))+1,key,2)+1;
                temp_e(i,key,inp,2)=bitxor(s_out_f,key-1);
                temp_h_e_f(sum(de2bi(temp_e(i,key,inp,2))>0)+1,key,1)=temp_h_e_f(sum(de2bi(temp_e(i,key,inp,2))>0)+1,key,1)+1;
                temp_h_e_f(sum(de2bi(bitxor(faultvalue,inp-1)))+1,key,2)=temp_h_e_f(sum(de2bi(bitxor(faultvalue,inp-1)))+1,key,2)+1;
                temp_h_ef(key,sum(de2bi(temp_e(i,key,inp,1))>0)+1,sum(de2bi(temp_e(i,key,inp,2))>0)+1)= temp_h_ef(key,sum(de2bi(temp_e(i,key,inp,1))>0)+1,sum(de2bi(temp_e(i,key,inp,2))>0)+1)+1;
            end
        end
    end
end
re=temp_h_e_f(:,:,1);
re_check=temp_h_e(:,:,1);
ri=temp_h_i(:,:,1);
re_t=temp_h_ef;
for i=1:256
re(:,i)=re(:,i)/sum(re(:,i));
ri(:,i)=ri(:,i)/sum(ri(:,i));
re_check(:,i)=re_check(:,i)/sum(re_check(:,i));
% re_t(:,i)=re(:,i);
end
temp_i=reshape(ri',256,9);
temp_e=reshape(re_check',256,9);
temp_e_f=reshape(re',256,9);
temp_joint=re_t;

end

```
**Key_col** is hexadecimal, it will change to integer value and denoted as **key_base**

correct ciphertext and fault ciphertext should change to HW. 
```
HW_c=zeros(16,sample_t,key_t);
HW_f=zeros(16,sample_t,key_t);
for i=1:key_t
    for j=1:sample_t
        for z=1:16
            HW_c(z,j,i)=sum(de2bi(cipherc(z,j,i)));
            HW_f(z,j,i)=sum(de2bi(cipherf(z,j,i)));
        end
    end
end
```
Then the HW should be change to noisy one. 

**L7** is denoted as noisy HW of correct ciphertext with variance of 0.7

**L15** is denoted as noisy HW of correct ciphertext with variance of 1.5

**Lf7** is denoted as noisy HW of correct ciphertext with variance of 0.7

**Lf15** is denoted as noisy HW of faulty ciphertext with variance of 1.5

All of the above paramteres have three indexes. For example **L7(i,s,t)**


```
variance=1.5;
W = sqrt(variance).*randn(size(HW_c)); %Gaussian white noise W
L15 =round(HW_c + W);
W = sqrt(variance).*randn(size(HW_f));
Lf15 =round(HW_f + W);
variance=0.7;
W = sqrt(variance).*randn(size(HW_c)); %Gaussian white noise W
L7 = round(HW_c + W);
W = sqrt(variance).*randn(size(HW_f));
Lf7 = round(HW_f + W);
Lf15(Lf15<0)=0;
Lf15(Lf15>8)=8;
L7(L7<0)=0;
L7(L7>8)=8;
L15(L15<0)=0;
L15(L15>8)=8;
Lf7(Lf7<0)=0;
Lf7(Lf7>8)=8;
```
now here we are going to use maximum likelihood to recover the key. 

First initialize
```
key_t=100;%number of different keys that might attacker consider
sample_t=35000;
PR_i_t15=zeros(key_t,sample_t);
PR_i_t7=zeros(key_t,sample_t);
PR_e_t15=zeros(key_t,sample_t);
PR_e_t7=zeros(key_t,sample_t);
PR_e_t_f15=zeros(key_t,sample_t);
PR_e_t_f7=zeros(key_t,sample_t);
PR_t=zeros(key_t,sample_t);
PR_t7=zeros(key_t,sample_t);
PR_t15=zeros(key_t,sample_t);
MA_joint_T15=zeros(key_t,sample_t);
MA_joint_T7=zeros(key_t,sample_t);
samp=100;

```
How we can calculate side channel profiling? we use the mean of output for the side channel attack


![image](https://user-images.githubusercontent.com/30938963/156991164-21bd900b-4677-4e7b-a6dc-0c2f72777084.png)


```

mean_h15=zeros(1,9);
mean_h7=zeros(1,9);
mean_hf15=zeros(1,9);
mean_hf7=zeros(1,9);

for i=1:1
    for j=1:20000
        for key=2:40
            mean_h15(1,HW_c(i,j,key)+1)=L15(i,j,key)+mean_h15(HW_c(i,j,key)+1);
            mean_h7(1,HW_c(i,j,key)+1)=L7(i,j,key)+mean_h7(HW_c(i,j,key)+1);
            mean_hf15(1,HW_f(i,j,key)+1)=Lf15(i,j,key)+mean_hf15(HW_f(i,j,key)+1);
            mean_hf7(1,HW_f(i,j,key)+1)=Lf7(i,j,key)+mean_hf7(HW_f(i,j,key)+1);
        end
    end
end
 for h=1:9
      mean_h15(1,h)=mean_h15(1,h)/sum(sum(HW_c(i,1:j,2:key)+1==h));
      mean_h7(1,h)=mean_h7(1,h)/sum(sum((HW_c(i,1:j,2:key)+1==h)));
      mean_hf15(1,h)=mean_hf15(1,h)/sum(sum(HW_f(i,1:j,2:key)+1==h));
      mean_hf7(1,h)=mean_hf7(1,h)/sum(sum(HW_f(i,1:j,2:key)+1==h));
 end

```
Then here We can calculate the average of key's rank for desired number!! here I calculate for 10 keys out of key_t number.


```
for key_n=1:10
    ni=0;
    ne=0;
    sample=1000;
    Ma_inef15=zeros(256,sample_t);
    Ma_inef7=zeros(256,sample_t);
    Ma_ef15=zeros(256,sample_t);
    Ma_ef7=zeros(256,sample_t);
    Ma_eff15=zeros(256,sample_t);
    Ma_eff7=zeros(256,sample_t);
    Ma_joint15=zeros(256,sample_t);
    Ma_joint7=zeros(256,sample_t);
    for j=1:sample_t
        S_i7=zeros(1,256);
        S_i15=zeros(1,256);
        S_e15=zeros(1,256);
        S_e7=zeros(1,256);
        S_F7=zeros(1,256);
        S_F15=zeros(1,256);
        S_joint7=zeros(256,1);
        S_joint15=zeros(256,1);
        if isequal(cipherc(:,j,key_n),cipherf(:,j,key_n))
            ni=ni+1;
            for i=1:256
                for h=1:9
                        S_i7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(temp_ineff(i,h)*1/256)+S_i7(i);
                        S_i15(i)=normpdf(((L15(1,j,key_n))),h-1,1.5)*(temp_ineff(i,h)*1/256)+S_i15(i);
                end
            end
            if ni==1
                Ma_inef15(:,ni)=log(S_i15);
                Ma_inef7(:,ni)=log(S_i7);
            else
                Ma_inef15(:,ni)=log(S_i15(:))+Ma_inef15(:,ni-1);
                Ma_inef7(:,ni)=log(S_i7(:))+Ma_inef7(:,ni-1);
            end
        else
            ne=ne+1;
            for i=1:256
                for h=1:9
                    S_e7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(temp_eff(i,h)*1/256)+S_e7(i);
                    S_e15(i)=normpdf(((L15(1,j,key_n))),h-1,1.5)*(temp_eff(i,h)*1/256)+S_e15(i);
                    S_F7(i)=normpdf(((Lf7(1,j,key_n))),h-1,0.7)*(temp_F(i,h)*1/256)+S_F7(i);
                    S_F15(i)=normpdf(((Lf15(1,j,key_n))),h-1,1.5)*(temp_F(i,h)*1/256)+S_F15(i);
                    for h2=1:9
                          S_joint7(i,1)=normpdf(((L7(1,j,key_n))),h2-1,0.7)*normpdf(((Lf7(1,j,key_n))),h-1,0.7)*(temp_joint(i,h2,h)*1/256*1/256)+S_joint7(i,1);
                          S_joint15(i,1)=normpdf(((L15(1,j,key_n))),h2-1,1.5)*normpdf(((Lf15(1,j,key_n))),h-1,1.5)*(temp_joint(i,h2,h)*1/256*1/256)+S_joint15(i,1);
                    end
                end
            end
            if ne==1
                Ma_ef15(:,ne)=log(S_e15);
                Ma_ef7(:,ne)=log(S_e7);
                Ma_eff15(:,ne)=log(S_F15);
                Ma_eff7(:,ne)=log(S_F7);
                Ma_joint15(:,ne)=log(S_joint15);
                Ma_joint7(:,ne)=log(S_joint7);
            else
                Ma_ef15(:,ne)=log(S_e15(:))+Ma_ef15(:,ne-1);
                Ma_ef7(:,ne)=log(S_e7(:))+Ma_ef7(:,ne-1);
                Ma_eff15(:,ne)=log(S_F15(:))+Ma_eff15(:,ne-1);
                Ma_eff7(:,ne)=log(S_F7(:))+Ma_eff7(:,ne-1);
                Ma_joint7(:,ne)=log(S_joint7)+Ma_joint7(:,ne-1);
                Ma_joint15(:,ne)=log(S_joint15)+Ma_joint15(:,ne-1);
            end
        end
        if sample==j
           sample=sample+samp;
            Pr_inef15(key_n,j/samp)=sum(Ma_inef15(key_base(key_n)+1,ni)<=Ma_inef15(:,ni));
            Pr_ef15(key_n,j/samp)=sum(Ma_ef15(key_base(key_n)+1,ne)<=Ma_ef15(:,ne));
            Pr_F15(key_n,j/samp)=sum(Ma_eff15(key_base(key_n)+1,ne)<=Ma_eff15(:,ne));
            Pr_inef7(key_n,j/samp)=sum(Ma_inef7(key_base(key_n)+1,ni)<=Ma_inef7(:,ni));
            Pr_ef7(key_n,j/samp)=sum(Ma_ef7(key_base(key_n)+1,ne)<=Ma_ef7(:,ne));
            Pr_F7(key_n,j/samp)=sum(Ma_eff7(key_base(key_n)+1,ne)<=Ma_eff7(:,ne));
            Pr_joint15(key_n,j/samp)=sum(Ma_joint15(key_base(key_n)+1,ne)<=Ma_joint15(:,ne));
            Pr_joint7(key_n,j/samp)=sum(Ma_joint7(key_base(key_n)+1,ne)<=Ma_joint7(:,ne));
        end
    end
end

x=zeros(1,sample_t/samp);
for i=1:sample_t/samp
    x(1,i)=i*samp;
end
key_s=2;
% key_n=10;
hold on
% plot(mean(Pr_e_t(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(mean(Pr_e_t_f(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(mean(Pr_i_t(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(mean(Pr_t(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)

% 
% plot(x,mean(Pr_ef15(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(x,mean(Pr_F15(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(x,mean(Pr_inef15(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(x,mean(Pr_joint15(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)


plot(x,mean(Pr_ef7(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
plot(x,mean(Pr_F7(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
plot(x,mean(Pr_inef7(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
plot(x,mean(Pr_joint7(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)


hold off
legend('sefa','sefaf','sifa','MaximumLintsec')
```

For the sake of clarity we consider a part of code here.

```
if isequal(cipherc(:,j,key_n),cipherf(:,j,key_n))
            ni=ni+1;
            for i=1:256
                for h=1:9
                        S_i7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(temp_i(i,h)*1/256)+S_i7(i);
                        S_i15(i)=normpdf(((L15(1,j,key_n))),h-1,1.5)*(temp_i(i,h)*1/256)+S_i15(i);
                end
            end
            if ni==1
                Ma_inef15(:,ni)=log(S_i15);
                Ma_inef7(:,ni)=log(S_i7);
            else
                Ma_inef15(:,ni)=log(S_i15(:))+Ma_inef15(:,ni-1);
                Ma_inef7(:,ni)=log(S_i7(:))+Ma_inef7(:,ni-1);
            end

```
How we consider profiling a noisy side channel?

![image](https://user-images.githubusercontent.com/30938963/155497712-c1a2999f-ff04-45a1-8ac9-5fe2b7946123.png)

Here is the key:
```
PR_i7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(temp_i(i,h)*1/256)+PR_i7(i);
```
y = normpdf(x,mu,sigma) returns the pdf of the normal distribution with mean mu and standard deviation sigma, evaluated at the values in x.

![image](https://user-images.githubusercontent.com/30938963/154865321-b368e5b5-23b2-4fa8-9723-f01a30f62c19.png)

So We can recover the Key when The HW  is noiseless for **fb=2** Here intersect is using for joint distribution SEFAF is using for HW of Faulty


![image](https://user-images.githubusercontent.com/30938963/154866149-65c12a0e-3fc6-49db-b492-5039f51280f3.png)





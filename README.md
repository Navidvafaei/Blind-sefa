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
Then here We can calculate the average of key's rank for desired number!! here I calculate for 10 keys out of key_t number.


```
for key_n=1:10
    ni=0;
    ne=0;
    sample=100;
    Ma_i_t15=zeros(256,sample_t);
    Ma_i_t7=zeros(256,sample_t);
    Ma_e_t15=zeros(256,sample_t);
    Ma_e_t7=zeros(256,sample_t);
    Ma_e_tf15=zeros(256,sample_t);
    Ma_e_tf7=zeros(256,sample_t);
    Ma_joint15=zeros(256,sample_t);
    joint_e15=zeros(256,1);
    Ma_joint7=zeros(256,sample_t);
    joint_e7=zeros(256,1);
    for j=1:35000
        PR_s=zeros(9,256);
        PR_s15=zeros(9,256);
        PR_s7=zeros(1,256);
        PR_i7=zeros(1,256);
        PR_i15=zeros(1,256);
        PR_e15=zeros(1,256);
        PR_e7=zeros(1,256);
        PR_e_f7=zeros(1,256);
        PR_e_f15=zeros(1,256);
        if isequal(cipherc(:,j,key_n),cipherf(:,j,key_n))
            ni=ni+1;
            for i=1:256
                for h=1:9
                    if  L7(1,j,key_n)>=0 && L7(1,j,key_n)<9
                        PR_i7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(temp_i(i,h)*1/256)+PR_i7(i);
                    elseif h==9
                        PR_i7(i)=1;
                    end
                    if L15(1,j,key_n)>=0 && L15(1,j,key_n)<9
                        PR_i15(i)=normpdf(((L15(1,j,key_n))),h-1,1.5)*(temp_i(i,h)*1/256)+PR_i15(i);
                    elseif h==9
                        PR_i15(i)=1;
                    end
                end
            end
            if ni==1
                Ma_i_t15(:,ni)=log(PR_i15);
                Ma_i_t7(:,ni)=log(PR_i7);
            else
                Ma_i_t15(:,ni)=log(PR_i15(:))+Ma_i_t15(:,ni-1);
                Ma_i_t7(:,ni)=log(PR_i7(:))+Ma_i_t7(:,ni-1);
            end
        else
            ne=ne+1;
            for i=1:256
                for h=1:9
                    if  L7(1,j,key_n)>=0 && L7(1,j,key_n)<9
                        PR_e7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(temp_e(i,h)*1/256)+PR_e7(i);
                    elseif h==9
                        PR_e7(i)=1;
                    end
                    if L15(1,j,key_n)>=0 && L15(1,j,key_n)<9
                        PR_e15(i)=normpdf(((L15(1,j,key_n))),h-1,1.5)*(temp_e(i,h)*1/256)+PR_e15(i);
                     elseif h==9
                        PR_e15(i)=1;
                    end
                    if  Lf7(1,j,key_n)>=0 && Lf7(1,j,key_n)<9
                        PR_e_f7(i)=normpdf(((Lf7(1,j,key_n))),h-1,0.7)*(temp_e_f(i,h)*1/256)+PR_e_f7(i);
                        if L7(1,j,key_n)>=0 && L7(1,j,key_n)<9
                            for h2=1:9
                                  joint_e7(i,1)=normpdf(((L7(1,j,key_n))),h2-1,0.7)*normpdf(((Lf7(1,j,key_n))),h-1,0.7)*(temp_joint(i,h2,h)*1/256*1/256)+joint_e7(i,1);
                            end
                        elseif h==9
                            joint_e7(i,1)=1;
                        end
                    elseif h==9
                        PR_e_f7(i)=1;
                        joint_e7(i,1)=1;
                    end
                    if Lf15(1,j,key_n)>=0 && Lf15(1,j,key_n)<9
                        PR_e_f15(i)=normpdf(((Lf15(1,j,key_n))),h-1,1.5)*(temp_e_f(i,h)*1/256)+PR_e_f15(i);
                        if L15(1,j,key_n)>=0 && L15(1,j,key_n)<9
                            for h2=1:9
                                  joint_e15(i,1)=normpdf(((L15(1,j,key_n))),h2-1,1.5)*normpdf(((Lf15(1,j,key_n))),h-1,1.5)*(temp(i,h2,h)*1/256*1/256)+joint_e15(i,1);
                            end
                        elseif h==9
                            joint_e15(i,1)=1;
                        end
                    elseif h==9
                        PR_e_f15(i)=1;
                        joint_e15(i,1)=1;
                    end
                end
            end
            if ne==1
                Ma_e_t15(:,ne)=log(PR_e15);
                Ma_e_t7(:,ne)=log(PR_e7);
                Ma_e_tf15(:,ne)=log(PR_e_f15);
                Ma_e_tf7(:,ne)=log(PR_e_f7);
                Ma_joint15(:,ne)=log(joint_e15);
                Ma_joint7(:,ne)=log(joint_e7);
            else
                Ma_e_t15(:,ne)=log(PR_e15(:))+Ma_e_t15(:,ne-1);
                Ma_e_t7(:,ne)=log(PR_e7(:))+Ma_e_t7(:,ne-1);
                Ma_e_tf15(:,ne)=log(PR_e_f15(:))+Ma_e_tf15(:,ne-1);
                Ma_e_tf7(:,ne)=log(PR_e_f7(:))+Ma_e_tf7(:,ne-1);
                Ma_joint7(:,ne)=log(joint_e7)+Ma_joint7(:,ne-1);
                Ma_joint15(:,ne)=log(joint_e15)+Ma_joint15(:,ne-1);
            end
        end
        if sample==j
           sample=sample+samp;
            PR_i_t15(key_n,j/samp)=sum(Ma_i_t15(key_base(key_n)+1,ni)<=Ma_i_t15(:,ni));
            PR_e_t15(key_n,j/samp)=sum(Ma_e_t15(key_base(key_n)+1,ne)<=Ma_e_t15(:,ne));
            PR_e_t_f15(key_n,j/samp)=sum(Ma_e_tf15(key_base(key_n)+1,ne)<=Ma_e_tf15(:,ne));
            PR_i_t7(key_n,j/samp)=sum(Ma_i_t7(key_base(key_n)+1,ni)<=Ma_i_t7(:,ni));
            PR_e_t7(key_n,j/samp)=sum(Ma_e_t7(key_base(key_n)+1,ne)<=Ma_e_t7(:,ne));
            PR_e_t_f7(key_n,j/samp)=sum(Ma_e_tf7(key_base(key_n)+1,ne)<=Ma_e_tf7(:,ne));
            MA_joint_T15(key_n,j/samp)=sum(Ma_joint15(key_base(key_n)+1,ne)<=Ma_joint15(:,ne));
            MA_joint_T7(key_n,j/samp)=sum(Ma_joint7(key_base(key_n)+1,ne)<=Ma_joint7(:,ne));
        end
    end
end
x=zeros(1,sample_t/samp);
for i=1:sample_t/samp
    x(1,i)=i*samp;
end
x=zeros(1,sample_t/samp);
for i=1:sample_t/samp
    x(1,i)=i*samp;
end
key_s=1;
key_n=10;
hold on
% plot(mean(PR_e_t(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(mean(PR_e_t_f(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(mean(PR_i_t(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(mean(PR_t(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)

% 
plot(x,mean(PR_e_t15(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
plot(x,mean(PR_e_t_f15(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
plot(x,mean(PR_i_t15(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
plot(x,mean(MA_joint_T15(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)


% plot(x,mean(PR_e_t7(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(x,mean(PR_e_t_f7(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(x,mean(PR_i_t7(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)
% plot(x,mean(MA_joint_T7(key_s:key_n,1:(sample_t/samp))),'LineWidth',1)


hold off
legend('sefa','sefaf','sifa','MaximumLintsec')

```

For the sake of clarity we consider a part of code here.

```
if isequal(cipherc(:,j,key_n),cipherf(:,j,key_n))
            ni=ni+1;
            for i=1:256
                for h=1:9
                    if  L7(1,j,key_n)>=0 && L7(1,j,key_n)<9
                        PR_i7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(temp_i(i,h)*1/256)+PR_i7(i);
                    elseif h==9
                        PR_i7(i)=1;
                    end
                    if L15(1,j,key_n)>=0 && L15(1,j,key_n)<9
                        PR_i15(i)=normpdf(((L15(1,j,key_n))),h-1,1.5)*(temp_i(i,h)*1/256)+PR_i15(i);
                    elseif h==9
                        PR_i15(i)=1;
                    end
                end
            end
            if ni==1
                Ma_i_t15(:,ni)=log(PR_i15);
                Ma_i_t7(:,ni)=log(PR_i7);
            else
                Ma_i_t15(:,ni)=log(PR_i15(:))+Ma_i_t15(:,ni-1);
                Ma_i_t7(:,ni)=log(PR_i7(:))+Ma_i_t7(:,ni-1);
            end


```
How we consider profiling a noisy side channel?
![image](https://user-images.githubusercontent.com/30938963/154936943-ee3ef54a-f24d-40a4-bffd-ab00a31fc09b.png)
Here is the key:
```
PR_i7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(temp_i(i,h)*1/256)+PR_i7(i);
```
y = normpdf(x,mu,sigma) returns the pdf of the normal distribution with mean mu and standard deviation sigma, evaluated at the values in x.

![image](https://user-images.githubusercontent.com/30938963/154865321-b368e5b5-23b2-4fa8-9723-f01a30f62c19.png)

So We can recover the Key when The HW  is noiseless for **fb=2** Here intersect is using for joint distribution SEFAF is using for HW of Faulty


![image](https://user-images.githubusercontent.com/30938963/154866149-65c12a0e-3fc6-49db-b492-5039f51280f3.png)






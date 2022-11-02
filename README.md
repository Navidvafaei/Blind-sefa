# Blind-SEFA


This repository contains the documentation of our simulations and Expriments on Blind-SEFA.We apply **Blind-SEFA** to AES cipher and retrive key without any direct access of plaintext and ciphertext.

<img src="https://user-images.githubusercontent.com/30938963/199470100-81a3244d-20be-4b32-af66-840ce330aaf9.png" alt="Your image title" width="400"/>


---

* [Simulation](https://github.com/Navidvafaei/Blind-sefa#simulation)
  * [Initialization](https://github.com/Navidvafaei/Blind-sefa#initialization)
  * [Fault Profiling](https://github.com/Navidvafaei/Blind-sefa#fault-profiling)
  * [Key Recovery Process](https://github.com/Navidvafaei/Blind-sefa#key-recovery-process)
* [Practical_Attacks](https://github.com/Navidvafaei/Blind-sefa#practical_attack)
  * [Template_Profiling](https://github.com/Navidvafaei/Blind-sefa#template-profiling)
  * [Online_Attack](https://github.com/Navidvafaei/Blind-sefa#online-attack)
  * [Simulated-HW](https://github.com/Navidvafaei/Blind-sefa#simulated-HW)


## Simulation
In this section, the simulation code is explained. 

[AES code](https://nevonprojects.com/aes-source-code-inmatlab/)  is written by J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de.
 We also used some part of [Improved blind side channel code](kpcrypto.net)

### Initialization

``` Sigma``` shows the standard deviation of gaussian noise ```fb``` shows number of bit taht the attacker wants to inject fault.

```matlab
close all;
load('s_box.mat');

no_bits =8;
fb=8;%no.faulty bits
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
```
### Fault Profiling

Fault Distribution during the **Fault profiling** is generated by following code:

```matlab

Distribution_SIFA = zeros(cardB, cardK);
Distribution_SEFA1 = zeros(cardB, cardK);
Distribution_SEFA2 = zeros(cardA, cardB, cardK);

for k=K
    for a=0:255
        for f=0:2^(fb)-1
            % compute the discrete part of the leakage function 
            faultvalue=bitxor((256-(2^fb)),f);
            b=bitand(a,faultvalue);
            hc=Hamming(bitxor(s_box(a+1),k)); %**Hamming Weight of corrosponding Correct Output**
            hf =Hamming(bitxor(s_box(b+1),k));%**Hamming Weight of corrosponding Faulty Output**
            hd=Hamming(bitxor(s_box(b+1),s_box(a+1)));%**Hamming Weight of corrosponding Xored Correct and Faulty Output**
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

```
For example , distribution of SEFA2 for key=57 is depicted as following:


<img src="https://user-images.githubusercontent.com/30938963/199102411-c1b7bb50-c791-44b0-b475-91bbf6be8725.png" alt="Your image title" width="400"/>



### Key Recovery Process

We have simulated by two scenarios:

**noiseless scenario**: All notations include ```sifa```,```sefa1``` and ```sefa2```.

**noiseless scenario**: All notations include ```SIFA_n```,```SEFA1_n``` and ```SEFA2_n```.

Adding gausain noise to HW is shown as follows:
```matlab
L = Hamming(c) + normrnd(0,sigma,no_traces,1);
```
Where ```L``` is a Hamming weight of Correct output which is added by gaussian noise.

In the following code, attacker uses Maximum likeliood to recover the key in a noisy enviroment:
It should be mentioned that in this formula hf is calculated once for SEFA1 and also it is used as distance of HW(c+c') for SEFA2
```matlab

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

```
The general code of the key recovery process for noiseless and nosiy scnarios for SIFA,SEFA1 and SEFA2 is as following:

```matlab

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
    est_sifa=ones(256,1);%estimation of SIFA for noiseless secnario
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
                    for hf = B %in this formula hf is calculated once for SEFA1 and also can be  distance of HW(c+c')
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
```
We can plot the mean of Key_rank as following:

```matlab

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

```



--- 

## Practical_Attack

we used [chipwhisperer](https://rtfm.newae.com) For the practical attack with the following setup (required library can be found at the [Chipwhisperer Github](https://github.com/newaetech/chipwhisperer)):

 <img src="https://user-images.githubusercontent.com/30938963/199026183-dd10d4a7-6fd3-4711-8d65-10e594688304.png" alt="Your image title" width="700"/>

Here, we inject faults in sbox by using trigger high and low.
```C

#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h> // CBC mode, for memset
#include "aes.h"


AES_CONST_VAR uint8_t sbox[256] =   {
  //0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

static uint8_t getSBoxValue(uint8_t num)
{
  return sbox[num];
}
uint8_t get_key(uint8_t *k, uint8_t len)
{
    return 0x00;
}

uint8_t SubBytes(uint8_t* state, uint8_t len)
{
  uint8_t i;
        for(i = 0; i < 16; ++i)
        { 
          trigger_high();
          state[i] = getSBoxValue((state[i]));
          trigger_low();
        }

  simpleserial_put('r', 16, state);
  return 0x00;
}

uint8_t reset(uint8_t *x, uint8_t len)
{
    return 0x00;
}

int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();


    simpleserial_init();
    simpleserial_addcmd('k', 16, get_key);
    simpleserial_addcmd('p', 16, SubBytes);
    simpleserial_addcmd('x', 0, reset);

    while (1)
        simpleserial_get();
}
```
It is worth mentioning that, we can inject faults at the input of AES at 10th round, but for the sake of simplicity we just inject fault in the input of 10th round.


### Template-Profiling
Here, for each input of S-Box, we inject fault to generate fault distribution.

**Glitch Setup**
---
```python
import time
import chipwhisperer.common.results.glitch as glitch

def reboot_flush():            
    scope.io.pdic = False
    time.sleep(0.1)
    scope.io.pdic = \high_z\
    time.sleep(0.1)
    target.flush()
gc = glitch.GlitchController(groups=[\success\, \reset\, \normal\], parameters=[\width\, \offset\, \ext_offset\])

gc.set_range(\width\, 5, 5)
gc.set_range(\offset\, -10, -10)
gc.set_range(\ext_offset\, 355, 355)#10362#10562
scope.glitch.clk_src = \clkgen\ # set glitch input clock
scope.glitch.output = \clock_xor\ # glitch_out = clk ^ glitch
scope.glitch.trigger_src = \ext_single\ # glitch only after scope.arm() called
scope.io.hs2 = \glitch\  # output glitch_out on the clock line
scope.adc.timeout = 0.1
```
**Profiling**
---
```python
RUNS = 400,
data1 = np.zeros(((RUNS, 16, 256)))  
data2= np.zeros(((RUNS, 16, 256)))  
data3= np.zeros(((RUNS, 16, 256)))  
data4= np.zeros(((RUNS, 16, 256)))  
S_out=np.zeros(16)  
coll_n=np.zeros((256 2))  
gc.display_stats()  
for glitch_settings in gc.glitch_values():  
    scope.glitch.width = glitch_settings[0]  
    scope.glitch.offset = glitch_settings[1]  
    scope.glitch.ext_offset = glitch_settings[2]  
    reset_cnt = 0  
    valid_run = 1  
    key  text = ktp.next()  
    for j in range (256):  
        ne=0  
        ni=0  
        text[0]=j  
        for r in range(RUNS):  
            for i in range (16):  
                S_out[i]=((sbox[text[i]]))  
            S_out=S_out.astype(int)  
            trace = cw.capture_trace(scope  target  text  key)  
            if scope.adc.state:  
                gc.add(\reset\  (scope.glitch.width  scope.glitch.offset  scope.glitch.ext_offset))  
                reboot_flush()  
                reset_cnt += 1  
                if reset_cnt >= reset_num:  
                    valid_run = 0  
                    break  
                continue  
            try:  
                trace = cw.capture_trace(scope  target  text  key)  
            except:  
                gc.add(\reset\  (scope.glitch.width  scope.glitch.offset  scope.glitch.ext_offset))  
                reboot_flush()  
                reset_cnt += 1  
                if reset_cnt >= reset_num:  
                    valid_run = 0  
                    break  
                continue  
            if trace is None:  
                gc.add(\reset\  (scope.glitch.width  scope.glitch.offset  scope.glitch.ext_offset))  
                reboot_flush()  
                reset_cnt += 1  
                if reset_cnt >= reset_num:  
                    valid_run = 0  
                    break  
                continue  
            ciphertext = trace.textout  
            if ciphertext is None:  
                gc.add(\reset\  (scope.glitch.width  scope.glitch.offset  scope.glitch.ext_offset))  
                reboot_flush()  
                reset_cnt += 1  
                if reset_cnt >= reset_num:  
                    valid_run = 0  
                    break  
                continue  
            sbox_in = text  
            sbox_out = ciphertext   
        if sbox_out[0] != S_out[0]:  
                gc.add(\success\  (scope.glitch.width  scope.glitch.offset  scope.glitch.ext_offset))  
                if sbox_out[1] !=[] :#and sbox_out[1] !=0  
                    data2[ne, :, j] =[sbox_out[0] S_out[0]]  
                    ne=ne+1  
                    coll_n[j 0]=ne  
            else:  
                gc.add(\normal\  (scope.glitch.width  scope.glitch.offset  scope.glitch.ext_offset))  
                data1[ni ,: ,j] =[sbox_out[0] S_out[0]]  
                ni=ni+1  
                coll_n[j 1]=ni
```
### Online-Attack
Same as profiling phase, we inject fault with the aformentioned setup to generate same distribution for fault. Then by using Simulated HW, we can recover the key.

### Simulated-HW

Key-recovery for **noiseless scenario** is depicted as following:

 <img src="https://user-images.githubusercontent.com/30938963/199211768-56ac72a0-6439-45fd-b414-bfd662c28f68.png" alt="Your image title" width="400"/>

Here, by using faulty and correct output, HW can be simulated.


```matlab
load('inv_s_box.mat');
data =  csvread('data.csv');
data1 =  csvread('data1.csv');
ne=79384;
ni=20616;
input_text=csvread('input_text4.csv');
correct_c=csvread('correct_c4.csv');
faulty_c=csvread('faulty_c4.csv');
no_tempi=26435;
no_temp=75900;
input_ineff4=csvread('input_ineff4.csv');
Distribution_Ineff =zeros(9,256);
Distribution_Faulty = zeros(9,256);
Distribution_Eff = zeros(9,256);
Distribution_joint=zeros(9,9,256);
Distribution_joint_HD=zeros(9,9,256);

for key=1:256
 for s=1:no_temp
     Distribution_Eff(Hamming(correct_c(s,key))+1,key)=Distribution_Eff(Hamming(correct_c(s,key))+1,key)+1;
     Distribution_Faulty(Hamming(faulty_c(s,key))+1,key)=Distribution_Faulty(Hamming(faulty_c(s,key))+1,key)+1;
     Distribution_joint(Hamming(correct_c(s,key))+1,Hamming(faulty_c(s,key))+1,key)=Distribution_joint(Hamming(correct_c(s,key))+1,Hamming(faulty_c(s,key))+1,key)+1;
     Xor=(bitxor(correct_c(s,key),faulty_c(s,key)));
     Distribution_joint_HD(Hamming(correct_c(s,key))+1,Hamming(Xor)+1,key)=Distribution_joint_HD(Hamming(correct_c(s,key))+1,Hamming(Xor)+1,key)+1;
 end
end
correct_ci=zeros(no_tempi,256);
for key=0:255
for si=1:no_tempi
    correct_ci(si,key+1)=bitxor(s_box(input_ineff4(si,1)+1),key);
    Distribution_Ineff(Hamming(correct_ci(si,key+1))+1,key+1)=Distribution_Ineff(Hamming(correct_ci(si,key+1))+1,key+1)+1;
end
end
data_n=zeros(100000,2);
s=1;
i=1;
n=1;
while n<100000
 data_n(n:n+3,:)=inv_s_box(data(s:s+3,:)+1);
 s=s+4;
 data_n(n+4,:)=inv_s_box(data1(i,:)+1);
 i=i+1;
 n=n+5;
end
```















































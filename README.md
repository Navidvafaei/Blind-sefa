# Blind-SEFA


This repository contains the documentation of our simulations and Expriments on Blind-SEFA

---

* [Simulation](https://github.com/Navidvafaei/Blind-sefa#simulation)
* [Practical_Attacks](https://github.com/Navidvafaei/Blind-sefa#practical_attack)
  * [Simulated_HW](https://github.com/Navidvafaei/Blind-sefa#simulated_hw)


## Simulation
**Random Fault Value**
First of all we should create cipher with a faulty value:
```
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

% % PART1: JOINT DISTRIBUTION ATTACK ON UNPROTECTED WITH ML DISTINGUISHER
% 
% % PART1A: Compute the theoretical distributions


Distribution_Ineff = zeros(cardB, cardK);
Distribution_Faulty = zeros(cardB, cardK);
Distribution_Eff = zeros(cardB, cardK);
Distribution_joint = zeros(cardA, cardB, cardK);
Distribution_joint_HD = zeros(cardA, cardB, cardK);

for k=K
    for a=0:255
        for f=0:2^(fb)-1
            % compute the discrete part of the leakage function 
            faultvalue=bitxor((256-(2^fb)),f);
            b=bitand(a,faultvalue);
            hc=Hamming(bitxor(s_box(a+1),k));
            hf =Hamming(bitxor(s_box(b+1),k));
            hz=Hamming(bitxor(s_box(b+1),s_box(a+1)));
            % compute the joint distribution for every k
            if (b)==(a)
                Distribution_Ineff(hf+1, k+1) = Distribution_Ineff(hf+1, k+1) + pr_uni;
            else
                Distribution_Faulty(hf+1, k+1) = Distribution_Faulty(hf+1, k+1) + pr_uni;
                Distribution_Eff(hc+1, k+1) = Distribution_Eff(hc+1, k+1) + pr_uni;
                Distribution_joint( hc+1, hf+1, k+1) = Distribution_joint( hc+1, hf+1, k+1) + (pr_uni);
                Distribution_joint_HD( hc+1, hz+1, k+1) = Distribution_joint_HD( hc+1, hz+1, k+1) + (pr_uni);
            end
        end
    end
end
```
This fault value will be bitand in the last round of encryption with the following code:

 


## Practical_Attack

 <img src="https://user-images.githubusercontent.com/30938963/199026183-dd10d4a7-6fd3-4711-8d65-10e594688304.png" alt="Your image title" width="700"/>



### Simulated_HW







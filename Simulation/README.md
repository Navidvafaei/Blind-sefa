# Blind-sefa

[AES code](https://nevonprojects.com/aes-source-code-inmatlab/)  is written by J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de


 We also used some part of [Improved blind side channel code](kpcrypto.net)


First of all we should create cipher with a faulty value:
```
 faultvalue=bitxor((256-(2^fb)),round(rand*((2^fb)-1)));
```
This fault value will be bitand in the last round of encryption with the following code:
```







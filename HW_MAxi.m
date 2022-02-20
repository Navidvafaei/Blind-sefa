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
% data_d_h_i=zeros(key_t,9);
% data_d_h_i5=zeros(key_t,9);
% data_d_h_i7=zeros(key_t,9);
% data_d_h_e=zeros(key_t,9);
% data_d_h_e5=zeros(key_t,9);
% data_d_h_e7=zeros(key_t,9);
% data_d_h_e_f=zeros(key_t,9);
% data_d_h_ef5=zeros(key_t,9);
% data_d_h_ef7=zeros(key_t,9);



samp=100;



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
                        PR_i7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(s_i(i,h)*1/256)+PR_i7(i);%*(1/sum(s_i(:,h))))
                    elseif h==9
                        PR_i7(i)=1;
                    end
                    if L15(1,j,key_n)>=0 && L15(1,j,key_n)<9
                        PR_i15(i)=normpdf(((L15(1,j,key_n))),h-1,1.5)*(s_i(i,h)*1/256)+PR_i15(i);
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
                        PR_e7(i)=normpdf(((L7(1,j,key_n))),h-1,0.7)*(s_e(i,h)*1/256)+PR_e7(i);
                    elseif h==9
                        PR_e7(i)=1;
                    end
                    if L15(1,j,key_n)>=0 && L15(1,j,key_n)<9
                        PR_e15(i)=normpdf(((L15(1,j,key_n))),h-1,1.5)*(s_e(i,h)*1/256)+PR_e15(i);
                     elseif h==9
                        PR_e15(i)=1;
                    end
                    if  HW_f7(1,j,key_n)>=0 && HW_f7(1,j,key_n)<9
                        PR_e_f7(i)=normpdf(((HW_f7(1,j,key_n))),h-1,0.7)*(s_e_f(i,h)*1/256)+PR_e_f7(i);
                        if L7(1,j,key_n)>=0 && L7(1,j,key_n)<9
                            for h2=1:9
                                  joint_e7(i,1)=normpdf(((L7(1,j,key_n))),h2-1,0.7)*normpdf(((HW_f7(1,j,key_n))),h-1,0.7)*(s_e_tf(i,h2,h)*1/256*1/256)+joint_e7(i,1);
                            end
                        elseif h==9
                            joint_e7(i,1)=1;
                        end
                    elseif h==9
                        PR_e_f7(i)=1;
                        joint_e7(i,1)=1;
                    end
                    if HW_f15(1,j,key_n)>=0 && HW_f15(1,j,key_n)<9
                        PR_e_f15(i)=normpdf(((HW_f15(1,j,key_n))),h-1,1.5)*(s_e_f(i,h)*1/256)+PR_e_f15(i);
                        if L15(1,j,key_n)>=0 && L15(1,j,key_n)<9
                            for h2=1:9
                                  joint_e15(i,1)=normpdf(((L15(1,j,key_n))),h2-1,1.5)*normpdf(((HW_f15(1,j,key_n))),h-1,1.5)*(s_e_tf(i,h2,h)*1/256*1/256)+joint_e15(i,1);
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
%             PR_s=PR_i7.*PR_e_f; 
%             PR_s15=PR_i5.*PR_e_f5;
%             PR_s7=PR_i7.*PR_e_f7;
            PR_i_t15(key_n,j/samp)=sum(Ma_i_t15(key_base(key_n)+1,ni)<=Ma_i_t15(:,ni));
            PR_e_t15(key_n,j/samp)=sum(Ma_e_t15(key_base(key_n)+1,ne)<=Ma_e_t15(:,ne));
            PR_e_t_f15(key_n,j/samp)=sum(Ma_e_tf15(key_base(key_n)+1,ne)<=Ma_e_tf15(:,ne));
            PR_i_t7(key_n,j/samp)=sum(Ma_i_t7(key_base(key_n)+1,ni)<=Ma_i_t7(:,ni));
            PR_e_t7(key_n,j/samp)=sum(Ma_e_t7(key_base(key_n)+1,ne)<=Ma_e_t7(:,ne));
            PR_e_t_f7(key_n,j/samp)=sum(Ma_e_tf7(key_base(key_n)+1,ne)<=Ma_e_tf7(:,ne));
            MA_joint_T15(key_n,j/samp)=sum(Ma_joint15(key_base(key_n)+1,ne)<=Ma_joint15(:,ne));
            MA_joint_T7(key_n,j/samp)=sum(Ma_joint7(key_base(key_n)+1,ne)<=Ma_joint7(:,ne));
%             PR_t(key_n,j/samp)=sum(Ma_e_t15(key_base(key_n)+1,j)<=Ma_i_t15(:,j));
%             PR_t7(key_n,j/samp)=sum(Ma_e_t15(key_base(key_n)+1,j)<=Ma_i_t15(:,j));
%             PR_t15(key_n,j/samp)=sum(Ma_e_t15(key_base(key_n)+1,j)<=Ma_i_t15(:,j));
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
% legend('sefa','sefaf','sifa','sefa15','sefa15f','sifa15','sefa7','sefa7f','sifa7')












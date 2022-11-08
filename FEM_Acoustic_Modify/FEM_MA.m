function [Ax,Az,M] = FEM_MA(p_loc_all,element_all)

len2 = size(element_all,1);
len = size(p_loc_all,1);
M=zeros(len,len);
Ax=zeros(len,len);
Az=zeros(len,len);
for ee=1:len2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% T Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ement_num=element_all(ee,:);
    ement_loc=p_loc_all(ement_num,:);
    T1=ones(3,3);
    T1(:,1)=ement_loc([1 2 3],1);
    T1(:,2)=ement_loc([1 2 3],2);
    T2=inv(T1);
    T_ksi=T2*[0;1;0];
    T_eta=T2*[0;0;1];
    T_ksi_eta=[T_ksi,T_eta];
    T=inv(T_ksi_eta(1:2,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% N Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N=[2 1 1;1 2 1;1 1 2]/24;
    N1=[1 -1 0;-1 1 0;0 0 0]/2;
    N2=[1 0 -1;-1 0 1;0 0 0]/2;
    N3=[1 -1 0;0 0 0;-1 1 0]/2;
    N4=[1 0 -1;0 0 0;-1 0 1]/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%% m ax az Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=abs(det(T))*N;
    
    N_ax1=T_ksi_eta(1,1)*T_ksi_eta(1,1)*N1;
    N_ax2=T_ksi_eta(1,1)*T_ksi_eta(1,2)*N2;
    N_ax3=T_ksi_eta(1,1)*T_ksi_eta(1,2)*N3;
    N_ax4=T_ksi_eta(1,2)*T_ksi_eta(1,2)*N4;
    ax=abs(det(T))*(N_ax1+N_ax2+N_ax3+N_ax4);
    
    N_az1=T_ksi_eta(2,1)*T_ksi_eta(2,1)*N1;
    N_az2=T_ksi_eta(2,1)*T_ksi_eta(2,2)*N2;
    N_az3=T_ksi_eta(2,1)*T_ksi_eta(2,2)*N3;
    N_az4=T_ksi_eta(2,2)*T_ksi_eta(2,2)*N4;
    az=abs(det(T))*(N_az1+N_az2+N_az3+N_az4);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% M Ax Az Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M(ement_num,ement_num)=M(ement_num,ement_num)+m;
    Ax(ement_num,ement_num)=Ax(ement_num,ement_num)+ax;
    Az(ement_num,ement_num)=Az(ement_num,ement_num)+az;
end

% MM=inv(M);
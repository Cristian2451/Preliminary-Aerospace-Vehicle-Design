function [K] = K_interp(h_b, ts_t,M)

tRatio2 = M(1:25,1:2); %M(1:15,1:2);
tRatio1_5 = M(1:25, 3:4); %M(1:16, 3:4);
tRatio1_25 = M(1:25, 5:6); %M(1:23, 5:6);
tRatio1 = M(1:25, 7:8);%M(1:19, 7:8);
tRatio0_9 = M(1:25, 9:10);%M(1:19, 9:10);
tRatio0_8 = M(1:25, 11:12);%M(1:21, 11:12);
tRatio0_7 = M(1:25, 13:14);%M(1:22, 13:14);
tRatio0_6 = M(1:25, 15:16);%M(1:19, 15:16);
tRatio0_5 = M(1:25, 17:18);%M(1:25, 17:18);

r(:,:,1)= tRatio2;
r(:,:,2) = tRatio1_5;
r(:,:,3)= tRatio1_25;
r(:,:,4)= tRatio1;
r(:,:,5) = tRatio0_9;
r(:,:,6) = tRatio0_8;
r(:,:,7) = tRatio0_7;
r(:,:,8) = tRatio0_6;
r(:,:,9) = tRatio0_5;

arr = [2, 1.5, 1.25,1,0.9,0.8,0.7,0.6,0.5];
value = ts_t;
if ismember(value,arr)
    index = find(arr==value);
    x1 = r(:,1,index);
    x2 = r(:,2,index);
    x1(x1==0)=[];
    x2(x2==0)=[];
    K = interp1(x1,x2 , h_b);
elseif value > 2
    x1 = r(:,1,1);
    x2 = r(:,2,1);
    x1(x1==0)=[];
    x2(x2==0)=[];
    K = interp1(x1,x2 , h_b);
    disp('Your ts/t factor is > 2')
elseif value < 0.5
    x1 = r(:,1,9);
    x2 = r(:,2,9);
    x1(x1==0)=[];
    x2(x2==0)=[];
    K = interp1(x1,x2 , h_b);
    disp('Your ts/t factor is < 0.5')
else
    %disp('interpolating')
    arr_dif = abs(arr - value);
    [a,idx] = min(arr_dif);
    a = arr(idx);
    arr_dif(idx) = 100;
    [b,idx2] = min(arr_dif);
    b=arr(idx2);
    x1 = r(:,1,idx);
    x2 = r(:,2,idx);
    x1(x1==0)=[];
    x2(x2==0)=[];
    
    x3=r(:,1,idx2);
    x3(x3==0)=[];
    x4=r(:,2,idx2);
    x4(x4==0)=[];

    int1 = interp1(x1, x2, h_b);
    int2 = interp1(x3, x4, h_b);
    K = interp1([a,b],[int1,int2],value);
end
end
